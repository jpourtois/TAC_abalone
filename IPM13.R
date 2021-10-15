# IPM implementation of abalone model

# Script for "Modeling the effect of habitat and fishing heterogeneity on the performance of a TAC-regulated fishery" (Pourtois et al., In Review)
# Adapted from script for "Catastrophic Mortality, Allee Effects, and Marine Protected Areas" (Aalto et al., 2019)

# General model
# use this for analysis

library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)
library(tidyverse)
library(ggpubr)

# Set directory

setwd("~/Documents/Stanford/Research /First Year/TAC Summer Project/IPM/abalone code/common")

# library of utility functions

source("myUtil.r")
source("myData.r")
source('myGraphics.r')

# Naming conventions:
# ALL_CAPS: hard-coded constants (includes defaults, type indicators, like PREY_LARGE, PREY_SMALL, etc.)
# lowerCamelCase: variables, functions
# myXYZ: designates global utility functions from my general library (with some historical exceptions)

# Code structure:
#   CONSTANTS
#   Convenience functions
#   Model:
#     runModel():  Top-level function
#     Initialization functions
#     Presimulation function
#     Main spatial simulation function
#     Plotting

########################################
# ABALONE CONSTANTS 
########################################
###?????????????????????  Spatial Arrangement  ???????????????????????????????????????????????

# Block and coast layout
BLOCK_WIDTH  = 500                              # Width of the block in m (away from coastline)
BLOCK_LENGTH = 100                              # Length of the block in m (along coastline)
COAST_LENGTH = 15000                	          # Length of the coastline in m
NUM_BLOCKS   = round(COAST_LENGTH/BLOCK_LENGTH) # Number of blocks
BLOCK_AREA   = BLOCK_WIDTH*BLOCK_LENGTH/10000   # Area of each block in ha
TOTAL_AREA   = BLOCK_AREA*NUM_BLOCKS            # Total area in ha

# Dispersal
DISP_PERC    = 99                               # Define percentage of larvae that are retained between +- dispersal distance
DISP_QUANT   = (1-DISP_PERC/100)/2

# Set Class Limits [25 mm width adjusted to 155 mm]
# Actually, the first size class should include 1 year individuals
# Initial size of larvae in Bardos et al. 2006 is 0.3 mm
# Viana et al. measures 2000 postlarvae (unknown age): 909 +- 130 um
# Leighton et al.1981 uses postlarvae (unknown age) of 1-1.5 mm
# Leighton 1974: postlarvae of 30-40 days are 1.7-2 mm SL
# In collectors at Isla Natividad, larvae at settlment were ~500 um
# Using Bardos(2005) and starting with larvae all of 500um, >90% is between 3.5 and 27 mm
# We could simplify and say that recruitment goes all to 1st size class

###?????????????????????  Demographic Settings  ???????????????????????????????????????????????
MIN_SIZE      = 1     # in mm
MAX_SIZE      = 306   # in mm
MIN_LAND_SIZE = 155   # set Minimum Landing Size
CLASS_WIDTH   = 5     # default to 5mm size classes
NUM_CLASSES   = (MAX_SIZE-MIN_SIZE)/CLASS_WIDTH

LARVAL_SURVIVAL  = 0.00309*1.846
# larval survival adjusted by approx. x2 to account for 
# Allee effect in Rosetto calculation (her observed densities were ~0.015/m2)
# This is somewhat deprecated
LARVAL_ALLEE_ADJ = 1 #1.846;
# approximation from original model
LARVAL_CONST_SIZE = 17.5

############################
# CONVENIENCE FUNCTIONS
############################

# Weight [g] at Length relationship da Shepherd 1998 [with shell]
# to be multiplied by 0.4 to obtain weight [g] without shell [From Natividad Data]
# W = a*L^b
WEIGHT_A  = 2.24 * 10^-5  
WEIGHT_B  = 3.36   
W_NOSHELL = 0.4
wForL <- function(L) { return(WEIGHT_A*(L^(WEIGHT_B))) }     

# Determine proportion dispersing to block at given position
propDispToBlock <- function(position, stdDev) { 
  return(pnorm(BLOCK_LENGTH*(position+1/2),0,stdDev) - 
           pnorm(BLOCK_LENGTH*(position-1/2),0,stdDev))
}

# Calculate dispersal SD for a given random value
getDispSD <- function(randVal) { 
  return(uniroot(function(x)qnorm(DISP_QUANT,mean=0,sd=x)+randVal,lower=0,upper=8000)$root) 
}

#####################
### Constant global vectors
#####################
# Mean class length vector
# Precompute for easy vector operations later
gMeanSizeV = rep(0, NUM_CLASSES)      # initialize below

# Vector giving distribution of new recruits across size classes
gRecruitDistribV = rep(0, NUM_CLASSES) # initialize below

# Size at sexual maturity, Rossetto et al.(2013), not accounting for 1:1 Sex Ratio
# Determines proportion mature at each size
# P(mature) = a / 1+e^(-(size-matSize)/b)

MAT_A    = 1 #0.5  DON't CUT IN HALF, modeling all individuals for Allee effect
MAT_B    = 30.20
MAT_SIZE = 135.99
gPropMatureV = c() # initialize below

# Eggs per individual for each size
# Convert size (length) to weight using wForL
# Eggs only from mature individuals (use gPropMatureV)
# Eggs = a*weight*prop

EGGS_A = 3772      # from Tutschulte 1976
gEggsPerIndV = c() # initialize below

DEF_CAT_MORT = 0.75

# Default **recruit** K param in each block
#DEF_K = 5.43*10^7  # Tuned to unfished density of ~0.4/m2
#DEF_K = 2.59*10^7  # Tuned to unfished density of ~0.2/m2 XXX readj for fixed Allee
DEF_K = 1.29*10^7  # Tuned to unfished density of ~0.2/m2

K_SD = 2*0.6*10^7 

DEF_HIGH_K = 1 # Default is 2 high K zones 
HET_K = rep(rep(c(DEF_K - K_SD, DEF_K + K_SD), each = NUM_BLOCKS/(2*DEF_HIGH_K)), DEF_HIGH_K)
if (length(HET_K) < NUM_BLOCKS) HET_K = c(DEF_K - K_SD,HET_K, DEF_K + K_SD)

HOM_K = rep(DEF_K, NUM_BLOCKS)

# To test dependence of fMSY and catch on K
HIGH_K = rep(DEF_K + K_SD, NUM_BLOCKS)
LOW_K = rep(DEF_K - K_SD, NUM_BLOCKS)

# Parameters for the arima.sim function
AR = 0.90 #0.6
MA = NULL

set.seed(2020) # 2020 for all analysis

# Build coast from parameters AR and MA
corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = AR, ma = MA), n = NUM_BLOCKS)
#Bring to positive
corr2 = corr1 - range(corr1)[1]
corr2[corr2 == 0] = 0.1
# Bring to right scale
corrK = corr2*DEF_K/mean(corr2)
# Get mean to DEF_K
corrK = corrK-mean(corrK)+DEF_K

# Dispersal constants:
DISP_GAMMA_SHAPE = 3
DISP_GAMMA_RATE  = 0.006

# Larval survival constants:
# Gamma and LN tuned to produce ~15x ratio of largest:smallest (after trimming outliers)
# Shepherd 1990
REC_GAMMA_SHAPE  = 2.5
REC_GAMMA_RATE   = 0.5
REC_GB_THRESHOLD = 0.8
REC_GB_GOOD      = 5
REC_GB_BAD       = 0.5
REC_LN_MEANLOG   = 0
REC_LN_SDLOG     = 0.7

# uses the probability of mixed-gender aggregations, based on local density (#/m2)
# Good shorthand
# From Button 2008 thesis, Fig 2.4 and Fig. 5.8
# Changed from 11.33 to 11.6 to set 80% threshold at 0.2/m^2 (well within CI on A value)
AGG_A = 11.6
AGG_B = 1  #1.959  # linear fit gives 1.959, but density~0 should have agg=1

# Fishing aggregation constants
BIOMASS_FACTOR = 0 #0.5 to 2
DIST_FACTOR = 0
CATCHABILITY = 0.1 
TOTAL_EFFORT = -3*log(1 - 0.09)/CATCHABILITY*NUM_BLOCKS #3x number of boats needed for MSY

###????????????????????? INITIALIZATION???????????????????????????????????????????????

##### INITIALIZATION CONSTANTS
NUM_REPLICATIONS    = 5    # number of replicates 
MAX_REPLICATE_GROUP = 999999   # break apart into groups
T_MAX               = 250   # simulation length
T_PRESIM            = 100   # Pre-simulation time
H_PRESIM            = 0.0  # No Pre-simulation harvest right now
INIT_ABUN_PRESIM    = 100  # Default initial density for pre-simulation in each size class (#/ha)
INIT_ABUN_NO_PRESIM = 100  # Default initial density if pre-simulation is skipped (#/ha)

#####  Allee types
ALLEE_PROB = 1
ALLEE_LIN  = 2
ALLEE_EXP  = 3

# Recreate using same randomness
RAND_SEED = 123456
# Collapse criterion: <10% of pre-catastrophe equilibrium
COLLAPSE_THRESHOLD = 0.1

# Fraction of fMSY for non-spatial
MSY_PROP = 0.67
DISCOUNT_RATE = 0.95

# Parameters for the dynamic management strategy
DYNAMIC_DELAY  = -999
DYNAMIC_THRESH = 0.5

# Default MPA parameters
MPA_SIZE = 5
MPA_PROP = 0

# Parameter labels
PARAM_CAT   = "Cat"
PARAM_AGG   = "Agg"
PARAM_K     = "K"
PARAM_DISP  = "Disp"
PARAM_MPA_W = "MPASize"
PARAM_MPA_P = "MPAProp"
PARAM_MSY_P = "MSYProp"
PARAM_CLOSE = "CloseThresh"
PARAM_DISC  = "DiscRate"
PARAM_MULT  = "MultiCat"

initializeGlobalValues <- function(classWidth=CLASS_WIDTH, larvaeConst=FALSE) {
  # initialize NUM_CLASSES
  if (classWidth<=1) MIN_SIZE = 1
  else MIN_SIZE = 5;
  assign("MIN_SIZE", MIN_SIZE, envir=.GlobalEnv);
  assign("CLASS_WIDTH", classWidth, envir=.GlobalEnv);
  NUM_CLASSES   = (MAX_SIZE-MIN_SIZE)/classWidth;
  assign("NUM_CLASSES", NUM_CLASSES, envir=.GlobalEnv);
  # initialize gMeanSizeV
  gMeanSizeV = rep(0, NUM_CLASSES);
  for (i in 1:NUM_CLASSES) { gMeanSizeV[i] = MIN_SIZE + (i-1/2)*CLASS_WIDTH; } # size is average between i and i+1  
  assign("gMeanSizeV", gMeanSizeV, envir=.GlobalEnv);
  # initialize gRecruitDistribV
  gRecruitDistribV = getRecruitDistributionVector(larvaeConst);
  assign("gRecruitDistribV", gRecruitDistribV, envir=.GlobalEnv);
  # initialize gPropMatureV
  gPropMatureV = round(MAT_A/(1+exp(-(gMeanSizeV-MAT_SIZE)/MAT_B)),1)        
  assign("gPropMatureV", gPropMatureV, envir=.GlobalEnv)
  # re-initialize gEggsPerIndV
  # add 0.5 to adjust for 1:1 sex ratio
  gEggsPerIndV = 0.5*EGGS_A*wForL(gMeanSizeV)*gPropMatureV
  assign("gEggsPerIndV", gEggsPerIndV, envir=.GlobalEnv)
  
  # use the same seed for now
  if (RAND_SEED>0) set.seed(RAND_SEED)
}

#################################################################
####                    SETS TO RUN
#################################################################

plotRunModel <- function(resultRunModel){
  
  spatialData = resultRunModel$res
  
  spatialMassA = spatialData$spatialMassA                       # [rep, maxT, NUM_BLOCK]
  spatialSSBA = spatialData$spatialSSBA                         # [rep, maxT, NUM_BLOCK]
  spatialCatchA = spatialData$spatialCatchA                     # [rep, maxT, NUM_BLOCK]
  spatialCatchShellA = spatialData$spatialCatchShellA           # [rep, maxT, NUM_BLOCK]
  spatialDensA = spatialData$spatialDensA                       # [rep, maxT, NUM_BLOCK]
  spatialSSBDensA = spatialData$spatialSSBDensA                 # [rep, maxT, NUM_BLOCK]
  actualEffort = spatialData$actualEffort                       # [rep, maxT, NUM_BLOCK]
  
  spatialMassM    = apply(spatialMassA, c(2,3), mean)           # [maxT, NUM_BLOCK]
  spatialSSBM     = apply(spatialSSBA, c(2,3), mean)            # [maxT, NUM_BLOCK]
  spatialCatchM   = apply(spatialCatchA, c(2,3), mean)          # [maxT, NUM_BLOCK]
  spatialCatchShellM = apply(spatialCatchShellA, c(2,3), mean)  # [maxT, NUM_BLOCK]
  spatialDensM    = apply(spatialDensA, c(2,3), mean)           # [maxT, NUM_BLOCK]
  spatialSSBDensM = apply(spatialSSBDensA, c(2,3), mean)        # [maxT, NUM_BLOCK]
  spatialEffort   = apply(actualEffort, c(2,3), mean)           # [maxT, NUM_BLOCK]
  
  totalSSB = apply(spatialSSBM,1,sum) # [maxT]
  totalCatch = apply(spatialCatchShellM,1,sum) #[maxT]
  
  propHarvested = totalCatch/totalSSB
  
  meansData = resultRunModel$means
  
  # Organize data for spatial plotting
  t_selected = c(1,10,T_MAX-10)
  t_length = length(t_selected)
  
  rm(spatialSummary)
  rm(spatialBind)
  
  for (i in 1:t_length) {
    
    t = t_selected[i]
    spatialBind = data.frame(time = rep(t, NUM_BLOCKS))
    spatialBind$Coast = 1:NUM_BLOCKS
    spatialBind$MassM = spatialMassM[t,]
    spatialBind$SSBM = spatialSSBM[t,]
    spatialBind$CatchM = spatialCatchM[t,]
    spatialBind$DensM = spatialDensM[t,]
    spatialBind$SSBDensM = spatialSSBDensM[t,]
    spatialBind$Effort = spatialEffort[t,]
    
    if (i == 1) { spatialSummary = spatialBind 
    } else { spatialSummary = rbind(spatialSummary,spatialBind) }
    
  } 
  
  
  # Plot spatial evolution 
  
  # Adult density
  thePlot = ggplot(spatialSummary, aes(x = Coast)) + 
    geom_line(aes(y=SSBDensM/10000,group=as.factor(time), color=as.factor(time)))+
    labs(title="", x="Coast (block #)", y="Adult density (#/m^2)", color = "Time (years)") +
    scale_y_continuous(limits=c(0,0.24))
  
  myPPlot(thePlot) 
  
  # Total density
  thePlot = ggplot(spatialSummary, aes(x = Coast)) + 
    geom_line(aes(y=DensM/10000,group=as.factor(time), color=as.factor(time)))+
    labs(title="", x="Coast (block #)", y="Total density (#/m^2)", color = "Time (years)")
  
  myPPlot(thePlot) 
  
  # Catch
  thePlot = ggplot(spatialSummary, aes(x = Coast)) + 
    geom_line(aes(y=CatchM,group=as.factor(time), color=as.factor(time)))+
    labs(title="", x="Coast (block #)", y="Catch (kg)", color = "Time (years)")
  
  myPPlot(thePlot) 
  
  # Effort
  thePlot = ggplot(spatialSummary, aes(x = Coast)) + 
    geom_line(aes(y=Effort,group=as.factor(time), color=as.factor(time)))+
    labs(title="", x="Coast (block #)", y="Effort", color = "Time (years)") + 
    scale_y_continuous(limits=c(0,1.3))
  
  
  myPPlot(thePlot) 
  
  # Plot heat map
  
  plotHeatMap(spatialCatchM, 'Catch (kg)')
  plotHeatMap(spatialSSBDensM/10000, 'Adult density (#/m^2)')
  plotHeatMap(spatialEffort, 'Effort')
  
  
}

function_analysis <- function(h = 0.09, Allee = TRUE, randomK = corrK, rep = 20) {
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, hetK = randomK,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=Allee, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, hetK = randomK,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=Allee, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=h, habitatHOM = TRUE, hetK = randomK,
                                maxT=T_MAX, numReplic=rep, biomass_factor = 0, dist_factor = 0,
                                isAllee=Allee, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST0 = runModel(hProp=h, habitatHOM = FALSE, hetK = randomK,
                            maxT=T_MAX, numReplic=rep, biomass_factor = 0, dist_factor = 0,
                            isAllee=Allee, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST0$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 0
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST0 = runModel(hProp=h, habitatHOM = FALSE, hetK = randomK,
                              maxT=T_MAX, numReplic=rep, biomass_factor = 100, dist_factor = 0,
                              isAllee=Allee, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST0$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 0
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST0 = runModel(hProp=h, habitatHOM = FALSE, hetK = randomK,
                              maxT=T_MAX, numReplic=rep, biomass_factor = 200, dist_factor = 0,
                              isAllee=Allee, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST0$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 0
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapse = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapse$bFactor = as.factor(dt_patch_collapse$bFactor)
  dt_patch_collapse$dFactor = as.factor(dt_patch_collapse$dFactor)
  dt_patch_collapse$habitat = as.factor(dt_patch_collapse$habitat)
  
  return(dt_patch_collapse)
  
}

analysis_all <- function(){
  
  # MSY =  0.09, Catch = 24015.55
  
  # Calculate baseline (no harvest) to determine collapse
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, hetK = corrK,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, hetK = corrK,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=0.09, habitatHOM = TRUE, hetK = corrK,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, hetK = corrK,
                            maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                            isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST0$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 0
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, hetK = corrK,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST0$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 0
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, hetK = corrK,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST0$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 0
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapse = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapse$bFactor = as.factor(dt_patch_collapse$bFactor)
  dt_patch_collapse$dFactor = as.factor(dt_patch_collapse$dFactor)
  dt_patch_collapse$habitat = as.factor(dt_patch_collapse$habitat)
  
  saveMatrixToCSVFile('figure2A_2021', dt_patch_collapse)
  
  # Above MSY
  
  # Case 4
  
  OVER_BIO0_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, hetK = corrK,
                             maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                             isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                             hSDProp=0.1, classicSurv=FALSE,
                             runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                             showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                             limitOutput = TRUE)
  
  spatialAdultMass_4 = OVER_BIO0_DIST0$adultMass
  patch_collapse_4 = matrix(0,length(spatialAdultMass_4[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_4[,1])) {
    patch_collapse_4[i,] = spatialAdultMass_4[i,] < collapse_mass 
  }
  
  collapse_per_rep_4 = rowMeans(1*patch_collapse_4)
  dt_4 = data.frame("Patch collapse" = collapse_per_rep_4)
  dt_4$bFactor = 0
  dt_4$dFactor = 0
  dt_4$AdultMass = apply(spatialAdultMass_4, 1, sum)
  dt_4$catch = OVER_BIO0_DIST0$catch
  dt_4$habitat = 'HET'
  
  
  # Case 5
  
  OVER_BIO100_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, hetK = corrK,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_5 = OVER_BIO100_DIST0$adultMass
  patch_collapse_5 = matrix(0,length(spatialAdultMass_5[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_5[,1])) {
    patch_collapse_5[i,] = spatialAdultMass_5[i,] < collapse_mass 
  }
  
  collapse_per_rep_5 = rowMeans(1*patch_collapse_5)
  dt_5 = data.frame("Patch collapse" = collapse_per_rep_5)
  dt_5$bFactor = 100
  dt_5$dFactor = 0
  dt_5$AdultMass = apply(spatialAdultMass_5, 1, sum)
  dt_5$catch = OVER_BIO100_DIST0$catch
  dt_5$habitat = 'HET'
  
  
  #Case 6
  
  OVER_BIO200_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, hetK = corrK,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_6 = OVER_BIO200_DIST0$adultMass
  patch_collapse_6 = matrix(0,length(spatialAdultMass_6[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_6[,1])) {
    patch_collapse_6[i,] = spatialAdultMass_6[i,] < collapse_mass 
  }
  
  collapse_per_rep_6 = rowMeans(1*patch_collapse_6)
  dt_6 = data.frame("Patch collapse" = collapse_per_rep_6)
  dt_6$bFactor = 200
  dt_6$dFactor = 0
  dt_6$AdultMass = apply(spatialAdultMass_6, 1, sum)
  dt_6$catch = OVER_BIO200_DIST0$catch
  dt_6$habitat = 'HET'
  
  #Case 7
  
  OVER_BIO0_DIST0_HOM = runModel(hProp=0.12, habitatHOM = TRUE, hetK = corrK,
                                 maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                 isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                 hSDProp=0.1, classicSurv=FALSE,
                                 runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                 showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                 limitOutput = TRUE)
  
  
  spatialAdultMass_7 = OVER_BIO0_DIST0_HOM$adultMass
  patch_collapse_7 = matrix(0,length(spatialAdultMass_7[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_7[,1])) {
    patch_collapse_7[i,] = spatialAdultMass_7[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_7 = rowMeans(1*patch_collapse_7)
  dt_7 = data.frame("Patch collapse" = collapse_per_rep_7)
  dt_7$bFactor = 'HOM'
  dt_7$dFactor = 0
  dt_7$AdultMass = apply(spatialAdultMass_7, 1, sum)
  dt_7$catch = OVER_BIO0_DIST0_HOM$catch
  dt_7$habitat = 'HOM'
  
  dt_patch_collapse_overMSY = rbind(dt_4,dt_5,dt_6, dt_7)
  dt_patch_collapse_overMSY$bFactor = as.factor(dt_patch_collapse_overMSY$bFactor)
  dt_patch_collapse_overMSY$dFactor = as.factor(dt_patch_collapse_overMSY$dFactor)
  dt_patch_collapse_overMSY$habitat = as.factor(dt_patch_collapse_overMSY$habitat)
  
  saveMatrixToCSVFile('figure2B_2021', dt_patch_collapse_overMSY)
  
  # Analysis for Figure 1
  
  fig1a = data_frame()
  fig1b = data_frame()
  fig1c = data_frame()
  fig1d = data_frame()
  fig1e = data_frame()
  fig1f = data_frame()
  df = data_frame(coast = 1:150)
  
  for (i in 1:length(spatialAdultMass_1[,1])){
    
    df$coast = 1:150
    df$rep = i
    df$rep = as.factor(df$rep)
    
    df$mass = spatialAdultMass_1[i,]
    fig1a = rbind(fig1a,df)
    
    df$mass = spatialAdultMass_2[i,]
    fig1b = rbind(fig1b,df)
    
    df$mass = spatialAdultMass_3[i,]
    fig1c = rbind(fig1c,df)
    
    df$mass = spatialAdultMass_4[i,]
    fig1d = rbind(fig1d,df)
    
    df$mass = spatialAdultMass_5[i,]
    fig1e = rbind(fig1e,df)
    
    df$mass = spatialAdultMass_6[i,]
    fig1f = rbind(fig1f,df)
    
  }
  
  collapse_thresh_df = tibble(spatialAdultMassEqui)
  collapse_thresh_df = cbind(collapse_thresh_df, block = 1:NUM_BLOCKS)
  
  write.csv(fig1a, 'output/figure1a_2021.csv' )
  write.csv(fig1b, 'output/figure1b_2021.csv' )
  write.csv(fig1c, 'output/figure1c_2021.csv' )
  write.csv(fig1d, 'output/figure1d_2021.csv' )
  write.csv(fig1e, 'output/figure1e_2021.csv' )
  write.csv(fig1f, 'output/figure1f_2021.csv' )
  
  write.csv(collapse_thresh_df, 'output/collapse_thresh_df_2021.csv' )
  
}

# Run analysis for figure 1 and figure 2 
analysis_fig1_2 <- function(){
  
  # MSY =  0.09, Catch = 24015.55
  
  # Calculate baseline (no harvest) to determine collapse
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, hetK = corrK,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  
  #saveMatrixToCSVFile('output/test',spatialAdultMassEqui)
  
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, hetK = corrK,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=0.09, habitatHOM = TRUE, hetK = corrK,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, hetK = corrK,
                            maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                            isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST0$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 0
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, hetK = corrK,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST0$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 0
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, hetK = corrK,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST0$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 0
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapse = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapse$bFactor = as.factor(dt_patch_collapse$bFactor)
  dt_patch_collapse$dFactor = as.factor(dt_patch_collapse$dFactor)
  dt_patch_collapse$habitat = as.factor(dt_patch_collapse$habitat)
  
  #saveMatrixToCSVFile('output/figure2A_2021', dt_patch_collapse)
  write.csv(dt_patch_collapse, 'output/figure2A_2021.csv')
  # Above MSY
  
  # Case 4
  
  OVER_BIO0_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, hetK = corrK,
                             maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                             isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                             hSDProp=0.1, classicSurv=FALSE,
                             runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                             showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                             limitOutput = TRUE)
  
  spatialAdultMass_4 = OVER_BIO0_DIST0$adultMass
  patch_collapse_4 = matrix(0,length(spatialAdultMass_4[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_4[,1])) {
    patch_collapse_4[i,] = spatialAdultMass_4[i,] < collapse_mass 
  }
  
  collapse_per_rep_4 = rowMeans(1*patch_collapse_4)
  dt_4 = data.frame("Patch collapse" = collapse_per_rep_4)
  dt_4$bFactor = 0
  dt_4$dFactor = 0
  dt_4$AdultMass = apply(spatialAdultMass_4, 1, sum)
  dt_4$catch = OVER_BIO0_DIST0$catch
  dt_4$habitat = 'HET'
  
  
  # Case 5
  
  OVER_BIO100_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, hetK = corrK,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_5 = OVER_BIO100_DIST0$adultMass
  patch_collapse_5 = matrix(0,length(spatialAdultMass_5[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_5[,1])) {
    patch_collapse_5[i,] = spatialAdultMass_5[i,] < collapse_mass 
  }
  
  collapse_per_rep_5 = rowMeans(1*patch_collapse_5)
  dt_5 = data.frame("Patch collapse" = collapse_per_rep_5)
  dt_5$bFactor = 100
  dt_5$dFactor = 0
  dt_5$AdultMass = apply(spatialAdultMass_5, 1, sum)
  dt_5$catch = OVER_BIO100_DIST0$catch
  dt_5$habitat = 'HET'
  
  
  #Case 6
  
  OVER_BIO200_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, hetK = corrK,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_6 = OVER_BIO200_DIST0$adultMass
  patch_collapse_6 = matrix(0,length(spatialAdultMass_6[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_6[,1])) {
    patch_collapse_6[i,] = spatialAdultMass_6[i,] < collapse_mass 
  }
  
  collapse_per_rep_6 = rowMeans(1*patch_collapse_6)
  dt_6 = data.frame("Patch collapse" = collapse_per_rep_6)
  dt_6$bFactor = 200
  dt_6$dFactor = 0
  dt_6$AdultMass = apply(spatialAdultMass_6, 1, sum)
  dt_6$catch = OVER_BIO200_DIST0$catch
  dt_6$habitat = 'HET'
  
  #Case 7
  
  OVER_BIO0_DIST0_HOM = runModel(hProp=0.12, habitatHOM = TRUE, hetK = corrK,
                                 maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                 isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                 hSDProp=0.1, classicSurv=FALSE,
                                 runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                 showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                 limitOutput = TRUE)
  
  
  spatialAdultMass_7 = OVER_BIO0_DIST0_HOM$adultMass
  patch_collapse_7 = matrix(0,length(spatialAdultMass_7[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_7[,1])) {
    patch_collapse_7[i,] = spatialAdultMass_7[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_7 = rowMeans(1*patch_collapse_7)
  dt_7 = data.frame("Patch collapse" = collapse_per_rep_7)
  dt_7$bFactor = 'HOM'
  dt_7$dFactor = 0
  dt_7$AdultMass = apply(spatialAdultMass_7, 1, sum)
  dt_7$catch = OVER_BIO0_DIST0_HOM$catch
  dt_7$habitat = 'HOM'
  
  dt_patch_collapse_overMSY = rbind(dt_4,dt_5,dt_6, dt_7)
  dt_patch_collapse_overMSY$bFactor = as.factor(dt_patch_collapse_overMSY$bFactor)
  dt_patch_collapse_overMSY$dFactor = as.factor(dt_patch_collapse_overMSY$dFactor)
  dt_patch_collapse_overMSY$habitat = as.factor(dt_patch_collapse_overMSY$habitat)
  
  #saveMatrixToCSVFile('output/figure2B_2021', dt_patch_collapse_overMSY)
  write.csv(dt_patch_collapse_overMSY, 'output/figure2B_2021.csv')
  
  # Analysis for Figure 1
  
  fig1a = data_frame()
  fig1b = data_frame()
  fig1c = data_frame()
  fig1d = data_frame()
  fig1e = data_frame()
  fig1f = data_frame()
  df = data_frame(coast = 1:150)
  
  for (i in 1:length(spatialAdultMass_1[,1])){
    
    df$coast = 1:150
    df$rep = i
    df$rep = as.factor(df$rep)
    
    df$mass = spatialAdultMass_1[i,]
    fig1a = rbind(fig1a,df)
    
    df$mass = spatialAdultMass_2[i,]
    fig1b = rbind(fig1b,df)
    
    df$mass = spatialAdultMass_3[i,]
    fig1c = rbind(fig1c,df)
    
    df$mass = spatialAdultMass_4[i,]
    fig1d = rbind(fig1d,df)
    
    df$mass = spatialAdultMass_5[i,]
    fig1e = rbind(fig1e,df)
    
    df$mass = spatialAdultMass_6[i,]
    fig1f = rbind(fig1f,df)
    
  }
  
  collapse_thresh_df = tibble(spatialAdultMassEqui)
  collapse_thresh_df = cbind(collapse_thresh_df, block = 1:NUM_BLOCKS)
  
  write.csv(fig1a, 'output/figure1a_2021.csv' )
  write.csv(fig1b, 'output/figure1b_2021.csv' )
  write.csv(fig1c, 'output/figure1c_2021.csv' )
  write.csv(fig1d, 'output/figure1d_2021.csv' )
  write.csv(fig1e, 'output/figure1e_2021.csv' )
  write.csv(fig1f, 'output/figure1f_2021.csv' )
  
  write.csv(collapse_thresh_df, 'output/collapse_thresh_df_2021.csv' )
}

makeFigure0 <- function(){
  
  set.seed(2020)
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.6, ma = MA), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK06 = corrK-mean(corrK)+DEF_K
  
  set.seed(2020)
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.9, ma = MA), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK09 = corrK-mean(corrK)+DEF_K
  
  set.seed(2020)
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.95, ma = MA), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK095 = corrK-mean(corrK)+DEF_K
  
  
  df_06 = data_frame(k = as.numeric(corrK06))
  df_06$corr = 0.6
  df_06$coast = 1:150
  
  df_09 = data_frame(k = as.numeric(corrK09))
  df_09$corr = 0.9
  df_09$coast = 1:150
  
  df_095 = data_frame(k = as.numeric(corrK095))
  df_095$corr = 0.95
  df_095$coast = 1:150
  
  df_corr = rbind(df_06, df_09, df_095)
  
  df_corr$corr = as.factor(df_corr$corr)
  
  levels(df_corr$corr) <- c("AR = 0.6", "AR = 0.9", "AR = 0.95")
  fig0 <- ggplot(df_corr, aes(x = coast, y = k, group = corr, color = corr)) +
    geom_line() +
    scale_y_continuous(trans='log10') +
    labs(x="Coast (block #)", y="Carrying capacity (# recruits)", group = expression(atop("Auto-correlation", "coefficient"))) +
    facet_grid(cols = vars(corr)) +
    theme_minimal() +
    theme(legend.position = 'none') 
  
  myPPlot(fig0)
  
}

makeFigure1 <- function(){
  
  fig1a <- read.csv('output/figure1a_2021.csv')
  fig1a$rep <- as.factor(fig1a$rep)
  fig1b <- read.csv('output/figure1b_2021.csv')
  fig1b$rep <- as.factor(fig1b$rep)
  fig1c <- read.csv('output/figure1c_2021.csv')
  fig1c$rep <- as.factor(fig1c$rep)
  fig1d <- read.csv('output/figure1d_2021.csv')
  fig1d$rep <- as.factor(fig1d$rep)
  fig1e <- read.csv('output/figure1e_2021.csv')
  fig1e$rep <- as.factor(fig1e$rep)
  fig1f <- read.csv('output/figure1f_2021.csv')
  fig1f$rep <- as.factor(fig1f$rep)
  
  collapse_thresh_df <- read.csv('output/collapse_thresh_df_2021.csv')
  
  p1a <- ggplot() +
    geom_line(aes(x = block, y = spatialAdultMassEqui/10000),size = 0.20,collapse_thresh_df) +
    geom_line(aes(x=coast, y=mass/1000, group = rep, colour = rep),size = 0.20,alpha = 1, fig1a) +
    annotate("text",label = "Agg. level = none",x = 115, y = 2.25, color = "black",size = 2.5) +
    ylim(0, 2.3) +
    labs(x="Coast (block #)", y="Adult biomass (tons)") +
    theme_minimal() +
    theme(legend.position = "none",axis.title = element_text(size = 8), axis.text = element_text(size = 7)) 
  
  #myPPlot(p1a)
  
  p1b <- ggplot() +
    geom_line(aes(x=coast, y=mass/1000, group = rep, colour = rep),size = 0.20,alpha = 1, fig1b) +
    geom_line(aes(x = block, y = spatialAdultMassEqui/10000),size = 0.20,collapse_thresh_df) +
    annotate("text",label = "Agg. level = low",x = 115, y = 2.25, color = "black",size = 2.5) +
    ylim(0, 2.3) + 
    labs(x="Coast (block #)", y="") +
    theme_minimal() +
    theme(legend.position = "none",axis.title.x = element_text(size = 8),axis.text = element_text(size = 7))
  #myPPlot(p1b)
  
  p1c <- ggplot() +
    geom_line(aes(x=coast, y=mass/1000, group = rep, colour = rep),size = 0.20,alpha = 1, fig1c) +
    geom_line(aes(x = block, y = spatialAdultMassEqui/10000),size = 0.20,collapse_thresh_df) +
    annotate("text",label = "Agg. level = high",x = 115, y = 2.25, color = "black",size = 2.5) +
    ylim(0, 2.3) + 
    labs(x="Coast (block #)", y="") +
    theme_minimal() +
    theme(legend.position = "none",axis.title.x = element_text(size = 8),axis.text = element_text(size = 7))
  #myPPlot(p1c)
  
  p1d <- ggplot() +
    geom_line(aes(x=coast, y=mass/1000, group = rep, colour = rep),size = 0.20,alpha = 1, fig1d) +
    geom_line(aes(x = block, y = spatialAdultMassEqui/10000),size = 0.20,collapse_thresh_df) +
    annotate("text",label = "Agg. level = none",x = 115, y = 2.25, color = "black",size = 2.5) +
    ylim(0, 2.3) +
    labs(x="Coast (block #)", y="Adult biomass (tons)") +
    theme_minimal() +
    theme(legend.position = "none",axis.title = element_text(size = 8),axis.text = element_text(size = 7))
  #myPPlot(p1d)   
  
  p1e <- ggplot() +
    geom_line(aes(x=coast, y=mass/1000, group = rep, colour = rep),size = 0.20,alpha = 1, fig1e) +
    geom_line(aes(x = block, y = spatialAdultMassEqui/10000),size = 0.20,collapse_thresh_df) +
    annotate("text",label = "Agg. level = low",x = 115, y = 2.25, color = "black",size = 2.5) +
    ylim(0, 2.3) +
    labs(x="Coast (block #)", y="") +
    theme_minimal() +
    theme(legend.position = "none",axis.title.x = element_text(size = 8),axis.text = element_text(size = 7))
  #myPPlot(p1e)
  
  p1f <- ggplot() +
    geom_line(aes(x=coast, y=mass/1000, group = rep, colour = rep),size = 0.20,alpha = 1, fig1f) +
    geom_line(aes(x = block, y = spatialAdultMassEqui/10000),size = 0.20,collapse_thresh_df) +
    annotate("text",label = "Agg. level = high",x = 115, y = 2.25, color = "black",size = 2.5) +
    ylim(0, 2.3) +
    labs(x="Coast (block #)", y="") +
    theme_minimal() +
    theme(legend.position = "none", axis.title.x = element_text(size = 8),axis.text = element_text(size = 7)) 
  #myPPlot(p1f)
  
  #fig1 <- ggarrange(p1a,p1d,p1b,p1e,p1c,p1f, labels = c("A", "B","C","D","E","F"), 
  #ncol = 2, nrow = 3, vjust = 1.3,font.label = list(size = 12))
  
  fig1 <- ggarrange(p1a,p1b,p1c,p1d,p1e,p1f, labels = c("A", "","","B","",""), 
                    ncol = 3, nrow = 2, vjust = 1.3,font.label = list(size = 10))
  
  
  ggsave("fig1.tiff",fig1, units="mm", width=170, height= 110, dpi=300)
  
  myPPlot(fig1)
}

makeFigure2 <- function(){
  
  #setwd("~/Documents/Stanford/First Year/TAC Summer Project/IPM/abalone code/common/output")
  dt_patch_collapse = read.csv('output/figure2A.csv')
  dt_patch_collapse_overMSY = read.csv('output/figure2B.csv')
  #setwd("~/Documents/Stanford/First Year/TAC Summer Project/IPM/abalone code/common")
  
  pa <- ggplot(dt_patch_collapse, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "",title="Exploitation at MSY", y="Extinct patches (%)", fill = 'Habitat') +
    scale_y_continuous(limits = c(0,118), breaks = c(0,25,50,75,100)) + 
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = 0, ymax = 115, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = 0, ymax = 115, alpha = .2) + 
    annotate("text", x = 3.1, y=107, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=107, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, vjust = 2),axis.title.x=element_blank())
  #myPPlot(pa)
  
  pb <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "", y = "",title="Overexploitation", fill = "Habitat") +
    scale_y_continuous(limits = c(0,118), breaks = c(0,25,50,75,100)) + 
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = 0, ymax = 115, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = 0, ymax = 115, alpha = .2) + 
    annotate("text", x = 3.1, y=107, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=107, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, vjust = 2),axis.title.x=element_blank())
  #myPPlot(p)
  
  pc <- ggplot(dt_patch_collapse, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "", y="Adult biomass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x=element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 175, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 175, alpha = .2) + 
    annotate("text", x = 3.1, y=165, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=165, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10),axis.title.x=element_blank())
  #myPPlot(p)
  
  pd <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "", y = "",fill = "Habitat") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -5, ymax = 175, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -5, ymax = 175, alpha = .2) + 
    annotate("text", x = 3.1, y=165, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=165, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10),axis.title.x=element_blank()) +
    ylim(-5, 175)
  #myPPlot(p)
  
  pe <- ggplot(dt_patch_collapse, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs( x="Fishing aggregation level", y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 35, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 35, alpha = .2) + 
    annotate("text", x = 3.1, y=32, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=32, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10))
  
  #myPPlot(p)
  
  
  pf <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=catch/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x="Fishing aggregation level", y = "", fill = "Habitat") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8))+
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 35, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 35, alpha = .2) + 
    annotate("text", x = 3.1, y=32, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=32, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none", axis.title = element_text(size = 10))
  #myPPlot(p)
  
  fig2 <- ggarrange(pa,pb,pc,pd,pe,pf, labels = c("A", "D","B","E","C","F"), 
                    ncol = 2, nrow = 3, vjust = 1.3,font.label = list(size = 12))
  myPPlot(fig2, w = 14, h = 7)
  
  ggsave("fig2.tiff",fig2, units="mm", width=170, height= 200, dpi=300)
  
  
}

statsAnalysis2 <- function(){
  
  # Check normality
  data_HOM <- dt_patch_collapse[dt_patch_collapse$habitat == "HOM",]
  
  shapiro.test(data_HOM$AdultMass) # p: 0.2283
  shapiro.test(data_HOM$catch) # p: 0.2249
  
  data_HET_0 <- dt_patch_collapse[dt_patch_collapse$bFactor == "0",]
  
  
  shapiro.test(data_HET_0$AdultMass) # p:0.3987
  shapiro.test(data_HET_0$catch) # p: 0.3005
  
  t.test(data_HOM$AdultMass, data_HET_0$AdultMass, paired = TRUE) #p: 1.708e-12
  t.test(data_HOM$catch, data_HET_0$catch, paired = TRUE)# p: 2.562e-12
  
  tapply(dt_patch_collapse$AdultMass, dt_patch_collapse$bFactor,sd)
  
  pre_anova = dt_patch_collapse[dt_patch_collapse$habitat == "HET",]
  pre_anova$reps = rep(1:20,3)
  #pre_anova$bFactor = as.factor(pre_anova$bFactor)
  #pre_anova$reps = factor(pre_anova$reps)
  
  pre_anova %>% anova_test(dv = AdultMass, wid = reps, within = bFactor)
  
  pairwise.t.test(data = pre_anova, x = pre_anova$AdultMass, g = pre_anova$bFactor, paired = TRUE,p.adjust.method = "bonferroni")
  
  relevel(pre_anova$bFactor, ref="0")
  lm_model <- lm(data = pre_anova, AdultMass/1000 ~ bFactor)
  summary(lm_model)
  lme_model <- lme(data = pre_anova, AdultMass/1000 ~ factor(bFactor), random = ~1|reps)
  summary(lme_model)
  ranef(lme_model)
  anova(lme_model)
  intervals(lme_model)
  
}

figure3 <- function(){
  
  
  # Figure 3a
  
  dt_patch_collapse = read.csv('output/figure3A.csv')
  
  p3a <- ggplot(dt_patch_collapse, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x=expression(atop("Fishing heterogeneity level", "")), y="Adult biomass (tons)",title = "Exploitation at MSY", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = 50, ymax = 125, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = 50, ymax = 125, alpha = .2) + 
    annotate("text", x = 3.1, y=120, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=120, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10),plot.title = element_text(hjust = 0.5, vjust = 2)) +
    ylim(50, 125)
  
  #myPPlot(p)
  
  # Figure 3b
  
  dt_patch_collapse_overMSY = read.csv('output/figure3B.csv')
  
  p3b <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x=expression(atop("Fishing heterogeneity level", "")), y="Adult biomass (tons)",title = "Overexploitation", fill = "Habitat") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = 50, ymax = 125, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = 50, ymax = 125, alpha = .2) + 
    annotate("text", x = 3.1, y=120, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=120, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10),plot.title = element_text(hjust = 0.5, vjust = 2)) +
    ylim(50, 125)
  
  #myPPlot(p)
  
  # Figure 3c
  
  MSY_Allee = read.csv('output/figure2A.csv')
  MSY_NoAllee = read.csv('output/figure3A.csv')
  
  # Create MSY data set
  
  # Part 1
  dt_MSY_Allee = data.frame("type" = c("None (hom. habitat)","None (het. habitat)","Low","High"))
  
  dt_MSY_Allee$Allee = 1
  
  dt_MSY_Allee$medianCollapse = 0
  dt_MSY_Allee$medianCollapse[1] = median(MSY_Allee$Patch.collapse[MSY_Allee$bFactor == "HOM"])
  dt_MSY_Allee$medianCollapse[2] = median(MSY_Allee$Patch.collapse[MSY_Allee$bFactor == "0"])
  dt_MSY_Allee$medianCollapse[3] = median(MSY_Allee$Patch.collapse[MSY_Allee$bFactor == "100"])
  dt_MSY_Allee$medianCollapse[4] = median(MSY_Allee$Patch.collapse[MSY_Allee$bFactor == "200"])
  
  dt_MSY_Allee$medianMass = 0
  dt_MSY_Allee$medianMass[1] = median(MSY_Allee$AdultMass[MSY_Allee$bFactor == "HOM"])
  dt_MSY_Allee$medianMass[2] = median(MSY_Allee$AdultMass[MSY_Allee$bFactor == "0"])
  dt_MSY_Allee$medianMass[3] = median(MSY_Allee$AdultMass[MSY_Allee$bFactor == "100"])
  dt_MSY_Allee$medianMass[4] = median(MSY_Allee$AdultMass[MSY_Allee$bFactor == "200"])
  
  #Part 2
  dt_MSY_NoAllee = data.frame("type" = c("None (hom. habitat)","None (het. habitat)","Low","High"))
  
  dt_MSY_NoAllee$Allee = 0
  
  dt_MSY_NoAllee$medianCollapse = 0
  dt_MSY_NoAllee$medianCollapse[1] = median(MSY_NoAllee$Patch.collapse[MSY_NoAllee$bFactor == "HOM"])
  dt_MSY_NoAllee$medianCollapse[2] = median(MSY_NoAllee$Patch.collapse[MSY_NoAllee$bFactor == "0"])
  dt_MSY_NoAllee$medianCollapse[3] = median(MSY_NoAllee$Patch.collapse[MSY_NoAllee$bFactor == "100"])
  dt_MSY_NoAllee$medianCollapse[4] = median(MSY_NoAllee$Patch.collapse[MSY_NoAllee$bFactor == "200"])
  
  dt_MSY_NoAllee$medianMass = 0
  dt_MSY_NoAllee$medianMass[1] = median(MSY_NoAllee$AdultMass[MSY_NoAllee$bFactor == "HOM"])
  dt_MSY_NoAllee$medianMass[2] = median(MSY_NoAllee$AdultMass[MSY_NoAllee$bFactor == "0"])
  dt_MSY_NoAllee$medianMass[3] = median(MSY_NoAllee$AdultMass[MSY_NoAllee$bFactor == "100"])
  dt_MSY_NoAllee$medianMass[4] = median(MSY_NoAllee$AdultMass[MSY_NoAllee$bFactor == "200"])
  
  dt_MSY = rbind(dt_MSY_Allee, dt_MSY_NoAllee)
  dt_MSY$h = "MSY"
  dt_MSY$Allee = as.factor(dt_MSY$Allee)
  
  # Create OVER data set
  
  OVER_Allee = read.csv('output/figure2B.csv')
  OVER_NoAllee = read.csv('output/figure3B.csv')
  
  # Part 1
  dt_OVER_Allee = data.frame("type" = c("None (hom. habitat)","None (het. habitat)","Low","High"))
  
  dt_OVER_Allee$Allee = 1
  
  dt_OVER_Allee$medianCollapse = 0
  dt_OVER_Allee$medianCollapse[1] = median(OVER_Allee$Patch.collapse[OVER_Allee$bFactor == "HOM"])
  dt_OVER_Allee$medianCollapse[2] = median(OVER_Allee$Patch.collapse[OVER_Allee$bFactor == "0"])
  dt_OVER_Allee$medianCollapse[3] = median(OVER_Allee$Patch.collapse[OVER_Allee$bFactor == "100"])
  dt_OVER_Allee$medianCollapse[4] = median(OVER_Allee$Patch.collapse[OVER_Allee$bFactor == "200"])
  
  dt_OVER_Allee$medianMass = 0
  dt_OVER_Allee$medianMass[1] = median(OVER_Allee$AdultMass[OVER_Allee$bFactor == "HOM"])
  dt_OVER_Allee$medianMass[2] = median(OVER_Allee$AdultMass[OVER_Allee$bFactor == "0"])
  dt_OVER_Allee$medianMass[3] = median(OVER_Allee$AdultMass[OVER_Allee$bFactor == "100"])
  dt_OVER_Allee$medianMass[4] = median(OVER_Allee$AdultMass[OVER_Allee$bFactor == "200"])
  
  #Part 2
  dt_OVER_NoAllee = data.frame("type" = c("None (hom. habitat)","None (het. habitat)","Low","High"))
  
  dt_OVER_NoAllee$Allee = 0
  
  dt_OVER_NoAllee$medianCollapse = 0
  dt_OVER_NoAllee$medianCollapse[1] = median(OVER_NoAllee$Patch.collapse[OVER_NoAllee$bFactor == "HOM"])
  dt_OVER_NoAllee$medianCollapse[2] = median(OVER_NoAllee$Patch.collapse[OVER_NoAllee$bFactor == "0"])
  dt_OVER_NoAllee$medianCollapse[3] = median(OVER_NoAllee$Patch.collapse[OVER_NoAllee$bFactor == "100"])
  dt_OVER_NoAllee$medianCollapse[4] = median(OVER_NoAllee$Patch.collapse[OVER_NoAllee$bFactor == "200"])
  
  dt_OVER_NoAllee$medianMass = 0
  dt_OVER_NoAllee$medianMass[1] = median(OVER_NoAllee$AdultMass[OVER_NoAllee$bFactor == "HOM"])
  dt_OVER_NoAllee$medianMass[2] = median(OVER_NoAllee$AdultMass[OVER_NoAllee$bFactor == "0"])
  dt_OVER_NoAllee$medianMass[3] = median(OVER_NoAllee$AdultMass[OVER_NoAllee$bFactor == "100"])
  dt_OVER_NoAllee$medianMass[4] = median(OVER_NoAllee$AdultMass[OVER_NoAllee$bFactor == "200"])
  
  dt_OVER = rbind(dt_OVER_Allee, dt_OVER_NoAllee)
  dt_OVER$h = "OVER"
  dt_OVER$Allee = as.factor(dt_OVER$Allee)
  
  dt = rbind(dt_MSY, dt_OVER)
  
  dt_allee = dt[dt$Allee == 1,]
  dt_noAllee = dt[dt$Allee == 0,]
  
  
  dt$fish = '0'
  
  for(i in 1:length(dt$type)){
    dt$fish[i] = paste(as.character(dt$Allee[i]),as.character(dt$h[i]))
  }
  
  dt$fish = as.factor(dt$fish)
  
  p3c <- ggplot(dt, aes(x=medianCollapse, y=medianMass/1000, color=fish, shape = type)) +
    geom_point(size=5, alpha = 0.8) + 
    labs(x="Median proportion of extinct patches", y="Median adult biomass (tons)", color = "", shape = "Fishing Heterogeneity") +
    scale_color_discrete(labels=c("No Allee effect: MSY", "No Allee effect: Over MSY", "Allee effect: MSY", "Allee effect: Over MSY")) +
    theme_minimal() +
    theme(axis.text.x =element_text(size=8),legend.title = element_text(size = 9), 
          legend.text = element_text(size = 7))
  
  fig3 <- ggarrange(ggarrange(p3a,p3b,ncol = 2, labels = c("A","B"),font.label = list(size = 12)),p3c, labels = c("","C"), 
                    ncol = 1, nrow = 2, vjust = 0.5,font.label = list(size = 12))
  myPPlot(fig3)
  
  ggsave("fig3.tiff",fig3, units="mm", width=170, height= 150, dpi=300)
  
  
}

statsAnalysis3 <- function(){
  
  dt_patch_collapse = read.csv('output/figure3A.csv')
  
  # Check normality
  data_HOM <- dt_patch_collapse[dt_patch_collapse$habitat == "HOM",]
  
  shapiro.test(data_HOM$AdultMass) # p: 0.6001
  shapiro.test(data_HOM$catch) # p: 0.3332
  
  data_HET_0 <- dt_patch_collapse[dt_patch_collapse$bFactor == "0",]
  
  shapiro.test(data_HET_0$AdultMass) # p:0.6921
  shapiro.test(data_HET_0$catch) # p: 0.4695
  
  t.test(data_HOM$AdultMass, data_HET_0$AdultMass, paired = TRUE) #p: 1.467e-12
  t.test(data_HOM$catch, data_HET_0$catch, paired = TRUE)# p: 4.388e-12
  
  pre_anova = dt_patch_collapse[dt_patch_collapse$habitat == "HET",]
  pre_anova$reps = rep(1:20,3)
  pre_anova$bFactor = factor(pre_anova$bFactor)
  pre_anova$reps = factor(pre_anova$reps)
  
  #pairwise.t.test(data = pre_anova, x = pre_anova$AdultMass, g = pre_anova$bFactor, paired = TRUE, p.adjust.method = "bonferroni")
  
  lme_model <- lme(data = pre_anova, AdultMass/1000 ~ factor(bFactor), random = ~1|reps)
  summary(lme_model)
  ranef(lme_model)
  anova(lme_model)
  intervals(lme_model)
  
  
  ## OVER MSY
  dt_patch_collapse_overMSY = read.csv('output/figure3B.csv')
  
  # Check normality
  data_HOM <- dt_patch_collapse_overMSY[dt_patch_collapse_overMSY$habitat == "HOM",]
  
  shapiro.test(data_HOM$AdultMass) # p: 0.4901
  shapiro.test(data_HOM$catch) # p: 0.6756
  
  data_HET_0 <- dt_patch_collapse_overMSY[dt_patch_collapse_overMSY$bFactor == "0",]
  
  shapiro.test(data_HET_0$AdultMass) # p:0.5911
  shapiro.test(data_HET_0$catch) # p: 0.6999
  
  t.test(data_HOM$AdultMass, data_HET_0$AdultMass, paired = TRUE) #p: 6.67e-13
  t.test(data_HOM$catch, data_HET_0$catch, paired = TRUE)# p: 1.47e-12
  
  mean(data_HOM$AdultMass -  data_HET_0$AdultMass)/mean(data_HOM$AdultMass)
  
  pre_anova = dt_patch_collapse_overMSY[dt_patch_collapse_overMSY$habitat == "HET",]
  pre_anova$reps = rep(1:20,3)
  pre_anova$bFactor = factor(pre_anova$bFactor)
  pre_anova$reps = factor(pre_anova$reps)
  
  lme_model <- lme(data = pre_anova, AdultMass/1000 ~ factor(bFactor), random = ~1|reps)
  summary(lme_model)
  ranef(lme_model)
  anova(lme_model)
  intervals(lme_model)
  
  pairwise.t.test(data = pre_anova, x = pre_anova$AdultMass, g = pre_anova$bFactor, paired = TRUE, p.adjust.method = "bonferroni")
  
  
}

analysis_figure4 <- function(){
  
  ### AR = 0.95 ###
  
  set.seed(2020) # 2020 for all analysis
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.95, ma = MA), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK = corrK-mean(corrK)+DEF_K
  
  # Calculate baseline (no harvest) to determine collapse
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, ar = AR, ma = MA,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, ar = AR, ma = MA,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=0.09, habitatHOM = TRUE, ar = AR, ma = MA,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                            maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                            isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST0$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 0
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST0$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 0
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST0$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 0
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapse06 = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapse06$bFactor = as.factor(dt_patch_collapse06$bFactor)
  dt_patch_collapse06$dFactor = as.factor(dt_patch_collapse06$dFactor)
  dt_patch_collapse06$habitat = as.factor(dt_patch_collapse06$habitat)
  
  saveMatrixToCSVFile('figure4_095', dt_patch_collapse06)
  
  p <- ggplot(dt_patch_collapse06, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=10))
  
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse06, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Total adult mass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse06, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  
  # Above MSY
  
  # Case 4
  
  OVER_BIO0_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                             maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                             isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                             hSDProp=0.1, classicSurv=FALSE,
                             runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                             showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                             limitOutput = TRUE)
  
  spatialAdultMass_4 = OVER_BIO0_DIST0$adultMass
  patch_collapse_4 = matrix(0,length(spatialAdultMass_4[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_4[,1])) {
    patch_collapse_4[i,] = spatialAdultMass_4[i,] < collapse_mass 
  }
  
  collapse_per_rep_4 = rowMeans(1*patch_collapse_4)
  dt_4 = data.frame("Patch collapse" = collapse_per_rep_4)
  dt_4$bFactor = 0
  dt_4$dFactor = 0
  dt_4$AdultMass = apply(spatialAdultMass_4, 1, sum)
  dt_4$catch = OVER_BIO0_DIST0$catch
  dt_4$habitat = 'HET'
  
  
  # Case 5
  
  OVER_BIO100_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_5 = OVER_BIO100_DIST0$adultMass
  patch_collapse_5 = matrix(0,length(spatialAdultMass_5[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_5[,1])) {
    patch_collapse_5[i,] = spatialAdultMass_5[i,] < collapse_mass 
  }
  
  collapse_per_rep_5 = rowMeans(1*patch_collapse_5)
  dt_5 = data.frame("Patch collapse" = collapse_per_rep_5)
  dt_5$bFactor = 100
  dt_5$dFactor = 0
  dt_5$AdultMass = apply(spatialAdultMass_5, 1, sum)
  dt_5$catch = OVER_BIO100_DIST0$catch
  dt_5$habitat = 'HET'
  
  
  #Case 6
  
  OVER_BIO200_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_6 = OVER_BIO200_DIST0$adultMass
  patch_collapse_6 = matrix(0,length(spatialAdultMass_6[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_6[,1])) {
    patch_collapse_6[i,] = spatialAdultMass_6[i,] < collapse_mass 
  }
  
  collapse_per_rep_6 = rowMeans(1*patch_collapse_6)
  dt_6 = data.frame("Patch collapse" = collapse_per_rep_6)
  dt_6$bFactor = 200
  dt_6$dFactor = 0
  dt_6$AdultMass = apply(spatialAdultMass_6, 1, sum)
  dt_6$catch = OVER_BIO200_DIST0$catch
  dt_6$habitat = 'HET'
  
  #Case 7
  
  OVER_BIO0_DIST0_HOM = runModel(hProp=0.12, habitatHOM = TRUE, ar = AR, ma = MA,
                                 maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                 isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                 hSDProp=0.1, classicSurv=FALSE,
                                 runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                 showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                 limitOutput = TRUE)
  
  
  spatialAdultMass_7 = OVER_BIO0_DIST0_HOM$adultMass
  patch_collapse_7 = matrix(0,length(spatialAdultMass_7[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_7[,1])) {
    patch_collapse_7[i,] = spatialAdultMass_7[i,] < collapse_mass 
  }
  
  collapse_per_rep_7 = rowMeans(1*patch_collapse_7)
  dt_7 = data.frame("Patch collapse" = collapse_per_rep_7)
  dt_7$bFactor = 'HOM'
  dt_7$dFactor = 0
  dt_7$AdultMass = apply(spatialAdultMass_7, 1, sum)
  dt_7$catch = OVER_BIO0_DIST0_HOM$catch
  dt_7$habitat = 'HOM'
  
  
  
  
  dt_patch_collapse_overMSY06 = rbind(dt_4,dt_5,dt_6, dt_7)
  dt_patch_collapse_overMSY06$bFactor = as.factor(dt_patch_collapse_overMSY06$bFactor)
  dt_patch_collapse_overMSY06$dFactor = as.factor(dt_patch_collapse_overMSY06$dFactor)
  dt_patch_collapse_overMSY06$habitat = as.factor(dt_patch_collapse_overMSY06$habitat)
  
  saveMatrixToCSVFile('figure4_095_OVER', dt_patch_collapse_overMSY06)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="Over MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Total Adult Mass (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=catch/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))+
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=5, color="black") 
  myPPlot(p)
  
  
  ## AR = 0.6
  # MSY =  0.09, Catch = 24015.55
  
  # Calculate baseline (no harvest) to determine collapse
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, ar = AR, ma = MA,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, ar = AR, ma = MA,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=0.09, habitatHOM = TRUE, ar = AR, ma = MA,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                            maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                            isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST0$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 0
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST0$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 0
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST0$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 0
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapse06 = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapse06$bFactor = as.factor(dt_patch_collapse06$bFactor)
  dt_patch_collapse06$dFactor = as.factor(dt_patch_collapse06$dFactor)
  dt_patch_collapse06$habitat = as.factor(dt_patch_collapse06$habitat)
  
  saveMatrixToCSVFile('figure4A', dt_patch_collapse06)
  
  p <- ggplot(dt_patch_collapse06, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=10))
  
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse06, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Total adult mass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse06, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  
  # Above MSY
  
  # Case 4
  
  OVER_BIO0_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                             maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                             isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                             hSDProp=0.1, classicSurv=FALSE,
                             runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                             showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                             limitOutput = TRUE)
  
  spatialAdultMass_4 = OVER_BIO0_DIST0$adultMass
  patch_collapse_4 = matrix(0,length(spatialAdultMass_4[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_4[,1])) {
    patch_collapse_4[i,] = spatialAdultMass_4[i,] < collapse_mass 
  }
  
  collapse_per_rep_4 = rowMeans(1*patch_collapse_4)
  dt_4 = data.frame("Patch collapse" = collapse_per_rep_4)
  dt_4$bFactor = 0
  dt_4$dFactor = 0
  dt_4$AdultMass = apply(spatialAdultMass_4, 1, sum)
  dt_4$catch = OVER_BIO0_DIST0$catch
  dt_4$habitat = 'HET'
  
  
  # Case 5
  
  OVER_BIO100_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_5 = OVER_BIO100_DIST0$adultMass
  patch_collapse_5 = matrix(0,length(spatialAdultMass_5[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_5[,1])) {
    patch_collapse_5[i,] = spatialAdultMass_5[i,] < collapse_mass 
  }
  
  collapse_per_rep_5 = rowMeans(1*patch_collapse_5)
  dt_5 = data.frame("Patch collapse" = collapse_per_rep_5)
  dt_5$bFactor = 100
  dt_5$dFactor = 0
  dt_5$AdultMass = apply(spatialAdultMass_5, 1, sum)
  dt_5$catch = OVER_BIO100_DIST0$catch
  dt_5$habitat = 'HET'
  
  
  #Case 6
  
  OVER_BIO200_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_6 = OVER_BIO200_DIST0$adultMass
  patch_collapse_6 = matrix(0,length(spatialAdultMass_6[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_6[,1])) {
    patch_collapse_6[i,] = spatialAdultMass_6[i,] < collapse_mass 
  }
  
  collapse_per_rep_6 = rowMeans(1*patch_collapse_6)
  dt_6 = data.frame("Patch collapse" = collapse_per_rep_6)
  dt_6$bFactor = 200
  dt_6$dFactor = 0
  dt_6$AdultMass = apply(spatialAdultMass_6, 1, sum)
  dt_6$catch = OVER_BIO200_DIST0$catch
  dt_6$habitat = 'HET'
  
  #Case 7
  
  OVER_BIO0_DIST0_HOM = runModel(hProp=0.12, habitatHOM = TRUE, ar = AR, ma = MA,
                                 maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                 isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                 hSDProp=0.1, classicSurv=FALSE,
                                 runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                 showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                 limitOutput = TRUE)
  
  
  spatialAdultMass_7 = OVER_BIO0_DIST0_HOM$adultMass
  patch_collapse_7 = matrix(0,length(spatialAdultMass_7[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_7[,1])) {
    patch_collapse_7[i,] = spatialAdultMass_7[i,] < collapse_mass 
  }
  
  collapse_per_rep_7 = rowMeans(1*patch_collapse_7)
  dt_7 = data.frame("Patch collapse" = collapse_per_rep_7)
  dt_7$bFactor = 'HOM'
  dt_7$dFactor = 0
  dt_7$AdultMass = apply(spatialAdultMass_7, 1, sum)
  dt_7$catch = OVER_BIO0_DIST0_HOM$catch
  dt_7$habitat = 'HOM'
  
  
  
  
  dt_patch_collapse_overMSY06 = rbind(dt_4,dt_5,dt_6, dt_7)
  dt_patch_collapse_overMSY06$bFactor = as.factor(dt_patch_collapse_overMSY06$bFactor)
  dt_patch_collapse_overMSY06$dFactor = as.factor(dt_patch_collapse_overMSY06$dFactor)
  dt_patch_collapse_overMSY06$habitat = as.factor(dt_patch_collapse_overMSY06$habitat)
  
  saveMatrixToCSVFile('figure4B', dt_patch_collapse_overMSY06)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="Over MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Total Adult Mass (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=catch/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))+
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=5, color="black") 
  myPPlot(p)
  
  ## AR = 0.3
  
  # MSY =  0.09, Catch = 24015.55
  
  set.seed(2020)
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.3, ma = MA), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK = corrK-mean(corrK)+DEF_K
  
  # Calculate baseline (no harvest) to determine collapse
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, ar = AR, ma = MA,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, ar = AR, ma = MA,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=0.09, habitatHOM = TRUE, ar = AR, ma = MA,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                            maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                            isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST0$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 0
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST0$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 0
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST0 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST0$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 0
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapse03 = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapse03$bFactor = as.factor(dt_patch_collapse03$bFactor)
  dt_patch_collapse03$dFactor = as.factor(dt_patch_collapse03$dFactor)
  dt_patch_collapse03$habitat = as.factor(dt_patch_collapse03$habitat)
  
  saveMatrixToCSVFile('figure4C', dt_patch_collapse03)
  
  p <- ggplot(dt_patch_collapse03, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=10))
  
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse03, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Total adult mass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse06, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  
  # Above MSY
  
  # Case 4
  
  OVER_BIO0_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                             maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                             isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                             hSDProp=0.1, classicSurv=FALSE,
                             runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                             showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                             limitOutput = TRUE)
  
  spatialAdultMass_4 = OVER_BIO0_DIST0$adultMass
  patch_collapse_4 = matrix(0,length(spatialAdultMass_4[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_4[,1])) {
    patch_collapse_4[i,] = spatialAdultMass_4[i,] < collapse_mass 
  }
  
  collapse_per_rep_4 = rowMeans(1*patch_collapse_4)
  dt_4 = data.frame("Patch collapse" = collapse_per_rep_4)
  dt_4$bFactor = 0
  dt_4$dFactor = 0
  dt_4$AdultMass = apply(spatialAdultMass_4, 1, sum)
  dt_4$catch = OVER_BIO0_DIST0$catch
  dt_4$habitat = 'HET'
  
  
  # Case 5
  
  OVER_BIO100_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_5 = OVER_BIO100_DIST0$adultMass
  patch_collapse_5 = matrix(0,length(spatialAdultMass_5[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_5[,1])) {
    patch_collapse_5[i,] = spatialAdultMass_5[i,] < collapse_mass 
  }
  
  collapse_per_rep_5 = rowMeans(1*patch_collapse_5)
  dt_5 = data.frame("Patch collapse" = collapse_per_rep_5)
  dt_5$bFactor = 100
  dt_5$dFactor = 0
  dt_5$AdultMass = apply(spatialAdultMass_5, 1, sum)
  dt_5$catch = OVER_BIO100_DIST0$catch
  dt_5$habitat = 'HET'
  
  
  #Case 6
  
  OVER_BIO200_DIST0 = runModel(hProp=0.12, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                               isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_6 = OVER_BIO200_DIST0$adultMass
  patch_collapse_6 = matrix(0,length(spatialAdultMass_6[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_6[,1])) {
    patch_collapse_6[i,] = spatialAdultMass_6[i,] < collapse_mass 
  }
  
  collapse_per_rep_6 = rowMeans(1*patch_collapse_6)
  dt_6 = data.frame("Patch collapse" = collapse_per_rep_6)
  dt_6$bFactor = 200
  dt_6$dFactor = 0
  dt_6$AdultMass = apply(spatialAdultMass_6, 1, sum)
  dt_6$catch = OVER_BIO200_DIST0$catch
  dt_6$habitat = 'HET'
  
  #Case 7
  
  OVER_BIO0_DIST0_HOM = runModel(hProp=0.12, habitatHOM = TRUE, ar = AR, ma = MA,
                                 maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                 isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                 hSDProp=0.1, classicSurv=FALSE,
                                 runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                 showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                 limitOutput = TRUE)
  
  
  spatialAdultMass_7 = OVER_BIO0_DIST0_HOM$adultMass
  patch_collapse_7 = matrix(0,length(spatialAdultMass_7[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_7[,1])) {
    patch_collapse_7[i,] = spatialAdultMass_7[i,] < collapse_mass 
  }
  
  collapse_per_rep_7 = rowMeans(1*patch_collapse_7)
  dt_7 = data.frame("Patch collapse" = collapse_per_rep_7)
  dt_7$bFactor = 'HOM'
  dt_7$dFactor = 0
  dt_7$AdultMass = apply(spatialAdultMass_7, 1, sum)
  dt_7$catch = OVER_BIO0_DIST0_HOM$catch
  dt_7$habitat = 'HOM'
  
  dt_patch_collapse_overMSY03 = rbind(dt_4,dt_5,dt_6, dt_7)
  dt_patch_collapse_overMSY03$bFactor = as.factor(dt_patch_collapse_overMSY03$bFactor)
  dt_patch_collapse_overMSY03$dFactor = as.factor(dt_patch_collapse_overMSY03$dFactor)
  dt_patch_collapse_overMSY03$habitat = as.factor(dt_patch_collapse_overMSY03$habitat)
  
  saveMatrixToCSVFile('figure4D', dt_patch_collapse_overMSY03)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="Over MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Total Adult Mass (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY06, aes(x=bFactor, y=catch/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))+
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=5, color="black") 
  myPPlot(p)
  
  
}

makeFigure4 <- function(){
  
  MSY_03 <- read.csv('output/figure4C.csv')
  OVERMSY_03 <- read.csv('output/figure4D.csv')
  
  MSY_06 <- read.csv('output/figure4A.csv')
  OVERMSY_06 <- read.csv('output/figure4B.csv')
  
  MSY_095 <- read.csv('output/figure4_095.csv')
  OVERMSY_095 <- read.csv('output/figure4_095_OVER.csv')
  
  
  MSY_03$corr = 0.3
  MSY_03$h = 'MSY'
  
  OVERMSY_03$corr = 0.3
  OVERMSY_03$h = 'OVER'
  
  MSY_06$corr = 0.6
  MSY_06$h = 'MSY'
  
  OVERMSY_06$corr = 0.6
  OVERMSY_06$h = 'OVER'
  
  
  MSY_095$corr = 0.95
  MSY_095$h = 'MSY'
  
  OVERMSY_095$corr = 0.95
  OVERMSY_095$h = 'OVER'
  
  
  fig4df = rbind(MSY_06,OVERMSY_06, MSY_095, OVERMSY_095)
  fig4df$corr = as.factor(fig4df$corr)
  levels(fig4df$corr) <- c("AR = 0.6", "AR = 0.95")
  
  fig4a <- ggplot(fig4df, aes(x=bFactor, y=AdultMass/1000, fill = habitat, color = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x="Fishing heterogeneity level", y="Adult biomass (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=2, color="black") +
    guides(color = FALSE) +  
    facet_grid(rows = vars(h), cols = vars(corr))+ 
    theme_light() 
  
  #myPPlot(fig4a)
  
  fig4b <- ggplot(fig4df, aes(x=bFactor, y=Patch.collapse*100, fill = habitat, color = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x="Fishing heterogeneity level", y="Blocks collapsed (%)", fill = "Habitat", color = "") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=2, color="black") +
    guides(color = FALSE) +  
    facet_grid(rows = vars(h), cols = vars(corr)) +
    theme_light() 
  
  #myPPlot(fig4b)
  
  ## Assembling the figure
  
  fig4 = ggarrange(fig4b,fig4a, ncol = 1, nrow = 2,common.legend = TRUE, labels = c("A","B"), legend = 'right')
  
  #myPPlot(fig4,  w = 4, h = 6)
  
  ggsave("fig4.tiff",fig4, units="mm", width=170, height= 130, dpi=300)
  
}

statsAnalysis4 <- function(){
  
  MSY_06 <- read.csv('output/figure4A.csv')
  OVERMSY_06 <- read.csv('output/figure4B.csv')
  
  MSY06_HOM <- MSY_06[MSY_06$habitat == "HOM",]
  
  shapiro.test(MSY06_HOM$AdultMass) # p: 0.2283
  shapiro.test(MSY06_HOM$catch) # p: 0.2249
  
  MSY06_HET_0 <- MSY_06[MSY_06$bFactor == "0",]
  
  shapiro.test(MSY06_HET_0$AdultMass) # p:0.2584
  shapiro.test(MSY06_HET_0$catch) # p: 0.4186
  
  t.test(MSY06_HOM$AdultMass,MSY06_HET_0$AdultMass, paired = TRUE) #p: 1.905e-06
  t.test(MSY06_HOM$catch, MSY06_HET_0$catch, paired = TRUE)# p: 1.754e-06
  
  mean(MSY06_HOM$AdultMass -  MSY06_HET_0$AdultMass)/mean(MSY06_HOM$AdultMass)
  
  pre_anova = MSY_06[MSY_06$habitat == "HET",]
  pre_anova$reps = rep(1:20,3)
  pre_anova$bFactor = factor(pre_anova$bFactor)
  pre_anova$reps = factor(pre_anova$reps)
  
  pairwise.t.test(data = pre_anova, x = pre_anova$AdultMass, g = pre_anova$bFactor, paired = TRUE, p.adjust.method = "bonferroni")
  
  lme_model <- lme(data = pre_anova, AdultMass/1000 ~ factor(bFactor), random = ~1|reps)
  summary(lme_model)
  ranef(lme_model)
  anova(lme_model)
  intervals(lme_model)
  
  
  MSY_095 <- read.csv('output/figure4_095.csv')
  OVERMSY_095 <- read.csv('output/figure4_095_OVER.csv')
  
  MSY095_HOM <- MSY_095[MSY_095$habitat == "HOM",]
  
  
  shapiro.test(MSY095_HOM$AdultMass) # p: 0.2283
  shapiro.test(MSY095_HOM$catch) # p: 0.2249
  
  MSY095_HET_0 <- MSY_095[MSY_095$bFactor == "0",]
  
  shapiro.test(MSY095_HET_0$AdultMass) # p:0.2584
  shapiro.test(MSY095_HET_0$catch) # p: 0.4186
  
  t.test(MSY095_HOM$AdultMass,MSY095_HET_0$AdultMass, paired = TRUE) #p: 1.723e-09
  t.test(MSY095_HOM$catch, MSY095_HET_0$catch, paired = TRUE)# p: 2.264e-09
  
  mean(MSY095_HOM$AdultMass -  MSY095_HET_0$AdultMass)/mean(MSY095_HOM$AdultMass)
  
  pre_anova = MSY_095[MSY_095$habitat == "HET",]
  pre_anova$reps = rep(1:20,3)
  pre_anova$bFactor = factor(pre_anova$bFactor)
  pre_anova$reps = factor(pre_anova$reps)
  
  lme_model <- lme(data = pre_anova, AdultMass/1000 ~ factor(bFactor), random = ~1|reps)
  summary(lme_model)
  ranef(lme_model)
  anova(lme_model)
  intervals(lme_model)
  
  
  pairwise.t.test(data = pre_anova, x = pre_anova$AdultMass, g = pre_anova$bFactor, paired = TRUE, p.adjust.method = "bonferroni")
  
  
  
  
}

figure5 <- function() {
  
  set.seed(2020)
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.9, ma = MA), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK = corrK-mean(corrK)+DEF_K
  
  # Calculate baseline (no harvest) to determine collapse
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, ar = AR, ma = MA,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, ar = AR, ma = MA,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=0.09, habitatHOM = TRUE, ar = AR, ma = MA,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST2 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                            maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 2,
                            isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST2$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 2
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST2 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 2,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST2$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 2
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST2 = runModel(hProp=0.09, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 2,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST2$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 2
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapseDist = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapseDist$bFactor = as.factor(dt_patch_collapseDist$bFactor)
  dt_patch_collapseDist$dFactor = as.factor(dt_patch_collapseDist$dFactor)
  dt_patch_collapseDist$habitat = as.factor(dt_patch_collapseDist$habitat)
  
  saveMatrixToCSVFile('figure5', dt_patch_collapseDist)
  
  p5d <- ggplot(dt_patch_collapseDist, aes(x=bFactor, y=Patch.collapse*100, fill = habitat, color = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x="Fishing heterogeneity", y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    guides(color = FALSE) + 
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=3, color="black") + 
    theme_minimal() + 
    theme(axis.text.x =element_text(size=10))
  
  myPPlot(p5d)
  
  p5e <- ggplot(dt_patch_collapseDist, aes(x=bFactor, y=AdultMass/1000, fill = habitat, color = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x="Fishing heterogeneity", y="Total adult mass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    guides(color = FALSE) + 
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=3, color="black") +
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  
  myPPlot(p5e)
  
  p <- ggplot(dt_patch_collapseDist, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p5de <- ggarrange(p5d,p5e,labels = c("B","C"),common.legend = TRUE, legend = "bottom")
  
  myPPlot(p5de)
  
  saveMatrixToCSVFile('figure5A', spatialAdultMass_1)
  saveMatrixToCSVFile('figure5B', spatialAdultMass_2)
  saveMatrixToCSVFile('figure5C', spatialAdultMass_3)
  
  # Figure 5a
  
  fig5a = data_frame()
  fig5b = data_frame()
  fig5c = data_frame()
  df = data_frame(coast = 1:150)
  
  for (i in 1:length(spatialAdultMass_1[,1])){
    
    df$coast = 1:150
    df$rep = i
    df$rep = as.factor(df$rep)
    
    df$mass = spatialAdultMass_1[i,]
    fig5a = rbind(fig5a,df)
    
    df$mass = spatialAdultMass_2[i,]
    fig5b = rbind(fig5b,df)
    
    df$mass = spatialAdultMass_3[i,]
    fig5c = rbind(fig5c,df)
    
  }
  
  
  
  
  p5a <- ggplot(fig5a, aes(x=coast, y=mass/1000, group = rep, colour = rep)) +
    geom_line(size=0.25, alpha = 1) +
    ylim(0, 2.3) +
    labs(x="Coast", y="Adult mass (tons)") +
    theme(legend.position = "none",axis.title.x=element_blank())
  
  myPPlot(p5a)
  
  p5b <- ggplot(fig5b, aes(x=coast, y=mass/1000, group = rep, colour = rep)) +
    geom_line(size=0.25, alpha = 1) +
    ylim(0, 2.3) +
    labs(x="", y="Adult mass (tons)") +
    theme(legend.position = "none",axis.title.x=element_blank())
  
  myPPlot(p5b)
  
  p5c <- ggplot(fig5c, aes(x=coast, y=mass/1000, group = rep, colour = rep)) +
    geom_line(size=0.25, alpha = 1) +
    ylim(0, 2.3) +
    labs(x="", y="Adult mass (tons)") +
    theme(legend.position = "none",axis.title.x=element_blank())
  
  myPPlot(p5c)
  
  fig5a$fishing = "Null"
  fig5b$fishing = "Int."
  fig5c$fishing = "High"
  
  fig5abc = rbind(fig5a,fig5b,fig5c)
  fig5abc$fishing = as.factor(fig5abc$fishing,levels = c("Null", "Int.", "High"))
  fig5abc$fishing = factor(fig5abc$fishing, levels=rev(levels(fig5abc$fishing)))
  
  p5abc <- ggplot(fig5abc, aes(x=coast, y=mass/1000, group = rep, colour = rep)) +
    geom_line(size=0.25, alpha = 1) +
    ylim(0, 2.3) +
    labs(x="Coast (block #)", y="Adult mass (tons)") +
    theme(legend.position = "none") +
    facet_grid(cols = vars(fishing))
  
  myPPlot(p5abc)
  
  p5 = ggarrange(p5abc,p5de, nrow = 2, ncol = 1,labels = "A")
  myPPlot(p5)
  
}  

makeFigure5 <- function(){
  
  
  dt_patch_collapseDist = read.csv('output/figure5.csv')
  
  spatialAdultMass_1 = read.csv('output/figure5A.csv')
  spatialAdultMass_1$X = NULL
  
  spatialAdultMass_2 = read.csv('output/figure5B.csv')
  spatialAdultMass_2$X = NULL
  spatialAdultMass_3 = read.csv('output/figure5C.csv')
  spatialAdultMass_3$X = NULL
  
  fig5a = data_frame()
  fig5b = data_frame()
  fig5c = data_frame()
  df = data_frame(coast = 1:150)
  
  for (i in 1:length(spatialAdultMass_1[,1])){
    
    df$coast = 1:150
    df$rep = i
    df$rep = as.factor(df$rep)
    
    df$mass = t(spatialAdultMass_1[i,])
    fig5a = rbind(fig5a,df)
    
    df$mass = t(spatialAdultMass_2[i,])
    fig5b = rbind(fig5b,df)
    
    df$mass = t(spatialAdultMass_3[i,])
    fig5c = rbind(fig5c,df)
    
  }
  
  p5a <- ggplot(fig5a, aes(x=coast, y=mass/1000, group = rep, colour = rep)) +
    geom_line(size=0.25, alpha = 1) +
    ylim(0, 2.3) +
    labs(x="Coast", y="Adult mass (tons)") +
    theme(legend.position = "none",axis.title = element_text(size = 9))
  
  #myPPlot(p5a)
  
  p5b <- ggplot(fig5b, aes(x=coast, y=mass/1000, group = rep, colour = rep)) +
    geom_line(size=0.25, alpha = 1) +
    ylim(0, 2.3) +
    labs(x="", y="Adult mass (tons)") +
    theme(legend.position = "none",axis.title = element_text(size = 9))
  
  #myPPlot(p5b)
  
  p5c <- ggplot(fig5c, aes(x=coast, y=mass/1000, group = rep, colour = rep)) +
    geom_line(size=0.25, alpha = 1) +
    ylim(0, 2.3) +
    labs(x="", y="Adult mass (tons)") +
    theme(legend.position = "none",axis.title = element_text(size = 9))
  
  #myPPlot(p5c)
  
  fig5a$fishing = "None"
  fig5b$fishing = "Low"
  fig5c$fishing = "High"
  
  fig5abc = rbind(fig5a,fig5b,fig5c)
  fig5abc$fishing = as.factor(fig5abc$fishing)
  fig5abc$fishing = factor(fig5abc$fishing, levels=rev(levels(fig5abc$fishing)))
  
  p5abc <- ggplot(fig5abc, aes(x=coast, y=mass/1000, group = rep, colour = rep)) +
    geom_line(size=0.25, alpha = 1) +
    ylim(0, 2.3) +
    labs(x="Coast (block #)", y="Adult biomass (tons)") +
    theme_light() + 
    theme(legend.position = "none",axis.title = element_text(size = 9)) +
    facet_grid(cols = vars(fishing))
  
  #myPPlot(p5abc)
  
  p5d <- ggplot(dt_patch_collapseDist, aes(x=bFactor, y=Patch.collapse*100, fill = habitat, color = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x="Fishing heterogeneity", y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    guides(color = FALSE) + 
    theme_minimal() + 
    theme(axis.text.x =element_text(size=8)) + 
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.45, ymin = -1, ymax = 80, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 80, alpha = .2) + 
    annotate("text", x = 3.1, y=75, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=75, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 9))
  
  #myPPlot(p5d)
  
  p5e <- ggplot(dt_patch_collapseDist, aes(x=bFactor, y=AdultMass/1000, fill = habitat, color = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x="Fishing heterogeneity", y="Adult biomass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    guides(color = FALSE) + 
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) + 
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.45, ymin = 50, ymax = 165, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = 50, ymax = 165, alpha = .2) + 
    annotate("text", x = 3.1, y=158, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=158, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 9))
  
  #myPPlot(p5e)
  
  p <- ggplot(dt_patch_collapseDist, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8))
  #myPPlot(p)
  
  p5de <- ggarrange(p5d,p5e,labels = c("B","C"),font.label = list(size = 12))
  
  p5 = ggarrange(p5abc,p5de, nrow = 2, ncol = 1,labels = "A",font.label = list(size = 12))
  #myPPlot(p5, h = 5, w = 3)
  
  ggsave("fig5.tiff",p5, units="mm", width=170, height= 130, dpi=300)
  
}

figureS1_analysis <- function(){
  
  set.seed(2021) # 2020 for all analysis
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.9, ma = NULL), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK = corrK-mean(corrK)+DEF_K
  
  env1 <- function_analysis(h = 0.09, Allee = TRUE, randomK = corrK, rep = 10) # Start at 12:45
  env1$rand <- '1'
  
  env1b <- function_analysis(h = 0.12, Allee = TRUE, randomK = corrK, rep = 10)
  env1b$rand <- '1b'
  
  set.seed(2022) # 2020 for all analysis
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.9, ma = NULL), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK = corrK-mean(corrK)+DEF_K
  
  env2 <- function_analysis(h = 0.09, Allee = TRUE, randomK = corrK, rep = 10)
  env2$rand <- '2'
  
  env2b <- function_analysis(h = 0.12, Allee = TRUE, randomK = corrK, rep = 10)
  env2b$rand <- '2b'
  
  set.seed(2023) # 2020 for all analysis
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.9, ma = NULL), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK = corrK-mean(corrK)+DEF_K
  
  env3 <- function_analysis(h = 0.09, Allee = TRUE, randomK = corrK, rep = 10)
  env3$rand <- '3'
  
  env3b <- function_analysis(h = 0.12, Allee = TRUE, randomK = corrK, rep = 10)
  env3b$rand <- '3b'
  
  rand_ds <- rbind(env1,env2,env3)
  
  rand_ds_over <- rbind(env1b, env2b, env3b)
  
  tapply(rand_ds$AdultMass,rand_ds$rand,mean)
  
  write.csv(rand_ds,'FigureS1.csv')
  write.csv(rand_ds_over,'FigureS1b.csv')
  
  
}

makefigureS1 <- function(){
  
  rand_ds <- read.csv('FigureS1.csv')
  rand_ds_over <- read.csv('FigureS1b.csv')
  
  p1 <- ggplot(rand_ds, aes(x=bFactor, y=Patch.collapse*100, fill = rand))+
    geom_boxplot() +
    labs(x = "Fishing aggregation level",title="Exploitation at MSY", y="Extinct patches (%)") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=8)) +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 50, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 50, alpha = .2) + 
    annotate("text", x = 3.1, y=45, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=45, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, vjust = 2))
  
  
  p2 <- ggplot(rand_ds, aes(x=bFactor, y=AdultMass/1000, fill = rand))+
    geom_boxplot() +
    labs(x = "Fishing aggregation level", y="Adult biomass (tons)") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -5, ymax = 175, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -5, ymax = 175, alpha = .2) + 
    annotate("text", x = 3.1, y=165, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=165, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10)) +
    ylim(-5, 175)
  
  p3 <- ggplot(rand_ds_over, aes(x=bFactor, y=Patch.collapse*100, fill = rand))+
    geom_boxplot() +
    labs(x = "Fishing aggregation level",title="Exploitation above MSY", y="Extinct patches (%)") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=8)) +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 110, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 110, alpha = .2) + 
    annotate("text", x = 3.1, y=105, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=105, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, vjust = 2))
  
  p4 <- ggplot(rand_ds_over, aes(x=bFactor, y=AdultMass/1000, fill = rand))+
    geom_boxplot() +
    labs(x = "Fishing aggregation level", y="Adult biomass (tons)") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -5, ymax = 175, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -5, ymax = 175, alpha = .2) + 
    annotate("text", x = 3.1, y=165, label = "Het. habitat",size = 2) +
    annotate("text", x = 4, y=165, label = "Hom. habitat",size = 2) +
    theme(legend.position = "none",axis.title = element_text(size = 10)) +
    ylim(-5, 175)
  p4
  
  p <- ggarrange(p1,p3,p2,p4, nrow = 2, ncol = 2)
  
  myPPlot(p)
  
  
}

figureS2_analysis <- function(){
  
  ## 2/3 MSY
  
  # Case 8 
  
  BELOW_BIO0_DIST0 = runModel(hProp=0.055, habitatHOM = FALSE, hetK = corrK,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                              isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_8 = BELOW_BIO0_DIST0$adultMass
  patch_collapse_8 = matrix(0,length(spatialAdultMass_8[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_8[,1])) {
    patch_collapse_8[i,] = spatialAdultMass_8[i,] < collapse_mass 
  }
  
  collapse_per_rep_8 = rowMeans(1*patch_collapse_8)
  dt_8 = data.frame("Patch collapse" = collapse_per_rep_8)
  dt_8$bFactor = 0
  dt_8$dFactor = 0
  dt_8$AdultMass = apply(spatialAdultMass_8, 1, sum)
  dt_8$catch = BELOW_BIO0_DIST0$catch
  dt_8$habitat = 'HET'
  
  
  # Case 9
  
  BELOW_BIO100_DIST0 = runModel(hProp=0.055, habitatHOM = FALSE, hetK = corrK,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                                isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_9 = BELOW_BIO100_DIST0$adultMass
  patch_collapse_9 = matrix(0,length(spatialAdultMass_9[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_9[,1])) {
    patch_collapse_9[i,] = spatialAdultMass_9[i,] < collapse_mass 
  }
  
  collapse_per_rep_9 = rowMeans(1*patch_collapse_9)
  dt_9 = data.frame("Patch collapse" = collapse_per_rep_9)
  dt_9$bFactor = 100
  dt_9$dFactor = 0
  dt_9$AdultMass = apply(spatialAdultMass_9, 1, sum)
  dt_9$catch = BELOW_BIO100_DIST0$catch
  dt_9$habitat = 'HET'
  
  #Case 10
  
  BELOW_BIO200_DIST0 = runModel(hProp=0.055, habitatHOM = FALSE, hetK = corrK,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                                isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_10 = BELOW_BIO200_DIST0$adultMass
  patch_collapse_10 = matrix(0,length(spatialAdultMass_10[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_10[,1])) {
    patch_collapse_10[i,] = spatialAdultMass_10[i,] < collapse_mass 
  }
  
  collapse_per_rep_10 = rowMeans(1*patch_collapse_10)
  dt_10 = data.frame("Patch collapse" = collapse_per_rep_10)
  dt_10$bFactor = 200
  dt_10$dFactor = 0
  dt_10$AdultMass = apply(spatialAdultMass_10, 1, sum)
  dt_10$catch = BELOW_BIO200_DIST0$catch
  dt_10$habitat = 'HET'
  
  #Case 11
  
  BELOW_BIO0_DIST0_HOM = runModel(hProp=0.055, habitatHOM = TRUE, hetK = corrK,
                                  maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                  isAllee=TRUE, randDisp = TRUE, randRec= TRUE,
                                  hSDProp=0.1, classicSurv=FALSE,
                                  runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                  showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                  limitOutput = TRUE)
  
  
  spatialAdultMass_11 = BELOW_BIO0_DIST0_HOM$adultMass
  patch_collapse_11 = matrix(0,length(spatialAdultMass_11[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_11[,1])) {
    patch_collapse_11[i,] = spatialAdultMass_11[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_11 = rowMeans(1*patch_collapse_11)
  dt_11 = data.frame("Patch collapse" = collapse_per_rep_11)
  dt_11$bFactor = 'HOM'
  dt_11$dFactor = 0
  dt_11$AdultMass = apply(spatialAdultMass_11, 1, sum)
  dt_11$catch = BELOW_BIO0_DIST0_HOM$catch
  dt_11$habitat = 'HOM'
  
  dt_patch_collapse_belowMSY = rbind(dt_8,dt_9,dt_10, dt_11)
  dt_patch_collapse_belowMSY$bFactor = as.factor(dt_patch_collapse_belowMSY$bFactor)
  dt_patch_collapse_belowMSY$dFactor = as.factor(dt_patch_collapse_belowMSY$dFactor)
  dt_patch_collapse_belowMSY$habitat = as.factor(dt_patch_collapse_belowMSY$habitat)
  
  write.csv(dt_patch_collapse_belowMSY, 'figureS2.csv')
  
}

makefigureS2 <- function(){
  
  dt_patch_collapse_belowMSY <- read.csv('figureS2.csv')
  
  pa <- ggplot(dt_patch_collapse_belowMSY, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "", y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 22.5, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 22.5, alpha = .2) + 
    annotate("text", x = 3.1, y=21, label = "Het. habitat",size = 2.5) +
    annotate("text", x = 4, y=21, label = "Hom. habitat",size = 2.5) +
    theme(legend.position = "none",axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, vjust = 2),axis.title.x=element_blank())
  
  
  pb <- ggplot(dt_patch_collapse_belowMSY, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "", y="Adult biomass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x=element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 220, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 220, alpha = .2) + 
    annotate("text", x = 3.1, y=210, label = "Het. habitat",size = 2.5) +
    annotate("text", x = 4, y=210, label = "Hom. habitat",size = 2.5) +
    theme(legend.position = "none",axis.title = element_text(size = 10),axis.title.x=element_blank())
  
  pc <- ggplot(dt_patch_collapse_belowMSY, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs( x="Fishing aggregation level", y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 35, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 35, alpha = .2) + 
    annotate("text", x = 3.1, y=32, label = "Het. habitat",size = 2.5) +
    annotate("text", x = 4, y=32, label = "Hom. habitat",size = 2.5) +
    theme(legend.position = "none",axis.title = element_text(size = 10))
  
  p <- ggarrange(pa,pb,pc, nrow = 2, ncol = 2, vjust = 1.3, labels = c("A","B","C"),font.label = list(size = 12))
  
  myPPlot(p)
  
}

figureS3 <- function(){
  
  REC_A = 0.01
  settlers <- seq(0,3*DEF_K,10000)
  recruits = REC_A*settlers*(exp(-settlers/DEF_K))
  
  recruits_top = REC_A*DEF_K*(exp(-1))
  recruits_half = REC_A*DEF_K/2*(exp(-1/2))
  
  recruits_half/recruits_top
  
  appA <- data.frame(set = settlers, rec = recruits)
  
  p <- ggplot(appA, aes(x = set, y = rec))+
    geom_line() +
    labs(x = "Number of Settlers", y = "Number of Recruits") +
    ylim(0,50000)+
    annotate("segment", x = DEF_K, xend = DEF_K, y = 0, yend = 50000,
             colour = "red", linetype = 2)
  
  myPPlot(p)
  
  
}

appendixC_03 <- function(){
  
  MSY_03 <- read.csv('output/figure4C.csv')
  OVERMSY_03 <- read.csv('output/figure4D.csv')
  
  MSY_03$corr = 0.3
  MSY_03$h = 'MSY'
  
  OVERMSY_03$corr = 0.3
  OVERMSY_03$h = 'OVER'
  
  fig03 = rbind(MSY_03,OVERMSY_03)
  fig03$corr = as.factor(fig03$corr)
  levels(fig03$corr) <- c("AR = 0.3")
  
  p_MSY_coll <- ggplot(MSY_03, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "Fishing aggregation level",title="Exploitation at MSY", y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 22.5, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.5, ymin = -1, ymax = 22.5, alpha = .2) + 
    annotate("text", x = 3.1, y=21, label = "Het. habitat",size = 2.5) +
    annotate("text", x = 4, y=21, label = "Hom. habitat",size = 2.5) +
    theme(legend.position = "none",axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, vjust = 2))
  
  p_OVER_coll <- ggplot(OVERMSY_03, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "Fishing aggregation level",title="Exploitation over MSY", y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 110, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 110, alpha = .2) + 
    annotate("text", x = 3.1, y=105, label = "Het. habitat",size = 2.5) +
    annotate("text", x = 4, y=105, label = "Hom. habitat",size = 2.5) +
    theme(legend.position = "none",axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, vjust = 2))
  
  p_MSY_mass <- ggplot(MSY_03, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "Fishing aggregation level", y="Adult biomass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x=element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -1, ymax = 175, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -1, ymax = 175, alpha = .2) + 
    annotate("text", x = 3.1, y=165, label = "Het. habitat",size = 2.5) +
    annotate("text", x = 4, y=165, label = "Hom. habitat",size = 2.5) +
    theme(legend.position = "none",axis.title = element_text(size = 10))
  
  p_OVER_mass <- ggplot(OVERMSY_03, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(x = "Fishing aggregation level", y = "Adult biomass (tons)",fill = "Habitat") +
    scale_x_discrete(labels=c("None", "Low", "High", "None")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=8)) +
    stat_summary(fun.y = mean, geom = "point",show.legend = FALSE, shape=23,size=1.5, fill = 'white', color='black') +
    annotate("rect", xmin = 3.55, xmax = 4.4, ymin = -5, ymax = 175, alpha = .2) +
    annotate("rect", xmin = 0.55, xmax = 3.45, ymin = -5, ymax = 175, alpha = .2) + 
    annotate("text", x = 3.1, y=165, label = "Het. habitat",size = 2.5) +
    annotate("text", x = 4, y=165, label = "Hom. habitat",size = 2.5) +
    theme(legend.position = "none",axis.title = element_text(size = 10)) +
    ylim(-5, 175)
  
  p <- ggarrange(p_MSY_coll,p_OVER_coll,p_MSY_mass,p_OVER_mass,nrow = 2, ncol = 2,
                labels = c("A","B","C","D"),font.label = list(size = 12))
  myPPlot(p)
  
}

analysis_NoAllee <- function(){
  
  # MSY =  0.09, Catch = 24015.55
  
  # Calculate baseline (no harvest) to determine collapse
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, ar = AR, ma = MA,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, ar = AR, ma = MA,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=0.18, habitatHOM = TRUE, ar = AR, ma = MA,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST0 = runModel(hProp=0.18, habitatHOM = FALSE, ar = AR, ma = MA,
                            maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                            isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST0$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 0
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST0 = runModel(hProp=0.18, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                              isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST0$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 0
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST0 = runModel(hProp=0.18, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                              isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST0$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 0
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapse = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapse$bFactor = as.factor(dt_patch_collapse$bFactor)
  dt_patch_collapse$dFactor = as.factor(dt_patch_collapse$dFactor)
  dt_patch_collapse$habitat = as.factor(dt_patch_collapse$habitat)
  
  p <- ggplot(dt_patch_collapse, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=10))
  
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Total adult mass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  saveMatrixToCSVFile('figure3A', dt_patch_collapse)
  
  # Above MSY
  
  # Case 4
  
  OVER_BIO0_DIST0 = runModel(hProp=0.24, habitatHOM = FALSE, ar = AR, ma = MA,
                             maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                             isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                             hSDProp=0.1, classicSurv=FALSE,
                             runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                             showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                             limitOutput = TRUE)
  
  spatialAdultMass_4 = OVER_BIO0_DIST0$adultMass
  patch_collapse_4 = matrix(0,length(spatialAdultMass_4[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_4[,1])) {
    patch_collapse_4[i,] = spatialAdultMass_4[i,] < collapse_mass 
  }
  
  collapse_per_rep_4 = rowMeans(1*patch_collapse_4)
  dt_4 = data.frame("Patch collapse" = collapse_per_rep_4)
  dt_4$bFactor = 0
  dt_4$dFactor = 0
  dt_4$AdultMass = apply(spatialAdultMass_4, 1, sum)
  dt_4$catch = OVER_BIO0_DIST0$catch
  dt_4$habitat = 'HET'
  
  
  # Case 5
  
  OVER_BIO100_DIST0 = runModel(hProp=0.24, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                               isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_5 = OVER_BIO100_DIST0$adultMass
  patch_collapse_5 = matrix(0,length(spatialAdultMass_5[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_5[,1])) {
    patch_collapse_5[i,] = spatialAdultMass_5[i,] < collapse_mass 
  }
  
  collapse_per_rep_5 = rowMeans(1*patch_collapse_5)
  dt_5 = data.frame("Patch collapse" = collapse_per_rep_5)
  dt_5$bFactor = 100
  dt_5$dFactor = 0
  dt_5$AdultMass = apply(spatialAdultMass_5, 1, sum)
  dt_5$catch = OVER_BIO100_DIST0$catch
  dt_5$habitat = 'HET'
  
  
  #Case 6
  
  OVER_BIO200_DIST0 = runModel(hProp=0.24, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                               isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_6 = OVER_BIO200_DIST0$adultMass
  patch_collapse_6 = matrix(0,length(spatialAdultMass_6[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_6[,1])) {
    patch_collapse_6[i,] = spatialAdultMass_6[i,] < collapse_mass 
  }
  
  collapse_per_rep_6 = rowMeans(1*patch_collapse_6)
  dt_6 = data.frame("Patch collapse" = collapse_per_rep_6)
  dt_6$bFactor = 200
  dt_6$dFactor = 0
  dt_6$AdultMass = apply(spatialAdultMass_6, 1, sum)
  dt_6$catch = OVER_BIO200_DIST0$catch
  dt_6$habitat = 'HET'
  
  #Case 7
  
  OVER_BIO0_DIST0_HOM = runModel(hProp=0.24, habitatHOM = TRUE, ar = AR, ma = MA,
                                 maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                 isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                                 hSDProp=0.1, classicSurv=FALSE,
                                 runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                 showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                 limitOutput = TRUE)
  
  
  spatialAdultMass_7 = OVER_BIO0_DIST0_HOM$adultMass
  patch_collapse_7 = matrix(0,length(spatialAdultMass_7[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_7[,1])) {
    patch_collapse_7[i,] = spatialAdultMass_7[i,] < collapse_mass 
  }
  
  collapse_per_rep_7 = rowMeans(1*patch_collapse_7)
  dt_7 = data.frame("Patch collapse" = collapse_per_rep_7)
  dt_7$bFactor = 'HOM'
  dt_7$dFactor = 0
  dt_7$AdultMass = apply(spatialAdultMass_7, 1, sum)
  dt_7$catch = OVER_BIO0_DIST0_HOM$catch
  dt_7$habitat = 'HOM'
  
  dt_patch_collapse_overMSY = rbind(dt_4,dt_5,dt_6, dt_7)
  dt_patch_collapse_overMSY$bFactor = as.factor(dt_patch_collapse_overMSY$bFactor)
  dt_patch_collapse_overMSY$dFactor = as.factor(dt_patch_collapse_overMSY$dFactor)
  dt_patch_collapse_overMSY$habitat = as.factor(dt_patch_collapse_overMSY$habitat)
  
  dt_patch_collapse_overMSY = read.csv('output/figure3B.csv')
  
  p <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="Over MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Total Adult Mass (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=catch/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))+
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=5, color="black") 
  myPPlot(p)
  
  saveMatrixToCSVFile('figure3B', dt_patch_collapse_overMSY)
  
}

analysis_NoAllee_AR095 <- function(){
  
  set.seed(2020)
  # Build coast from parameters AR and MA
  corr1 = arima.sim(model = list(order = c(1, 0, 0),ar = 0.95, ma = MA), n = NUM_BLOCKS)
  #Bring to positive
  corr2 = corr1 - range(corr1)[1]
  corr2[corr2 == 0] = 0.1
  # Bring to right scale
  corrK = corr2*DEF_K/mean(corr2)
  # Get mean to corrK
  corrK = corrK-mean(corrK)+DEF_K
  
  
  baseline = runModel(hProp=0, habitatHOM = FALSE, ar = AR, ma = MA,
                      maxT=T_MAX, numReplic=5, 
                      isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                      hSDProp=0, classicSurv=FALSE,
                      runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                      showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                      limitOutput = TRUE)
  
  
  spatialAdultMassEqui = colMeans(baseline$adultMass)
  collapse_mass = COLLAPSE_THRESHOLD*spatialAdultMassEqui
  
  baseline_HOM = runModel(hProp=0, habitatHOM = TRUE, ar = AR, ma = MA,
                          maxT=T_MAX, numReplic=5, 
                          isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                          hSDProp=0, classicSurv=FALSE,
                          runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                          showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                          limitOutput = TRUE)
  
  
  spatialAdultMassEqui_HOM = colMeans(baseline_HOM$adultMass)
  collapse_mass_HOM = COLLAPSE_THRESHOLD*spatialAdultMassEqui_HOM
  
  #plot(spatialAdultMassEqui,type = 'n')
  #lines(1:150,spatialAdultMassEqui)
  
  # At MSY
  
  NUM_REP = 20
  
  # Case 0 
  
  MSY_BIO0_DIST0_HOM = runModel(hProp=0.18, habitatHOM = TRUE, ar = AR, ma = MA,
                                maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                                hSDProp=0.1, classicSurv=FALSE,
                                runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                limitOutput = TRUE)
  
  spatialAdultMass_0 = MSY_BIO0_DIST0_HOM$adultMass
  patch_collapse_0 = matrix(0,length(spatialAdultMass_0[,1]),NUM_BLOCKS)
  
  
  for (i in 1:length(spatialAdultMass_0[,1])) {
    patch_collapse_0[i,] = spatialAdultMass_0[i,] < collapse_mass_HOM 
  }
  
  collapse_per_rep_0 = rowMeans(1*patch_collapse_0)
  dt_0 = data.frame("Patch collapse" = collapse_per_rep_0)
  dt_0$bFactor = 'HOM'
  dt_0$dFactor = 0
  dt_0$AdultMass = apply(spatialAdultMass_0, 1, sum)
  dt_0$catch = MSY_BIO0_DIST0_HOM$catch
  dt_0$habitat = 'HOM'
  
  
  # Case 1
  
  MSY_BIO0_DIST0 = runModel(hProp=0.18, habitatHOM = FALSE, ar = AR, ma = MA,
                            maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                            isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                            hSDProp=0.1, classicSurv=FALSE,
                            runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                            showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                            limitOutput = TRUE)
  
  
  spatialAdultMass_1 = MSY_BIO0_DIST0$adultMass
  patch_collapse_1 = matrix(0,length(spatialAdultMass_1[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_1[,1])) {
    patch_collapse_1[i,] = spatialAdultMass_1[i,] < collapse_mass 
  }
  
  collapse_per_rep_1 = rowMeans(1*patch_collapse_1)
  dt_1 = data.frame("Patch collapse" = collapse_per_rep_1)
  dt_1$bFactor = 0
  dt_1$dFactor = 0
  dt_1$AdultMass = apply(spatialAdultMass_1, 1, sum)
  dt_1$catch = MSY_BIO0_DIST0$catch
  dt_1$habitat = 'HET'
  
  
  # Case 2
  
  MSY_BIO100_DIST0 = runModel(hProp=0.18, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                              isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_2 = MSY_BIO100_DIST0$adultMass
  patch_collapse_2 = matrix(0,length(spatialAdultMass_2[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_2[,1])) {
    patch_collapse_2[i,] = spatialAdultMass_2[i,] < collapse_mass 
  }
  
  collapse_per_rep_2 = rowMeans(1*patch_collapse_2)
  dt_2 = data.frame("Patch collapse" = collapse_per_rep_2)
  dt_2$bFactor = 100
  dt_2$dFactor = 0
  dt_2$AdultMass = apply(spatialAdultMass_2, 1, sum)
  dt_2$catch = MSY_BIO100_DIST0$catch
  dt_2$habitat = 'HET'
  
  #Case 3
  
  MSY_BIO200_DIST0 = runModel(hProp=0.18, habitatHOM = FALSE, ar = AR, ma = MA,
                              maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                              isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                              hSDProp=0.1, classicSurv=FALSE,
                              runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                              showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                              limitOutput = TRUE)
  
  spatialAdultMass_3 = MSY_BIO200_DIST0$adultMass
  patch_collapse_3 = matrix(0,length(spatialAdultMass_3[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_3[,1])) {
    patch_collapse_3[i,] = spatialAdultMass_3[i,] < collapse_mass 
  }
  
  collapse_per_rep_3 = rowMeans(1*patch_collapse_3)
  dt_3 = data.frame("Patch collapse" = collapse_per_rep_3)
  dt_3$bFactor = 200
  dt_3$dFactor = 0
  dt_3$AdultMass = apply(spatialAdultMass_3, 1, sum)
  dt_3$catch = MSY_BIO200_DIST0$catch
  dt_3$habitat = 'HET'
  
  dt_patch_collapse = rbind(dt_0,dt_1,dt_2,dt_3)
  dt_patch_collapse$bFactor = as.factor(dt_patch_collapse$bFactor)
  dt_patch_collapse$dFactor = as.factor(dt_patch_collapse$dFactor)
  dt_patch_collapse$habitat = as.factor(dt_patch_collapse$habitat)
  
  p <- ggplot(dt_patch_collapse, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() + 
    theme(axis.text.x =element_text(size=10))
  
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Total adult mass (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse, aes(x=bFactor, y=catch/1000, fill = habitat)) + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="MSY", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = 'Habitat') +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  saveMatrixToCSVFile('NoAllee_AR095_MSY', dt_patch_collapse)
  
  # Above MSY
  
  # Case 4
  
  OVER_BIO0_DIST0 = runModel(hProp=0.24, habitatHOM = FALSE, ar = AR, ma = MA,
                             maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                             isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                             hSDProp=0.1, classicSurv=FALSE,
                             runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                             showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                             limitOutput = TRUE)
  
  spatialAdultMass_4 = OVER_BIO0_DIST0$adultMass
  patch_collapse_4 = matrix(0,length(spatialAdultMass_4[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_4[,1])) {
    patch_collapse_4[i,] = spatialAdultMass_4[i,] < collapse_mass 
  }
  
  collapse_per_rep_4 = rowMeans(1*patch_collapse_4)
  dt_4 = data.frame("Patch collapse" = collapse_per_rep_4)
  dt_4$bFactor = 0
  dt_4$dFactor = 0
  dt_4$AdultMass = apply(spatialAdultMass_4, 1, sum)
  dt_4$catch = OVER_BIO0_DIST0$catch
  dt_4$habitat = 'HET'
  
  
  # Case 5
  
  OVER_BIO100_DIST0 = runModel(hProp=0.24, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 100, dist_factor = 0,
                               isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_5 = OVER_BIO100_DIST0$adultMass
  patch_collapse_5 = matrix(0,length(spatialAdultMass_5[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_5[,1])) {
    patch_collapse_5[i,] = spatialAdultMass_5[i,] < collapse_mass 
  }
  
  collapse_per_rep_5 = rowMeans(1*patch_collapse_5)
  dt_5 = data.frame("Patch collapse" = collapse_per_rep_5)
  dt_5$bFactor = 100
  dt_5$dFactor = 0
  dt_5$AdultMass = apply(spatialAdultMass_5, 1, sum)
  dt_5$catch = OVER_BIO100_DIST0$catch
  dt_5$habitat = 'HET'
  
  
  #Case 6
  
  OVER_BIO200_DIST0 = runModel(hProp=0.24, habitatHOM = FALSE, ar = AR, ma = MA,
                               maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 200, dist_factor = 0,
                               isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                               hSDProp=0.1, classicSurv=FALSE,
                               runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                               showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                               limitOutput = TRUE)
  
  spatialAdultMass_6 = OVER_BIO200_DIST0$adultMass
  patch_collapse_6 = matrix(0,length(spatialAdultMass_6[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_6[,1])) {
    patch_collapse_6[i,] = spatialAdultMass_6[i,] < collapse_mass 
  }
  
  collapse_per_rep_6 = rowMeans(1*patch_collapse_6)
  dt_6 = data.frame("Patch collapse" = collapse_per_rep_6)
  dt_6$bFactor = 200
  dt_6$dFactor = 0
  dt_6$AdultMass = apply(spatialAdultMass_6, 1, sum)
  dt_6$catch = OVER_BIO200_DIST0$catch
  dt_6$habitat = 'HET'
  
  #Case 7
  
  OVER_BIO0_DIST0_HOM = runModel(hProp=0.24, habitatHOM = TRUE, ar = AR, ma = MA,
                                 maxT=T_MAX, numReplic=NUM_REP, biomass_factor = 0, dist_factor = 0,
                                 isAllee=FALSE, randDisp = TRUE, randRec= TRUE,
                                 hSDProp=0.1, classicSurv=FALSE,
                                 runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                                 showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                                 limitOutput = TRUE)
  
  
  spatialAdultMass_7 = OVER_BIO0_DIST0_HOM$adultMass
  patch_collapse_7 = matrix(0,length(spatialAdultMass_7[,1]),NUM_BLOCKS)
  
  for (i in 1:length(spatialAdultMass_7[,1])) {
    patch_collapse_7[i,] = spatialAdultMass_7[i,] < collapse_mass 
  }
  
  collapse_per_rep_7 = rowMeans(1*patch_collapse_7)
  dt_7 = data.frame("Patch collapse" = collapse_per_rep_7)
  dt_7$bFactor = 'HOM'
  dt_7$dFactor = 0
  dt_7$AdultMass = apply(spatialAdultMass_7, 1, sum)
  dt_7$catch = OVER_BIO0_DIST0_HOM$catch
  dt_7$habitat = 'HOM'
  
  dt_patch_collapse_overMSY = rbind(dt_4,dt_5,dt_6, dt_7)
  dt_patch_collapse_overMSY$bFactor = as.factor(dt_patch_collapse_overMSY$bFactor)
  dt_patch_collapse_overMSY$dFactor = as.factor(dt_patch_collapse_overMSY$dFactor)
  dt_patch_collapse_overMSY$habitat = as.factor(dt_patch_collapse_overMSY$habitat)
  
  
  p <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=Patch.collapse*100, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="Over MSY", x=expression(atop("", "Fishing heterogeneity")), y="Extinct patches (%)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=AdultMass/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Total Adult Mass (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))
  myPPlot(p)
  
  p <- ggplot(dt_patch_collapse_overMSY, aes(x=bFactor, y=catch/1000, fill = habitat)) +
    geom_dotplot(binaxis='y', stackdir='center', stackratio=0.6, dotsize=0.8) + 
    labs(title="", x=expression(atop("", "Fishing heterogeneity")), y="Catch (tons)", fill = "Habitat") +
    scale_x_discrete(labels=c("Null", "Intermediate", "High", "Null")) + 
    scale_fill_discrete(labels=c("Heterogeneous", "Homogeneous"))+
    theme_minimal() +
    theme(axis.text.x =element_text(size=10))+
    stat_summary(fun.y = median, geom = "point",show.legend = FALSE, shape=18,size=5, color="black") 
  myPPlot(p)
  
  saveMatrixToCSVFile('figure3B', dt_patch_collapse_overMSY)
  
  
}

tuneUnfishedDensity <- function(targD=0.4, numReplic=3, maxT=50) {
  dispSDM = getRandomDispSDMatrix(maxT, numReplic, rType="gamma")
  # loop to find K value
  theK = DEF_K
  pL   = getParamList(hProp=0, mpaProp=0, blockK=theK, maxT=maxT, numReplic=numReplic, addCat=FALSE)
  totalRuns = 0
  aveFrom   = maxT-20
  adjBy     = 2
  wasLarger = FALSE
  while (TRUE) {
    myP("Running with K=", theK)
    pL$blockK = theK
    noHarvL   = runModelFromList(pL, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM)
    meanD     = noHarvL$means
    meanDens  = mean(meanD[meanD$t>aveFrom,]$meanSSB) / (TOTAL_AREA*10000)
    myP("Density:", meanDens)
    if (my.eq(targD, meanDens, minDiff=targD*0.01)) break
    else if (totalRuns>50) {
      myP("Max runs reached.")
      break
    } else if (meanDens>targD) {
      if (!wasLarger) {
        adjBy = adjBy * 0.5
        wasLarger = TRUE
      }
      theK = theK / (1+adjBy)
    } else {
      if (wasLarger) {
        adjBy = adjBy * 0.5
        wasLarger = FALSE
      }
      theK = theK * (1+adjBy)
    }
  }
  
  myP("Final K:", theK, ", density=", meanDens)
}

########################
##  MSY Calculation
########################

recurseFindMSY <- function(paramL, oldH, hInc=0.01, oldCatch=-1, dispSDM=NULL, randRecM=NULL, aveOver=40) {
  #if (is.null(dispSDM)) dispSDM = getRandomDispSDMatrix(paramL$maxT, paramL$numReplic, rType="gamma");
  #if (is.null(randRecM)) randRecM = getRandomRecruitMatrix(paramL$maxT, paramL$numReplic, rType="log-normal", shape=REC_LN_SDLOG);
  
  lastT     = paramL$maxT;
  aveV      = (lastT-aveOver+1):lastT;
  firstTest = (oldCatch==-1)
  
  if (firstTest) {
    myP("Generating starting value for h=", oldH)
    hParamL         = paramL;
    hParamL$hProp   = oldH;
    hParamL$hPreSim = oldH;
    meanD = runModelFromList(hParamL, showPlots=FALSE, saveToFile=FALSE, dispSDM=NULL, randRecM=NULL, prefix=combine("test",oldH))$means;
    oldCatch = mean(meanD$meanCatch[aveV]);
  } 
  
  newH = oldH + hInc;
  if ((newH<0) || (newH>1)) return(list(msyH=oldH, msyCatch=oldCatch))
  
  hParamL         = paramL;
  hParamL$hProp   = newH;
  hParamL$hPreSim = newH;
  myP("Running test catch for h=", newH)
  meanD = runModelFromList(hParamL, showPlots=FALSE, saveToFile=FALSE, dispSDM=NULL, randRecM=NULL, prefix=combine("test",newH))$means;
  newCatch = mean(meanD$meanCatch[aveV]);
  
  if (newCatch>oldCatch) {
    # recurse on next value (either up 1 or down 1)
    # make sure to pass on hInc since it could be negative now
    myP("New catch greater. Recursing...")
    return(recurseFindMSY(paramL, oldH=newH, hInc=hInc, oldCatch=newCatch, dispSDM=NULL, randRecM=NULL))
  } else if (newCatch<=oldCatch) {
    # don't increase any more
    if (firstTest) {
      # if first, then go down.  Otherwise, we came from down
      myP("New catch smaller. Decrementing...")
      return(recurseFindMSY(paramL, oldH=oldH, hInc=-0.01, oldCatch=oldCatch, dispSDM=NULL, randRecM=NULL))
    } else {
      # not the first, so we've already looked above and come from or looked at below
      myP("New catch smaller. Stopping...")
      return(list(msyH=oldH, msyCatch=oldCatch))
    }
  }
}

#####################################
#           RUN MODEL
#####################################

getParamList <- function(mpaSize=5, mpaProp=0, hProp=0.09, hPresim=H_PRESIM, preK = 0.1, habitatHOM = TRUE, hetK = corrK,
                         biomass_factor = BIOMASS_FACTOR, dist_factor = DIST_FACTOR, maxT=T_MAX, numReplic=3, classWidth=CLASS_WIDTH, 
                         isAllee=FALSE, aggA=AGG_A, randDisp = TRUE, randRec=TRUE,
                         hSDProp=0.1, classicSurv=FALSE) {
  return(list(mpaSize=mpaSize, mpaProp=mpaProp, hProp=hProp, hPreSim=hPresim, preK = preK, habitatHOM = habitatHOM, hetK = hetK,
              biomass_factor = biomass_factor, dist_factor = dist_factor,maxT=maxT, numReplic=numReplic, classWidth=classWidth, isAllee=isAllee, aggA=aggA, randDisp = randDisp, randRec=randRec,
              hSDProp=hSDProp, classicSurv=classicSurv))  
} 

runModelFromList <- function(paramL, dispSDM=NULL, randRecM=NULL, runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE, showPlots=FALSE, saveToFile=TRUE, createMovie=FALSE, prefix="") {
  runModel(mpaSize=paramL$mpaSize, mpaProp=paramL$mpaProp, hProp=paramL$hProp, hPreSim=paramL$hPreSim, preK = paramL$preK, habitatHOM = paramL$habitatHOM, hetK=paramL$hetK, 
           biomass_factor = paramL$biomass_factor, dist_factor = paramL$dist_factor, 
           maxT=paramL$maxT, numReplic=paramL$numReplic, classWidth=paramL$classWidth, isAllee=paramL$isAllee, aggA=paramL$aggA, randDisp = paramL$randDisp, randRec=paramL$randRec,
           dispSDM=dispSDM, randRecM=randRecM, hSDProp=paramL$hSDProp, classicSurv=paramL$classicSurv, 
           runPreSim=runPreSim, showPreSim=showPreSim, runPreSimOnly=runPreSimOnly, showPlots=showPlots, saveToFile=saveToFile, createMovie=createMovie, prefix=prefix)
} 

runModel <- function(mpaSize=5, mpaProp=0, hProp=0.09, hPreSim=H_PRESIM, preK = 1, 
                     habitatHOM = FALSE, hetK = corrK, biomass_factor = BIOMASS_FACTOR, dist_factor = DIST_FACTOR,
                     maxT=T_MAX, numReplic=NUM_REPLICATIONS, classWidth=CLASS_WIDTH,
                     isAllee=TRUE, aggA=AGG_A, randDisp = TRUE, randRec= TRUE,
                     hSDProp=0, dispSDM=NULL, randRecM=NULL, classicSurv=FALSE,
                     runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
                     showPlots=FALSE, saveToFile=FALSE, createMovie= FALSE, prefix="",
                     limitOutput = FALSE) {
  
  initializeGlobalValues(classWidth)
  
  if (runPreSim) {
    myP("Running preliminary simulation...")
    myP("Replicates:", numReplic, "Time:", T_PRESIM)
    myP("Harvest rate (after 2/3 time):", hPreSim)
    preSimAbunV = runPreSimulationSpatial(habitatHOM = habitatHOM, hetK = corrK, isAllee, aggA)
  } else { 
    myP("Skipping preliminary simulation...")
    myP("Setting initial abundance to", INIT_ABUN_NO_PRESIM,"per size...")
    preSimAbunV = rep(INIT_ABUN_NO_PRESIM, NUM_CLASSES)
  } 
  
  if (!runPreSimOnly) {
    myP("Running simulation...")
    myP("Replicates:", numReplic, "Time:", maxT)
    myP("Proportion of coastline protected:", mpaProp)
    myP("Reserve block size:", mpaSize)
    myP("Harvest proportion:", hProp, "with SD proportion", hSDProp)
    
    modelResL = runSimulation(mpaSize = mpaSize, mpaProp = mpaProp, hProp = hProp, initAbunV = preSimAbunV, maxT = maxT, numReplic = numReplic, 
                              biomass_factor = biomass_factor, dist_factor = dist_factor, habitatHOM = habitatHOM, hetK = hetK, 
                              isAllee = isAllee, aggA = aggA, randDisp = randDisp, randRec = randRec, hSDProp = hSDProp, dispSDM = dispSDM, randRecM = randRecM, classicSurv = classicSurv)
    
    # Output from runSimulation: totalAbunM, totalAbunSSBM, totalCatchM, spatialMassA, spatialSSBA, spatialCatchA, spatialSSBDensA,
    # compSizeDataA, actualEffort
    
    spatialAdultMass = apply(modelResL$spatialSSBA[,(T_MAX - 9):T_MAX,],c(1,3),mean)
    totalCatch = apply(modelResL$totalCatchM[,(T_MAX - 9):T_MAX], 1, mean)
    
    
    effort = apply(modelResL$actualEffort, c(2,3), mean)
    effortT = apply(effort,1,sum)
    
    profit = getProfit(modelResL,maxT)
    
    # save size structure matrix
    #    compSizeDataA   = array(0, c(numReplic,NUM_CLASSES,maxT,NUM_BLOCKS)) # complete size/space data
    # sum across blocks
    compSizeRepA = apply(modelResL$compSizeDataA, c(1,2,3), sum)
    # ave across replicates	
    compSizeM    = apply(compSizeRepA, c(2,3), mean)
    saveMatrixToCSVFile(combine(prefix, "size"), compSizeM)
    
    # plot all graphs and save the matrices
    if (showPlots) {plotResults(modelResL, showPlots, saveToFile, createMovie, prefix)}
    
    if (limitOutput){
      
      return(list(adultMass = spatialAdultMass, catch = totalCatch))
      
    } else{
      
      return(list(means=getMeanResults(modelResL, isAllee, prefix), res=modelResL, profit = profit))
      
    }
    
  }
}

getProfit <- function(modelResL, maxT = maxT) {
  
  FIXED_COST = 5000 # Dollars per boat and per year (license, wages, maintenance, insurance)
  HARVEST_COST = 1000 # Dollars per boat and per year (fuel,ice,supplies)
  PRICE = 50 # dollars/kg 
  DISCOUNT_RATE = 0.95
  
  meanCatch = colMeans(modelResL$totalCatchM)
  timeEffort = apply(modelResL$actualEffort, c(1,2), sum)
  meanEffort = colMeans(timeEffort)
  
  totalProfit = 0
  
  
  for (i in (maxT/2):maxT) {
    
    costs = HARVEST_COST*meanEffort[i] + FIXED_COST*TOTAL_EFFORT
    profit = meanCatch[i]*PRICE - costs 
    
    totalProfit  =  totalProfit + profit*DISCOUNT_RATE^i
  }
  
  return(totalProfit)
  
}

# Run replicated pre-simulation with spatial structure
# Returns a vector of mean abundance for each size
runPreSimulationSpatial <- function(habitatHOM = habitatHOM, hetK=corrK, isAllee=TRUE, aggA=AGG_A) {
  
  densPreA   = array(0, c(NUM_CLASSES,T_PRESIM,NUM_BLOCKS))   # Pre-simulation density
  eggsProdV = rep(0, NUM_BLOCKS)                              # Eggs produced vector
  eggsArrV  = rep(0, NUM_BLOCKS)                              # Eggs arrived vector
  eggDispM  = matrix(0, nrow=NUM_BLOCKS, ncol=NUM_BLOCKS)     # Egg dispersal matrix
  
  # Set initial abundance across replicates
  densPreA[,1,] = matrix(rep(INIT_ABUN_PRESIM, floor(NUM_CLASSES)*NUM_BLOCKS),ncol = NUM_BLOCKS);           
  
  # Adjust by 1/2 class_width to get lower limit
  lowerSizeV = gMeanSizeV - CLASS_WIDTH/2
  P = getTransitionMatrix(lowerSizeV)                 # Growth matrix
  S = getSurvivalMatrix(lowerSizeV, P)   # Survival matrix
  PSMatrix = P*S                                      # in modern version, survival is based on both initial and final length
  
  # put these in global to compare with other P,S
  assign("newP", P, envir=.GlobalEnv)
  assign("newS", S, envir=.GlobalEnv)
  
  if (habitatHOM) recK = HOM_K
  else recK = hetK
  
  # loop through presim time
  for (t in 1:(T_PRESIM-1)) {
    
    if (t%%10==0) myP("Time step t =", t, "...")
    
    dispSD = 0
    
    # Calculate egg dispersal matrix for this time step			
    for (j in 1:(NUM_BLOCKS/2)) { 
      eggDispM[,j] = propDispToBlock(c((1-j):(NUM_BLOCKS-j)), dispSD) + 
        rev(propDispToBlock(c(j:(NUM_BLOCKS+(j-1))), dispSD))
    } 
    
    # adds circular nature to coastline
    eggDispM[1:NUM_BLOCKS,(NUM_BLOCKS/2+1):NUM_BLOCKS] = eggDispM[NUM_BLOCKS:1,(NUM_BLOCKS/2):1]
    
    # apply Allee effect to each block, if any
    if (isAllee) effDensM  = getEffectiveDensityMatrix(densPreA[,t,], aggA)
    else effDensM = densPreA[,t,]
    
    eggsProdV = apply(sweep(effDensM,gEggsPerIndV,MARGIN=1,'*'), 2, sum) # vectorized 
    eggsArrV  = apply(sweep(eggDispM,MARGIN=2,eggsProdV,'*'), 1, sum) 
    
    for (block in 1:NUM_BLOCKS) {
    
      # apply larval survival
      numSettlers = LARVAL_SURVIVAL*eggsArrV[block]
  
      recruitV = getRecruitV(numSettlers,recK[block])
      
      densPreA[,t+1,block] = PSMatrix%*%densPreA[,t,block] + recruitV
      
    }
    
  } 
  
  # take mean for each size class across replicates
  preSimAbunV = densPreA[,T_PRESIM,]
  
  return(preSimAbunV)
  
}

###????????????????????? SIMULATION MODEL???????????????????????????????????????????????

# function that returns matrix of totalAbunM, totalCatchM, totalCatchMassM                                                          
# given the spatial arrangement
# arguments are:
# I. harvest proportion (from 0 to 1)
# Returns list with three elements:
#  -totalAbunM			Total numerical abundance
#  -totalAbunSSBM		Total number of mature individuals
#  -totalCatchM         Total catch biomass
#  ALSO RETURNS SPATIAL VALUES
#  Dimensions: row=replicate, col=time
runSimulation <- function(mpaSize, mpaProp, hProp, initAbunV, maxT=T_MAX, numReplic=NUM_REPLICATIONS, 
                          biomass_factor = BIOMASS_FACTOR, dist_factor = DIST_FACTOR,   
                          habitatHOM = TRUE, hetK = corrK, isAllee=TRUE, aggA=AGG_A, randDisp=FALSE, randRec=FALSE,
                          hSDProp=0, dispSDM=NULL, randRecM=NULL, classicSurv=FALSE) { 
  
  # Arrays for model dynamics
  # IMPORTANT: density is #/ha, not #/m^2
  densA     = array(0, c(NUM_CLASSES,maxT,NUM_BLOCKS))    # Initialize Dens_SEXPL 
  catchA    = array(0, c(NUM_CLASSES,maxT,NUM_BLOCKS))    # Initialize Catch_SEXPL 
  TACCatch  = array(0, c(NUM_CLASSES,maxT,NUM_BLOCKS))    
  potCatch  = array(0, c(NUM_CLASSES,maxT,NUM_BLOCKS))
  eggsProdV = rep(0, NUM_BLOCKS)                          # Eggs produced vector
  eggsArrV  = rep(0, NUM_BLOCKS)                          # Eggs arrived vector
  eggDispM  = matrix(0, nrow=NUM_BLOCKS, ncol=NUM_BLOCKS) # Egg dispersal matrix
  distance = rep(0, NUM_BLOCKS)
  
  if (habitatHOM == TRUE) recK = HOM_K
  else {
    recK = hetK
  }
  
  # Adjust by 1/2 class_width to get lower limit
  lowerSizeV = gMeanSizeV - CLASS_WIDTH/2
  
  P = getTransitionMatrix(lowerSizeV)     # Growth matrix
  S = getSurvivalMatrix(lowerSizeV, P, classicSurv)
  
  PSMatrix = P*S  # survival matrix based on initial and final length
  
  # get random dispersal SD value matrix
  # shape: time x replicates
  if (!is.null(dispSDM)) randDisp = TRUE
  else dispSDM = getRandomDispSDMatrix(maxT, numReplic, rType="gamma")
  
  if (!is.null(randRecM)) randRec = TRUE
  else randRecM = getRandomRecruitMatrix(maxT, numReplic, rType="log-normal")
  
  # Matrices for summed results
  # Row = replicate
  # Col = time
  totalAbunM    = matrix(0, nrow=numReplic, ncol=maxT)  # Total abundance (individuals)
  totalAbunSSBM = matrix(0, nrow=numReplic, ncol=maxT)  # Total abundance (breeders)
  totalCatchM   = matrix(0, nrow=numReplic, ncol=maxT)  # Total catch (weight in tons)
  
  # Arrays for spatial results (mass)
  # dim = rep, time, blocks
  spatialMassA    = array(0, c(numReplic,maxT,NUM_BLOCKS))  # Spatial mass (individuals)
  spatialSSBA     = array(0, c(numReplic,maxT,NUM_BLOCKS))  # Spatial mass (breeders)
  spatialCatchA   = array(0, c(numReplic,maxT,NUM_BLOCKS))  # Spatial catch (weight in tons)
  spatialCatchShellA = array(0, c(numReplic,maxT,NUM_BLOCKS)) # Spatial catch (weight in tons) 
  spatialCatchInd = array(0, c(numReplic,maxT,NUM_BLOCKS))  # Spatial catch (individuals)
  spatialDensA = array(0, c(numReplic,maxT,NUM_BLOCKS))     # Spatial density (total)
  spatialSSBDensA = array(0, c(numReplic,maxT,NUM_BLOCKS))  # Spatial density (breeders)
  actualEffort = array(0, c(numReplic, maxT, NUM_BLOCKS))   # Spatial fishing effort
  compSizeDataA   = array(0, c(numReplic,NUM_CLASSES,maxT,NUM_BLOCKS)) # complete size/space data
  densA_eggs = array(0, c(numReplic,maxT,NUM_BLOCKS)) # Number of eggs getting to each block
  
  # Run simulation
  
  harvestableV   = gMeanSizeV>MIN_LAND_SIZE
  
  # Get effort vector across blocks (h=0 for reserves)
  # Determines reserve structure
  blockEffortV   = getBlockEffortVector(mpaSize, mpaProp, hProp) # num[1:NUM_BLOCKS], 0 or hprop
  
  for (rep in 1:numReplic) { # Start loop on replicates
    
    myP("Starting replicate", rep, "...")
    densA[,1,] = initAbunV
    
    for (t in 1:(maxT-1)) {  # Start loop on time
      
      if (t%%10==0) myP("Time step t =", t, "...")
      
      # Get random dispersal SD value
      if (randDisp) dispSD = getDispSD(dispSDM[t,rep]) 
      else dispSD = 0
      
      # Calculate egg dispersal matrix for this time step									 
      for (j in 1:(NUM_BLOCKS/2)) { 
        eggDispM[,j] = propDispToBlock(c((1-j):(NUM_BLOCKS-j)), dispSD) + 
          rev(propDispToBlock(c(j:(NUM_BLOCKS+(j-1))), dispSD))
      } 
      
      # adds circular nature to coastline
      eggDispM[1:NUM_BLOCKS,(NUM_BLOCKS/2+1):NUM_BLOCKS] = eggDispM[NUM_BLOCKS:1,(NUM_BLOCKS/2):1]
      
      # apply Allee effect to each block, if any
      
      if (isAllee) effDensM  = getEffectiveDensityMatrix(densA[,t,], aggA)
      else effDensM = densA[,t,]
      eggsProdV = apply(sweep(effDensM,gEggsPerIndV,MARGIN=1,'*'), 2, sum) # vectorized 
      eggsArrV  = apply(sweep(eggDispM,MARGIN=2,eggsProdV,'*'), 1, sum)    # vectorized
      if (isAllee) eggsArrV = eggsArrV*LARVAL_ALLEE_ADJ
      
      # Store number of eggs
      
      densA_eggs[rep,t,] = eggsArrV
      
      # Calculate mature biomass
      biomass = sweep(densA[,t,]*gPropMatureV, wForL(gMeanSizeV)/1000, MARGIN=1, '*') # num[1:NUM_CLASSES, 1:NUM_BLOCK]
      biomass.block = apply(biomass, 2, sum) # num[1:NUM_BLOCK]
      biomass.total = sum(biomass.block) # num 
      
      for (block in 1:NUM_BLOCKS) {
        distance[block] = abs(NUM_BLOCKS/2 - block)/(NUM_BLOCKS/2)
      }
      
      biomass.weight = exp(-biomass_factor*(1-biomass.block/sum(biomass.block)) - dist_factor*(distance)) # num[1:NUM_BLOCK]
      effort.block = TOTAL_EFFORT*biomass.weight/sum(biomass.weight) # num[1:NUM_BLOCK]
      
      potentialH = 1 - exp(-CATCHABILITY*effort.block) # num[1:NUM_BLOCK]
      
      
      # Calculate TAC and potential catch value for aggregated fishing effort
      
      blockH = numeric(NUM_BLOCKS)
      potBlockH = numeric(NUM_BLOCKS)
      
      for (block in 1:NUM_BLOCKS) {
        
        blockH[block] = blockEffortV[block]
        
        if (blockEffortV[block] == 0) {potBlockH[block] = 0 } # num[1:NUM_BLOCK]
        else {potBlockH[block] = potentialH[block]}
        
        TACCatch[,t,block]  = (diag(1,NUM_CLASSES) - getHarvestMatrix(blockH[block],verbose=FALSE))%*%densA[,t,block]
        potCatch[,t,block] = (diag(1,NUM_CLASSES) - getHarvestMatrix(potBlockH[block],verbose=FALSE))%*%densA[,t,block]
        
      } 
      
      # Calculate TAC if h is applied uniformly
      TAC.class = sweep(TACCatch[,t,], wForL(gMeanSizeV)/1000, MARGIN=1, '*') # num[1:NUM_CLASSES, 1:NUM_BLOCK]
      TAC.block = apply(TAC.class, 2, sum) #num[1:NUM_BLOCK]
      TAC.total = sum(TAC.block) #num
      
      # Calculate total catch for aggregated fishing effort 
      potCatch.class = sweep(potCatch[,t,], wForL(gMeanSizeV)/1000, MARGIN=1, '*') # num[1:NUM_CLASSES, 1:NUM_BLOCK]
      potCatch.block = apply(potCatch.class, 2, sum) #num[1:NUM_BLOCK]
      potCatch.total = sum(potCatch.block) #num
      
      # If potential catch exceeds TAC, reduce harvest season
      if (potCatch.total <= TAC.total) {
        actualH = potBlockH
        #print('Below TAC')
      }
      else {
        harvSeason = TAC.total/potCatch.total
        actualH = potBlockH*harvSeason
        #print('Above TAC')
      }  
      
      actualEffort[rep,t,] = -log(1-actualH)/CATCHABILITY
      
      # Apply growth, survival, and harvest and add in new recruits
      for (block in 1:NUM_BLOCKS) {
        
        # pull H from unique vector
        if (actualH[block]==0) blockH = diag(NUM_CLASSES)
        else {
          blockH = getHarvestMatrix(actualH[block], hSDProp, verbose=FALSE)
        } 
        
        # apply larval survival
        numSettlers = LARVAL_SURVIVAL*eggsArrV[block]
        
        # add random recruitment pre-DD
        if (randRec) numSettlers = numSettlers*randRecM[t,rep]
        recruitV = getRecruitV(numSettlers, recK[block])
        
        #		  if (randRec) recruitV = recruitV*randRecM[t,rep]
        # adjust if random recruitment survival
        densA[,t+1,block] = ((PSMatrix%*%blockH))%*%densA[,t,block] + recruitV
        # include all sizes, even non-harvested
        catchA[,t,block]  = (diag(1,NUM_CLASSES) - blockH)%*%densA[,t,block]
      }
      
      
    } # End loop on time
    
    for (block in 1:NUM_BLOCKS) {
      blockH = getHarvestMatrix(actualH[block], hSDProp, verbose=FALSE)
      catchA[,maxT,block] = (diag(1,NUM_CLASSES)-blockH)%*%densA[,maxT,block]
    } 
    
    # Weight of catch [kg]
    abunWeightA       = sweep(densA, wForL(gMeanSizeV)/1000, MARGIN=1, '*') 
    ssbWeightA        = sweep(densA[,,]*gPropMatureV, wForL(gMeanSizeV)/1000, MARGIN=1, '*') 
    catchWeightShellA = sweep(catchA, wForL(gMeanSizeV)/1000, MARGIN=1, '*')
    catchWeightA      = sweep(catchA, wForL(gMeanSizeV)*W_NOSHELL/1000, MARGIN=1, '*') 
    
    
    # mass
    # row, col now time, block
    spatialMassA[rep,,]  = apply(abunWeightA, c(2,3), sum)
    spatialSSBA[rep,,]   = apply(ssbWeightA, c(2,3), sum)
    blockCatchMassM      = apply(catchWeightA, c(2,3), sum)
    spatialCatchA[rep,,] = blockCatchMassM
    spatialCatchShellA[rep,,] = apply(catchWeightShellA, c(2,3), sum)
    
    # still density here
    # row, col now time, block
    blockDensM      = apply(densA, c(2,3), sum)
    blockDensSSBM   = apply(densA[,,]*gPropMatureV, c(2,3), sum)
    spatialCatchInd[rep,,] = apply(catchA, c(2,3), sum) 
    
    totalAbunM[rep,]    = apply(blockDensM, 1, sum)*BLOCK_AREA
    totalAbunSSBM[rep,] = apply(blockDensSSBM, 1, sum)*BLOCK_AREA
    totalCatchM[rep,]   = apply(blockCatchMassM, 1, sum)*BLOCK_AREA
    
    
    spatialDensA[rep,,] = blockDensM # not summed across blocks
    spatialSSBDensA[rep,,] = blockDensSSBM # not summed across blocks
    compSizeDataA[rep,,,]  = densA # not summed across blocks or sizes
  } # End loop on replicates 
  
  
  return(list(totalAbunM=totalAbunM, totalAbunSSBM=totalAbunSSBM, totalCatchM=totalCatchM,
              spatialMassA=spatialMassA, spatialSSBA=spatialSSBA, spatialCatchInd = spatialCatchInd, spatialCatchA=spatialCatchA, spatialCatchShellA = spatialCatchShellA, 
              spatialDensA = spatialDensA,spatialSSBDensA=spatialSSBDensA, compSizeDataA=compSizeDataA, actualEffort = actualEffort, densA_eggs = densA_eggs))
}  

plotResults <- function(modelResL, showPlots=TRUE, saveToFile=FALSE, createMovie=TRUE, prefix="") {
  # Dimensions: row=replicate, col=time
  totalAbunM    = modelResL$totalAbunM
  totalAbunSSBM = modelResL$totalAbunSSBM
  totalCatchM   = modelResL$totalCatchM
  totalEffort = modelResL$actualEffort
  
  numReplic = nrow(totalAbunM)
  # create mean/sd/95% rows for each
  myP("Summarizing replicate data...")
  myP("Plotting results...")
  fileStr = "totalAbun"
  if (prefix!="") fileStr = combine(prefix, fileStr)
  titleStr = paste(prefix, paste("Total abundance, Rep:", numReplic))
  plotAndSaveMatrix(totalAbunM, title=titleStr, xLab="Time", yLab="Count", fileStr=fileStr, 
                    showPlot=showPlots, savePlotToFile=saveToFile)
  
  fileStr = "totalAbunSSB"
  if (prefix!="") fileStr = combine(prefix, fileStr)
  titleStr = paste(prefix, paste("Total SSB, Rep:", numReplic))
  plotAndSaveMatrix(totalAbunSSBM, title=titleStr, xLab="Time", yLab="Count", fileStr=fileStr, 
                    showPlot=showPlots, savePlotToFile=saveToFile)
  
  fileStr = "totalCatchMass";
  if (prefix!="") fileStr = combine(prefix, fileStr)
  titleStr = paste(prefix, paste("Total catch, Rep:", numReplic))
  plotAndSaveMatrix(totalCatchM, title=titleStr, xLab="Time", yLab="Kg", fileStr=fileStr, 
                    showPlot=showPlots, savePlotToFile=saveToFile)
  
  spatialMassA    = modelResL$spatialMassA
  spatialSSBA     = modelResL$spatialSSBA
  spatialCatchA   = modelResL$spatialCatchA
  spatialDensA    = modelResL$spatialDensA
  spatialSSBDensA = modelResL$spatialSSBDensA
  
  # average across replicates
  # row,col time, block
  spatialMassM    = apply(spatialMassA, c(2,3), mean)
  spatialSSBM     = apply(spatialSSBA, c(2,3), mean)
  spatialCatchM   = apply(spatialCatchA, c(2,3), mean)
  spatialDensM    = apply(spatialDensA, c(2,3), mean)
  spatialSSBDensM = apply(spatialSSBDensA, c(2,3), mean)
  spatialEffortM = apply(totalEffort, c(2,3), mean)
  
  plotHeatMap(spatialDensM/10000)
  plotHeatMap(spatialSSBDensM/10000)
  plotHeatMap(spatialCatchM)
  
  # Mean of last 10 time steps
  spatialMassE = apply(spatialMassM[(length(spatialMassM[,1])-10):length(spatialMassM[,1]),], 2, mean)
  spatialSSBE = apply(spatialSSBM[(length(spatialSSBM[,1])-10):length(spatialSSBM[,1]),], 2, mean)
  spatialCatchE = apply(spatialCatchM[(length(spatialCatchM[,1])-10):length(spatialCatchM[,1]),], 2, mean)
  spatialSSBDensE = apply(spatialSSBDensM[(length(spatialSSBDensM[,1])-10):length(spatialSSBDensM[,1]),], 2, mean)
  spatialEffortE = apply(spatialEffortM[(length(spatialEffortM[,1])-10):length(spatialEffortM[,1]),], 2, mean)
  
  plotSpatial(spatialMassE, title = "Mass") 
  plotSpatial(spatialSSBE, title = "SSB")
  plotSpatial(spatialCatchE, title = "Catch")
  plotSpatial(spatialSSBDensE, title = "SSBDens")
  plotSpatial(spatialEffortE, title = "Effort")
  
  # leave now if we're not saving a movie
  if (!createMovie) return()
  
  # create movies of spatial mass patterns
  # reform matrices as melted data
  massD  = melt(spatialMassM, varnames=c("t", "block"), value.name="mass")
  ssbD   = melt(spatialSSBM, varnames=c("t", "block"), value.name="mass")
  catchD = melt(spatialCatchM, varnames=c("t", "block"), value.name="mass")
  massD$Type  = "All"
  ssbD$Type   = "SSB"
  catchD$Type = "Catch"
  
  # merge together
  allMassD = rbind(massD, ssbD, catchD)
  
  # keep density separate
  # (keeping name as "mass" to simplify shared movie code
  ssbDensD = melt(spatialSSBDensM, varnames=c("t", "block"), value.name="mass")
  ssbDensD$Type = "SSB"
  # convert from /ha to /m2
  ssbDensD$mass = ssbDensD$mass / 10000
  
  # assemble total
  totalD       = data.frame(t=1:nrow(spatialMassM))
  totalD$mass  = apply(spatialMassM, 1, sum)
  totalD$ssb   = apply(spatialSSBM, 1, sum)
  totalD$catch = apply(spatialCatchM, 1, sum)
  
  fileStr = "spatialMass"
  if (prefix!="") fileStr = combine(prefix, fileStr)
  createSpatialMovie(allMassD, totalD, fileStr=fileStr, prefix=prefix)
  
  fileStr = "spatialDens"
  if (prefix!="") fileStr = combine(prefix, fileStr)
  createSpatialMovie(ssbDensD, totalD, yLab="#/m2", fileStr=fileStr, prefix=prefix)
  
}

# Summarizes the result matrix, plots summary, and saves both matrices to file.
plotAndSaveMatrix <- function(resM, title="", xV=NULL, xLab="", yLab="", fileStr, 
                              sepExtinct=TRUE, showPlot=TRUE, savePlotToFile=FALSE) {
  if (is.null(xV)) xV = 1:ncol(resM)
  resD       = data.frame(t=xV)
  resD$mean  = apply(resM, 2, mean, na.rm=TRUE)
  resD$sd    = apply(resM, 2, sd, na.rm=TRUE)
  resD$q_025 = apply(resM, 2, quantile, 0.025, na.rm=TRUE)
  resD$q_975 = apply(resM, 2, quantile, 0.975, na.rm=TRUE)
  
  if (sepExtinct) {
    extinctV   = resM[,ncol(resM)] < 10
    numExtinct = sum(extinctV)
    if (numExtinct<nrow(resM)) {
      # survivors
      if (numExtinct==(nrow(resM)-1)) {
        resD$surv_mean  = resM[!extinctV,]
        resD$surv_sd    = NA
        resD$surv_q_025 = NA
        resD$surv_q_975 = NA
      } else { 
        resD$surv_mean  = apply(resM[!extinctV,], 2, mean, na.rm=TRUE)
        resD$surv_sd    = apply(resM[!extinctV,], 2, sd, na.rm=TRUE)
        resD$surv_q_025 = apply(resM[!extinctV,], 2, quantile, 0.025, na.rm=TRUE)
        resD$surv_q_975 = apply(resM[!extinctV,], 2, quantile, 0.975, na.rm=TRUE)
      }
    }
    
    if (numExtinct>0) {
      # extinct
      if (numExtinct==1) {
        resD$ext_mean  = resM[extinctV,]
        resD$ext_sd    = NA
        resD$ext_q_025 = NA
        resD$ext_q_975 = NA
      } else { 
        resD$ext_mean  = apply(resM[extinctV,], 2, mean, na.rm=TRUE)
        resD$ext_sd    = apply(resM[extinctV,], 2, sd, na.rm=TRUE)
        resD$ext_q_025 = apply(resM[extinctV,], 2, quantile, 0.025, na.rm=TRUE)
        resD$ext_q_975 = apply(resM[extinctV,], 2, quantile, 0.975, na.rm=TRUE)
      } 
    }
  }
  
  # draw total abundance plot
  thePlot = ggplot(resD, aes(x=t)) + 
    geom_line(aes(y=mean)) +
    geom_line(aes(y=q_025), linetype="dashed") +
    geom_line(aes(y=q_975), linetype="dashed")
  thePlot = thePlot + labs(title=title, x=xLab, y=yLab)
  if (savePlotToFile) {
    myP("Saving plot to file...")
    fileName = paste("output/", paste(fileStr, ".png", sep=""), sep="")
    myPTiff(thePlot, fileName, h=5, w=7)
  } else if (showPlot) myPPlot(thePlot)
  
  # save the raw results and the summaries  
  saveMatrixToCSVFile(fileStr, resM)
  fileName = combine(fileStr, "summary")
  saveMatrixToCSVFile(fileName, resD)
}

plotSpatial <- function(spatialRes, title=titleStr, xLab = "", yLab = "", fileStr=fileStr, 
                        showPlot=showPlots, savePlotToFile=saveToFile) {
  
  spatialDF <- data.frame(patch = 1:length(spatialRes), spatialdata = spatialRes)
  thePlot = ggplot(spatialDF, aes(patch, spatialdata)) + 
    geom_line(aes(y=spatialdata)) +
    labs(title=title, x=xLab, y=yLab)
  
  myPPlot(thePlot) 
  
}

plotHeatMap <- function(spatialData, colorName = 'variable'){
  
  colnames(spatialData) = 1:length(spatialData[1,])
  
  spatialData %>%
    
    as_tibble() %>%
    rowid_to_column(var="X") %>%
    gather(key="Y", value="Z", -1) %>%
    
    # Change Y to numeric
    mutate(Y=as.numeric(gsub("V","",Y))) %>%
    
    # Viz
    ggplot(aes(X, Y, fill= Z)) + 
    geom_tile() +
    labs(title="", x="Time (years)", y="Coast", fill = colorName) + 
    theme(legend.position="right") -> thePlot
  
  myPPlot(thePlot)
  
}

# Summarizes the result matrix, plots summary, and saves both matrices to file.
createSpatialMovie <- function(allMassD, totalD, xLab="Block", yLab="Kg", fileStr, prefix="", delFrames=TRUE) {
  myP("Making movie...");
  maxVal     = max(allMassD$mass);
  maxT       = max(allMassD$t);
  maxTot     = max(totalD$mass);
  
  # function to draw each frame and save as PNG
  drawTimeStep <- function(ts) {
    if (ts%%10==0) myP("Drawing frame t=", ts, "...");
    tsStr = as.character(ts);
    if (ts<10) tsStr = paste("0", tsStr, sep="");
    if (ts<100) tsStr = paste("0", tsStr, sep="");
    fileName = paste(combine(prefix, combine("timestep", tsStr)), ".png", sep="");
    png(file=fileName, width=1200, height=800);
    tsD    = allMassD[allMassD$t==ts,];
    tsPlot = ggplot(data=tsD, aes(x=block, y=mass));
    tsPlot = tsPlot + geom_line(aes(group=Type, color=Type), size=2);
    tsPlot = tsPlot + scale_y_continuous(limits=c(0, maxVal)) +
      scale_color_manual(values=c("All"="blue", "SSB"="green", "Catch"="red"), breaks=c("All", "SSB", "Catch"));
    tsPlot = tsPlot + labs(title=combine(prefix, paste("t=",ts,sep="")), x=xLab, y=yLab) +
      theme(plot.title=element_text(size=rel(5)),
            legend.text=element_text(size=rel(4)), 
            legend.title=element_text(size=rel(4)),
            axis.text.x=element_text(size=rel(5)), 
            axis.text.y=element_text(size=rel(5)),
            axis.title.x=element_text(size=rel(5)),
            axis.title.y=element_text(size=rel(5)));
    # draw sub plot
    print(tsPlot);
    vp = viewport(width=0.15, height=0.15, x=0.9, y=0.75);
    # create inset temp plot with cats
    totTSD = totalD[totalD$t<=ts,];
    totPlot = ggplot(data=totTSD, aes(x=t)) + 
      geom_line(aes(y=mass), size=1.5, color="blue") + # mass
      geom_line(aes(y=ssb), size=1.5, color="green") + # ssb
      geom_line(aes(y=catch), size=1.5, color="red")  # catch
    totPlot = totPlot + labs(title="Total", x="time", y="kg") +
      scale_x_continuous(limits=c(0, maxT)) +
      scale_y_continuous(limits=c(0, maxTot)) +
      theme(axis.title.x=element_text(size=rel(1)),
            axis.title.y=element_text(size=rel(1)));
    
    print(totPlot, vp=vp);
    dev.off();
  }
  
  # loop function for use in saveGIF()
  loopTimeSteps <- function() {
    lapply(1:maxT, function(t) { drawTimeStep(t); });
  }
  
  # saveGIF PATH info not working quite right, do by hand
  loopTimeSteps();
  # hardcode PATH info
  # NOTE: incorporates ALL PNG files
  #       Must delete earlier files before running again
  system(paste('"C:\\Program Files\\ImageMagick-6.9.1-Q16\\convert.exe" -delay 15 *.png', 
               paste("C:\\Research\\Code\\abalone\\output\\", paste(fileStr, ".gif", sep=""), sep="")));
  if (delFrames) file.remove(list.files(pattern=".png"));
}

createSizeMovie <- function(allSizeD, totalD, catV, tV=NULL, xLab="Size", yLab="Density (#/m2)", fileStr, prefix="", delFrames=TRUE) {
  myP("Making movie...");
  maxVal = max(allSizeD$density);
  maxT   = max(allSizeD$t);
  if (is.null(tV)) tV = c(1, maxT);
  maxTot = max(totalD$resDens, totalD$fishDens);
  # quick loop to drop all 1 points and set others=0
  for (i in 1:maxT) catV[i] = ifelse(catV[i]<1, 0, -100);
  catD = data.frame(t=1:maxT, cat=catV);
  
  # function to draw each frame and save as PNG
  drawTimeStep <- function(ts) {
    if (ts%%10==0) myP("Drawing frame t=", ts, "...");
    tsStr = as.character(ts);
    if (ts<10) tsStr = paste("0", tsStr, sep="");
    if (ts<100) tsStr = paste("0", tsStr, sep="");
    fileName = paste(combine(prefix, combine("timestep", tsStr)), ".png", sep="");
    png(file=fileName, width=1200, height=800);
    tsD    = allSizeD[allSizeD$t==ts,];
    tsPlot = ggplot(data=tsD, aes(x=size, y=density));
    tsPlot = tsPlot + geom_line(aes(group=Type, color=Type), size=2);
    tsPlot = tsPlot + scale_y_continuous(limits=c(0, maxVal)) +
      scale_color_manual(values=c("Reserve"="blue", "Fishery"="red"), breaks=c("Reserve", "Fishery"));
    tsPlot = tsPlot + labs(title=combine(prefix, paste("t=",ts,sep="")), x=xLab, y=yLab) +
      theme(plot.title=element_text(size=rel(5)),
            legend.text=element_text(size=rel(4)), 
            legend.title=element_text(size=rel(4)),
            axis.text.x=element_text(size=rel(5)), 
            axis.text.y=element_text(size=rel(5)),
            axis.title.x=element_text(size=rel(5)),
            axis.title.y=element_text(size=rel(5)));
    # draw sub plot
    print(tsPlot);
    vp = viewport(width=0.15, height=0.15, x=0.9, y=0.75);
    # create inset temp plot with cats
    totTSD = totalD[totalD$t<=ts,];
    totPlot = ggplot(data=totTSD, aes(x=t)) + 
      geom_line(aes(y=resDens), size=1.5, color="blue") + # mass
      geom_line(aes(y=fishDens), size=1.5, color="red") + # catch
      geom_point(data=catD, aes(x=t, y=cat), size=3, color="yellow"); # catastrophes
    totPlot = totPlot + labs(title="Total", x="time", y="kg") +
      scale_x_continuous(limits=c(0, maxT)) +
      scale_y_continuous(limits=c(0, maxTot)) +
      theme(axis.title.x=element_text(size=rel(1)),
            axis.title.y=element_text(size=rel(1)));
    
    print(totPlot, vp=vp);
    dev.off();
  }
  
  # loop function for use in saveGIF()
  loopTimeSteps <- function(startT=1, endT=maxT) {
    lapply(startT:endT, function(t) { drawTimeStep(t); });
  }
  
  # saveGIF PATH info not working quite right, do by hand
  loopTimeSteps(startT=tV[1], endT=tV[2]);
  # hardcode PATH info
  # NOTE: incorporates ALL PNG files
  #       Must delete earlier files before running again
  system(paste('"C:\\Program Files\\ImageMagick-6.9.1-Q16\\convert.exe" -delay 30 *.png', 
               paste("C:\\Research\\Code\\abalone\\output\\", paste(fileStr, ".gif", sep=""), sep="")));
  if (delFrames) file.remove(list.files(pattern=".png"));
}

# returns a data frame with the mean time series
getMeanResults <- function(modelResL, isAllee, prefix) {
  # Dimensions: row=replicate, col=time
  totalAbunM    = modelResL$totalAbunM
  totalAbunSSBM = modelResL$totalAbunSSBM
  totalCatchM   = modelResL$totalCatchM
  
  
  maxT          = ncol(totalAbunM)
  meanResD      = data.frame(t=1:maxT)
  
  # Calculate mean over replicates
  meanResD$meanAbun  = apply(totalAbunM, 2, mean, na.rm=TRUE)
  meanResD$meanSSB   = apply(totalAbunSSBM, 2, mean, na.rm=TRUE)
  meanResD$meanCatch = apply(totalCatchM, 2, mean, na.rm=TRUE)
  meanResD$cumCatch  = 0
  
  # add cumulative catch
  for (cT in (maxT/2):maxT) {  
    lastCumCatch = meanResD[meanResD$t==(cT-1),]$cumCatch 
    meanResD[meanResD$t==cT,]$cumCatch = lastCumCatch + meanResD[meanResD$t==cT,]$meanCatch
  }
  
  # function to find mean size
  meanSizeF <- function(xV, ssbOnly=FALSE) {
    if (ssbOnly) xV = xV*gPropMatureV
    
    totalSizeV = xV*gMeanSizeV
    return(sum(totalSizeV)/sum(xV))
  }
  # function to find fecundity
  fecF <- function(xV, total=FALSE) {
    eggsV = xV*gEggsPerIndV
    
    if (total) return(sum(eggsV))
    else return(sum(eggsV)/sum(xV))
  }
  
  # determine which have collapsed
  extinctV = totalCatchM[,maxT] < (COLLAPSE_THRESHOLD*totalCatchM[,(maxT/2)-1])
  
  # calculate mean size, mean SSB size, mean fecundity, total fecundity
  # array(c(numReplic,NUM_CLASSES,maxT,NUM_BLOCKS)) # complete size/space data
  repSizeDataA         = apply(modelResL$compSizeDataA, c(1,2,3), sum, na.rm=TRUE) # sum across all blocks
  
  # replicate-specific values
  # for each rep, metrics at equil: abun, abun SSB, catch, size, size SSB, mean fec, total fec
  repEquilD = data.frame(rep=1:nrow(totalCatchM))
  AVE_OVER  = 10
  aveV      = (maxT-AVE_OVER+1):maxT
  repEquilD$collapse = extinctV
  repEquilD$abun     = apply(totalAbunM[,aveV], 1, mean, na.rm=TRUE)   
  repEquilD$abunSSB  = apply(totalAbunSSBM[,aveV], 1, mean, na.rm=TRUE)   
  repEquilD$catch    = apply(totalCatchM[,aveV], 1, mean, na.rm=TRUE)
  equilRepSizeM      = apply(repSizeDataA[,,aveV], c(1,2), mean, na.rm=TRUE)
  repEquilD$size     = apply(equilRepSizeM, 1, meanSizeF)  
  repEquilD$sizeSSB  = apply(equilRepSizeM, 1, meanSizeF, ssbOnly=TRUE)  
  repEquilD$meanFec  = apply(equilRepSizeM, 1, fecF)  
  repEquilD$totalFec = apply(equilRepSizeM, 1, fecF, total=TRUE)  
  fileName = combine(prefix, "repData")
  saveMatrixToCSVFile(fileName, repEquilD)
  
  # metrics across replicates
  # note that the averages are slightly different for size and fec. than when averaging the by-replicate values
  # average of mean sizes is different from mean size of average
  if (sum(extinctV)>0) {
    # average across extinct replicates
    if (sum(extinctV)==1) meanSizeDataExtM = repSizeDataA[extinctV,,] # don't average
    else meanSizeDataExtM = apply(repSizeDataA[extinctV,,], c(2,3), mean, na.rm=TRUE) 
    
    meanResD$meanSizeExt    = apply(meanSizeDataExtM, 2, meanSizeF)
    meanResD$meanSizeSSBExt = apply(meanSizeDataExtM, 2, meanSizeF, ssbOnly=TRUE)
    meanResD$meanFecExt     = apply(meanSizeDataExtM, 2, fecF)
    meanResD$totalFecExt    = apply(meanSizeDataExtM, 2, fecF, total=TRUE)
  } else {
    meanResD$meanSizeExt    = rep(-1, maxT)
    meanResD$meanSizeSSBExt = rep(-1, maxT)
    meanResD$meanFecExt     = rep(-1, maxT)
    meanResD$totalFecExt    = rep(-1, maxT)
  }
  if (sum(!extinctV)>0) {
    # average across non-extinct replicates
    if (sum(!extinctV)==1) meanSizeDataSurvM = repSizeDataA[!extinctV,,] # don't average
    else meanSizeDataSurvM = apply(repSizeDataA[!extinctV,,], c(2,3), mean, na.rm=TRUE) 
    
    meanResD$meanSizeSurv    = apply(meanSizeDataSurvM, 2, meanSizeF)
    meanResD$meanSizeSSBSurv = apply(meanSizeDataSurvM, 2, meanSizeF, ssbOnly=TRUE)
    meanResD$meanFecSurv     = apply(meanSizeDataSurvM, 2, fecF)
    meanResD$totalFecSurv    = apply(meanSizeDataSurvM, 2, fecF, total=TRUE)
  } else {
    meanResD$meanSizeSurv    = rep(-1, maxT)
    meanResD$meanSizeSSBSurv = rep(-1, maxT)
    meanResD$meanFecSurv     = rep(-1, maxT)
    meanResD$totalFecSurv    = rep(-1, maxT)
  }
  
  meanResD$name      = prefix
  meanResD$isAllee   = isAllee
  meanResD$numExtinct = sum(extinctV)
  
  return(meanResD)
}

# adds mean spatial density in reserve and outside
# creates separate means for time series which go extinct
calcSpatialDensity <- function(theD, resL, paramL, sepExtinct=TRUE) {
  numReplic       = paramL$numReplic
  # get location of reserves
  blockEffortV    = getBlockEffortVector(paramL$mpaSize, paramL$mpaProp, paramL$hProp)
  fishV           = blockEffortV>0
  resV            = !fishV
  spatialSSBDensA = resL$spatialSSBDensA
  totalCatchM     = resL$totalCatchM
  # sum across relevant blocks and replicates
  if (numReplic>1) sumCol = 2  
  else sumCol = 1 # deal with dimension loss for numReplic=1
  
  # Sum for reserves and fishable blocks
  theD$resDens  = apply(spatialSSBDensA[,,resV], sumCol, sum) / (sum(resV)*numReplic*10000)
  theD$fishDens = apply(spatialSSBDensA[,,fishV], sumCol, sum) / (sum(fishV)*numReplic*10000)
  
  maxT         = ncol(spatialSSBDensA)
  # Collapse criterion: <90% pre-catastrophe catch (Worm et al. 2006)
  # XXX assumes 1st catastrophe at maxT/2
  extinctV = totalCatchM[,maxT] < (COLLAPSE_THRESHOLD*totalCatchM[,(maxT/2)-1])
  
  # Abundance
  #  replicFinalV = apply(spatialSSBDensA[,maxT,], 1, sum);
  #  initialV     = apply(spatialSSBDensA[,(maxT/2)-1,], 1, sum);
  #  extinctV     = replicFinalV<(COLLAPSE_THRESHOLD*initialV);
  propExtinct  = sum(extinctV) / numReplic
  myP("Proportion extinct:", propExtinct)
  
  if (sepExtinct) {
    # separate out extinct time series from non-
    if (propExtinct>0) {
      if (sum(extinctV)==1) sumCol = 1
      else sumCol = 2
      theD$resDensExt      = apply(spatialSSBDensA[extinctV,,resV], sumCol, sum) / (sum(resV)*numReplic*10000*propExtinct)
      theD$fishDensExt     = apply(spatialSSBDensA[extinctV,,fishV], sumCol, sum) / (sum(fishV)*numReplic*10000*propExtinct)
      theD$resDensExtSD    = apply(spatialSSBDensA[extinctV,,resV], sumCol, sd) / (10000)
      theD$fishDensExtSD   = apply(spatialSSBDensA[extinctV,,fishV], sumCol, sd) / (10000)
      theD$resDensExtQ025  = apply(spatialSSBDensA[extinctV,,resV], sumCol, quantile, 0.025) / (10000)
      theD$fishDensExtQ025 = apply(spatialSSBDensA[extinctV,,fishV], sumCol, quantile, 0.025) / (10000)
      theD$resDensExtQ975  = apply(spatialSSBDensA[extinctV,,resV], sumCol, quantile, 0.975) / (10000)
      theD$fishDensExtQ975 = apply(spatialSSBDensA[extinctV,,fishV], sumCol, quantile, 0.975) / (10000)
      if (sumCol==1) theD$catchExt = totalCatchM[extinctV,]
      else theD$catchExt = apply(totalCatchM[extinctV,], sumCol, sum) / (numReplic*propExtinct)
    } else {
      myP("No extinct replicates.")
      theD$resDensExt      = rep(-1, maxT)
      theD$fishDensExt     = rep(-1, maxT)
      theD$resDensExtSD    = rep(-1, maxT)
      theD$fishDensExtSD   = rep(-1, maxT)
      theD$resDensExtQ025  = rep(-1, maxT)
      theD$fishDensExtQ025 = rep(-1, maxT)
      theD$resDensExtQ975  = rep(-1, maxT)
      theD$fishDensExtQ975 = rep(-1, maxT)
      theD$catchExt        = rep(-1, maxT)
    }
    if (propExtinct<1) {
      if (sum(!extinctV)==1) sumCol = 1
      else sumCol = 2;
      theD$resDensSurv      = apply(spatialSSBDensA[!extinctV,,resV], sumCol, sum) / (sum(resV)*numReplic*10000*(1-propExtinct))
      theD$fishDensSurv     = apply(spatialSSBDensA[!extinctV,,fishV], sumCol, sum) / (sum(fishV)*numReplic*10000*(1-propExtinct))
      theD$resDensSurvSD    = apply(spatialSSBDensA[!extinctV,,resV], sumCol, sd) / (10000)
      theD$fishDensSurvSD   = apply(spatialSSBDensA[!extinctV,,fishV], sumCol, sd) / (10000)
      theD$resDensSurvQ025  = apply(spatialSSBDensA[!extinctV,,resV], sumCol, quantile, 0.025) / (10000)
      theD$fishDensSurvQ025 = apply(spatialSSBDensA[!extinctV,,fishV], sumCol, quantile, 0.025) / (10000)
      theD$resDensSurvQ975  = apply(spatialSSBDensA[!extinctV,,resV], sumCol, quantile, 0.975) / (10000)
      theD$fishDensSurvQ975 = apply(spatialSSBDensA[!extinctV,,fishV], sumCol, quantile, 0.975) / (10000)
      if (sumCol==1) theD$catchSurv = totalCatchM[!extinctV,]
      else theD$catchSurv = apply(totalCatchM[!extinctV,], sumCol, sum) / (numReplic*(1-propExtinct))
    } else {
      myP("No non-extinct replicates.")
      theD$resDensSurv      = rep(-1, maxT)
      theD$fishDensSurv     = rep(-1, maxT)
      theD$resDensSurvSD    = rep(-1, maxT)
      theD$fishDensSurvSD   = rep(-1, maxT)
      theD$resDensSurvQ025  = rep(-1, maxT)
      theD$fishDensSurvQ025 = rep(-1, maxT)
      theD$resDensSurvQ975  = rep(-1, maxT)
      theD$fishDensSurvQ975 = rep(-1, maxT)
      theD$catchSurv        = rep(-1, maxT)
    }
    
    cumCatchSurv = rep(0, (maxT/2)-1)
    cumCatchExt  = rep(0, (maxT/2)-1)
    for (t in (maxT/2):maxT) {
      cumCatchSurv = c(cumCatchSurv, cumCatchSurv[t-1] + theD$catchSurv[t])
      cumCatchExt  = c(cumCatchExt, cumCatchExt[t-1] + theD$catchExt[t])
    }
    theD$cumCatchSurv = cumCatchSurv
    theD$cumCatchExt  = cumCatchExt
  }
  
  # add extinction prop to data for ease of access
  theD$propExtinct = propExtinct
  return(list(densD=theD, propExtinct=propExtinct))
}

#============================================================================# 
#  KERNEL FUNCTIONS
#============================================================================# 
ab_survX <- function(x) { 
  # s = e^(a-b*log(weight))
  SURV_A = 0.6346126
  SURV_B = 0.317388
  mort   = exp(SURV_A-SURV_B*log(wForL(x))) ### ALL NAT DATA
  return(exp(-mort))
} 
ab_harvX <- function(x, hProp, minSize=MIN_LAND_SIZE, maxSize=9999) { 
  # XXX ignoring maxsize for now
  ifelse((x>=minSize), 1-hProp, 1); 
} 
#ab_growYX <- function(y,x,meanC=1, meanXC=0.7, sd=0.3) { dnorm(y, mean=meanC+meanXC*x, sd=sd); }
#ab_fecYX <- function(y,x, fec=0.75, meanC=0.8, meanXC=0.2, sd=0.35) { fec*dnorm(y, mean=meanC+meanXC*x, sd=sd); }
#ab_kernelYX <- function(y,x) { ab_survX(x)*(ab_growYX(y,x)+ab_fecYX(y,x)); }

getHarvestMatrix <- function(hProp, hSDProp=0, minSize=MIN_LAND_SIZE, maxSize=9999, verbose=TRUE) {
  if (verbose) {
    myP("Creating harvest matrix for h =", hProp);
    if (hSDProp>0) myP("Harvest is uncertain with SD proportion =", hSDProp);
  }
  
  # calculate uncertain hProp
  # no error if h=0
  if ((hProp==0) || (hSDProp==0)) hVal = hProp
  else {
    # draw from rand, using SD value as a proportion of the mean
    # draw from here so all size classes have same random harvest
    hVal = rnorm(1, mean=hProp, sd=hSDProp*hProp);
    if (hVal<0) hVal = 0
    else if (hVal>1) hVal = 1;
  }
  
  H = matrix(0,NUM_CLASSES,NUM_CLASSES);  # Harvest matrix (survival from harvest)
  diag(H) <- ab_harvX(gMeanSizeV, hVal);
  return(H);
}

getEventSurvivalArray <- function(maxT, numReplic=NUM_REPLICATIONS) {
  eventA = array(1, dim=c(NUM_CLASSES, NUM_BLOCKS, maxT, numReplic));
  return(eventA);
}

getAlleeEffect <- function(densHa, aggA=AGG_A, alleeType=ALLEE_PROB) {
  # IMPORTANT! Density comes in as #/ha, not #/m2
  densM2 = densHa / 10000;
  
  # age_agg = a*d + b
  aveAggSize = aggA*densM2 + AGG_B
  
  if (alleeType==ALLEE_PROB) {
    # Probabilistic Allee:
    # P_mixed(x) = 1 - 0.5^(x-1)
    alleeMult = 1 - 0.5^(aveAggSize-1)
  } else if (alleeType==ALLEE_LIN) {
    # Linear Allee:
    # Assumes 0.2 = 80% 
    # A = 4*dens if dens<0.25, else =1-1
    if (densM2<0.25) alleeMult = 4*densM2
    else alleeMult = 1
  } else if (alleeType==ALLEE_EXP) {
    # Exponential Allee:
    # Assumes 0.2 = 80% 
    # A = exp(-(3.38-agg)*2) if agg<3.38
    if (aveAggSize<3.38) alleeMult = exp(-2*(3.38-aveAggSize))
    else alleeMult = 1
  } else alleeMult = 1  # no Allee effect
  return(alleeMult);
}

# Calculates the effective density of each block utilizing the Allee effect
# Adjusts for possible OA effect on fertilization success
getEffectiveDensityMatrix <- function(localDensM, aggA=AGG_A, oaFertMult=1, oaDensityMult=1) {
  # assumes row,col = class,block
  for (i in 1:ncol(localDensM)) {
    # calculate total density of mature individuals
    # apply OA multiplier
    totalMature = sum(gPropMatureV*localDensM[,i])*oaDensityMult;
    # apply Allee effect, if any
    # IMPORTANT! Density here is in #/ha, not #/m2
    # will convert inside
    localDensM[,i] = localDensM[,i]*getAlleeEffect(totalMature, aggA)*oaFertMult;
  }
  return(localDensM);
}

# Creates a size distribution vector for year-1 recruits 
# Must include possible OA effect
getRecruitDistributionVector <- function(isConst=FALSE, oaMult=1) {
  gRecruitDistribV = rep(0, NUM_CLASSES);
  if (isConst) {
    recruitInd = 1;
    for (i in 2:NUM_CLASSES) {
      if (gMeanSizeV[i]>LARVAL_CONST_SIZE) break
      else recruitInd = i;
    }
    gRecruitDistribV[recruitInd] = 1;
  } else {
    # Using Bardos(2005) and starting with larvae all of 500um, determine share in each size class after 1 year
    lowerSizesV = gMeanSizeV - CLASS_WIDTH/2;
    # Piggy-backing on the transition matrix code here
    gRecruitDistribV = getTransitionMatrix(lowerSizesV, recruitVectorOnly=TRUE);
  }
  return(gRecruitDistribV);
}

# Making recruit vector
# Applies density dependence on incoming settlers
getRecruitV <- function(numSettlers, recK=DEF_K) { # totSettlers is the number of arrived settlers    
  
  REC_A = 0.01; 
  
  #  REC_B = 5*10^7;
  newRecruits = REC_A*numSettlers*(exp(-numSettlers/recK));
  
  # multiply by global constant distribution vector
  recruitV = newRecruits*gRecruitDistribV;
  return(recruitV);
} 

# Creates the per-block effort vector which "creates" the MPAs
getBlockEffortVector <- function(mpaSize, mpaProp, hProp) {
  # Determine block structure of reserve network
  numBlockProt = round(mpaProp*NUM_BLOCKS, 0);           # Total number of blocks protected
  numReserve   = round(numBlockProt/mpaSize, 0); # Total number of reserves
  numBlockFish = NUM_BLOCKS - numBlockProt;               # Total number of fished blocks
  
  x  = floor(numBlockFish/numReserve);
  y  = ceiling(numBlockFish/numReserve);
  nx = numBlockFish - x*numReserve;
  ny = numReserve - nx;
  
  # Reserves exist as h=0 blocks
  # generate spatial effort vector across blocks
  if (is.na(nx)|is.na(ny)) { 
    # No reserves
    blockEffortV = rep(hProp, NUM_BLOCKS); 
  } else if (nx==ny) { 
    blockEffortV = rep(c(rep(0,mpaSize),rep(hProp,x),rep(0,mpaSize),rep(hProp,y)), ny); 
  } else if (nx>ny) { 
    blockEffortV = c(rep(c(rep(0,mpaSize),rep(hProp,x),rep(0,mpaSize),rep(hProp,y)), ny),
                     rep(c(rep(0,mpaSize),rep(hProp,y)), (nx-ny))); 
  } else if (nx<ny) { 
    blockEffortV = c(rep(c(rep(0,mpaSize),rep(hProp,x),rep(0,mpaSize),rep(hProp,y)), nx), 
                     rep(c(rep(0,mpaSize),rep(hProp,x)), (ny-nx)));
  }
  
  return(blockEffortV);
}

# Random matrix of dispersal SD values across time and replicates
# Currently only supports gamma distribution
getRandomDispSDMatrix <- function(maxT, nReplic, rType="gamma", shape=DISP_GAMMA_SHAPE) {
  # generate random value vector
  if (rType=="gamma") randV = rgamma(maxT*nReplic, shape=shape, rate=DISP_GAMMA_RATE)
  else {
    myP("WARNING: Distribution", rType, "unsupported for getRandomDispSDMatrix.");
    randV = rep(0, maxT*nReplic);
  }
  
  # create and return matrix for time and replicates
  return(matrix(randV, nrow=maxT, ncol=nReplic));
}

# Random matrix of recruitment SD values across time and replicates
# Values vary from 
getRandomRecruitMatrix <- function(maxT, nReplic, rType="log-normal", shape=REC_LN_SDLOG) {
  # generate random value vector
  if (rType=="gamma") {
    myP("Creating random recruit matrix with gamma distribution...");
    randV = rgamma(maxT*nReplic, shape=shape, rate=REC_GAMMA_RATE);
  } else if (rType=="log-normal") {
    myP("Creating random recruit matrix with log-normal distribution...");
    randV   = rlnorm(maxT*nReplic, meanlog=REC_LN_MEANLOG, sdlog=shape);
  }	else if (rType=="good-bad") {
    myP("Creating random recruit matrix with good-bad years...");
    randV   = runif(maxT*nReplic);
    threshV = randV > REC_GB_THRESHOLD;
    randV[threshV]  = shape;
    randV[!threshV] = REC_GB_BAD;
  } else {
    myP("WARNING: Distribution", rType, "unsupported for getRandomRecruitMatrix.");
    randV = rep(0, maxT*nReplic);
  }
  
  # create and return matrix for time and replicates
  return(matrix(randV, nrow=maxT, ncol=nReplic));
}

# Creates transition matrix by integrating fitted gamma function across each size class
# Takes a size class vector of lower bounds, not midpoints
getTransitionMatrix <- function(sizesV, recruitVectorOnly=FALSE) {
  # growth parameters from "Fitting_Bardos" MLE fit
  # load("Fitting_Bardos"); estimate5_3;
  
  #XXX test cynthia's growth function (MEPS catton)
  B_GI    = 0.5635249;
  B_SIGMA = 55.9560109;
  B_LN    = 150.3854074;
  B_GAMMA = 1.7190492;
  B_BETA  = 1.4780344;
  
  m     = 3;
  n     = 3;
  
  # local functions to generate gamma values
  # calculate Linf (Eq. A2 in Appendix)
  LinfF <- function(dL, L1) { return(((L1+dL)*L1^(-exp(-B_GI)))^(1/(1-exp(-B_GI)))); }
  
  # nF and dF are numerator/denominator functions for lambdaF and roF
  numF <- function(L1) { return(B_LN/(1+(B_BETA*L1/B_LN)^m)); }
  denF <- function(L1) { return(B_SIGMA/(1+(B_GAMMA*L1/B_LN)^n)); }
  
  # gamma rate (Eq. A3 in Appendix)
  lambdaF <- function(L1) { return(numF(L1)/(denF(L1))^2); }
  
  # gamma shape (Eq. A4 in Appendix)
  roF <- function(L1) { return((numF(L1))^2/(denF(L1))^2); }
  
  # function to generate the matrix itself
  # NOTE: Needs sizes vector with one extra large size (i.e. beyond the limit)
  createTransMatrix <- function(sizesPlusV) {
    myP("Creating transition matrix...");
    numSizesM1 = length(sizesPlusV) - 1;
    P = matrix(0, nrow=numSizesM1, ncol=numSizesM1);
    
    # integrate k across size class
    for (i in 1:numSizesM1) {
      for (j in 1:numSizesM1) {
        # integrand function
        integrand <- function(L1) { 
          return(pgamma(LinfF(sizesPlusV[i+1]-L1,L1)-L1, rate=lambdaF(L1), shape=roF(L1)) - 
                   pgamma(LinfF(sizesPlusV[i]-L1,L1)-L1, rate=lambdaF(L1), shape=roF(L1)));
        }
        P[i, j] = round((integrate(integrand,lower=sizesPlusV[j],upper=sizesPlusV[j+1]))$value/(sizesPlusV[j+1]-sizesPlusV[j]), 4);
      }
    }
    P[numSizesM1, numSizesM1] = 1;
    
    # make sure each column adds to 1
    for (col in 1:ncol(P)) P[,col] = P[,col] / sum(P[,col]);
    
    return(P);
  }
  
  if (recruitVectorOnly) {
    # Assume recruits start at L0 and grow for one year
    L0 = 1;
    # Using Bardos(2005) and starting with larvae all of 500um, determine share in each size class after 1 year
    recruitV = rep(0, length(sizesV));
    # convert lower sizes vector to upper and add L0
    upperSizesPlusV = c(L0, sizesV + CLASS_WIDTH);
    # integrate growth across size class
    # note that we assume that anything below MIN_SIZE is bumped up to size class 1
    for (i in 1:length(sizesV)) {
      recruitV[i] = pgamma(LinfF(upperSizesPlusV[i+1]-L0, L0)-L0, rate=lambdaF(L0), shape=roF(L0)) - 
        pgamma(LinfF(upperSizesPlusV[i]-L0, L0)-L0, rate=lambdaF(L0), shape=roF(L0));
    }
    #	recruitV = recruitV / sum(recruitV);
    return(recruitV);
  } else {
    # create sizes vector with one extra large size
    sizesPlusV = c(sizesV, sizesV[length(sizesV)]+CLASS_WIDTH);
    # call and return
    return(createTransMatrix(sizesPlusV));
  }
}

# Creates transition matrix by integrating fitted gamma function across each size class
# Takes a size class vector of lower bounds, not midpoints
getSurvivalMatrix <- function(sizesV, pM, classicSurv=FALSE) {
  # This is more complicated now
  # Need to integrate across expected growth, or survival is way too low (due to low mean size)
  # Depends a lot on size-class
  
  if (classicSurv) {
    # generate simple length-based diagonal using class midpoint
    S = matrix(0, nrow=length(gMeanSizeV), ncol=length(gMeanSizeV));
    diag(S) <- ab_survX(gMeanSizeV);
    return(S);
  }
  
  # function to generate the matrix itself
  # NOTE: Needs sizes vector with one extra large size (i.e. beyond the limit)
  createSurvMatrix <- function(sizesPlusV) {
    myP("Creating survival matrix...");
    numSizesM1 = length(sizesPlusV) - 1;
    S = matrix(0, nrow=numSizesM1, ncol=numSizesM1);
    
    # integrate across size class
    for (i in 1:numSizesM1) {
      for (j in 1:numSizesM1) {
        if (i>=j) { # no negative growth
          # integrand function
          if (i==j) midPoint1 = sizesPlusV[j]
          else midPoint1 = sizesPlusV[j] + CLASS_WIDTH/2;
          integrand <- function(L1) { 
            intSum = 0;
            dL     = L1 - midPoint1;
            if (dL<0) return(log(ab_survX(midPoint1)));
            D_K    = 100;
            kdL    = dL/D_K;
            for (k in 1:D_K) intSum = intSum + log(ab_survX(midPoint1+k*kdL));
            return(intSum/D_K);
          }
          ans = integrate(integrand,lower=sizesPlusV[i],upper=sizesPlusV[i+1])$value/(sizesPlusV[i+1]-sizesPlusV[i]);
          #		  if (j==1) myP(ans, exp(ans));
          S[i, j] = exp(ans);
          #		  S[i,j] = exp((log(ab_survX(midPoint1)) + log(ab_survX(midPoint2)))/2);
        }
      }
    }
    
    return(S)
  }
  
  # create sizes vector with one extra large size
  sizesPlusV = c(sizesV, sizesV[length(sizesV)]+CLASS_WIDTH)
  # call and return
  return(createSurvMatrix(sizesPlusV))
}

initializeGlobalValues()

################################
#      PLOTTING RESULTS
################################
plotGrowthLikelihood <- function(pM, lV=NULL) {
  # pMatrix: row = Lt+1, col=Lt
  if (is.null(lV)) lV = gMeanSizeV;
  growthM = matrix(0, nrow=length(lV)^2, ncol=3);
  for (i in 1:nrow(pM)) {
    Lt_1 = lV[i];
    for (j in 1:ncol(pM)) {
      Lt = lV[j];
      dL = Lt_1 - Lt;
      prop = pM[i,j];
      prop = prop / max(pM[,j]);
      growthM[(i-1)*ncol(pM)+j,] = c(Lt, dL, prop);
    }
  }
  growthD = data.frame(l=growthM[,1], dL=growthM[,2], prop=growthM[,3]);
  growthD = growthD[growthD$prop>0.00001,];
  thePlot = ggplot(growthD, aes(x=l, y=dL)) + theme_bw() +
    geom_point(aes(color=prop), size=3) +
    labs(x="Length (mm)", y="Growth (mm)") + 
    theme(legend.text=element_text(size=rel(1.25)), 
          legend.title=element_text(size=0),
          legend.position=c(.90, .8)) +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)));
  
  myPPlot(thePlot);
}

plotCohort <- function(numInd, pM, sM, maxT=10, showHist=FALSE) {
  indM = matrix(0, nrow=numInd, ncol=maxT);
  indM[,1] = 8; # start at ave size class
  cumPM    = pM;
  totalRec = 0;
  for (i in 2:nrow(pM)) cumPM[i,] = cumPM[i,] + cumPM[i-1,];
  for (t in 2:maxT) {
    growRandV = runif(numInd);
    survRandV = runif(numInd);
    for (i in 1:numInd) {
      indL = indM[i,t-1];
      if (indL>0) { # alive
        if (survRandV[i]>sM[indL,indL]) indM[i,t] = 0 # dies
        else {
          # grows
          newIndL   = sum(cumPM[,indL]<growRandV[i])+1;
          if (t==maxT) totalRec  = totalRec + gEggsPerIndV[newIndL];
          #		  if (indL==1) myP(growRandV[i], indL, newIndL, gMeanSizeV[indL], gMeanSizeV[newIndL]);
          indM[i,t] = newIndL;
        }
      }
    }    
  }
  
  # count # of >0 values in each column
  survV = apply(apply(indM, 2, function(x) { ifelse(x>0,1,0);}), 2, sum) / numInd;
  survD = data.frame(t=1:maxT, surv=survV);
  
  indD  = melt(indM, c("ind", "t"), value.name="length");
  indD  = indD[indD$length>0,]; # drop all dead individuals
  # adjust from index
  indD$length = indD$length * CLASS_WIDTH + MIN_SIZE - CLASS_WIDTH/2;
  
  thePlot = ggplot(indD, aes(x=t, y=length)) + theme_bw() + # geom_line(aes(group=ind)) +
    geom_smooth(size=2, color="black") + labs(title="Growth", x="Time", y="Length") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)));
  if (showHist)	myPPlot(thePlot, plotDim=c(3,1))
  else myPPlot(thePlot, plotDim=c(2,1));
  thePlot = ggplot(survD, aes(x=t, y=surv)) + theme_bw() +
    geom_line(size=2, color="black") + labs(title="Survival", x="Time", y="Prop. alive") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)));
  
  myPPlot(thePlot, plotPos=c(2,1), newWindow=FALSE);
  if (showHist) {
    thePlot = ggplot(indD[indD$t==maxT,], aes(x=length)) + theme_bw() +
      geom_histogram() + labs(x="Length", y="Count");
    myPPlot(thePlot, plotPos=c(3,1), newWindow=FALSE);
  }
  myP("Total survivors:", survV[maxT]*numInd);
  myP("Total mean per-capita eggs output in final year:", myF(totalRec/numInd));  
  
  
}

plotRecruitmentOptions <- function(lnVal=REC_LN_SDLOG, gamVal=REC_GAMMA_SHAPE, gbVal=REC_GB_GOOD) {
  numVal = 100000;
  
  rand1V = getRandomRecruitMatrix(numVal, 1, "log-normal", lnVal);
  rand2V = getRandomRecruitMatrix(numVal, 1, "gamma", gamVal);
  rand3V = getRandomRecruitMatrix(numVal, 1, "good-bad", gbVal);
  
  # drop outliers
  rand1V = dropOutliers(rand1V, c(0.025, 0.975));
  rand2V = dropOutliers(rand2V, c(0.025, 0.975));
  #  rand3V = dropOutliers(rand3V, c(0.025, 0.975));
  
  myP("Ratio of largest to smallest:");
  myP("Log-normal:", max(rand1V)/min(rand1V));
  myP("Gamma:", max(rand2V)/min(rand2V));
  myP("Good-bad:", max(rand3V)/min(rand3V));
  
  theD    = data.frame(x=rand1V);
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram() + labs(title="Log-normal") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)));
  
  myPPlot(thePlot, plotDim=c(3,1), newWindow=TRUE);
  theD    = data.frame(x=rand2V);
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram() + labs(title="Gamma") 
  theme(axis.text.x=element_text(size=rel(1.25)), 
        axis.text.y=element_text(size=rel(1.25)),
        axis.title.x=element_text(size=rel(1.25)),
        axis.title.y=element_text(size=rel(1.25)));
  
  myPPlot(thePlot, plotPos=c(2,1), newWindow=FALSE);
  theD    = data.frame(x=rand3V);
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram() + labs(title="Good-bad") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)));
  myPPlot(thePlot, plotPos=c(3,1), newWindow=FALSE);
}

plotStochasticElements <- function(numVal=100000, lnVal=REC_LN_SDLOG, gamVal=DISP_GAMMA_SHAPE) {
  
  rand1V = getRandomRecruitMatrix(numVal, 1, "log-normal", lnVal)
  rand2V = getRandomDispSDMatrix(numVal, 1, "gamma", gamVal)
  
  # drop outliers
  rand1V = dropOutliers(rand1V, c(0.025, 0.975))
  rand2V = dropOutliers(rand2V, c(0.025, 0.975))
  
  myP("Ratio of largest to smallest:")
  myP("Log-normal:", max(rand1V)/min(rand1V))
  myP("Gamma:", max(rand2V)/min(rand2V))
  
  myMode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  myP("Mean/mode recruit:", mean(rand1V), myMode(rand1V))
  myP("Mean/mode disp:", mean(rand2V), myMode(rand2V))
  
  theD    = data.frame(x=rand1V)
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram(binwidth=(max(rand1V)-min(rand1V))/100) + theme_bw() +
    labs(x="Recruitment strength", y="Count") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)))
  
  myPPlot(thePlot)
  theD    = data.frame(x=rand2V)
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram(binwidth=(max(rand2V)-min(rand2V))/100)  + theme_bw() +
    labs(x="Mean dispersal distance (m)", y="Count") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)))
  
  myPPlot(thePlot)
}

plotConstVectors <- function() {
  theD    = data.frame(x=gMeanSizeV, eggs=gEggsPerIndV, pMat=gPropMatureV, rec=getRecruitDistributionVector())
  thePlot = ggplot(theD, aes(x=x)) + geom_line(aes(y=eggs)) +
    labs(x="Size", y="Eggs") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)))
  
  myPPlot(thePlot, plotDim=c(2,2))
  thePlot = ggplot(theD, aes(x=x)) + geom_line(aes(y=pMat)) +
    labs(x="Size", y="Prop. Maturity") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)))
  myPPlot(thePlot, newWindow=FALSE, plotPos=c(1,2))
  
  thePlot = ggplot(theD, aes(x=x)) + geom_line(aes(y=rec)) +
    labs(x="Size", y="Recruits") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)))
  myPPlot(thePlot, newWindow=FALSE, plotPos=c(2,1))
  
  densityV = seq(0,1,by=0.1)
  theD = data.frame(x=densityV, allee=getAlleeEffect(densityV*10000))
  thePlot = ggplot(theD, aes(x=x)) + geom_line(aes(y=allee)) +
    labs(x="Density", y="Allee effect") +
    theme(axis.text.x=element_text(size=rel(1.25)), 
          axis.text.y=element_text(size=rel(1.25)),
          axis.title.x=element_text(size=rel(1.25)),
          axis.title.y=element_text(size=rel(1.25)))
  myPPlot(thePlot, newWindow=FALSE, plotPos=c(2,2))
}
