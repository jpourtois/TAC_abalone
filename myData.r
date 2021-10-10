# My data functions
# Script from Aalto et al. 2019 ('Catastrophic Mortality, Allee Effects, and Marine Protected Areas), written by Emilius A. Aalto

# NOTE: Many of these functions make changes *in place*, without creating a new column
#library(plyr);
library(nlme);
#source("../common/borrowed.r");

# DATA MANIPULATION FUNCTIONS

# <myPrediction>
# Creates predicted values for data from model
# Adds "predicted" and "se" columns
myPrediction <- function(dataIn, theModel) {
	theData = as.data.frame(dataIn);
    pred    = predict(theModel, theData, se=T);
	theData$predicted = pred;
#	theData$se        = pred$se.fit
	return(theData);
}

# <myTrim>
# Trims data column in place, applying outlier calculations
# Options: trim vs. winsorize (replace with threshold value)
#          use outlier labeling calculations vs. 3 s.d. trimming
myTrim <- function(dataIn, colName, colOut="", doWinsor=FALSE, doLabeling=TRUE, gVal=2.2, removeNA=TRUE, opl=2) {
  # XXX currently ignores 'colOut'
#  if (colOut=="") {
#    colOut = colName;
  myP2(opl, "Trimming column '", colName, "' in place...");
#  } else myP2(opl, "Trimming column '", colName, "' -> '", colOut, "'...");
  trimmedData = as.data.frame(dataIn);
  namesV      = names(trimmedData)
  if (removeNA) {
    curCnt      = nrow(trimmedData);
	# added to protect against single column data errors
    trimmedData = as.data.frame(trimmedData[!is.na(trimmedData[[colName]]),]);
	names(trimmedData) <- namesV
	trimCnt     = nrow(trimmedData);
	if (curCnt>trimCnt) {
	  myP3(opl, "Removing", (curCnt-trimCnt), "NA rows...");
	} else myP3(opl, "No NA rows removed.");
  }
  if (doWinsor) trimStr = "Winsoring via"
  else trimStr = "Trimming via";
  if (doLabeling) {
    # use the +/- g*(Q3-Q1) labeling rule
    q25 = quantile(trimmedData[[colName]], 0.25, na.rm=TRUE);
    q75 = quantile(trimmedData[[colName]], 0.75, na.rm=TRUE);
	gPrime = gVal*(q75-q25);
	btm = q25 - gPrime;
	top = q75 + gPrime;
	trimStr = paste(trimStr, "outlier labeling: <", btm, " >", top);
  } else {
    # remove 3rd std dev
    btm = quantile(trimmedData[[colName]], 0.021, na.rm=TRUE);
    top = quantile(trimmedData[[colName]], 0.979, na.rm=TRUE);
	trimStr = paste(trimStr, "3rd std. dev.: <", btm, " >", top);
  }
  myP2(opl, trimStr);
  btmCnt = nrow(trimmedData[trimmedData[[colName]]<btm,]);
  topCnt = nrow(trimmedData[trimmedData[[colName]]>top,]);
  if (is.null(btmCnt)) btmCnt = 0
  if (is.null(topCnt)) topCnt = 0
  myP3(opl, btmCnt, ":", topCnt);
  if ((btmCnt+topCnt)>0) {
    if (doWinsor) {
      myP3(opl, "Pre, min:", min(trimmedData[[colName]], na.rm=TRUE), "max:", max(trimmedData[[colName]], na.rm=TRUE));
	  trimmedData[trimmedData[[colName]]<btm, colName] = btm;
	  trimmedData[trimmedData[[colName]]>top, colName] = top;
      myP3(opl, "Post, min:", min(trimmedData[[colName]], na.rm=TRUE), "max:", max(trimmedData[[colName]], na.rm=TRUE));
      myP2(opl, paste((topCnt+btmCnt), "rows replaced."));
    } else {
      myP3(opl, "Pre, min:", min(trimmedData[[colName]], na.rm=TRUE), "max:", max(trimmedData[[colName]], na.rm=TRUE));
      trimmedData = as.data.frame(trimmedData[trimmedData[[colName]]>=btm,]);  
      namesV      = names(trimmedData)
      trimmedData = as.data.frame(trimmedData[trimmedData[[colName]]<=top,]);  
      namesV      = names(trimmedData)
      myP3(opl, "Post, min:", min(trimmedData[[colName]], na.rm=TRUE), "max:", max(trimmedData[[colName]], na.rm=TRUE));
      myP2(opl, paste((topCnt+btmCnt), "rows dropped."));
    }
  } else myP2(opl, "No outliers found.")
  return(trimmedData);
}

# <myCenter>
# Centers column
myCenter <- function(dataIn, colName, colOut="", geoMean=FALSE, opl=2) {
  # XXX fix to use geomean when asked
  if (colOut=="") {
    colOut = colName;
    myP2(opl, "Centering column '", colName, "' in place...");
  } else myP2(opl, "Centering column '", colName, "' -> '", colOut, "'...");
  centeredData = as.data.frame(dataIn);
  meanVal      = mean(centeredData[[colName]], na.rm=TRUE);
  centeredData[[colOut]] = centeredData[[colName]] - meanVal;
  return(centeredData);
}

# <myScale>
# Scales column
myScale <- function(dataIn, colName, colOut="", scaleBy=1, opl=2) {
  if (colOut=="") {
    colOut = colName;
    myP2(opl, "Scaling column '", colName, "' by", scaleBy, "in place...");
  } else myP2(opl, "Scaling column '", colName, "' by", scaleBy, "-> '", colOut, "'...");
  normalizedData = as.data.frame(dataIn);
  normalizedData[[colOut]] = normalizedData[[colName]] / (1.0*scaleBy);
  return(normalizedData);
}

# <myNormalize>
# Normalizes column
myNormalize <- function(dataIn, colName, colOut="", opl=2) {
  myP2(opl, "Normalizing column", colName);
  maxVal = max(dataIn[[colName]], na.rm=TRUE);
  return(myScale(dataIn, colName, colOut, scaleBy=maxVal));
}

# <myTransform>
# Transforms column
# Options: log, logit
myTransform <- function(dataIn, colName, colOut="", type="log", opl=2) {
  if (colOut=="") {
    colOut = colName;
	colStr = "in place...";
  } else colStr = paste("-> '", colOut, "'...");
  myP2(opl, paste(type, "-transforming column '", sep=""), colName, "'", colStr);

  transData = as.data.frame(dataIn);
  if (type=="log") transData[[colOut]] = log(as.numeric(transData[[colName]]))
  else if (type=="logit") transData[[colOut]] = qlogis(as.numeric(transData[[colName]]));
  return(transData);
}

# <myAugmentZero>
# Adds a delta value to every zero value in the column
myAugmentZero <- function(dataIn, colName, colOut="", useDelta=0.01, useHalfMin=TRUE, opl=2) {
  transData = as.data.frame(dataIn);
  if (useHalfMin) {
    nonZeroR = transData[transData[[colName]]>0,];
    minDelta = min(nonZeroR[[colName]], na.rm=TRUE) / 2.0;
	delStr   = paste("d=", minDelta, "(half-min)");
  } else {
    minDelta = useDelta;
	delStr   = paste("d=", useDelta);
  }
  if (colOut=="") {
    colOut = colName;
    myP2(opl, "Augmenting zeroes in column '", colName, "' in place...");
  } else myP2(opl, "Augmenting zeroes in column '", colName, "' -> '", colOut, "'...");
  
  transData[[colOut]] = transData[[colName]];
  numZeroes = 0;
  for (i in 1:nrow(transData)) {
    iVal = transData[[colOut]][i];
	if ((!is.na(iVal))&&((iVal==0))) {
	  transData[[colOut]][i] = minDelta;
	  numZeroes = numZeroes + 1;
	}
  }  
  if (numZeroes>0) myP2(opl, numZeroes, "zeroes augmented by", delStr);
  return(transData);
}

# <myAddIndicator>
# Adds an indicator value
myAddIndicator <- function(dataIn, testCol, indCol="", minTrue, lowTrue=FALSE, opl=2) {
  if (indCol=="") indCol = "indicator";
  if (lowTrue) compStr = "<"
  else compStr = ">=";
  myP2(opl, "Adding indicator column '", indCol, "' which = 1 if ", testCol, compStr, minTrue, "...");
  indData = as.data.frame(dataIn);
  
  indTest <- function(x) {
    if (x<minTrue) {
	  if (lowTrue) indVal = 1.0
	  else indVal = 0.0;
	} else {
	  if (lowTrue) indVal = 0.0
	  else indVal = 1.0;
	}
	return(indVal);
  }
  
  indData[[indCol]] = sapply(indData[[testCol]], indTest);
  return(indData);
}

# <myFindDuplicates>
myFindDuplicates <- function(dataIn, testCol, grpCol=NULL, valueColV=c(), mergeDuplicates=FALSE) {
  dupData = as.data.frame(dataIn);
  testGroups = !is.null(grpCol);
  if (testGroups) {
    myP("Testing uniqueness of column '",testCol,"' with group '",grpCol,"'...");  
    grpV = unique(dupData[[grpCol]]);
  } else myP("Testing uniqueness of column '",testCol,"'...");
  testV = unique(dupData[[testCol]]);
  
  dupFound  = FALSE;
  noDupData = NULL;
  for (x in testV) {
    testRows = dupData[dupData[[testCol]]==x,];
    if (testGroups) {
	  newRows = NULL;
	  for (g in grpV) {
	    grpRows = testRows[testRows[[grpCol]]==g,];
		numRows = nrow(grpRows);
		if (numRows>1) {
		  dupFound = TRUE;
		  if (mergeDuplicates) {
		    if (length(valueColV)==0) myP("Cannot merge duplicates with no value columns.")
			else {
  		      myP("Merging", numRows,"duplicates for", testCol,"=",x, ",",grpCol,"=",g,"...");
			  # drop all duplicated rows
			  newGrpRow = grpRows[1,];
			  for (v in valueColV) newGrpRow[,v] = mean(grpRows[,v]);
			  grpRows   = newGrpRow;
			}
		  } else myP("",numRows,"duplicates found for", testCol,"=",x, ",",grpCol,"=",g,".");
		} 
		
	    if (is.null(newRows)) newRows = grpRows
     	else newRows = rbind(newRows, grpRows);
	  }
	  
	  testRows = newRows;
	  
	} else if (nrow(testRows)>1) {
	  numRows   = nrow(testRows);
      dupFound = TRUE;
	  if (mergeDuplicates) {
	    if (length(valueColV)==0) myP("Cannot merge duplicates with no value columns.")
		else {
  		  myP("Merging", numRows,"duplicates for", testCol,"=",x, "...");
		  # drop all duplicated rows
		  newRow = testRows[1,];
		  for (v in valueColV) newRow[,v] = mean(testRows[,v]);
		  testRows = newRow;
		}
	  } else myP("",numRows,"duplicates found for", testCol,"=",x,".");
	} 
	
	# add in the row(s) for testCol=x
	if (is.null(noDupData)) noDupData = testRows
	else noDupData = rbind(noDupData, testRows);
  }
  
  if (!dupFound) myP("No duplicates found.");
  return(noDupData);
}

# <trim start/end NA>
# Pretty inefficient and sloppy
trimNA <- function(tsV, opl=2) {
  index = 0;
  isTS = is.ts(tsV);
  if (isTS) startTime = start(tsV)[1];
  for (i in 1:length(tsV)) {
    if (is.na(tsV[i])) index = i
	else break;
  }
  if (index>0) {
    myP3(opl, "Trimmed", index, "NA values from front.");
    tsV = tsV[-1:-index];
	if (isTS) startTime = startTime + index;
  } 

  index = 0;
  for (i in length(tsV):1) {
    if (is.na(tsV[i])) index = i
	else break;
  }
  if (index>0) {
    myP3(opl, "Trimmed", (length(tsV) - index), "NA values from back.");
    tsV = tsV[1:(index-1)];
  }

  if (isTS) return(ts(tsV, start=startTime))
  else return(tsV);
}

# <longestRun>
# Returns the longest continuous run of non-NA data
# Returns as vector of row numbers
longestRun <- function(dataV) {
  firstNonNA = 0;
  longestRun = 0;
  
  tempNonNA  = 0;
  lengthCur  = 0;
  for (i in 1:length(dataV)) {
    if (is.na(dataV[i])) {
	  # found NA
	  if (lengthCur>0) {
	    # end current run and set if longest
		if (lengthCur>longestRun) {
		  firstNonNA = tempNonNA;
		  longestRun = lengthCur;
		}
		# reset current
		tempNonNA = 0;
		lengthCur = 0;
	  }
	} else {
	  # start new run?
	  if (lengthCur==0) tempNonNA = i;
	  lengthCur = lengthCur + 1;
	}
  }
  # check final run
  if (lengthCur>longestRun) {
    firstNonNA = tempNonNA;
	longestRun = lengthCur;
  }
  
  if (longestRun==0) return(c(0))
  else return(c(firstNonNA:(firstNonNA + longestRun - 1)));
    
}

# <dropOutliers>
# Drops values less than first quantile or greater than second quantile
dropOutliers <- function(dataV, quantiles=c(0.025, 0.975)) {
  qV = quantile(dataV, quantiles);
  dataV = dataV[dataV>qV[1]];
  dataV = dataV[dataV<qV[2]];
  return(dataV);
}

# READ/WRITE FUNCTIONS

# <saveMatrixToCSVFile>
# writes the matrices to a CSV file
# uses "as.matrix" not "my.matrix" because we *want* vectors transposed
saveMatrixToCSVFile <- function(fileName, theM, dirName="output", label=NULL, opl=1) {
  if (!is.null(label)) outputM = rbind(label, as.matrix(theM))
  else outputM = theM
  fileStr = paste(fileName, "csv", sep=".")
  if (dirName!="") fileStr = paste(dirName, fileStr, sep="/")
  write.csv(outputM, file=fileStr)
  myP1(opl, "File", fileName, "written.")
}

loadMatrixFromCSVFile <- function(fileName, dirName="output", hasHeader=TRUE, dropFirstCol=TRUE, opl=1, altType=NULL) {
  if (is.null(altType)) fileStr = paste(fileName, "csv", sep=".")
  else fileStr = paste(fileName, altType, sep=".")
  if (dirName!="") fileStr = paste(dirName, fileStr, sep="/")
  if (!file.exists(fileStr)) {
    myP2(opl, "File", fileName, "not found. Returning NULL...")
	return(NULL)
  }
  theM = read.csv(fileStr, header=hasHeader, as.is=TRUE)
  myP1(opl, "File", fileName, "loaded.")
  if (dropFirstCol) return(theM[,-1])
  else return(theM)
}

isCSVFile <- function(fileName, dirName="output") {
  fileStr = paste(fileName, "csv", sep=".")
  if (dirName!="") fileStr = paste(dirName, fileStr, sep="/")
  return(file.exists(fileStr))
}

# <myMeltAndSave>
# Reads a CSV, melts, and resaves
myMeltAndSave <- function(fileName, prefix="M", id, variable, value) {
  shortM = loadMatrixFromCSVFile(fileName)
  longM  = melt(shortM, id.vars=id, variable.name=variable, value.name=value)
  # need to add the l_id columns by hand again
  saveMatrixToCSVFile(combine(prefix, fileName), longM)
}

# TS FUNCTIONS

# <tsMovingAverage>
# Performs a moving average on the column
# Uses R cookbook code
tsMovingAverage <- function(theTS, n=1, centered=FALSE, label="", opl=2) {
  if (label=="") label = "time series"
  myP2(opl, "Performing", n, "time step moving average on", label, "...")

#  transData = as.data.frame(dataIn);
  # code from R cookbook
  if (n>1) aveV = movingAverage(theTS, n=n, centered=centered)
  else aveV = theTS
  return(ts(aveV, start=start(theTS)))
}

# <tsFillGaps>
# Fills in holes in a time series by extrapolating linearly across each
tsFillGaps <- function(theTS, label="", opl=2, fillEnds=TRUE) {
  if (label=="") label = "time series"
  myP2(opl, "Filling gaps in", label, "...")

  filledV   = theTS
  gapLength = 0
  for (i in 1:length(filledV)) {
    if (is.na(filledV[i])) {
	  if (gapLength==0) gapStart = i
	  gapLength = gapLength + 1
	} else {
	  if (gapLength==0) next
	  # value after gap
	  afterVal  = filledV[i]
	  # value before gap
	  if (gapStart==1) {
	    if (!fillEnds) {
		  gapLength = 0
		  next
		} else beforeVal = afterVal
	  } else beforeVal = filledV[gapStart - 1]
	  # step size within gap (+1 to incorporate start/end, so gap=1 is halfway)
	  stepSize  = (afterVal - beforeVal) / (gapLength + 1)
	  # fill in missing values linearly
	  for (j in 1:gapLength) filledV[gapStart+(j-1)] = beforeVal + j*stepSize
	  # all done
	  if (gapStart==1) myP3(opl, "Filled starting gap (length =", gapLength, ").")
	  else myP3(opl, "Filled gap (position =", gapStart, ", length =", gapLength, ", step =", stepSize, ").")
      gapLength = 0
	} 
  }
  
  if (fillEnds && (gapLength>0)) {
    # gap at the end of the TS
	lastVal = filledV[gapStart - 1]
	filledV[gapStart:length(filledV)] = lastVal
    myP3(opl, "Filled end gap (length =", gapLength, ").")
  }
  
  return(ts(filledV, start=start(theTS)))
}

# <tsOverlap>
# Returns TRUE if two time series overlap by at least minOver units
tsOverlap <- function(ts1, ts2, minOver=5) {
  time1 = time(trimNA(ts1));
  time2 = time(trimNA(ts2));
  return(length(intersect(time1, time2))>=minOver);
}

# <tsFirstDiff>
# Returns the first difference of a time series
tsFirstDiff <- function(theTS, label="", opl=2) {
  if (label=="") label = "time series"
  myP2(opl, "Performing first difference on", label, "...");

  diffV     = c();
  for (i in 1:(length(theTS)-1)) {
    diff  = theTS[i+1] - theTS[i];
    diffV = c(diffV, diff);
  }
  
  return(ts(diffV, start=start(theTS)));
}

# <tsFitTrend>
# Returns fit and residuals (detrended)
tsFitTrend <- function(theTS, centerYear=TRUE, label="", opl=2) {
  if (label=="") label = "time series"
  myP2(opl, "Fitting linear trend to", label, "...");

  if (centerYear) yr = time(theTS) - mean(time(theTS))
  else yr = time(theTS);
  tsFit   = lm(theTS~yr, na.action=na.exclude);
  tsCoef  = coef(tsFit)[2];
  tsR     = resid(tsFit);
  # detrended residuals
  dTS     = ts(resid(tsFit), start=start(theTS));
  return(list(fit=tsFit, coef=tsCoef, resid=tsR, detrendTS=dTS));
}

# <tsCUSUM>
# Applies CUSUM to time series
# Sp = sum(x) from 1-p  - pk
tsCUSUM <- function(theTS, label="", opl=2) {
  if (label=="") label = "time series"
  myP2(opl, "Performing CUSUM on", label, "...");

  k  = mean(theTS);
  cV = c();
  for (p in 1:length(theTS)) {
    Sp = sum(theTS[1:p]) - p*k;
	cV = c(cV, Sp);
  } 
  return(ts(cV, start=start(theTS)));  
}

# <examineTS>
# Shows the time series under various transformations
examineTS <- function(theTS, fillGaps=TRUE, diffNum=1, aveNum=5, incACF=TRUE, incPACF=TRUE, label="") {
  mainLabel = "Raw"
  if (fillGaps) theTS = tsFillGaps(theTS, label=label)
  if (diffNum==1) { 
    theTS = tsFirstDiff(theTS, label=label)
	mainLabel = "First-differenced"
  }
  if (aveNum>1) {
    theTS = tsMovingAverage(theTS, n=aveNum, centered=TRUE, label=label)
	mainLabel = paste("Moving Average =", aveNum)
  } 
  deTS  = tsFitTrend(theTS, label=label)$detrendTS
  csTS  = tsCUSUM(theTS, label=label)
        
  myWin();
  if (incPACF) par(mfrow=c(3,3))
  else par(mfrow=c(3,2))
  plot.ts(theTS, main=mainLabel)
  if (incACF) acf(theTS)
  if (incPACF) pacf(theTS)
  plot.ts(deTS, main="Detrended")
  if (incACF) acf(deTS)
  if (incPACF) pacf(deTS)
  plot.ts(csTS, main="CUSUM")
  if (incACF) acf(csTS)
  if (incPACF) pacf(csTS)
 
}

showSmoothedTS <- function(theTS, fillGaps=TRUE, diffNum=0, aveV=c(5,10), label="") {
  mainLabel = "Raw"
  if (fillGaps) theTS = tsFillGaps(theTS, label=label)
  if (diffNum==1) { 
    theTS = tsFirstDiff(theTS, label=label)
	mainLabel = "First-differenced"
  }
        
  myWin();
  numAve = length(aveV)
  par(mfrow=c(numAve+1,1))
  plot.ts(theTS, main=mainLabel)
  for (i in 1:numAve) {
    aveTS = tsMovingAverage(theTS, n=aveV[i], label=label)
	mainLabel = paste("Moving Average =", aveV[i])
    plot.ts(aveTS, main=mainLabel)
  }
}

