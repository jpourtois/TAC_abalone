# Script from Aalto et al. 2019 ('Catastrophic Mortality, Allee Effects, and Marine Protected Areas), written by Emilius A. Aalto

# My utility functions

# <makeGlobal>
# Assigns a variable to the global namespace
makeGlobal <- function(nameStr, value) { assign(nameStr, value, envir=.GlobalEnv) }

# <buildList>
# Makes a list from names and values
buildList <- function(nameV, valV) {
  bL = list()
  for (i in 1:length(nameV)) bL[[nameV[i]]] = valV[i]
  
  invisible(bL)
}

# <padVector>
# Pads out a vector to a given length
padVector <- function(theV, newLen, val=0) {
  oldLen = length(theV)
  if (oldLen<newLen) {
    for (i in 1:(newLen-oldLen)) {
	  theV[oldLen+i] = val
	}
  }
  invisible(theV)
}


# <combine>
# Combines two strings with a separator, unless one is ""
combine <- function(s1, s2, sep="_") {
  if (is.na(s1)||is.null(s1)||(s1=="")) {
    if (is.na(s2)||is.null(s2)||(s2=="")) return("_")
	else invisible(s2)
  } else if (is.na(s2)||is.null(s2)||(s2=="")) invisible(s1)
  else invisible(paste(s1, s2, sep=sep))
}

# <meanNZ>
# Calculates the mean of non-zero terms in a vector
meanNZ <- function(mV) {
  meanTot = 0
  meanN   = 0
  for (i in 1:length(mV)) {
    val = as.numeric(mV[i])
    if (val>0) {
	  meanTot = meanTot + val
	  meanN   = meanN + 1
	}
  }
  if (meanN==0) meanN = 1
  invisible(meanTot/meanN)
}

# <cz>
# Converts NA and NULL to zero when making a vector
cz <- function(...) {
  v = c()
  args = list(...)
  for (a in args) {
    if ((is.null(a))||(is.na(a))) v = c(v, 0)
    else v = c(v, a)
  }
  invisible(v)
}

# <contains>
# checks if a matrix contains the given value (handles NULL, NaN, and Inf also)
# can be passed a single value or a vector as well
contains <- function(mat, val) {
  mM = as.matrix(mat)
  contains = FALSE
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
	  if (is.nan(val)) contains = (contains||is.nan(mat[i,j]))
      else if (is.null(val)) contains = (contains||is.null(mat[i,j]))
      else if (is.infinite(val)) contains = (contains||is.infinite(mat[i,j]))
      else contains = (contains||my.equals(mat[i,j], val))
	}
  }
  
  invisible(contains)
}

# <my.unique>
# Extracts unique numeric values from a vector, with a minimum difference
my.unique <- function(mV, minDiff=0.001) {
  uV = c(mV[1]);
  for (i in 2:length(mV)) {
    mVal     = as.numeric(mV[i])
	isUnique = TRUE
    for (j in 1:length(uV)) {
      uVal = as.numeric(uV[j])
      if (my.equals(uVal, mVal, minDiff)) {
        isUnique = FALSE
        break
	  }		
	}
    if (isUnique) uV = c(uV, mVal)
  }
  invisible(uV)
}

# <my.gt>
# <my.gte>
# <my.lt>
# <my.lte>
# >, >=, <, <= with a min diff
# <my.equals>
# Returns true if the two numbers are within minDiff of each other
my.gt <- function(v1, v2, minDiff=0.001) { return(my.diffComp(v1,v2,gt=TRUE,minDiff=minDiff)) }
my.gte <- function(v1, v2, minDiff=0.001) { return(my.diffComp(v1,v2,gt=TRUE,eq=TRUE,minDiff=minDiff)) }
my.lt <- function(v1, v2, minDiff=0.001) { return(my.diffComp(v1,v2,lt=TRUE,minDiff=minDiff)) }
my.lte <- function(v1, v2, minDiff=0.001) { return(my.diffComp(v1,v2,lt=TRUE,eq=TRUE,minDiff=minDiff)) }
my.eq <- function(v1, v2, minDiff=0.001) { return(my.diffComp(v1,v2,eq=TRUE,minDiff=minDiff)) }
# <my.diffComp>
# Generic implementation
my.diffComp <- function(v1, v2, gt=FALSE, lt=FALSE, eq=FALSE, minDiff=0.001) { 
  isEqual = abs(as.numeric(v1) - as.numeric(v2)) < minDiff
  if (isEqual) return(eq) 
  else {
    if (v1>v2) return(gt)
	else return(lt)
  } 
}

# <my.equalsR> 

my.equalsR <- function(v1, v2, minDiff=0.001) {
  if (as.numeric(v2)!=0) {
    return(my.equals((as.numeric(v1)/as.numeric(v2)), 1.0, minDiff))
  } else {
    # default to diff, not ratio
    return(my.equals(v1, v2, minDiff))
  }
}


# <myF>
# shorthand for numerical output formatting
myF <- function(value, d=4) { format(value, digits=d) }

# <myP>
# shorthand for print(paste())
# additional versions for various levels of output
# Level 0: silent
# Level 1: standard output
# Level 2: verbose output
# Level 3: complete output
myPBase <- function(opl, bar, ...) { if (opl>=bar) print(paste(...)) }
myP <- function(...) { myPBase(1, 1, ...) }
myP1 <- function(opl, ...) { myPBase(opl, 1, ...) }
myP2 <- function(opl, ...) { myPBase(opl, 2, ...) }
myP3 <- function(opl, ...) { myPBase(opl, 3, ...) }

# Used for printing nice matrix
myPM <- function(theM) {
  for (i in 1:nrow(theM)) {
    rowStr = ""
    for (j in 1:ncol(theM)) {
	  theVal = theM[i,j]
	  if (is.na(as.numeric(theVal))) rowStr = paste(rowStr, theVal)
	  else rowStr = paste(rowStr, myF(as.numeric(theVal)))
	}
	print(rowStr)
  }
}

# Used for iterating towards a target
numAdjuster <- function(n, calcVal, target, wasLarger, adjBy, adjInv=FALSE, posN=FALSE, minDiff=0.001) {
  adj = adjBy
  if (adjInv) adj = -adjBy
  wasL = wasLarger
  if (my.equals(calcVal, target, minDiff)) {
    return(list(newN=n, equal=TRUE))
  } else {
	if (calcVal>target) {
	  if (wasL) newN = n - adj
	  else {
	    wasL = TRUE
	    adj  = adj / 2
		newN = n - adj
      }
	  if (posN&&(newN<=0)) {
	    wasL = FALSE
	    newN = adjBy / 2
	  }
	} else {
	  if (!wasL) newN = n + adj
	  else {
	    wasL = FALSE
	    adj = adj / 2
		newN = n + adj
	  }
	}
	# invert again to return
	if (adjInv) adj = -adj
  }
  return(list(newN=newN, wasLarger=wasL, adjBy=adj, equal=FALSE))
}

# <my.matrix>
# as.matrix transposes vectors into columns, so this makes sure single rows
# are cast as one row matrices
my.matrix <- function(datain) {
  if (!is.null(datain)) {
    dM = as.matrix(datain)
    if (ncol(dM)==1)
      dM = t(dM)
    dM
  } else NULL
}

# <se>
# Calculates the standard error of the sample  
se <- function(x) { sd(x)/sqrt(length(x)) }

# <Christmas calculator>
# Picks match-ups for Christmas gift giving
xmasCalc <- function(encrypt=FALSE) {
  pV = c("Emil", "Brie", "Karen", "Charley", "Cal", "Crystal")
  matchM = matrix(nrow=length(pV), ncol=2)
  matchM[,1] = pV
  
  totalRuns = 0
  while (TRUE) {
    # randomly sample from the giver vector
    totalRuns  = totalRuns + 1
    matchM[,2] = sample(pV)
	  valid      = TRUE
	  for (i in 1:length(pV)) {
	    giver     = matchM[i,1]
	    recipient = matchM[i,2]
	    if (giver==recipient) valid = FALSE
	    else if (giver==pV[1] && recipient==pV[2]) valid = FALSE
	    else if (giver==pV[2] && recipient==pV[1]) valid = FALSE
	    else if (giver==pV[3] && recipient==pV[4]) valid = FALSE
	    else if (giver==pV[4] && recipient==pV[3]) valid = FALSE
	    else if (giver==pV[5] && recipient==pV[6]) valid = FALSE
	    else if (giver==pV[6] && recipient==pV[5]) valid = FALSE
	    else {
	      for (j in i:length(pV)) if (matchM[j,1]==recipient && matchM[j,2]==giver) valid = FALSE
	    }
	  }
    if (valid) break
  }
  myP("Total runs:", totalRuns)

  # employ rot13 encryption plus generate some garbage
  if (encrypt) {
    recV = matchM[,2]
	encV = c();
	
	for (i in 1:length(recV)) {
	  prefL   = sample(2:4, 1)
	  sufL    = sample(2:4, 1)
	  prefStr = paste(sample(c(letters, LETTERS), prefL), collapse="")
	  sufStr  = paste(sample(c(letters, LETTERS), sufL), collapse="")
	  encStr  = paste(paste(prefStr, recV[i], sep=""), sufStr, sep="")
	  encStr  = rot13(encStr)
	  matchM[i,2] = encStr
	}
  }
  
  myPM(matchM)
}


# <rot13>
# ROT13 function  
rot13 <- function(x) {
  old <- paste(letters, LETTERS, collapse="", sep="")
  new <- paste(substr(old, 27, 52), substr(old, 1, 26), sep="")
  chartr(old, new, x)
}

# <myRprof>
# Wrapper to run Rprof around function and display results
myRprof <- function(theFunc, args=list()) {
  Rprof()
  print(system.time(invisible(do.call(theFunc, args))))
  Rprof(NULL)
  print(summaryRprof()$by.self)
}






