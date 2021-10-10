# Script from Aalto et al. 2019 ('Catastrophic Mortality, Allee Effects, and Marine Protected Areas), written by Emilius A. Aalto

library(ggplot2);
library(grid);
#source("/Users/julie/Documents/Stanford/First Year/TAC Summer Project/IPM/abalone code/common/myUtil.r");
#source("/Users/julie/Documents/Stanford/First Year/TAC Summer Project/IPM/abalone code/common/myData.r");
source('myUtil.r')
source('myData.r')

# My graphics functions

# <myBaseGraph>
# Base graphing function
# Creates a new window for each graph
myBaseGraph <- function(dataIn=NULL, xCol, yCol=NULL, grp=1, clr="", 
					    title=NULL, xLabel=NULL, yLabel=NULL,
						type="scatter", overlayX="", overlayY="", binWidth=-1,
					    addTrend=TRUE, overallTrend=FALSE, useLines=FALSE, custModel=NULL,
						returnPlot=FALSE, opl=1) {
  if (is.null(dataIn)) {
    # assume that x/y are vectors
	myP("No data frame given: printing plot as raw vectors...");
	print(qplot(xCol, yCol, geom="point"));
	return(invisible());
  } else {
    if (is.null(yCol)) thePlot = ggplot(dataIn, aes_string(x=xCol, group=grp, color=clr))
	else thePlot = ggplot(dataIn, aes_string(x=xCol, y=yCol, group=grp, color=clr));
  }
  if (type=="box") {
    myP1(opl, "Box plot for", yCol, "vs.", xCol, "with group =", grp);
    thePlot = thePlot + geom_boxplot();
  } else if (type=="hist") {
    myP1(opl, "Histogram for", xCol);
	if (binWidth>0) thePlot = thePlot + geom_histogram(binwidth=binWidth)
	else thePlot = thePlot + geom_histogram();
  } else if (type=="scatter") {
    myP1(opl, "Graphing", yCol, "vs.", xCol, "with group =", grp);
	if (useLines) thePlot = thePlot + geom_line()
	else thePlot = thePlot + geom_point();
	# add overall trend
	if (addTrend) {
	  if (overallTrend) grp = 1;
	  if (is.null(custModel)) {
	    myP1(opl, "Fitting glm model...");
		thePlot = thePlot + geom_smooth(aes_string(x=xCol, y=yCol, group=grp), method="glm");
	  } else {
	    myP1(opl, "Fitting custom model...");
		fitData = myPrediction(dataIn, custModel);
        thePlot = ggplot(fitData, aes_string(x=xCol, y=yCol, group=grp, color=clr));
		thePlot = thePlot + geom_point(aes_string(x=xCol, y=yCol, group=grp));
		thePlot = thePlot + geom_smooth(aes_string(x=xCol, y="predicted", group=grp), method="glm");
	  }
    }
  }
  
  if (overlayX!="") {
    if (clr=="red") altColor = "green"
	else altColor = "red";
    thePlot = thePlot + geom_line(aes_string(x=overlayX, y=overlayY), color=altColor);
  }
  
  if (is.null(title)) title = combine(xCol, yCol, " vs. ");
  thePlot = thePlot + labs(title=title);
  if (!is.null(xLabel)) thePlot = thePlot + xlab(xLabel);
  if (!is.null(yLabel)) thePlot = thePlot + ylab(yLabel);
  if (!returnPlot) myPPlot(thePlot)
  else return(thePlot);
}

# <myGraph>
# Convenience function for x-y graphs
myGraph <- function(dataIn, xCol, yCol, grp=1, clr=NULL, 
					title=NULL, xLabel=NULL, yLabel=NULL,
					addTrend=TRUE, overallTrend=FALSE, useLines=FALSE, custModel=NULL,
					returnPlot=FALSE, opl=1) {
  myBaseGraph(dataIn, xCol, yCol, grp, clr, title, xLabel, yLabel, type="scatter",
	          addTrend=addTrend, overallTrend=overallTrend, useLines=useLines, custModel=custModel, returnPlot=returnPlot, opl=opl);
}

# <myBox>
# Convenience function for boxplots
myBox <- function(dataIn, xCol, yCol, grp=1, clr=NULL,
				  title=NULL, xLabel=NULL, yLabel=NULL,
				  returnPlot=FALSE, opl=1) {
  myBaseGraph(dataIn, xCol, yCol, grp, clr, title, xLabel, yLabel, 
              type="box", addTrend=FALSE, returnPlot=returnPlot, opl=opl);
}

# <myHist>
# Convenience function for histograms
myHist <- function(dataIn, xCol, binWidth=-1,
			       title=NULL, xLabel=NULL, yLabel=NULL, overlayX="", overlayY="",
				   returnPlot=FALSE, opl=1) {
  myBaseGraph(dataIn, xCol, title=title, xLabel=xLabel, binWidth=binWidth,
              type="hist", addTrend=FALSE, overlayX=overlayX, overlayY=overlayY, returnPlot=returnPlot, opl=opl);
}

# <myQQNorm>
# Convenience function to display quant-quant graph for normal distribution
myQQNorm <- function(x, y=NULL, title="", addLine=TRUE) {
  myQQDistrib(x, y, distrib="normal", title, addLine);
}

# <myQQNorm>
# Displays quant-quant graph
myQQDistrib <- function(x, y=NULL, distrib="normal", title="", addLine=TRUE, normalizeX=TRUE) {
  myWin();

  if (normalizeX) {
    normX = (x - mean(x)) / sd(x);
  } else normX = x;
  
  if (is.null(y)) {
    if (distrib=="normal") {
      qqnorm(normX, main=title);
	  if (addLine) qqline(normX);
	}
  } else {
    qqplot(normX, y, main=title);
	abline(0,1);
  }
}

# <graphDataSet>
# Graphs a set of x-y-grp-color-title-label vectors
graphDataSet <- function(dataIn, xColV, yColV, grpV=NULL, colorV=NULL, titleV=NULL, xLabelV=NULL, yLabelV=NULL, useCustModel=FALSE, opl=1) {
  myP("Graphing data set...");
  for (i in 1:length(xColV)) {
    xCol = xColV[i];
	yCol = yColV[i];
    if ((is.null(grpV))||(is.na(grpV[i]))) grp = 1
	else grp = grpV[i];
    if ((is.null(colorV))||(is.na(colorV[i]))) color = NULL
	else color = colorV[i];
    if ((is.null(titleV))||(is.na(titleV[i]))) title = NULL
	else title = titleV[i];
    if ((is.null(xLabelV))||(is.na(xLabelV[i]))) xLabel = NULL
	else xLabel = xLabelV[i];
    if ((is.null(yLabelV))||(is.na(yLabelV[i]))) yLabel = NULL
	else yLabel = yLabelV[i];
    if (is.factor(dataIn[[xCol]])) {
	  grp = xCol;
	  myBox(dataIn, xCol, yCol, grp, color, title, xLabel, yLabel, opl=opl);
	} else {
	  # XXX currently pulling from eel/model.r because of problems passing the LME object
	  if (useCustModel) model = getBestModel(dataIn, yCol)
	  else model = NULL;
	  myGraph(dataIn, xCol, yCol, grp, color, title, xLabel, yLabel, custModel=model, opl=opl);
	}
  }
}

# <stripPlot>
# Strips elements off a plot, useful for sub-plotting.
stripPlot <- function(thePlot, stripX=FALSE, stripY=FALSE, stripTitle=FALSE, stripLegend=FALSE, stripGrid=FALSE) {
  if (stripX) thePlot = thePlot + theme(axis.title.x=element_blank());
  if (stripY) thePlot = thePlot + theme(axis.title.y=element_blank());
  if (stripTitle) thePlot = thePlot + theme(plot.title=element_blank());
  if (stripLegend) thePlot = thePlot + theme(legend.position="none");
  if (stripGrid) thePlot = thePlot + theme(panel.grid=element_blank());
  return(thePlot);
}

# <addStandardThemes>
# Convenience function to add standard sizing, etc.
addStandardThemes <- function(thePlot, axisTextSize=15, axisTitleSize=20, legendTextSize=20) {
  thePlot = thePlot + theme_bw() +
             theme(axis.text.x = element_text(size=axisTextSize), axis.text.y = element_text(size=axisTextSize),
				   axis.title.x = element_text(size=axisTitleSize), axis.title.y = element_text(size=axisTitleSize),
				   legend.text = element_text(size=legendTextSize),
				   panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(thePlot)
}

# <myPPlot>
# Convenience to open new window and print plot
# Allows printing subplots
myPPlot <- function(thePlot, w=16, h=12, plotPos=c(1,1), invertY=FALSE,
                    newWindow=TRUE, plotDim=c(1,1)) {
  if (newWindow) {
    myWin(w=w, h=h)
    pushViewport(viewport(layout=grid.layout(plotDim[1],plotDim[2])))
  }
  
  plotRow = plotPos[1]
  if (invertY) {
    # invert row order to match MatLab style (y=1 at the bottom, not the top)
    if (plotDim[1]>1) plotRow = plotDim[1] - plotRow + 1
  }
  print(thePlot, vp=viewport(layout.pos.row=plotRow, layout.pos.col=plotPos[2]))
}

getSubplotPosition <- function(plotNum, plotDimV) {
  totalPlots = plotDimV[1] * plotDimV[2]
  while (plotNum > totalPlots) plotNum = plotNum - totalPlots
  xPos = ((plotNum-1) %/% plotDimV[2]) + 1
  yPos = ((plotNum-1) %% plotDimV[2]) + 1 
  return(c(xPos, yPos))
}

myPPDF <- function(thePlot, filename="") {
  if (filename=="") filename = "image.pdf"
  else filename = paste(filename, ".pdf", sep="")
  pdf(filename)
  print(thePlot)
  dev.off()
}

# <myPTiff>
# Convenience to save as 600dpi TIFF
# Note that this has issues if h&w are too big for the res.
# Can take a list with members p1, p2, p3, etc. if ggSave is FALSE
# XXX list functionality not currently working
myPTiff <- function(thePlot, filename="", h, w, res=600, ggSave=TRUE, catchError=TRUE,
					plotDim=c(1,1)) {
  if (filename=="") filename = "image.tiff"
  else filename = paste(filename, ".tiff", sep="")
  
  if (ggSave) {
    if (catchError) {
	  myP("Saving TIFF with h=", h, ",w=", w)
      tryCatch(ggsave(filename, height=h, width=w, units='in', dpi=res),
	           error=function(e) { 
#			     dev.off();
             	 myP("Trying TIFF with h=", h*0.75, ",w=", w*0.75);
				 tryCatch(ggsave(filename, height=h*0.75, width=w*0.75, units='in', dpi=res),
	                      error=function(e) { 
#						    dev.off();
							myP("Trying TIFF with h=", h*0.5, ",w=", w*0.5);
						    tryCatch(ggsave(filename, height=h*0.5, width=w*0.5, units='in', dpi=res),
		                             error=function(e) { 
#									   dev.off();
									   myP("Could not save TIFF..."); 
									 }); 
						  }); 
			   });
	} else ggsave(filename, height=h, width=w, units='in', dpi=res);
  } else {
    tryTiff <- function(filename, res, h, w) {
	  myP("Saving TIFF with h=", h, ",w=", w);
      tiff(filename, res=res, compression = "lzw", height=h, width=w, units="in");
	  if (length(thePlot)==1) print(thePlot)
	  else {
        grid.newpage();
     	pushViewport(viewport(layout=grid.layout(plotDim[1],plotDim[2])));
        vplayout <- function(x,y) { viewport(layout.pos.row=x, layout.pos.col=y); }
	    
		pRow = 1;
		pCol = 1;
		for (i in 1:length(thePlot)) {
		  pName = paste(p, i, sep="");
		  print(thePlot[pName], vp=vplayout(pRow, pCol));
		  pCol = pCol + 1;
		  if (pCol>plotDim[2]) {
		    pRow = pRow + 1;
			pCol = 1;
		  }
		}
	  }
      dev.off();
	}
  
    if (catchError) {
      tryCatch(tryTiff(filename, res, h, w), 
	           error=function(e) { 
			     dev.off();
             	 myP("Trying TIFF with h=", h*0.75, ",w=", w*0.75);
				 tryCatch(tryTiff(filename, res, h*0.75, w*0.75), 
	                      error=function(e) { 
 			                dev.off();
							myP("Trying TIFF with h=", h*0.5, ",w=", w*0.5);
						    tryCatch(tryTiff(filename, res, h*0.5, w*0.5), 
		                             error=function(e) { 			     
									   dev.off();
									   myP("Could not save TIFF..."); 
									 }); 
						  }); 
			   });
	} else tryTiff(filename, res, h, w);
  }
}

# <myWin>
# Convenience to open new window
myWin <- function(w=16, h=12) {
  quartz(w,h)
  grid.newpage()
}
