##' Plot the Pareto Front with step functions.
##' 
##' @title Pareto Front visualization
##' @param nondominatedPoints points considered to plot the Pareto Front with segments, matrix with one point per column,
##' @param add boolean indicating whether a new graphic should be drawn.
##' @param ... additional values to be passed to the \code{\link[graphics]{lines}} function.
##' @export
##' @examples
##' #------------------------------------------------------------
##' # Simple example
##' #------------------------------------------------------------
##' 
##' x <- c(0.2, 0.4, 0.6, 0.8)
##' y <- c(0.8, 0.7, 0.5, 0.1)
##' 
##' plot(x, y, col = "green", pch = 20) 
##' 
##' plotParetoEmp(cbind(x, y), col = "green")
##' ## Alternative
##' plotParetoEmp(cbind(x, y), col = "red", add = FALSE)

plotParetoEmp <- function(nondominatedPoints, add = TRUE, ...){
  #     temp <- cbind(nondominatedPoints,c(min(nondominatedPoints[1,])-0.0001,10000))
  #    
  #     # Construction des step functions
  #     f1_pareto <- sort(temp[1,],index.return=TRUE)
  #     f2_pareto <- temp[2,f1_pareto$ix]
  #     f1_pareto <- f1_pareto$x
  #     
  #     steptest <- stepfun(f1_pareto,c(f2_pareto,f2_pareto[length(f2_pareto)]),right=TRUE)
  #     
  #     plot.stepfun(steptest,...)
  
  if(add == FALSE){
    plot(nondominatedPoints, pch = ".",...)
  }
  

  temp <- nondominatedPoints[order(nondominatedPoints[,1]),]
  lines(c(temp[1,1],temp[1,1]),c(abs(20*temp[1,2]),temp[1,2]),...) #first segment
  
  lines(c(temp[dim(temp)[1], 1], abs(20*temp[dim(temp)[1],1])),
        c(temp[dim(temp)[1], 2], temp[dim(temp)[1], 2]),...) #last segment
  
  # Segment in between
  for (i in 1:(dim(temp)[1]-1)){
    lines(c(temp[i,1], temp[i+1,1]),c(temp[i,2],temp[i,2]),...) #horizontal part
    lines(c(temp[i+1,1],temp[i+1,1]),c(temp[i,2],temp[i+1,2]),...) #vertical part
  }
  
}
