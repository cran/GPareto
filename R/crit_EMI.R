#' Expected Maximin Improvement with respect to the current Pareto front with Sample Average Approximation.
#' The semi-analytical formula is used in the bi-objective scale if the Pareto front is in [-2,2]^2,
#' for numerical stability reasons. 
#' To avoid numerical instabilities, the new point is penalized if it is too close to an existing observation.
#' 
#' @title Expected Maximin Improvement with m objectives
#' 
#' @param x a vector representing the input for which one wishes to calculate \code{EMI},
#' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
#' @param paretoFront (optional) matrix corresponding to the Pareto front of size \code{[n.pareto x n.obj]}, or any reference set of observations, 
#' @param critcontrol optional list with arguments (for more than 2 objectives only): 
#'   \itemize{
#'   \item \code{nb.samp} number of random samples from the posterior distribution,
#'        default to \code{50}, increasing gives more reliable results at the cost of longer computation time;
#'   \item  \code{seed} seed used for the random samples.
#'   }
#'        Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.
#' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation 
#'        has to be taken into account. 
#' @return The Expected Maximin Improvement at \code{x}.
#' @details It is recommanded to scale objectives, e.g. to \code{[0,1]}. 
#' If the Pareto front does not belong to [-2,2]^2, then SAA is used. 
#' @seealso \code{\link[DiceOptim]{EI}} from package DiceOptim, \code{\link[GPareto]{crit_EHI}}, \code{\link[GPareto]{crit_SUR}}, \code{\link[GPareto]{crit_SMS}}.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats dnorm integrate
#' @references J. D. Svenson & T. J. Santner (2010), Multiobjective Optimization of Expensive Black-Box
#' Functions via Expected Maximin Improvement, Technical Report.\cr
#' 
#' J. D. Svenson (2011), \emph{Computer Experiments: Multiobjective Optimization and Sensitivity Analysis}, Ohio State University, PhD thesis. 
#' 
#' @examples
#' #---------------------------------------------------------------------------
#' # Expected Maximin Improvement surface associated with the "P1" problem at a 15 points design
#' #---------------------------------------------------------------------------
#' set.seed(25468)
#' library(DiceDesign)
#' 
#' n_var <- 2 
#' f_name <- "P1" 
#' n.grid <- 21
#' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
#' n_appr <- 15 
#' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
#' response.grid <- t(apply(design.grid, 1, f_name))
#' Front_Pareto <- t(nondominated_points(t(response.grid)))
#' mf1 <- km(~., design = design.grid, response = response.grid[,1])
#' mf2 <- km(~., design = design.grid, response = response.grid[,2])
#' 
#' EMI_grid <- apply(test.grid, 1, crit_EMI, model = list(mf1, mf2), paretoFront = Front_Pareto,
#'                      critcontrol = list(nb_samp = 20))
#' 
#' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
#'                matrix(EMI_grid, nrow = n.grid), main = "Expected Maximin Improvement",
#'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
#'                plot.axes = {axis(1); axis(2);
#'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
#'                             }
#'               )
crit_EMI <- function(x, model, paretoFront = NULL, critcontrol = list(nb.samp = 50, seed = 42),
                     type = "UK"){
  nobj     <- length(model)
  if(nobj < 2){
    cat("Incorrect Number of objectives \n")
    return(NA)
  }
  
  if (is.null(critcontrol$nb.samp)) critcontrol$nb.samp <- 50
  if (is.null(critcontrol$seed)) critcontrol$seed <- 42
  
  if(nobj == 2)
    if(max(abs(model[[1]]@y)) <=2 && max(abs(model[[2]]@y)) <=2)
      return(EMI_2d(x = x, model = model, critcontrol = critcontrol, type = type, paretoFront = paretoFront))
  
  if(is.null(critcontrol)){
    critcontrol <- list()
  }
  
  critcontrol$type <- "maximin"
  return(SAA_mEI(x = x, model = model, critcontrol = critcontrol, type = "UK", paretoFront = paretoFront))
}






