##' Executes \code{nsteps} iterations of multi-objective EGO methods to objects of class \code{\link[DiceKriging]{km}}.
##' At each step, kriging models are re-estimated (including covariance parameters re-estimation)
##'  based on the initial design points plus the points visited during all previous iterations;
##'  then a new point is obtained by maximizing one of the four multi-objective Expected Improvement criteria available. 
##' @title Sequential multi-objective Expected Improvement maximization and model re-estimation,
##'  with a number of iterations fixed in advance by the user
##' @details Extension of the function \code{\link[DiceOptim]{EGO.nsteps}} for multi-objective optimization.\cr
##' Available infill criteria with \code{crit} are: \cr
##' \itemize{
##' \item Expected Hypervolume Improvement (\code{EHI}) \code{\link[GPareto]{crit_EHI}},
##' \item SMS criterion (\code{SMS}) \code{\link[GPareto]{crit_SMS}},
##' \item Expected Maximin Improvement (\code{EMI}) \code{\link[GPareto]{crit_EMI}},
##' \item Stepwise Uncertainty Reduction of the excursion volume (\code{SUR}) \code{\link[GPareto]{crit_SUR}}.
##' }
##' Depending on the selected criterion, parameters such as reference point for \code{SMS} and \code{EHI} or arguments for \code{\link[GPareto]{integration_design_optim}} with \code{SUR} can be given with \code{critcontrol}.
##' Also options for \code{\link[GPareto]{checkPredict}} are available.
##' More precisions are given in the corresponding help pages. 
##' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##' @param fun the multi-objective function to be minimized (vectorial output), found by a call to \code{\link[base]{match.fun}},
##' @param cheapfun optional additional fast-to-evaluate objective function (handled next with class \code{\link[GPareto]{fastfun}}), which does not need a kriging model, handled by a call to \code{\link[base]{match.fun}},
##' @param crit choice of multi-objective improvement function: "\code{SMS}", "\code{EHI}", "\code{EMI}" or "\code{SUR}",
##' see details below,
##' @param nsteps an integer representing the desired number of iterations,
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation has to be taken into account,
##' @param cov.reestim optional boolean specifying if the kriging hyperparameters should be re-estimated at each iteration,
##' @param critcontrol optional list of parameters for criterion \code{crit}, see details,
##' @param optimcontrol an optional list of control parameters for optimization of the selected infill criterion:
##' \itemize{
##' \item{"\code{method}" can be set to "\code{discrete}", "\code{pso}" or "\code{genoud}". For "\code{discrete}", a matrix \code{candidate.points} must be given,
##' For "\code{pso}" and "\code{genoud}", specific parameters to the chosen method can also be specified  (see \code{\link[rgenoud]{genoud}} and \code{\link[pso]{psoptim}}).}
##' \item{Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.}
##' }
##' Option \code{notrace} can be set to \code{TRUE} to suppress printing of the optimization progresses. 
##' 
##' @param ... additional parameters to be given to the objective \code{fun}.
##' @export
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{par}}{: a data frame representing the additional points visited during the algorithm,}
##' \item{\code{values}}{: a data frame representing the response values at the points given in \code{par},}
##' \item{\code{nsteps}}{: an integer representing the desired number of iterations (given in argument),}
##' \item{\code{lastmodel}}{: a list of objects of class \code{\link[DiceKriging]{km}} corresponding to the last kriging models fitted.}
##' If a problem occurs during either model updates or criterion maximization, the last working model and corresponding values are returned.
##' }
##' 
##' @references 
##' M. T. Emmerich, A. H. Deutz, J. W. Klinkenberg (2011), Hypervolume-based expected improvement: Monotonicity properties and exact computation,
##' \emph{Evolutionary Computation (CEC)}, 2147-2154. \cr \cr
##' V. Picheny (2014), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
##' \emph{Statistics and Computing}. \cr \cr
##' T. Wagner, M. Emmerich, A. Deutz, W. Ponweiser (2010), On expected-improvement criteria for model-based multi-objective optimization.   
##' \emph{Parallel Problem Solving from Nature}, 718-727, Springer, Berlin. \cr \cr
##' J. D. Svenson (2011), \emph{Computer Experiments: Multiobjective Optimization and Sensitivity Analysis}, Ohio State university, PhD thesis. 
##' 
##' @importFrom stats runif pnorm qnorm
##' 
##' @examples
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' d <- 2 
##' 
##' fname <- ZDT3
##' n.grid <- 21
##' test.grid = expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' nappr <- 15 
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname))
##' Front_Pareto <- t(nondominated_points(t(response.grid)))
##' 
##' mf1 <- km(~., design = design.grid, response = response.grid[, 1])
##' mf2 <- km(~., design = design.grid, response = response.grid[, 2])
##' model <- list(mf1, mf2)
##' 
##' nsteps <- 3
##' lower <- rep(0, d)
##' upper <- rep(1, d)
##' 
##' # Optimization 1: EHI with pso
##' optimcontrol <- list(method = "pso", maxit = 20)
##' critcontrol <- list(refPoint = c(1, 10))
##' omEGO1 <- GParetoptim(model = model, fun = fname, crit = "EHI", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO1$par)
##' print(omEGO1$values)
##' 
##' \dontrun{
##' # Optimization 2: SMS with discrete search
##' optimcontrol <- list(method = "discrete", candidate.points = test.grid)
##' critcontrol <- list(refPoint = c(1, 10))
##' omEGO2 <- GParetoptim(model = model, fun = fname, crit = "SMS", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO2$par)
##' print(omEGO2$values)
##' 
##' # Optimization 3: SUR with genoud
##' optimcontrol <- list(method = "genoud", pop.size = 20, max.generations = 10)
##' critcontrol <- list(SURcontrol = list(distrib = "SUR", n.points = 100))
##' omEGO3 <- GParetoptim(model = model, fun = fname, crit = "SUR", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO3$par)
##' print(omEGO3$values)
##' 
##' # Optimization 4: EMI with pso
##' optimcontrol <- list(method = "pso", maxit = 20)
##' critcontrol <- list(nbsamp = 200)
##' omEGO4 <- GParetoptim(model = model, fun = fname, crit = "EMI", nsteps = nsteps,
##'                      lower = lower, upper = upper, optimcontrol = optimcontrol)
##' print(omEGO4$par)
##' print(omEGO4$values)
##' 
##' # graphics
##' sol.grid <- apply(expand.grid(seq(0, 1, length.out = 100),
##'                               seq(0, 1, length.out = 100)), 1, fname)
##' plot(t(sol.grid), pch = 20, col = rgb(0, 0, 0, 0.05), xlim = c(0, 1),
##'      ylim = c(-2, 10), xlab = expression(f[1]), ylab = expression(f[2]))
##' points(response.grid[,1], response.grid[,2], col = "black", pch = 20)     
##' points(omEGO1$values, pch = 17, col = "blue")
##' text(omEGO1$values[,1], omEGO1$values[,2], labels = 1:nsteps, pos = 3, col = "blue")
##' points(omEGO2$values, pch = 17, col = "green")
##' text(omEGO2$values[,1], omEGO2$values[,2], labels = 1:nsteps, pos = 3, col = "green")
##' points(omEGO3$values, pch = 17, col = "red")
##' text(omEGO3$values[,1], omEGO3$values[,2], labels = 1:nsteps, pos = 3, col = "red")
##' points(omEGO4$values, pch = 17, col = "orange")
##' text(omEGO4$values[,1], omEGO4$values[,2], labels = 1:nsteps, pos = 3, col = "orange")
##' legend("topright", c("EHI", "SMS", "SUR", "EMI"), col = c("blue", "green", "red", "orange"),
##'  pch = rep(17,4))
##' }
GParetoptim <- function (model, fun, cheapfun=NULL, crit="SMS", nsteps, lower, upper, type="UK", cov.reestim=TRUE,
                         critcontrol = NULL,
                         optimcontrol = list(method="genoud", threshold = 1e-5, distance = "euclidean", notrace = FALSE), ...){
  ##########################################################################################
  # Inputs :
  # model: list of 2 or 3 models
  # fun: objective function, returns 2 or 3 objectives
  # nsteps: number of iterations
  # lower, upper: design region
  # optimcontrol, type, CovReEstimate: parameters as in DiceOptim & KrigInv
  ##########################################################################################
  n     <- nrow(model[[1]]@X)
  d     <- model[[1]]@d
  
  fun <- match.fun(fun)
  
  # Build fastfun if necessary
  if (!is.null(cheapfun)) {
    cheapfun <- match.fun(cheapfun)
    fastobs <- apply(model[[1]]@X, 1, cheapfun)
    fastmod <- fastfun(fn = cheapfun, design = model[[1]]@X, response = fastobs)
    model[[length(model)+1]] <- fastmod
  }else{
    Y.new.cheap = NULL
  }
  n.obj <- length(model)
  
  if (length(model) < 2 ){
    cat("Error in model definition: 'model' must be a list of 2 or more km models \n")
    return(NULL)
  }
  
  if (length(model) >= 3 && crit=="EHI"){
    cat("Analytical Hypervolume EI only works with 2 objectives; SAA approximation used. \n")
  }
  if(length(model) > 3 && crit == "SUR"){
    cat("SUR is available for 2 or 3 objectives \n")
  }
  
  # Regroup all observations
  observations <- c()
  for (i in 1:n.obj) observations <- cbind(observations, model[[i]]@y)
  
  if(is.null(optimcontrol$notrace)){
    notrace <- FALSE
  }else{
    notrace <- optimcontrol$notrace
  }
  
  if(!notrace){
    cat("----------------------------\n")
    cat("Starting optimization with : \n The criterion", crit, "\n The solver",  optimcontrol$method, "\n")
    cat("----------------------------\n")
    cat("Ite / Crit / New x / New y \n")
  }
  
  #### Main loop starts here ################
  for (i in 1:nsteps) {
    
    ## Compute current Pareto front
    paretoFront <- t(nondominated_points(t(observations)))
    n.pareto    <- nrow(paretoFront)
    
    ## Change the seeds for genoud to avoid selecting always the same initial values
    if(optimcontrol$method == "genoud" & is.null(optimcontrol$unif.seed)){
      optimcontrol$unif.seed <- runif(1)
    }
    
    # observations removed (could be reintegrated)
    sol <- try(crit_optimizer(crit=crit, model=model, lower=lower, upper=upper, 
                              optimcontrol=optimcontrol, type=type, paretoFront=paretoFront, 
                              critcontrol=critcontrol, nsteps.remaining=nsteps-i))
    
    if (typeof(sol) == "character") {
      if(!notrace){
        cat("Unable to maximize criterion at iteration ", i, "- optimization stopped \n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE]
        par <- model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE]
      }
      return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
    }
    
    ## Update
    X.new <- matrix(as.numeric(sol$par), nrow=1, ncol=d)
    Y.new <- try(fun(as.numeric(sol$par), ...))
    if (!is.null(cheapfun)) {
      Y.new.cheap <- try(cheapfun(as.numeric(sol$par)))
    }
    
    if (typeof(Y.new) == "character" || (!is.null(cheapfun) && typeof(Y.new.cheap) == "character")) {
      if(!notrace){
        cat("Unable to compute objective function at iteration ", i, "- optimization stopped \n")
        cat("Problem occured for the design: ", X.new, "\n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE]
        par <- model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE]
      }
      return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
    }
    Y.new <- c(Y.new, Y.new.cheap)
    if(!notrace) cat( i, signif(sol$val,3), signif(X.new,3), signif(Y.new,3), "\n")
    
    # Remove new observation from integration points if discrete case is used
    if (optimcontrol$method=="discrete") {
      optimcontrol$candidate.points <- optimcontrol$candidate.points[-sol$index,,drop=FALSE]
      if (crit=="SUR") { 
        critcontrol$integration.points <- critcontrol$integration.points[-sol$index,,drop=FALSE]
      }
    }
    #     print(model[[1]]@covariance)
    # Update models
    observations <- rbind(observations, Y.new)
    newmodel <- model
    for (j in 1:n.obj) {
      newmodel[[j]] <- try(update(object = model[[j]], newX = X.new, newy=Y.new[j], newX.alreadyExist=FALSE,
                                  cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE))), silent = TRUE)
      if (typeof(newmodel[[j]]) == "character" && cov.reestim) {
        cat("Error in hyperparameter estimation - old hyperparameter values used instead for model ", j, "\n")
        newmodel[[j]] <- try(update(object = model[[j]], newX = X.new, newy=Y.new[j], newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
      }
      if (typeof(newmodel[[j]]) == "character") {
        cat("Unable to udpate kriging model at iteration", i-1, "- optimization stopped \n")
        cat("lastmodel is the model at iteration", i-1, "\n")
        cat("par and values contain the ",i, "th observation \n \n")
        if (i > 1) allX.new <- rbind(model[[1]]@X[(n+1):(n+i-1),, drop=FALSE], X.new)
        return(list(
          par    = allX.new,
          values = observations[(n+1):(n+i),, drop=FALSE],
          nsteps = i, 
          lastmodel = model))
      } else {
        model[[j]] <- newmodel[[j]]
      }
    }
    
  }
  
  if(!notrace) cat("\n")
  #### End of main loop ################
  return(list(
    par=model[[1]]@X[(n+1):(n+nsteps),, drop=FALSE], 
    values=observations[(n+1):(n+nsteps),, drop=FALSE], 
    nsteps=nsteps, 
    lastmodel=model))
}
