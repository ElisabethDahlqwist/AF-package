############## Common functions ##############
is.binary <- function(v) {
  x <- unique(v)
  if((length(x) - sum(is.na(x)) == 2L) & (max(x) == 1 & min(x) == 0)) TRUE
  else FALSE
}

############## AF function for cross-sectional sampling design #####################
#' @title Attributable fraction cross-sectional sampling design
#' @description \code{AF.cs} estimate the model-based adjusted attributable fraction for data from a cross-sectional sampling design.
#' @param formula an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model used for adjusting for confounders. The independent variables should be specified as the exposure and confounders. The dependent variable should be specified as the outcome of interest. The formula is used to fit a logistic regression by \code{\link{glm}}.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which the function is called.
#' @param exposure the exposure variable. Must be binary (0/1) where 0 is coded as unexposed.
#' @param clusterid for clustered data specify the variable in the data frame which is the cluster id in order to calculate cluster robust standard errors.
#' @return \item{AF.est}{estimated attributable fraction}
#' @return \item{AF.var}{estimated variance of the AF estimate (\code{AF.est}). The variance is estimated by the Sandwich estimator and the delta method}
#' @return \item{P}{estimated factual proportion of cases}
#' @return \item{P.var}{estimated variance of the estimate \code{P}}
#' @return \item{P0}{estimated counterfactual proportion of cases if exposure would be eliminated}
#' @return \item{P0.var}{estimated variance of the estimate \code{P0}}
#' @details \code{Af.cs} estimate the attributable fraction for binary outcomes
#' under the scenario of an elimination of the exposure \code{X} in the population. 
#' The estimate is adjusted for confounders by logistic regression (\code{\link{glm}}).
#'  Let the AF be defined as 
#' \deqn{AF = 1 - \frac{Pr(Y_0=1)}{Pr(Y = 1)}}{AF = 1 - Pr(Y0 = 1) / Pr(Y = 1)}
#' where \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} denote the counterfactual probability of the outcome if
#' the exposure \code{X} would have been eliminated from the population.
#' By adjusting the attributable fraction for confounders the counterfactual probablity
#'  \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} can be estimated. 
#'  \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} is denoted as \code{P0} in the 
#'  output and \eqn{Pr(Y = 1)} is denoted as \code{P}.
#' @author Elisabeth Dahlqwist, Arvid Sjolander
#' @seealso \code{\link{glm}}.
#' @references Bruzzi, P., Green, S. B., Byar, D., Brinton, L. A., and Schairer, C. (1985). Estimating the population attributable risk for multiple risk factors using case-control data. \emph{American Journal of Epidemiology} \bold{122}, 904-914.
#' @references Greenland, S. and Drescher, K. (1993). Maximum Likelihood Estimation of the Attributable Fraction from logistic Models. \emph{Biometrics} \bold{49}, 865-872.
#' @export
AF.cs<- function(formula, data, exposure, clusterid){
  #### Preparation of dataset ####
  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
      data <- data[complete, ]
   outcome <- as.character(terms(formula)[[2]])
  n <- nrow(data)
  n.cases <- sum(data[, outcome])
  if(missing(clusterid)) n.cluster <- 0
  else {
    n.cluster <- length(unique(data[, clusterid]))
         data <- data[order(data[, clusterid]), ]
  }
  ## Checks ##
  if(!is.binary(data[, outcome]))
    stop("Only binary outcome (0/1) is accepted.", call. = FALSE)
  if(!is.binary(data[, exposure]))
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  ## Counterfactual dataset ##
  data0 <- data
  data0[, exposure] <- 0
  #### Estimate model ####
   fit <- glm(formula = formula, family = binomial, data = data)
  npar <- length(fit$coef)
  ## Design matrices ##
   design <- model.matrix(object = delete.response(terms(fit)), data = data)
  design0 <- model.matrix(object = delete.response(terms(fit)), data = data0)
  #### Meat: score equations ####
  ## Score equation 1 ##
  score.P <- data[, outcome]
   pred.Y <- predict(fit, newdata = data, type = "response")
  ## Score equation 2 ##
    score.P0 <- predict(fit, newdata = data0, type = "response")
  ## Score equation 3 ##
  score.beta <- design * (score.P - pred.Y)
  ### Meat ###
  score.equations <- cbind(score.P, score.P0, score.beta)
  if (!missing(clusterid)){
    score.equations <- aggregate(score.equations, list(data[, clusterid]), sum)[, - 1]
  }
  meat <- var(score.equations, na.rm = TRUE)
  #### Bread: hessian of score equations ####
  ## Hessian of score equation 1 ##
  hessian.P <- matrix(c(- 1, 0, rep(0,npar)), nrow = 1, ncol = 2 + npar)
  ## Hessian of score equation 2 ##
  g <- family(fit)$mu.eta
    dmu.deta <- g(predict(object = fit, newdata = data0))
  deta.dbeta <- design0
   dmu.dbeta <- dmu.deta * deta.dbeta
  hessian.P0 <- matrix(c(0, - 1, colMeans(dmu.dbeta)), nrow = 1, ncol = 2 + npar)
  ## Hessian of score equation 3 ##
  hessian.beta <- cbind(matrix(rep(0, npar * 2), nrow = npar, ncol = 2)
                        , - solve(vcov(object = fit)) / n)
  ### Bread ###
  bread <- rbind(hessian.P, hessian.P0, hessian.beta)
  #### Sandwich ####
  if (!missing(clusterid))
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster / n^2 ) [1:2, 1:2]
  else
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
  #### Point estimate of AF ####
       P <- mean(score.P, na.rm = TRUE)
      P0 <- mean(score.P0, na.rm = TRUE)
  AF.est <- 1 - P0 / P
  ## Delta method for variance estimate ##
  gradient <- as.matrix(c(P0 / P ^ 2, - 1 / P), nrow = 2, ncol = 1)
    AF.var <- t(gradient) %*% sandwich %*% gradient
     P.var <- sandwich[1, 1]
    P0.var <- sandwich[2, 2]
  #### Output ####
  out <- c(list(AF.est = AF.est, AF.var = AF.var, P = P, P0 = P0, P.var = P.var,
                P0.var = P0.var, call = fit$call, exposure = exposure, outcome = outcome,
                fit = fit, sandwich = sandwich, gradient = gradient, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster))
  class(out) <- "AF"
  return(out)
}

############## AF function for cohort time-to-event outcomes #####################
library(survival)
#' @title Attributable fraction function from cohort sampling design with time-to-event outcomes.
#' @description \code{AF.ch} estimate the model-based adjusted attributable fraction function for data from a cohort sampling design with time-to-event outcomes.
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the Surv function (\code{\link{Surv}}). A symbolic description of the model used for adjusting for confounders. The independent variables should be specified as the exposure and confounders. The dependent variable should be specified as the outcome of interest. The formula is used to fit a Cox Proportional Hazards model in the survival package. For details see documentation on coxph (\code{\link{coxph}}).
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which the function is called.
#' @param exposure the exposure variable. Must be binary (0/1) where 0 is coded as unexposed.
#' @param ties a character string specifying the method for tie handling. If there are no tied death times all the methods are equivalent. Uses the Breslow method by default.
#' @param time.sequence a scalar or vector of time points specified by the user for which the Attributable fraction function is estimated. If not specified the death times will be used.
#' @param clusterid for clustered data specify variable name for the cluster id in order to calculate cluster robust standard errors.
#' @return \item{AF.est}{estimated attributable fraction function for every time point specified by \code{time.sequence}}
#' @return \item{AF.var}{estimated variance of the AF estimate (\code{AF.est}) for every time point specified by \code{time.sequence}. The variance is estimated by the Sandwich estimator and the delta method}
#' @return \item{St.est}{estimated factual survival function}
#' @return \item{St.var}{estimated variance of the estimate \code{St.est}}
#' @return \item{St0.est}{estimated counterfactual survival function if exposure would be eliminated}
#' @return \item{St0.var}{estimated variance of the estimate \code{St.est}}
#' @details \code{Af.ch} estimate the attributable fraction for time-to-event outcomes
#' under the scenario of an elimination of the exposure \code{X} in the population. The estimate is adjusted for confounders
#' by the Cox Proportional Hazards model (\code{\link{coxph}}). Let the Attributable fraction function be defined as 
#' \deqn{AF = 1 -\frac{( 1 - St_0(t) )}{( 1 - St(t) )}}{AF = 1 - ( 1 - St0(t) ) / ( 1 - St(t) )}
#' where \eqn{St_0(t)}{St0(t)} denote the counterfactual survival function for the event if
#' the exposure would have been eliminated from the population at baseline and \eqn{St(t)} denote the factual survival function. 
#' \eqn{St_0(t)}{St0(t)} is estimated by adjusting for confounders at baseline using the Cox Proportional Hazards model.
#' \eqn{St_0(t)}{St0(t)} is denoted as \code{St0.est} in the output and \eqn{St(t)}
#' is denoted as \code{St.est}.
#' @author Elisabeth Dahlqwist, Arvid Sjolander
#' @seealso \code{\link{coxph}} and \code{\link{Surv}}.
#' @references Chen, L., Lin, D. Y., and Zeng, D. (2010). Attributable fraction functions for censored event times. \emph{Biometrika} \bold{97}, 713-726.
#' @references Sjolander, A. and Vansteelandt, S. (2014). Doubly robust estimation of attributable fractions in survival analysis. \emph{Statistical Methods in Medical Research}.
#' @export
AF.ch<- function(formula, data, exposure, ties="breslow",
                   time.sequence, clusterid){
  #### Preparation of dataset ####
  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
         m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
      data <- data[complete, ]
  ## If time.sequence is missing ##
  if(missing(time.sequence)){time.sequence <- fit.detail$time}
  ## Checks ##
  if(!is.binary(data[, exposure]))
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  if(missing(clusterid)) n.cluster <- 0
  else n.cluster <- length(unique(data[, clusterid]))
  ## Find names of end variable and event variable
  rr <- rownames(attr(terms(formula), "factors"))[1]
  temp <- gregexpr(", ", rr)[[1]]
  if(length(temp == 1)){
    endvar <- substr(rr, 6, temp[1] - 1)
    eventvar <- substr(rr, temp[1] + 2, nchar(rr) - 1)
  }
  if(length(temp) == 2){
    endvar <- substr(rr, temp[1] + 2, temp[2] - 1)
    eventvar <- substr(rr, temp[2] + 2, nchar(rr) - 1)
  }
  n <- nrow(data)
  n.cases <- sum(data[, eventvar])
  # Sort on "end-variable"
  data <- data[order(data[, endvar]), ]
  # Create dataset data0 for counterfactual X=0
  data0 <- data
  data0[, exposure] <- 0
  #### Estimate model ####
  ## Fit a Cox PH model ##
  environment(formula) <- new.env()
   fit <- coxph(formula = formula, data = data, ties = "breslow")
  npar <- length(fit$coef)
  fit.detail <- coxph.detail(object = fit)
  ## Design matrices ##
   design <- model.matrix(object = delete.response(terms(fit)), data = data)[, -1]
  design0 <- model.matrix(object = delete.response(terms(fit)), data = data0)[, -1]
  ### Estimate the survival functions ###
  ## Hazard increment ##
  dH0 <- fit.detail$hazard
   H0 <- cumsum(dH0)
  ## Baseline hazard function ##
      H0step <- stepfun(fit.detail$time, c(0, H0))
       H0res <- rep(0, n)
  dH0.untied <- rep(dH0, fit.detail$nevent) / rep(fit.detail$nevent, fit.detail$nevent)
  H0res[data[, eventvar] == 1] <- dH0.untied * n #handle ties
  #H0res[data[, eventvar] == 1] <- dH0 * n
  ## Predict based on the Cox PH model ##
   epred <- predict(object = fit, newdata = data, type = "risk")
  epred0 <- predict(object = fit, newdata = data0, type = "risk")
  ### Meat ###
  ## Score equation 4 ## for the Cox PH model (made outside of loop)
  score.beta <- residuals(object = fit, type = "score")
  ## Weighted mean of the variable at event for all at risk at that time ##
  E <- matrix(0, nrow = n, ncol = npar)
  means <- fit.detail$means
  means <- means[rep(1:nrow(means), fit.detail$nevent), ] #handle ties
  E[data[, eventvar] == 1, ] <- means
  #E[data[, eventvar] == 1, ] <- fit.detail$means
  ## One point and variance estimate for each time t in time.sequence ##
   St.est <- vector(length = length(time.sequence))
  St0.est <- vector(length = length(time.sequence))
   AF.var <- vector(length = length(time.sequence))
   St.var <- vector(length = length(time.sequence))
  St0.var <- vector(length = length(time.sequence))

  # Loop over all t in time.sequence
  for (i in 1:length(time.sequence)){
    t <- time.sequence[i]
    #### Meat: score equations ####
    ## Score equation 1 ## for the factual survival function
    score.S <- exp( - H0step(t) * epred)
    ## Score equation 2 ## for the counterfactual survival function
    score.S0 <- exp( - H0step(t) * epred0)
    ## Score equation 3 ##  for the Breslow estimator
    score.H0 <- H0res * (data[, endvar] <= t)
    ## Score equation 4 ## for the Cox PH model (made outside of loop)
    ### Meat ###
    score.equations <- cbind(score.S, score.S0, score.H0, score.beta)
    if (!missing(clusterid)){
      score.equations <- aggregate(score.equations, by = list(data[, clusterid]), sum)[, - 1]
    }
    meat <- var(score.equations, na.rm = TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
     hessian.S <- c(-1, 0, mean(epred * score.S), colMeans(design * H0step(t) * epred * score.S))
    ## Hessian of score equation 2 ##
    hessian.S0 <- c(0, -1, mean(epred0 * score.S0), colMeans(design0 * H0step(t) * epred0 * score.S0))
    ## Hessian of score equation 3 ##
    hessian.H0 <- c(rep(0,2), - 1, - colMeans(E * score.H0, na.rm = TRUE))
    ## Hessian of score equation 4 ##
    hessian.beta <- cbind(matrix(0, nrow = npar, ncol = 3), - solve(vcov(object = fit)) / n)
    ### Bread ###
    bread<-rbind(hessian.S, hessian.S0, hessian.H0, hessian.beta)
    ### Sandwich ###
    if (!missing(clusterid))
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster/ n^2 ) [1:2, 1:2]
    else
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    #### For point estimate ####
     St.est[i] <- mean(x = score.S, na.rm = TRUE)
    St0.est[i] <- mean(x = score.S0, na.rm = TRUE)
    #### Estimate of variance using the delta method ####
    gradient <- as.matrix(c( - (1 - St0.est[i]) / (1 - St.est[i]) ^ 2, 1 / (1 - St.est[i]))
                          , nrow = 2, ncol = 1)
     AF.var[i] <- t(gradient) %*% sandwich %*% gradient
     St.var[i] <- sandwich[1, 1]
    St0.var[i] <- sandwich[2, 2]
  }
  ### The AF function estimate ###
  AF.est <- 1 - (1 - St0.est) / (1 - St.est)
  #### Output ####
  out <- c(list(AF.est = AF.est, AF.var = AF.var, St.est = St.est,
                St0.est = St0.est, St.var = St.var, St0.var = St0.var,
                call = fit$call, exposure = exposure, outcome = eventvar, fit = fit,
                sandwich = sandwich, gradient = gradient, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster,  time.sequence = time.sequence))
  class(out) <- "AF"
  return(out)
}

############## AF function step 3 and 4 #####################
library(survival)
library(drgee)
#' @title Attributable fraction mached or non-matched case-control sampling design.
#' @description \code{AF.cc} estimate the attributable fraction for data from a mached or non-matched case-control sampling design.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model used for adjusting for confounders. The independent variables should be specified as the exposure and confounders. The dependent variable should be specified as the outcome of interest. The formula is used to fit a logistic regression by \code{\link{glm}} for non-matched case-control and conditional logistic regression by \code{\link{gee}} (in package \code{\link{drgee}}) for matched case-control.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment (formula), typically the environment from which the function is called.
#' @param exposure the exposure variable. Must be binary (0/1) where 0 is coded as unexposed.
#' @param sampling.design a string which specify if the sampling is matched or non-matched case-control. Default setting is non-matched.
#' @return \item{AF.est}{estimated attributable fraction}
#' @return \item{AF.var}{estimated variance of the AF estimate (\code{AF.est}). The variance is estimated by the Sandwich estimator and the delta method}
#' @return \item{log.or}{estimated log odds ratio adjusted for confounders. The sum of all parameters depending on the exposure variable \code{X} as specified by the user in the formula. If the model to be estimated is
#'  \deqn{y = \beta{x}+\gamma{z}}{y = \beta x + \gamma z} the \code{log.odds} is the estimate of \eqn{\beta}.
#'   If the model to be estimated is \deqn{y=\beta{x}+\gamma{z} +\psi{xz}}{y = \beta x +\gamma z +\psi xz} the \code{log.odds} is the estimate of \eqn{\beta + \psi}.}
#' @details \code{Af.cc} estimate the attributable fraction for binary outcomes
#' under the scenario of an elimination of the exposure \code{X} in the population. 
#' The estimate is adjusted for confounders by logistic regression for unmatched case-control (\code{\link{glm}}) and conditional logistic regression for matched case-control (\code{\link{gee}}). 
#' The estimation is an approximation based on the assumption that the outcome is rare so that the relative risk can be approximated by the odds ratio, for details see references \link{Bruzzi}.
#' Let the AF be defined as 
#' \deqn{AF = 1 - \frac{Pr(Y_0=1)}{Pr(Y = 1)}}{AF = 1 - Pr(Y0 = 1) / Pr(Y = 1)}
#' where \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} denote the counterfactual probability of the outcome if
#' the exposure \code{X} would have been eliminated from the population. Under the assumption that Z
#' is the only confounder the probability \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} can be expressed as
#' \deqn{Pr(Y_0=1)=E_z{Pr(Y=1\mid{X}=0,Z)}}{Pr(Y0=1) = E_z{Pr(Y = 1 | X = 0, Z)}}
#' Using Bayes' theorem this imply that the AF can be defined as the equality
#' \deqn{AF = 1-\frac{E_z\{Pr(Y=1\mid X=0,Z)\}}{Pr(Y=1)}=1-E{RR^{-X}(Z)|Y = 1}}{
#' AF = 1 - E_z{Pr( Y = 1 | X = 0, Z)} / Pr(Y = 1) = 1 - E{RR^{-X} (Z) | Y = 1}} 
#' where \eqn{E{RR^{-X}(Z) | Y = 1}} can be estimated by the relative risk.
#' Moreover, the relative risk can be approximated by the odds ratio if the outcome is rare.
#' The odds ratio is estimated by logistic regression or conditional logistic regression. Thus, 
#' \deqn{1 - E{RR^{-X}(Z) | Y = 1} \approx 1 - E{OR^{-X}(Z) | Y = 1}}{1 - E{RR^{-X}(Z) | Y = 1} is approximately 1 - E{OR^{-X}(Z) | Y = 1}}
#' where \eqn{E{OR^{-X}(Z) | Y = 1}} is the sum of the parameters depending on the exposure X.
#' @author Elisabeth Dahlqwist, Arvid Sjolander
#' @seealso \code{\link{glm}} and \code{\link{gee}} used for estimating the conditional logistic regression for mached case-control.
#' @references Bruzzi, P., Green, S. B., Byar, D., Brinton, L. A., and Schairer, C. (1985). Estimating the population attributable risk for multiple risk factors using case-control data. \emph{American Journal of Epidemiology} \bold{122}, 904-914.
#' @export
AF.cc<-function(formula, data, exposure, clusterid,
                  sampling.design = "non.matched"){
  #### Preparation of dataset ####
  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
         m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
      data <- data[complete, ]
   outcome <- as.character(terms(formula)[[2]])
  ## Checks ##
  if(is.binary(data[, outcome]) == FALSE)
    stop("Only binary outcome (0/1) is accepted.", call. = FALSE)
  if(is.binary(data[, exposure]) == FALSE)
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(!as.numeric(summary(all.vars(formula[[3]]) == exposure)[3]) == 1)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  if(missing(clusterid)) n.cluster <- 0
  else n.cluster <- length(unique(data[, clusterid]))
  #### Methods for non-matched or matched sampling designs ####
  #ARVID: I inserted these lines
  n <- nrow(data)
  n.cases <- sum(data[, outcome])
  if (!missing(clusterid))
   data <- data[order(data[,clusterid]), ]
  data0 <- data
  data0[, exposure] <- 0
  #### Estimate model ####
  if(sampling.design == "non.matched")
    fit <- glm(formula = formula, family = binomial, data = data)
  if(sampling.design == "matched")
    fit <- gee(formula, link = "logit", data, cond = TRUE, clusterid = clusterid)
  npar <- length(fit$coef)
  ## Design matrices ##
  if(sampling.design == "non.matched"){
     design <- model.matrix(object = delete.response(terms(fit)), data = data)
    design0 <- model.matrix(object = delete.response(terms(fit)), data = data0)
  }
  if(sampling.design == "matched"){
     design <- model.matrix(object = formula, data = data )[, - 1]
    design0 <- model.matrix(object = formula, data = data0)[, - 1]
  }
  ## Create linear predictors to estimate the log odds ratio ##
       diff.design <- design0 - design
   linearpredictor <- design  %*% coef(fit)
  linearpredictor0 <- design0 %*% coef(fit)
  #log odds ratio#
  log.or <- linearpredictor - linearpredictor0
  ## Estimate approximate AF ##
  AF.est <- 1 - sum(data[, outcome] * exp( - log.or)) / sum(data[, outcome])
  #### Meat: score equations ####
  ## Score equation 1 ## individual estimating equations of the estimate of AF
  score.AF <- data[, outcome] * (exp( - log.or) - AF.est)
  ## Score equation 2 ## individual estimating equations from conditional logistic reg.
  if(sampling.design == "non.matched")
    pred.diff <- data[, outcome] - predict(fit, newdata = data, type = "response")
  if(sampling.design == "matched")
    pred.diff <- fit$res

  #ARVID: if missing data, then fit$res and design have different dimensions!
  #Must talk to Johan about this!!
  score.beta <- design * pred.diff
  score.equations <- cbind(score.AF, score.beta)
  if (!missing(clusterid))
    score.equations <- aggregate(score.equations, list(data[, clusterid]), sum)[, - 1]
  meat <- var(score.equations, na.rm=TRUE)
  #### Bread: hessian of score equations ####
  ## Hessian of score equation 1 ##
  #### Estimating variance using Sandwich estimator ####
  hessian.AF1 <- - data[, outcome]
  hessian.AF2 <- (design0 - design) * as.vector(data[, outcome] * exp( - log.or))
  if (!missing(clusterid))
    hessian.AF <- cbind(mean(aggregate(hessian.AF1, list(data[, clusterid]), sum)[, - 1], na.rm=TRUE)
                        , t(colMeans(aggregate(hessian.AF2
                        , list(data[, clusterid]), sum)[, - 1], na.rm = TRUE)))
  else
    hessian.AF <- cbind(mean(hessian.AF1), t(colMeans(hessian.AF2, na.rm = TRUE)))
  hessian.beta <- cbind(matrix(rep(0, npar), nrow = npar, ncol = 1), - solve(vcov(object = fit)) / n)
  ### Bread ###
  bread <- rbind(hessian.AF, hessian.beta)
  #### Sandwich ####
  if (!missing(clusterid))
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster/ n ^ 2 ) [1:2, 1:2]
  else
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]

  if(sampling.design == "matched"){
    if (missing(clusterid))
      stop("Need to specify clusterid.", call. = FALSE)
    n.cluster <- length(unique(data[, clusterid]))
         data <- data[order(data[,clusterid]), ]
            n <- nrow(data)
      n.cases <- sum(data[, outcome])
        data0 <- data
    data0[, exposure] <- 0
    #### Estimate model ####
     fit <- gee(formula, link = "logit", data, cond = TRUE, clusterid = clusterid)
    npar <- length(fit$coef)
    ## Design matrices ##
     design <- model.matrix(object = formula, data = data)[, - 1]
    design0 <- model.matrix(object = formula, data = data0)[, - 1]
    ## Create linear predictors to estimate the log odds ratio ##
         diff.design <- design0 - design
     linearpredictor <- design %*% coef(fit)
    linearpredictor0 <- design0 %*% coef(fit)
    #log odds ratio#
    log.or <- linearpredictor - linearpredictor0
    ## Estimate approximate AF ##
    AF.est <- 1 - sum(data[, outcome] * exp( - log.or)) / sum(data[, outcome])
    #### Meat: score equations ####
    ## Score equation 1 ## individual estimating equations of the estimate of AF
    score.AF <- data[, outcome] * (exp( - log.or) - AF.est)
    ## Score equation 2 ## individual estimating equations from conditional logistic reg.

    #ARVID: if missing data, then fit$res and design have different dimensions!
    #Must talk to Johan about this!!

    score.beta <- fit$res * design
    ### Meat ###
    score.equations <- cbind(score.AF, score.beta)
    score.equations <- aggregate(score.equations, list(data[, clusterid]), sum)[, - 1]
    meat <- var(score.equations, na.rm=TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    #### Estimating variance using Sandwich estimator ####
    hessian.AF1 <- - data[, outcome]
    hessian.AF2 <- diff.design * as.vector(data[, outcome] * exp(-log.or))
    hessian.AF <- cbind(mean(aggregate(hessian.AF1, list(data[, clusterid]), sum)[, -1], na.rm=TRUE)
                        , t(colMeans(aggregate(hessian.AF2
                        , list(data[, clusterid]), sum)[, -1], na.rm = TRUE)))
    ## Hessian of score equation 2 ##

    #ARVID: if missing data, then design and fit$d.res have different dimensions!
    #Must talk to Johan about this!!
    #Can we simply replace crossprod(design, fit$d.res) / n with - solve(vcov(object = fit)) / n?
    hessian.beta <- cbind(matrix(rep(0,npar), nrow = npar, ncol = 1), crossprod(design, fit$d.res) / n)

    ### Bread ###
    bread <- rbind(hessian.AF, hessian.beta)
    #### Sandwich ####
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster / n ^ 2)
  }
  AF.var <- sandwich[1, 1]
  #### Output ####
  out <- c(list(AF.est = AF.est, AF.var = AF.var, log.or = log.or,
                call = fit$call, exposure = exposure, outcome = outcome, fit = fit,
                sandwich = sandwich, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster))
   class(out) <- "AF"
  return(out)
}

############## Summary and print functions ##############
print.AF<-function(object, digits = max(3L, getOption("digits") - 3L), ...){
  if(!object$n.cluster == 0) {
    Std.Error <- "Robust SE"
    se <- "cluster-robust standard error"
  }
  else {
    Std.Error <- "Std.Error"
    se <- "standard error"
  }
  cat("\nEstimated attributable fraction (AF) and", se, ":", "\n")
  cat("\n")
  table.est <- cbind(object$AF.est, sqrt(object$AF.var))
  colnames(table.est) <- c("AF", Std.Error)
  r <- rep("", , length(object$AF.est))
  rownames(table.est) <- c(r)
  modelcall <- as.character(object$fit$call[1])
  if(modelcall == "coxph") {
    table.est <- cbind(object$time.sequence, table.est)
    colnames(table.est) <- c("Time", "AF", Std.Error)
    print.default(table.est)
  }
  else {
    print.default(table.est)
  }
}

CI <- function(AF, Std.Error, confidence.level, CI.transform){
  if(CI.transform == "untransformed"){
    lower <- AF - abs(qnorm((1 - confidence.level) / 2)) * Std.Error
    upper <- AF + abs(qnorm((1 - confidence.level) / 2)) * Std.Error
  }
  if(CI.transform == "log"){
    lower <- AF * exp( - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / AF)
    upper <- AF * exp(abs(qnorm((1 - confidence.level) / 2)) * Std.Error / AF)
  }
  if(CI.transform == "logit"){
    logit <- function(x) log(x / (1 - x))
    lower <- exp(logit(AF) - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))) / (1 + exp(logit(AF) - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))))
    upper <- exp(logit(AF) + abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))) / (1 + exp(logit(AF) + abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))))
  }
  CI <- cbind(lower, upper)
  return(CI)
}
#' @title Summary function for objects of class "\code{AF}".
#' @description Gives a summary of the AF estimate(s) including zvalue, p-value and confidence interval(s).
#' @param object an object of class \code{AF} from \code{\link{AF.cs}}, \code{\link{AF.ch}} or \code{\link{AF.cc}} functions.
#' @param confidence.level user-specified confidence level for the confidence intervals. If not specified the default is 95 percent. Should be specified in decimals such as 0.95 for 95 percent.
#' @param CI.transform user-specified transformation of the Wald confidence interval(s). Options are \code{untransformed}, \code{log} and \code{logit}. If not specified untransformed will be calculated.
#' @author Elisabeth Dahlqwist, Arvid Sjolander
#' @export
summary.AF <- function(object, digits = max(3L, getOption("digits") - 3L),
                       confidence.level, CI.transform, ...){
  if(missing(confidence.level)) confidence.level <- 0.95
  if(missing(CI.transform)) CI.transform <- "untransformed"
  se <- sqrt(object$AF.var)
  zvalue <- object$AF.est / sqrt(object$AF.var)
  pvalue <- 2 * pnorm( - abs(zvalue))
  confidence.interval <- CI(AF = object$AF.est, Std.Error = sqrt(object$AF.var),
                            confidence.level = confidence.level,
                            CI.transform = CI.transform)
  colnames(confidence.interval) <- c("Lower limit", "Upper limit")

  if(!object$n.cluster == 0) Std.Error <- "Robust SE"
  else Std.Error <- "Std.Error"
  AF <- cbind(object$AF.est, se, zvalue, pvalue)
  colnames(AF) <- c("AF estimate", Std.Error, "z value", "Pr(>|z|)")

  modelcall <- as.character(object$fit$call[1])
  if(modelcall == "glm") method = "Logistic regression"
  if(modelcall == "coxph") method = "Cox Proportional Hazard model"
  if(modelcall == "gee") method = "Conditional logistic regression"

  fit <- summary(object$fit)

  if(modelcall == "coxph"){
    ans <- list(AF = AF, time.sequence = object$time.sequence,
                CI.transform = CI.transform, confidence.level = confidence.level,
                confidence.interval = confidence.interval, n.obs = object$n,
                n.cases = object$n.cases, n.cluster = object$n.cluster,
                modelcall = modelcall, method = method, formula = object$formula,
                exposure = object$exposure, outcome = object$outcome, fit = fit,
                sandwich = object$sandwich)
  }
  else{
    ans <- list(AF = AF, CI.transform = CI.transform, confidence.level = confidence.level,
                confidence.interval = confidence.interval, n.obs = object$n, n.cases = object$n.cases,
                n.cluster = object$n.cluster, modelcall = modelcall, method = method,
                formula = object$formula, exposure = object$exposure, outcome = object$outcome,
                fit = fit, sandwich = object$sandwich)
  }
  class(ans) <- "summary.AF"
  return(ans)
}

#ARVID: in print.summary, change this:
#delete blankspace before "%"

print.summary.AF <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  if(!x$n.cluster == 0) Std.Error <- "Robust SE"
  else Std.Error <- "Std.Error"
  if(x$CI.transform == "log") x$CI.transform <- "log transformed"
  if(x$CI.transform == "logit") x$CI.transform <- "logit transformed"
  cat("\nEstimated attributable fraction (AF) and", x$CI.transform, x$confidence.level * 100, "%",  "Wald CI:", "\n")
  cat("\n")
  table.est <- cbind(x$AF, x$confidence.interval)
  colnames(table.est) <- c("AF", Std.Error, "z value", "Pr(>|z|)",
                           "Lower limit", "Upper limit")
  r <- rep("", , nrow(x$AF))
  rownames(table.est) <- c(r)
  modelcall <- as.character(x$fit$call[1])
  if(x$modelcall == "coxph") {
    table.est <- cbind(x$time.sequence, table.est)
    colnames(table.est) <- c("Time", "AF", Std.Error, "z value", "Pr(>|z|)",
                             "Lower limit", "Upper limit")
    print.default(table.est)
  }
  else {
    print.default(table.est)
  }

  #cat("\n")
  #table.exposure.outcome <- cbind(x$exposure, x$outcome)
  #rownames(table.exposure.outcome) <- c("")
  #colnames(table.exposure.outcome) <- c("Exposure", "Outcome")
  #print.default(table.exposure.outcome)
  cat("\nExposure", ":", x$exposure, "\n")

  if(x$modelcall == "coxph") outcome <- "Event   "
  else outcome <- "Outcome "
  #cat("\n")
  cat(outcome, ":", x$outcome, "\n")

  cat("\n")
  table.nr <- cbind(x$n.obs, x$n.cases)
  rownames(table.nr) <- c("")
  if(x$modelcall == "coxph") number <- "Events"
  else number <- "Cases"
  colnames(table.nr) <- c("Observations", number)
  if (x$n.cluster == 0) print.default(table.nr)
  else{
    table.nr.cluster <- cbind(table.nr, x$n.cluster)
    colnames(table.nr.cluster) <- c("Observations", number, "Clusters")
    print.default(table.nr.cluster)
  }
  cat("\nMethod for confounder adjustment: ", x$method, "\n")
  formula <- as.character(x$formula)
  cat("\nFormula: ", formula[2], formula[1], formula[3], "\n")
}

#ARVID: in plot.AF, change this:
#delete blankspace before "%" in legend

#' @title Plot function for objects from the function \code{\link{AF.ch}}.
#' @description Creates a simple scatterplot for the Attributable fraction function with time sequence (defined by the user in the \code{\link{AF.ch}} function) on the x-axis and the Attributable fraction function estimate on the y-axis. 
#' @param object an object of class \code{AF} from the \code{\link{AF.ch}} function.
#' @param confidence.level user-specified confidence level for the confidence intervals. If not specified the default is 95 percent. Should be specified in decimals such as 0.95 for 95 percent.
#' @param CI.transform user-specified transformation of the Wald confidence interval(s). Options are \code{untransformed}, \code{log} and \code{logit}. If not specified untransformed will be calculated.
#' @param xlab label on the x-axis. If not specified \emph{"Time"} will be displayed.
#' @param main main title of the plot. If not specified  \emph{"Estimate of the attributable fraction function"} will be displayed.
#' @author Elisabeth Dahlqwist, Arvid Sjolander
#' @export
plot.AF <- function(object, digits = max(3L, getOption("digits") - 3L),
                    CI = FALSE, confidence.level, CI.transform, xlab, main, ...){
  modelcall <- as.character(object$fit$call[1])
  if(!modelcall == "coxph")
    stop("Plot function is only available for the attributable fraction function. That is objects from the AF.ch function", call. = FALSE)
  if(missing(confidence.level)) confidence.level <- 0.95
  if(missing(CI.transform)) CI.transform <- "untransformed"
  if(missing(xlab)) xlab <- "Time"
  if(missing(main)) main <- "Estimate of the attributable fraction function"
  plot(object$time.sequence, object$AF.est, main = main,
         ylab = "Attributable fraction function" , xlab = xlab, pch = 19, lty = 1, type = "o")

  if(CI == TRUE){
    confidence.interval <- CI(AF = object$AF.est, Std.Error = sqrt(object$AF.var),
                              confidence.level = confidence.level,
                              CI.transform = CI.transform)
  lines( object$time.sequence, confidence.interval[, 2], lty = 2)
  lines( object$time.sequence, confidence.interval[, 1], lty = 2)
  level <- confidence.level * 100
  CI <- sprintf("%i %s", level,  "% Conf. Interval")
  if(CI.transform == "log") transform <- "(log transformed)"
  if(CI.transform == "logit") transform <- "(logit transformed)"
  if(CI.transform == "untransformed") transform <- "(untransformed)"
  legend("topright", legend = c("AF estimate", CI, transform), pch = c(19, NA, NA),
         lty = c(1, 2, 0), bty = "n")
  }
}

#plot(est2, CI=T, confidence.level=0.97, time.sequence = c(1,2), CI.transform="logit",
#     xlab= "Time since diagnosis (months)", main = "AF function cancer survival")

#plot(est1)
#plot(est2, CI=T, confidence.level=0.97, CI.transform="logit")


#require(ggplot2)
#ggplot(df, aes(Time = Time, y = AF)) +
 # geom_point(size = 4) +
 # geom_errorbar(aes(ymax = U, ymin = L))
