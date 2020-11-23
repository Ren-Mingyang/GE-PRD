############################################################
# This document includes main functions of proposed methods 
# for joint and marginal analysis,
# which are used to support numerical simulation studies
# and real data analysis in the paper
# "Gene-environment Interaction Identification via Penalized Robust Divergence"
############################################################


#########################################################################################
# Functions for joint analysis:
# This section includes two proposed methods: 
# gamma divergence, density power divergence.
#########################################################################################


#---------------gamma divergence------------
BIC_sgMCP_gamma <- function(y, X, W, gamma = 0.5, b1, 
                             xi = 3, max_iter = 20, family, 
                             Matla, eps = 0.001, criterion="AIC") {
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: BIC_sgMCP_gamma
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Performing joint analysis based on gamma divergence with the sparse group penalty.
  ##            The regression coefficients under the given tuning parameters can be calculated, 
  ##            and the optimal regression coefficient is selected under the given criterion such as BIC.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: gammagaussian_loss()   gammabinomial_loss()
  ##            Rcpp functions: gammalinearsgMCP()  gammalogisticsgMCP()
  ##            R packages: Rcpp
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable
  ## @ X: n * q matrix, the design matrix corresponding environmental variables.
  ## @ W: n * ((q+1)*p) matrix, the design matrix corresponding gene variables and its interactive effects.
  ## @ gamma: A vector or a float value, the robust parameter in gamma divergence, the default setting is 0.5.
  ## @ b1: ((q+1)*(p+1)) * 1 vector, initial value of the regression coefficient (including intercept terms).
  ## @ xi: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ max_iter: int, Maximum number of cycles of the algorithm, the default setting is 20.
  ## @ family: types of response variables, which can be selected from "gaussian" and "binomial".
  ## @ Matla: s * 2 matrix, the tuning parameters, 
  ##         the s rows represent s choices of the tuning parameters,
  ##         the 2 columns represent lambda1, lambda2,
  ##         (in this paper, to simplify, lambda1 = lambda2 is set, but this function is written in such a way that the two can be not equal.)
  ## @ eps: a float value, algorithm termination threshold.
  ## @ criterion: the parameter selection criterion, which can be selected from "AIC" and "BIC".
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ beta: ( (q+1)*(p+1)-1 ) * 1 vector, the optimal regression coefficients selected according to set criterion (without intercept term).
  ## @ beta0: a float value, the intercept term in optimal regression coefficients.
  ## @ beta_all: ( (q+1)*(p+1)-1 ) * s vector, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ beta0_all: 1 * s vector, all intercept terms in regression coefficients corresponding s choices of given tuning parameters.
  ## @ loss: s * 1 vector, the values of the loss function (without penalty function) corresponding s choices of given tuning parameters.
  ## @ nonzero: s * 1 vector, the numbers of non-zero regression coefficients corresponding s choices of given tuning parameters.
  ## @ iter_all: Number of iterations.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  n <- length(y)
  p <- length(X[1,])+length(W[1,])
  gg <- length(gamma)
  L <- nrow(Matla)
  loss <- rep(0,L*gg)
  nonzero <- rep(0,L*gg)
  beta_all <- matrix(0,length(X[1,])+length(W[1,]),L*gg)
  beta0_all <- matrix(0,1,L*gg)
  iter_all <- rep(0,L*gg)
  Sigma <- rep(1,L*gg)
  for (l in 1:L) {
    matla <- Matla[l,]
    for (gi in 1:gg) {
      gammai <- gamma[gi]
      if (family == "gaussian"){
        ll <- 1
        repeat { 
          res = gammalinearsgMCP(y, X, W, gammai, b1, matla[1], matla[2], xi, max_iter, eps, 1)
          b <- res$beta
          ll <- ll + 1
          if(length(b[is.na(b)]) == 0 | ll>15) {
            break
          }
        }
        b0 <- res$beta0
        b[is.na(b)] <- 0
        b0[is.na(b0)] <- 0
        iter_all[(l-1)*gg+gi] <- res$iter
        beta_all[,(l-1)*gg+gi] <- b
        beta0_all[1,(l-1)*gg+gi] <- b0
        sigma <- res$sigma
        loss[(l-1)*gg+gi] <- gammagaussian_loss(y, cbind(X, W), gammai, b, b0, sigma)
        nonzero[(l-1)*gg+gi] <- length(which(b != 0))
        Sigma[(l-1)*gg+gi] <- sigma
      }else if (family == "binomial"){
        ll <- 1
        repeat { 
          res = gammalogisticsgMCP(y, X, W, gammai, b1, matla[1], matla[2], xi, max_iter, eps, 1)
          b <- res$beta
          ll <- ll + 1
          if(length(b[is.na(b)]) == 0 | ll>15) {
            break
          }
        }
        b0 <- res$beta0
        b[is.na(b)] <- 0
        b0[is.na(b0)] <- 0
        iter_all[(l-1)*gg+gi] <- res$iter
        beta_all[,(l-1)*gg+gi] <- b
        beta0_all[1,(l-1)*gg+gi] <- b0
        loss[l] <- gammabinomial_loss(y, cbind(X, W), gammai, b, b0)
        nonzero[(l-1)*gg+gi] <- length(which(b != 0))
      }
    }
  }
  if (criterion=="BIC"){
    bic <- 2*n*loss + nonzero*log(n)
  } else {
    bic <- 2*n*loss + nonzero/5
  }
  bic0 <- bic[!is.na(bic)]
  mm <- which(bic == min(bic0))
  beta <- as.matrix(beta_all[,mm[1]])
  beta0 <- as.matrix(beta0_all[mm[1]])
  re <- list(beta, beta_all, beta0_all, loss, nonzero, iter_all, beta0)
  return (re)
}

gammagaussian_loss <- function(y, x, gamma, beta, beta0, sigma){
  
  ## ----------------------------------------------------------------------------
  ## The name of the function: gammagaussian_loss
  ## ----------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the value of the loss function (without penalty) 
  ##            based on gamma divergence for linear regression.
  ## ----------------------------------------------------------------------------
  ## Required preceding functions or packages: no
  ## ----------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable.
  ## @ x: n-row matrix, the design matrix.
  ## @ gamma: a float value, the robust parameter in gamma divergence.
  ## @ beta: Estimated (or true) regression coefficients (without intercept term).
  ## @ beta0: Estimated (or true) intercept term in regression coefficients.
  ## @ sigma: Estimated (or true) sigma (Standard deviation of response variables)
  ## ----------------------------------------------------------------------------
  ## Output:
  ## @ val_loss: the value of the loss function (without penalty).
  ## ----------------------------------------------------------------------------
  
  ei <- exp( -gamma * (y - x%*%beta - beta0)^2 / (2*sigma^2) )
  val_loss <- -mean( ( (1 + gamma)/( 2*pi*sigma^2 ) )^( gamma/2/(gamma+1) ) * ei )
  return(val_loss)
}

gammabinomial_loss <- function(y, x, gamma, beta, beta0){
  
  ## ----------------------------------------------------------------------------
  ## The name of the function: gammabinomial_loss
  ## ----------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the value of the loss function (without penalty) 
  ##            based on gamma divergence for logistic regression.
  ## ----------------------------------------------------------------------------
  ## Required preceding functions or packages: no
  ## ----------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable.
  ## @ x: n-row matrix, the design matrix.
  ## @ gamma: a float value, the robust parameter in gamma divergence.
  ## @ beta: Estimated (or true) regression coefficients (without intercept term).
  ## @ beta0: Estimated (or true) intercept term in regression coefficients.
  ## ----------------------------------------------------------------------------
  ## Output:
  ## @ val_loss: the value of the loss function (without penalty).
  ## ----------------------------------------------------------------------------
  
  ei_gamma = exp((gamma+1)*(x%*%beta+beta0))
  wi = ( (ei_gamma^y)/(1+ei_gamma) )^(gamma/(gamma+1))
  val_loss = - mean( wi )
  return(val_loss)
}

gammabinomial_weight <- function(y, x, gamma, beta, beta0){
  
  ## ----------------------------------------------------------------------------
  ## The name of the function: gammabinomial_weight
  ## ----------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the weights of the samples
  ##            based on gamma divergence for logistic regression.
  ## ----------------------------------------------------------------------------
  ## Required preceding functions or packages: no
  ## ----------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable.
  ## @ x: n-row matrix, the design matrix.
  ## @ gamma: a float value, the robust parameter in gamma divergence.
  ## @ beta: Estimated (or true) regression coefficients (without intercept term).
  ## @ beta0: Estimated (or true) intercept term in regression coefficients.
  ## ----------------------------------------------------------------------------
  ## Output:
  ## @ wi: the weights of the samples.
  ## ----------------------------------------------------------------------------
  
  ei_gamma = exp((gamma+1)*(x%*%beta+beta0))
  wi = ( (ei_gamma^y)/(1+ei_gamma) )^(gamma/(gamma+1))
  return(wi)
}


#--------------- density power divergence ------------
BIC_sgMCP_mdpd <- function(y, X, W, gamma = 0.5, b1, 
                           xi = 3, max_iter = 20, family, 
                           Matla, eps = 0.001, criterion="AIC", 
                           Armijo=1, step.size=0.1, step.size.logi=0.001){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: BIC_sgMCP_mdpd
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Performing joint analysis based on density power divergence with the sparse group penalty.
  ##            The regression coefficients under the given tuning parameters can be calculated, 
  ##            and the optimal regression coefficient is selected under the given criterion such as BIC.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: mdpdgaussian_loss()   mdpdbinomial_loss()
  ##            Rcpp functions: mdpdlinearsgMCP()  mdpdlogisticsgMCP()
  ##            R packages: Rcpp
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable
  ## @ X: n * q matrix, the design matrix corresponding environmental variables.
  ## @ W: n * ((q+1)*p) matrix, the design matrix corresponding gene variables and its interactive effects.
  ## @ gamma: A vector or a float value, the robust parameter in density power divergence, the default setting is 0.5.
  ## @ b1: ((q+1)*(p+1)) * 1 vector, initial value of the regression coefficient (including intercept terms).
  ## @ xi: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ max_iter: int, Maximum number of cycles of the algorithm, the default setting is 20.
  ## @ family: types of response variables, which can be selected from "gaussian" and "binomial".
  ## @ Matla: s * 2 matrix, the tuning parameters, 
  ##         the s rows represent s choices of the tuning parameters,
  ##         the 2 columns represent lambda1, lambda2,
  ##         (in this paper, to simplify, lambda1 = lambda2 is set, but this function is written in such a way that the two can be not equal.)
  ## @ eps: a float value, algorithm termination threshold.
  ## @ criterion: the parameter selection criterion, which can be selected from "AIC" and "BIC".
  ## @ Armijo: when updating main E effects using the group coordinate descent algorithm, 
  ##           Whether the Armijo search is used when determining step size (the default setting is not used).
  ##           (This parameter is for main E effects only for preventing algorithms from falling into undesirable minimum point, 
  ##           when updating the genes and their interaction effects, Armijo search is both used).
  ## @ step.size/step.size.logi: the given step size for linear/logistic regression (generally a smaller value) 
  ##                             when the Armijo search is not used when updating main E effects.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ beta: ( (q+1)*(p+1)-1 ) * 1 vector, the optimal regression coefficients selected according to set criterion (without intercept term).
  ## @ beta0: a float value, the intercept term in optimal regression coefficients.
  ## @ beta_all: ( (q+1)*(p+1)-1 ) * s vector, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ beta0_all: 1 * s vector, all intercept terms in regression coefficients corresponding s choices of given tuning parameters.
  ## @ loss: s * 1 vector, the values of the loss function (without penalty function) corresponding s choices of given tuning parameters.
  ## @ nonzero: s * 1 vector, the numbers of non-zero regression coefficients corresponding s choices of given tuning parameters.
  ## @ iter_all: Number of iterations.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  n <- length(y)
  p <- length(X[1,])+length(W[1,])
  gg <- length(gamma)
  L <- nrow(Matla)
  loss <- rep(0,L*gg)
  nonzero <- rep(0,L*gg)
  beta_all <- matrix(0,length(X[1,])+length(W[1,]),L*gg)
  beta0_all <- matrix(0,1,L*gg)
  iter_all <- rep(0,L*gg)
  Sigma <- rep(1,L*gg)
  for (l in 1:L) {
    matla <- Matla[l,]
    for (gi in 1:gg) {
      gammai <- gamma[gi]
      if (family == "gaussian"){
        ll <- 1
        repeat { 
          res = mdpdlinearsgMCP(y, X, W, gammai, b1, matla[1], matla[2], xi, max_iter, eps, Armijo, step.size, 1, 1)
          b <- res$beta
          ll <- ll + 1
          if(length(b[is.na(b)]) == 0 | ll>15) {
            break
          }
        }
        b0 <- res$beta0
        b[is.na(b)] <- 0
        b0[is.na(b0)] <- 0
        iter_all[(l-1)*gg+gi] <- res$iter
        beta_all[,(l-1)*gg+gi] <- b
        beta0_all[1,(l-1)*gg+gi] <- b0
        sigma <- res$sigma
        loss[(l-1)*gg+gi] <- mdpdgaussian_loss(y, cbind(X, W), gamma, b, b0, sigma)
        nonzero[(l-1)*gg+gi] <- length(which(b != 0))
        Sigma[(l-1)*gg+gi] <- sigma
      }else if (family == "binomial"){
        ll <- 1
        repeat { 
          res = mdpdlogisticsgMCP(y, X, W, gammai, b1, matla[1], matla[2], xi, max_iter, eps, Armijo, step.size.logi)
          b <- res$beta
          ll <- ll + 1
          if(length(b[is.na(b)]) == 0 | ll>15) {
            break
          }
        }
        b0 <- res$beta0
        b[is.na(b)] <- 0
        b0[is.na(b0)] <- 0
        iter_all[(l-1)*gg+gi] <- res$iter
        beta_all[,(l-1)*gg+gi] <- b
        beta0_all[1,(l-1)*gg+gi] <- b0
        loss[(l-1)*gg+gi] <- mdpdbinomial_loss(y, cbind(X, W), gamma, b, b0)
        nonzero[(l-1)*gg+gi] <- length(which(b != 0))
      }
    }
  }
  if (criterion=="BIC"){
    bic <- 2*n*loss + nonzero*log(n)
  } else {
    if (family == "gaussian"){
      bic <- 10*n*loss + nonzero
    }
    if (family == "binomial"){
      bic <- 6*n*loss + nonzero
    }
  }
  bic0 <- bic[!is.na(bic)]
  mm <- which(bic == min(bic0))
  beta <- as.matrix(beta_all[,mm[1]])
  beta0 <- as.matrix(beta0_all[,mm[1]])
  re <- list(beta, beta_all, beta0_all, loss, nonzero, iter_all, beta0)
  return (re)
}

mdpdgaussian_loss <- function(y, x, gamma, beta, beta0, sigma){
  
  ## ----------------------------------------------------------------------------
  ## The name of the function: mdpdgaussian_loss
  ## ----------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the value of the loss function (without penalty) 
  ##            based on density power divergence for linear regression.
  ## ----------------------------------------------------------------------------
  ## Required preceding functions or packages: no
  ## ----------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable.
  ## @ x: n-row matrix, the design matrix.
  ## @ gamma: a float value, the robust parameter in density power divergence.
  ## @ beta: Estimated (or true) regression coefficients (without intercept term).
  ## @ beta0: Estimated (or true) intercept term in regression coefficients.
  ## @ sigma: Estimated (or true) sigma (Standard deviation of response variables)
  ## ----------------------------------------------------------------------------
  ## Output:
  ## @ val_loss: the value of the loss function (without penalty).
  ## ----------------------------------------------------------------------------
  
  ei = exp( -gamma * (y - x%*%beta - beta0)^2 / (2*sigma^2) )
  val_loss = mean( 1/( (2*pi)^(gamma/2) * sigma^gamma * sqrt( 1+gamma ) ) - (1+gamma)/gamma * 1/( (2*pi)^(gamma/2) * sigma^gamma ) * ei )
}

mdpdbinomial_loss <- function(y, x, gamma, beta, beta0){
  
  ## ----------------------------------------------------------------------------
  ## The name of the function: mdpdbinomial_loss
  ## ----------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the value of the loss function (without penalty) 
  ##            based on density power divergence for logistic regression.
  ## ----------------------------------------------------------------------------
  ## Required preceding functions or packages: no
  ## ----------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable.
  ## @ x: n-row matrix, the design matrix.
  ## @ gamma: a float value, the robust parameter in density power divergence.
  ## @ beta: Estimated (or true) regression coefficients (without intercept term).
  ## @ beta0: Estimated (or true) intercept term in regression coefficients.
  ## ----------------------------------------------------------------------------
  ## Output:
  ## @ val_loss: the value of the loss function (without penalty).
  ## ----------------------------------------------------------------------------
  
  ei = exp(x%*%beta+beta0)
  wi = ( (ei^y)/(1+ei) )^gamma
  val_loss = mean( ( 1 + ei^(1 + gamma) ) / (( 1 + ei)^(1 + gamma)) - (1 + 1/gamma) * wi )
  return(val_loss)
}

mdpdbinomial_weight <- function(y, x, gamma, beta, beta0){
  
  ## ----------------------------------------------------------------------------
  ## The name of the function: mdpdbinomial_weight
  ## ----------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the weights of the samples
  ##            based on density power divergence for logistic regression.
  ## ----------------------------------------------------------------------------
  ## Required preceding functions or packages: no
  ## ----------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable.
  ## @ x: n-row matrix, the design matrix.
  ## @ gamma: a float value, the robust parameter in density power divergence.
  ## @ beta: Estimated (or true) regression coefficients (without intercept term).
  ## @ beta0: Estimated (or true) intercept term in regression coefficients.
  ## ----------------------------------------------------------------------------
  ## Output:
  ## @ wi: the weights of the samples.
  ## ----------------------------------------------------------------------------
  
  ei = exp(x%*%beta+beta0)
  wi = ( (ei^y)/(1+ei) )^gamma / 0.5 - 1
  return(wi)
}



#########################################################################################


#########################################################################################
# Functions for marginal analysis:
# This section includes two proposed methods: 
# gamma divergence, density power divergence.
#########################################################################################

gamma_marginal <- function(y, X, W, gamma = 0.5, Matla, family, initial=F, b1){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: gamma_marginal
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Performing marginal analysis based on gamma divergence with the sparse group penalty.
  ##            A sequence of lambda is considered, and the regression coefficients are estimated at each value.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            Rcpp functions: gammalinearsgMCP()  gammalogisticsgMCP()
  ##            R packages: Rcpp
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable
  ## @ X: n * q matrix, the design matrix corresponding environmental variables.
  ## @ W: n * ((q+1)*p) matrix, the design matrix corresponding gene variables and its interactive effects.
  ## @ gamma: A vector or a float value, the robust parameter in gamma divergence, the default setting is 0.5.
  ## @ Matla: s * 2 matrix, the tuning parameters, 
  ##         the s rows represent s choices of the tuning parameters,
  ##         the 2 columns represent lambda1, lambda2,
  ##         (in this paper, to simplify, lambda1 = lambda2 is set, but this function is written in such a way that the two can be not equal.)
  ## @ family: types of response variables, which can be selected from "gaussian" and "binomial".
  ## @ initial: logical variable, the default setting is F, the initial value of the regression coefficient is 0s; 
  ##                              if initial = T, the initial value uses b1.
  ## @ b1: ((q+1)*(p+1)) * 1 vector, initial value of the regression coefficient (including intercept terms).
  ## 
  ## Remark: In this function, the following parameters are set as default values for functions gammalinearsgMCP() and gammalogisticsgMCP():
  ##                 @ xi(the regularization parameter in MCP) defaults to 3;
  ##                 @ max_iter(Maximum number of cycles of the algorithm) defaults to 20;
  ##                 @ eps(algorithm termination threshold) defaults to 0.001.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ b_all: ( (q+1)*(p+1)-1 ) * s vector, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ loss_all: p * s vector, the values of the loss function (without penalty function) corresponding 
  ##                           every marginal model(p-row) and s choices of given tuning parameters(s-column).
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  q <- dim(X)[2]
  p <- dim(W)[2]/(q+1)
  b00 <- rep(0,2*(q+1))
  gg <- length(gamma)
  nlambda1 <- length(Matla[,1])
  b_all <- matrix(0,q+(1+q)*p,nlambda1*gg)
  b_la <- rep(0,q+(1+q)*p)
  loss <- rep(0,p)
  loss_all <- matrix(0,p,nlambda1*gg)
  for (ll in 1:nlambda1) {
    lambda1 <- Matla[ll,1];lambda2 <- Matla[ll,2]
    for (gi in 1:gg) {
      gammai <- gamma[gi]
      for (j in 0:(p-1)) {
        Wj <- W[,((q+1)*j+1):((q+1)*j+q+1)]
        if(initial){b00 <- b1[c(1:(q+1),((j+1)*(q+1)+1) : ((j+2)*(q+1)))]}
        if (family == "gaussian"){
          res <- gammalinearsgMCP(y, X, Wj, gammai, b00, lambda1, lambda2, 3, 20, 0.001, 0)
        }else if (family == "binomial"){
          lll <- 1
          repeat { 
            res = gammalogisticsgMCP(y, X, Wj, gammai, b00, lambda1, lambda2, 3, 20, 0.001, 0)
            bj <- res$beta
            lll <- lll + 1
            if(length(bj[is.na(bj)]) == 0 | ll>5) {
              break
            }
          }
        }
        
        bj <- res$beta
        bj[is.na(bj)] <- 0
        b_la[c(1:q,q+((q+1)*j+1):((q+1)*j+q+1))] <- bj
        loss[j+1] <- res$obj
        
      }
      loss[is.na(loss)] <- 0
      b_all[,(ll-1)*gg+gi] <- b_la
      loss_all[,(ll-1)*gg+gi] <- loss
    }
  }
  return(list(b_all,loss_all))
}

mdpd_marginal <- function(y, X, W, gamma = 0.5, Matla, family, initial=F, b1, Armijo=1, step.size=0.001){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: mdpd_marginal
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Performing marginal analysis based on density power divergence with the sparse group penalty.
  ##            A sequence of lambda is considered, and the regression coefficients are estimated at each value.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            Rcpp functions: mdpdlinearsgMCP()  mdpdlogisticsgMCP()
  ##            R packages: Rcpp
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable
  ## @ X: n * q matrix, the design matrix corresponding environmental variables.
  ## @ W: n * ((q+1)*p) matrix, the design matrix corresponding gene variables and its interactive effects.
  ## @ gamma: A vector or a float value, the robust parameter in gamma divergence, the default setting is 0.5.
  ## @ Matla: s * 2 matrix, the tuning parameters, 
  ##         the s rows represent s choices of the tuning parameters,
  ##         the 2 columns represent lambda1, lambda2,
  ##         (in this paper, to simplify, lambda1 = lambda2 is set, but this function is written in such a way that the two can be not equal.)
  ## @ family: types of response variables, which can be selected from "gaussian" and "binomial".
  ## @ initial: logical variable, the default setting is F, the initial value of the regression coefficient is 0s; 
  ##                              if initial = T, the initial value uses b1.
  ## @ b1: ((q+1)*(p+1)) * 1 vector, initial value of the regression coefficient (including intercept terms).
  ## @ Armijo: when updating main E effects using the group coordinate descent algorithm, 
  ##           Whether the Armijo search is used when determining step size (the default setting is not used).
  ##           (This parameter is for main E effects only for preventing algorithms from falling into undesirable minimum point, 
  ##           when updating the genes and their interaction effects, Armijo search is both used).
  ## @ step.size: the given step size (generally a smaller value) 
  ##                             when the Armijo search is not used when updating main E effects.
  ## 
  ## Remark: In this function, the following parameters are set as default values for functions mdpdlinearsgMCP() and mdpdlogisticsgMCP():
  ##                 @ xi(the regularization parameter in MCP) defaults to 3;
  ##                 @ max_iter(Maximum number of cycles of the algorithm) defaults to 20;
  ##                 @ eps(algorithm termination threshold) defaults to 0.001.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ b_all: ( (q+1)*(p+1)-1 ) * s vector, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ loss_all: p * s vector, the values of the loss function (without penalty function) corresponding 
  ##                           every marginal model(p-row) and s choices of given tuning parameters(s-column).
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  q <- dim(X)[2]
  p <- dim(W)[2]/(q+1)
  b00 <- rep(0,2*(q+1))
  gg <- length(gamma)
  nlambda1 <- length(Matla[,1])
  b_all <- matrix(0,q+(1+q)*p,nlambda1*gg)
  b_la <- rep(0,q+(1+q)*p)
  loss <- rep(0,p)
  loss_all <- matrix(0,p,nlambda1*gg)
  for (ll in 1:nlambda1) {
    lambda1 <- Matla[ll,1];lambda2 <- Matla[ll,2]
    for (gi in 1:gg) {
      gammai <- gamma[gi]
      for (j in 0:(p-1)) {
        Wj <- W[,((q+1)*j+1):((q+1)*j+q+1)]
        if(initial){b00 <- b1[c(1:(q+1),((j+1)*(q+1)+1) : ((j+2)*(q+1)))]}
        if (family == "gaussian"){
          res <- mdpdlinearsgMCP(y, X, Wj, gammai, b00, lambda1, lambda2, 3, 20, 0.001, 1, 2, 0, 0)
        }else if (family == "binomial"){
          lll <- 1
          repeat { 
            res <- mdpdlogisticsgMCP(y, X, Wj, gammai, b00, lambda1, lambda2, 3, 20, 0.001, Armijo, step.size)
            bj <- res$beta
            lll <- lll + 1
            if(length(bj[is.na(bj)]) == 0 | ll>2) {
              break
            }
          }
        }
        bj <- res$beta
        bj[is.na(bj)] <- 0
        b_la[c(1:q,q+((q+1)*j+1):((q+1)*j+q+1))] <- bj
        loss[j+1] <- res$obj
      }
      b_la[is.na(b_la)] <- 0
      loss[is.na(loss)] <- 0
      b_all[,(ll-1)*gg+gi] <- b_la
      loss_all[,(ll-1)*gg+gi] <- loss
    }
  }
  return(list(b_all,loss_all))
}

mar_TFPR <- function(b_all,loss_all,coe=2, criterion="AIC"){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: mar_TFPR
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Performing selection of the optimal regression coefficient for marginal analysis
  ##            under the given criterion such as BIC.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ b_all: ( (q+1)*(p+1)-1 ) * s matrix, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ loss_all: p * s matrix, the values of the loss function (without penalty function) corresponding 
  ##                           every marginal model(p-row) and s choices of given tuning parameters(s-column).
  ## @ coe: The weight coefficient of the likelihood(loss function) in the selection criterion, the default setting is 2.
  ## @ criterion: the parameter selection criterion, which can be selected from "AIC" and "BIC".
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ beta: the optimal regression coefficient under the given criterion such as BIC.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  
  nonzero <- apply(b_all,2,function(a){sum(as.numeric(a != 0))})
  loss <- apply(loss_all,2,sum)
  if (criterion=="BIC"){
    bic <- 2*n*loss + nonzero*log(n)
  } else {
    bic <- coe*n*loss + nonzero
  }
  bic0 <- bic[!is.na(bic)]
  mm <- which(bic == min(bic0))
  beta <- as.matrix(b_all[,mm[1]])
  return(beta)
}




############################# The function for generating the lambda ############################
gen_lambda <- function(min,max,nlambda=20){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: gen_lambda
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating a sequence of the tuning parameters (lambda1 and lambda2).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ min: The minimum lambda.
  ## @ max: The maximum lambda.
  ## @ nlambda: The number of lambda.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ Matla: nlambda * 2 matrix,   
  ##          the nlambda rows represent nlambda choices of the tuning parameters,
  ##          the 2 columns represent lambda1 and lambda2 (lambda1 = lambda2 is defaulted).
  ## -----------------------------------------------------------------------------------------------------------------
  
  lambda1 = exp(seq(log(min),log(max),len= nlambda))
  Matla <- cbind(lambda1,lambda1)
  return(Matla)
}

############################# Functions for generate the initial beta ############################
initial <- function(y, X, W, family='gaussian', methods, b00=F, n_alter,S_diver="S1", realdata=F, G, s=0.05){
  if (!realdata){
    q <- dim(X)[2]
    p <- dim(W)[2]/(q+1)
    if(b00==T){
      b1 <- rep(0,(q+1)*(p+1))
    } else {
      if(methods=='null'){
        lasso <- glmnet(cbind(X,W),y,family='gaussian')
        coefficients <- coef(lasso,s=0.5)
        non <- which(coefficients!=0)[-1]
        b_grl <- rep(0,(q+1)*(p+1)-1)
        b_grl[non-1] <- coefficients[non]
        b1 <- c(0,b_grl)
      } else if(methods=='alter'){
        if(n_alter==1){
          lasso <- glmnet(cbind(X,W),y,family='gaussian')
          coefficients <- coef(lasso,s=0.1)
          non <- which(coefficients!=0)[-1]
          b_grl <- rep(0,(q+1)*(p+1)-1)
          b_grl[non-1] <- coefficients[non]
          b1 <- c(0,b_grl)
        }else if(n_alter==2){
          lasso <- glmnet(cbind(X,W),y,family='gaussian')
          coefficients <- coef(lasso,s=0.5)
          non <- which(coefficients!=0)[-1]
          b_grl <- rep(0,(q+1)*(p+1)-1)
          b_grl[non-1] <- coefficients[non]
          b1 <- c(0,b_grl)
        }
      } else if(methods=='divergence'){
        if(S_diver=="S0"){
          grlasso <- grpreg(cbind(X,W), y, group=group,penalty=c("grLasso"),family=c("gaussian"))
          fit <- select(grlasso, "BIC")
          b_grl <- fit$beta
          b_grl <- b_grl[-1]
          b1 <- c(0,b_grl)
        } else if (S_diver=="S00"){
          grlasso <- grpreg(cbind(X,W), y, group=group,penalty=c("grLasso"),family=c("gaussian"))
          fit <- select(grlasso, "AIC")
          b_grl <- fit$beta
          b_grl <- b_grl[-1]
          b1 <- c(0,b_grl)
        } else {
          grlasso <- cv.grpreg(cbind(X,W), y, group=group,penalty=c("grLasso"),family=c("gaussian"))
          betaall <- grlasso$fit[[1]]
          b_grl <- betaall[-1,match(grlasso$lambda.min,grlasso$lambda)]
          b1 <- c(0,b_grl)
        }
      }
    }
    return(b1)
  } else {
    lasso <- glmnet(as.matrix(cbind(G)),y,family='gaussian')
    coefficients <- coef(lasso,s=s)
    coefficients <- coefficients[-1]
    non <- which(coefficients!=0)
    b_grl <- rep(0,(q+1)*(p+1)-1)
    b_grl[non*(q+1)] <- coefficients[non]
    b1 <- c(0,b_grl)
    return(b1)
  }
  
}

initial_binomial <- function(y, X, W, family='binomial', methods, n_alter, b00=F,S_diver="S1"){
  q <- dim(X)[2]
  p <- dim(W)[2]/(q+1)
  if(b00==T){
    b1 <- rep(0,(q+1)*(p+1))
  } else {
    if(methods=='null'){
      grlasso <- grpreg(cbind(X,W), y, group=group,penalty=c("grLasso"),family=c("binomial"))
      fit <- select(grlasso, "AIC")
      b_grl <- fit$beta
      b_grl <- b_grl[-1]
      b1 <- c(0,b_grl)
    } else if(methods=='alter'){
      if(n_alter==1){
        lasso <- glmnet(cbind(X,W),y,family='binomial')
        coefficients <- coef(lasso,s=0.001) 
        non <- which(coefficients!=0)[-1]
        b_grl <- rep(0,(q+1)*(p+1)-1)
        b_grl[non-1] <- coefficients[non]
        b1 <- c(0,b_grl)
      }else if(n_alter==2){
        lasso <- glmnet(cbind(X,W),y,family='binomial')
        coefficients <- coef(lasso,s=0.5)
        non <- which(coefficients!=0)[-1]
        b_grl <- rep(0,(q+1)*(p+1)-1)
        b_grl[non-1] <- coefficients[non]
        b1 <- c(0,b_grl)
      }
    } else if(methods=='divergence'){
      if(S_diver=="S0"){
        grlasso <- grpreg(cbind(X,W), y, group=group,penalty=c("grLasso"),family=c("binomial"))
        fit <- select(grlasso, "BIC")
        b_grl <- fit$beta
        b_grl <- b_grl[-1]
        b1 <- c(0,b_grl)
      } else{
        grlasso <- grpreg(cbind(X,W), y, group=group, lambda = 0.01,penalty=c("grLasso"),family=c("binomial"))
        b_grl <- grlasso$beta
        b_grl <- b_grl[-1,]
        b1 <- c(0,b_grl)
      }
    }
  }
  return(b1)
}
