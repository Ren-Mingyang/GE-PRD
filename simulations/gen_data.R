#####################################################################################
# Codes for simulation studies: 
# The functions for generating simulated data & evaluating performances of competitors,
# and some general parameter settings,
# which are used to support numerical simulation studies.
#####################################################################################


#####################################################################################
# Functions for generation of simulated data
# This section includes three functions: GenGE()  CorMat()
#####################################################################################

GenGE <- function(n, p, q, Gcorr="En", Ecorr="En", mis,
               alpha0=0, a, b, family="gaussian", is.SNP = FALSE, 
               maf= 0.3, E.cat = 0, eta0=0.05, eta1=0.2){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: GenGE
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the simulated data.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: CorMat()
  ##            R packages: mvtnorm
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ n: int, the number of subjects.
  ## @ p: int, the dim of gene variables.
  ## @ q: int, the dim of environment variables.
  ## @ Gcorr: p * p matrix, the correlation matrix of gene variables.
  ## @ Ecorr: q * q matrix, the correlation matrix of environment variables.
  ## @ mis: the error distribution of the residual or mislabel, which can be selected from "S0","S1", and "S2".
  ## @ alpha0: a float value, the intercept term in true regression coefficients, the default setting is 0. 
  ## @ a: q * 1 vector, the true regression coefficients corresponding environmental variables. 
  ## @ b: p*(q+1) matrix, the true regression coefficients corresponding gene variables and its interactive effects. 
  ## @ family, types of response variables, which can be selected from "gaussian" and "binomial".
  ## @ E.cat, int, the number of the discrete environment variables.
  ## @ is.SNP, logical variable, if True, the gene variables are a SNP type.
  ## @ maf, float, maf \in (0,0.5), which needs to be set only when is.SNP = T.
  ## @ eta0: mislabeled probability (0 -- 1) for mislabel logistic regression, the default setting is 0.05. 
  ## @ eta1: mislabeled probability (1 -- 0) for mislabel logistic regression, the default setting is 0.2.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list:
  ## @ y: n * 1 vector, the response variable
  ## @ X: n * q matrix, the design matrix corresponding environmental variables.
  ## @ W: n * ((q+1)*p) matrix, the design matrix corresponding gene variables and its interactive effects.
  ## ---------------------------------------------------------------------------------------------------------------
  
  Zsigmat <- CorMat(p, Gcorr)
  Xsigmat <- CorMat(q, Ecorr)
  Z <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=Zsigmat) 
  X <- mvtnorm::rmvnorm(n, mean=rep(0,q), sigma=Xsigmat)

  if (is.SNP){
    MAF.vec <- rep(maf, p)
    cum.fre <- rbind((1-MAF.vec)^2, 1-MAF.vec^2)
    cutPoint <- qnorm(cum.fre)
    Z <- sapply(1:p,  function(j){
      sapply(Z[,j], function(x){ sum(x > cutPoint[,j])})
    })

  }

  if (E.cat > 0){
      X[, 1:E.cat] <- as.numeric(X[,1:E.cat] > 0 )
  }
  
  W <- matrix(NA, nrow=n, ncol=p*(q+1))
  X1 <- cbind(1,X)
  for (j in 1:p){
    W[,(1:(q+1))+(q+1)*(j-1)] <- Z[,j]*X1
  }

  ymean = alpha0+ X%*%a + W%*%b

  ### continue response
  if (family == "gaussian"){
     ##error distribution 
      ep = rep(0,n)
      if (mis == "M_NC"){
        ep =  ifelse(runif(n,0,1)< 0.6,rnorm(n),rcauchy(n))
      }else if(mis == "S1"){
        ep  = ifelse(runif(n,0,1)< 0.6,rnorm(n),exp(rnorm(n)))
      }else if(mis== "S2"){
        ep=rt(n,5)
      }else if(mis == "S0"){
        ep=rnorm(n)
      }
      y = ymean+ep
      data = list(Y=y,Z=Z, X=X, W=W, mis =ep, Ecor=Xsigmat, Gcor=Zsigmat)
  }else if(family == "binomial"){
      prob = 1/(1+exp(-ymean))
      ytrue <- rbinom(n,1,prob)
      if (mis == "S0"){
        y =  ytrue
      }else if(mis == "S1"){
        y = ytrue
        jjytrue1 = which(ytrue == 1)
        jjytrue0 = which(ytrue == 0)
        eta00 = rep(eta0,length(jjytrue0))
        eta11 = rep(eta1,length(jjytrue1))
        y[jjytrue1[rbinom(length(eta11),1,eta11) == 1]] = 0
        y[jjytrue0[rbinom(length(eta00),1,eta00) == 1]] = 1
      }else if(mis == "S2"){
        y = ytrue
        jjytrue1 = which(ytrue == 1)
        jjytrue0 = which(ytrue == 0)
        eta00 = rep(eta0,length(jjytrue0))
        eta11 = eta0 + (eta1 - eta0) * prob
        y[jjytrue1[rbinom(length(jjytrue1),1,eta11[jjytrue1]) == 1]] = 0
        y[jjytrue0[rbinom(length(eta00),1,eta00) == 1]] = 1
      }
      data = list(Y=y,Z=Z, X=X, W=W, Ytrue=ytrue, Ecor=Xsigmat, Gcor=Zsigmat)
    }  
  return(data)
}

CorMat = function(m, sig.index = "En"){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: CorMat
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the correlation structure of genes
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding data: No
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ m: Dimensions of the correlation matrix.
  ## @ sig.index: The selection of the correlation structure of the design matrix
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ sig.G: The correlation matrix.
  ## ---------------------------------------------------------------------------------------------------------------
  
  sig.G= matrix(0,m,m)
  if (sig.index != "En"){
    for (i in 1:m){
      for(j in i:m){
        if(sig.index=="AR1"){
          sig.G[i,j]=0.25^abs(i-j)
        }else if(sig.index=="AR2"){
          sig.G[i,j]=0.75^abs(i-j)
        }else if(sig.index=="B2"){
          if(abs(i-j) ==1){
            sig.G[i,j]= 0.6
          }else if (abs(i-j)==2){
            sig.G[i,j] = 0.3
          }
        }else if(sig.index=="B1"){
          if(abs(i-j) ==1){
            sig.G[i,j]= 0.3
          }
        }else if(sig.index=="AR_E"){
          sig.G[i,j]=0.2^abs(i-j)
        }
        sig.G[j,i] = sig.G[i,j]
      }
    }
  }
  diag(sig.G)= 1
  return(sig.G)
}


################################################################################
# Functions for evaluating performances of proposed methods:
# This section includes certain evaluate functions to evaluate 
# the identification performance for marginal and joint analysis
################################################################################

#------------evaluation function------------
index_over <- function (b_set,  k_main_nonzero, k_Interaction_nonzero, q, p){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: index_over
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the true positive rate (TPR) and  
  ##            false positive rate (FPR) of identifying important main and interaction effects.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: no
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ b_set: ( (q+1)*(p+1)-1 ) * s matrix, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ k_main_nonzero: (signumberG+signumberE) * 1 vector, the location of significant main effects.
  ## @ k_Interaction_nonzero: (signumberGE) * 1 vector, the location of significant interaction effects.
  ## @ p: int, the dim of gene variables.
  ## @ q: int, the dim of environment variables.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ index_all: s * 8 matrix, the number of correct/incorrect identifications for important main and interaction effects,
  ##                            the true positive rate (TPR), and false positive rate (FPR)
  ##                            corresponding s choices of given tuning parameters.
  ## -----------------------------------------------------------------------------------------------------------------
  
  k_all = 1:(q + p*(q+1))
  k_Main_effect = c(1:q,k_all[k_all%%(q+1) == 0])
  k_Interaction = setdiff(k_all, k_Main_effect)
  k_main_zero = setdiff(k_Main_effect, k_main_nonzero)
  k_main_zero = setdiff(k_main_zero, 1:q)
  k_Interaction_zero = setdiff(k_Interaction, k_Interaction_nonzero)
  b_set = as.matrix(b_set)
  pi = dim(b_set)[2]
  index_all = as.data.frame(matrix(rep(0,8*pi),ncol = 8))
  for(i in 1:pi){
    beta_index = b_set[,i]
    # Main_effect nonzero
    beta_main_nonzero = beta_index[k_main_nonzero]
    beta_main_zero = beta_index[k_main_zero]
    C_main = mean(sum(beta_main_nonzero != 0))
    IC_main = mean(sum(beta_main_zero != 0))
    
    # Interaction_effect nonzero
    beta_Interaction_nonzero = beta_index[k_Interaction_nonzero]
    beta_Interaction_zero = beta_index[k_Interaction_zero]
    C_Interaction = mean(sum(beta_Interaction_nonzero != 0))
    IC_Interaction = mean(sum(beta_Interaction_zero != 0))
    index_all[i,1:4] = c(C_main,IC_main,C_Interaction,IC_Interaction)
    
    TPR_main = C_main/length(k_main_nonzero)
    FPR_main = IC_main/(p - length(k_main_nonzero) )
    TPR_inte = C_Interaction/length(k_Interaction_nonzero)
    FPR_inte = IC_Interaction/(p*q - length(k_Interaction_nonzero) )
    index_all[i,5:8] = c(TPR_main,FPR_main,TPR_inte,FPR_inte)
  }
  names(index_all)[5:8] <- rep(c("TPR","FPR"),2)
  return(index_all)
}

index_AUC <- function(TFPR,r1=0.3,r2=0.5,n1=20,n2=40){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: index_AUC
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the AUC, pAUC, and top 20/40 of identifying important main and interaction effects
  ##            based on a sequence of lambda for marginal analysis.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: AUC()   pAUC()   top()
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ TFPR: s * 8 matrix, the output table of R function index_over().
  ## @ r1: Upper bound of FPR needed to calculate pAUC1(the AUC with 0 and r1 marking the range of FPR).
  ## @ r2: Upper bound of FPR needed to calculate pAUC2(the AUC with 0 and r2 marking the range of FPR).
  ## @ n1: top-n1(defined as the number of true positives when n1 variables are identified).
  ## @ n2: top-n2(defined as the number of true positives when n2 variables are identified).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ AUCtop: s * 10 matrix, the AUC, pAUC1, pAUC2, and top 20/40 of identifying important main and interaction effects
  ##                            and interaction effects corresponding s choices of given tuning parameters.
  ## -----------------------------------------------------------------------------------------------------------------
  
  AUC_main = AUC(TFPR[,5:6])
  AUC_inte = AUC(TFPR[,7:8])
  pAUC1_main = pAUC(TFPR[,5:6],r1)
  pAUC1_inte = pAUC(TFPR[,7:8],r1)
  pAUC2_main = pAUC(TFPR[,5:6],r2)
  pAUC2_inte = pAUC(TFPR[,7:8],r2)
  top1_main = top(TFPR[,1:2],n1)
  top1_inte = top(TFPR[,3:4],n1)
  top2_main = top(TFPR[,1:2],n2)
  top2_inte = top(TFPR[,3:4],n2)
  AUCtop = as.data.frame( t(as.matrix(c(AUC_main, AUC_inte, pAUC1_main, pAUC1_inte,
                                        pAUC2_main, pAUC2_inte,top1_main, top1_inte,top2_main, top2_inte))) )
  names(AUCtop) <- c("AUC_main", "AUC_inte", "pAUC0.3_main", "pAUC0.3_inte", 
                     "pAUC0.5_main", "pAUC0.5_inte", "top20_main", "top20_inte","top40_main", "top40_inte")
  return(AUCtop)
}

AUC <- function(TFPR){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: AUC
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the approximate estimated AUC of identifying important main and interaction effects 
  ##            based on a sequence of lambda for marginal analysis.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ TFPR: s * 8 matrix, the output table of R function index_over().
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ AUC: a float value, the approximate estimated AUC of identifying important main and interaction effects
  ##                            and interaction effects corresponding s choices of given tuning parameters.
  ## -----------------------------------------------------------------------------------------------------------------
  
  TFPR = TFPR[order(TFPR[,2]),]
  ss = dim(TFPR)[1]
  gg = dim(TFPR)[2]
  ss = max(ss,gg)
  
  TPRy = c(0,TFPR[,1],1)
  FPRx = c(0,TFPR[,2],1)
  ss = min(which(TPRy==1))
  TPRy = TPRy[1:ss]
  FPRx = FPRx[1:ss]
  TPRy = sort(c(0,TPRy,1))
  FPRx = sort(c(0,FPRx,1))
  AUC = 0
  for (s in 1:(ss+1)){
    AUC = AUC + (TPRy[s]+TPRy[s+1]) * abs(FPRx[s+1] - FPRx[s])/2
  }
  return(AUC)
}

pAUC <- function(TFPR,r){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: pAUC
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the approximate estimated pAUC of identifying important main and interaction effects 
  ##            based on a sequence of lambda for marginal analysis.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ TFPR: s * 8 matrix, the output table of R function index_over().
  ## @ r: Upper bound of FPR needed to calculate pAUC(the AUC with 0 and r marking the range of FPR).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ AUC: a float value, the approximate estimated pAUC of identifying important main and interaction effects
  ##                            and interaction effects corresponding s choices of given tuning parameters.
  ## -----------------------------------------------------------------------------------------------------------------
  
  TFPR = TFPR[order(TFPR[,2]),]
  FPRx = sort(TFPR[,2])
  a = FPRx[FPRx < r]
  ss  = length(a)
  TPRy = c(0,TFPR[1:length(a),1],TFPR[length(a),1])
  FPRx = c(0,TFPR[1:length(a),2],r)
  AUC = 0
  for (s in 1:(ss+1)){
    AUC = AUC + min(TPRy[s],TPRy[s+1]) * abs(FPRx[s+1] - FPRx[s])
  }
  AUC = AUC/r
  return(AUC)
}

top <- function(TFPR,n){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: top
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the estimated top 20/40 of identifying important main and interaction effects 
  ##            based on a sequence of lambda for marginal analysis.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ TFPR: s * 8 matrix, the output table of R function index_over().
  ## @ n: top-n(defined as the number of true positives when n variables are identified).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ nn: int, the estimated top 20/40 of identifying important main and interaction effects
  ##            and interaction effects corresponding s choices of given tuning parameters.
  ## -----------------------------------------------------------------------------------------------------------------
  
  a = apply(TFPR,1,sum)
  b = abs(n - a)
  c = which(b == min(b))
  c = max(c)
  nn = TFPR[c,1]
  return(nn)
}

# Integration of results into a standardized tabular format
trans_table_joint <- function(index, family, gaussian_rownames, binomial_rownames, or=4){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: trans_table_joint
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Integrating the output results for joint analysis into a standardized tabular format in the paper.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ index: 4 * 8 matrix, simulation results (including mean and standard deviation) after o replicates.
  ## @ family: types of response variables, which can be selected from "gaussian" and "binomial".
  ## @ gaussian_rownames: c("LS","LAD","dpd","$gamma$-divergence").
  ## @ binomial_rownames: c("logistic","constant","dpd","$gamma$-divergence").
  ## @ or: Numbers of columns in the final output table, the default value is 4.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ index0: s * 4 matrix, the final output standardized table in the paper.
  ## -----------------------------------------------------------------------------------------------------------------
  
  index0 <- index
  for (j in c(1:or)) {
    index0[,j] <- paste(round(index0[,j],3),"(",round(index0[,or+j],3),")",sep = "")
  }
  index0 <- index0[,1:or]
  names(index0) <- joint_colnames
  if(family == "gaussian"){
    rownames <- gaussian_rownames
  } 
  if(family == "binomial"){
    rownames <- binomial_rownames
  }
  row.names(index0) <- rownames
  return(index0)
}
trans_table_marginal <- function(index, family, gaussian_rownames, binomial_rownames, or=10){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: trans_table_marginal
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Integrating the output results for marginal analysis into a standardized tabular format in the paper.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ index: 4 * 8 matrix, simulation results (including mean and standard deviation) after o replicates.
  ## @ family: types of response variables, which can be selected from "gaussian" and "binomial".
  ## @ gaussian_rownames: c("LS","LAD","dpd","$gamma$-divergence").
  ## @ binomial_rownames: c("logistic","constant","dpd","$gamma$-divergence").
  ## @ or: Numbers of columns in the final output table, the default value is 10.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ index0: s * 4 matrix, the final output standardized table in the paper.
  ## -----------------------------------------------------------------------------------------------------------------
  
  index0 <- index
  for (j in c(1:or)) {
    index0[,j] <- paste(round(index0[,j],3),"(",round(index0[,or+j],3),")",sep = "")
  }
  index0 <- index0[,1:or]
  names(index0) <- names(index)[1:or]
  index0 <- index0[,c(c(1:(or/2))*2-1,c(1:(or/2))*2)]
  
  if(family == "gaussian"){
    rownames <- gaussian_rownames
  } 
  if(family == "binomial"){
    rownames <- binomial_rownames
  }
  row.names(index0) <- rownames
  return(index0)
}

#####################################################################################
# Some general parameter settings:
# This section includes some general parameter settings for numerical simulation studies
#####################################################################################
n.seq <- c(200, 400)                             # sample sizes
family.seq <- c("gaussian","binomial")           # types of response variables
corr.seq <- c("AR1","AR2","B1","B2")             # The correlation structure of the design matrix
mis.seq <- c("S0","S1","S2")                     # types of contamination or mislabeled mechanisms
Ecorr <- "AR_E"                                  # correlation structure of environmental variables
is.SNP <- FALSE                                  # the gene variables are not SNP type
E.cat <- 2                                       # the number of the discrete environment variables: 2.
eta1 <- 0.2                                      # mislabeled probability (1 -- 0) for mislabel logistic regression.
####---------Generation of coefficients with hierarchical structure-----
setG <- sort(sample(1:p, size = signumberG, replace = FALSE))
setE <- sort(sample(1:q, size = signumberE, replace = FALSE))
se = cbind(sort(rep(setG,signumberE)),rep(setE,signumberG))
setGE = matrix(se[sample(1:nrow(se),size = signumberGE,replace = F),],nrow= signumberGE,ncol=2)
alpha0 = 0
alpha = rep(0,q)
beta = rep(0,p)
ggamma = matrix(0,p,q)
alpha[setE] = runif(signumberE, 0.2,0.8)
beta[setG] = runif(signumberG, 0.2,0.8)
ggamma[setGE] = runif(signumberGE, 0.2,0.8)
a = alpha
b = as.vector(t(cbind(beta,ggamma)))
vgamma = as.vector(ggamma)
beta0 = c(a,b)
k_main_nonzero = setG*(q+1)                               # location of significant main effects
k_Interaction_nonzero = sort(setGE[,1]*(q+1)+setGE[,2])   # location of significant interaction effects
m_gamma <-c(0.3)
gamma <-c(0.5)
#  Group coding
group <- c()
group[1:q] <- 1
for (j in 1:p) {
  group[((q+1)*j):((q+1)*j+q)] <- j+1
}




