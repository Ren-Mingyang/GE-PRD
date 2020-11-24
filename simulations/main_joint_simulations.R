##################################################################
# Codes for simulation studies: 
# joint analysis for logistic regression and linear regression
##################################################################
rm(list = ls(all = TRUE))
ls()
##################################################################

library("Rcpp")
library("grpreg")
library("glmnet")
sourceCpp("Rcpp_fun_mdpdlogistic.cpp")
sourceCpp("Rcpp_fun_gammalogistic.cpp")
sourceCpp("Rcpp_fun_mdpdlinear.cpp")
sourceCpp("Rcpp_fun_gammalinear.cpp")
source("function.R")



# The column and row names of the output table
joint_colnames <- c("main-TPR","main-FPR","interaction-TPR","interaction-FPR")
binomial_rownames <- rep(c("dpd","$gamma$-divergence"),1)
gaussian_rownames <- rep(c("dpd","$gamma$-divergence"),1)

# the number of repetitions
o=2


# ----------------- Parameter setting --------------------
p <- 1000               # the dimension of genetic variables
q <- 5                  # the dimension of environmental variables
signumberG <- 10        # the dimension of important genetic variables
signumberGE <- 20       # the dimension of important interaction effects
signumberE <- q         # the dimension of important environmental variables
source("gen_data.R")    # Generating regression coefficients


############### logistic regression ###############
family <- family.seq[2]
############### S1 & AR1 & n1 ###############
mis <- mis.seq[2]
Gcorr <- corr.seq[1]
n <- 300


Matla_mdpd <- gen_lambda(0.02,0.005,10)
Matla_gamma <- gen_lambda(0.01,0.001,10)
filename <- paste(paste("joint_binomial",mis,paste("p",p,sep=""),Gcorr,paste("n",n,sep=""),sep = "_"),"RData",sep = ".")
filename <- paste(path,filename,sep = "")
mdpd_b <- as.data.frame(matrix(0,q + p*(q+1),o))
mdpd_beta_MSE <- rep(0,o)
mdpd_index <- as.data.frame(matrix(0,o,8))
gamma_b <- as.data.frame(matrix(0,q + p*(q+1),o))
gamma_beta_MSE <- rep(0,o)
gamma_index <- as.data.frame(matrix(0,o,8))
ptm <- proc.time()
for (i in 1:o) {
  ##generate data
  dat <- GenGE(n, p, q, Gcorr=Gcorr, Ecorr=Ecorr,mis=mis,
               alpha0=0, a=a, b=b, family=family, E.cat=E.cat)
  y <- dat$Y;X <- dat$X;W <- dat$W
  
  b1 <- initial_binomial(y, X, W, family=family, methods='divergence',S_diver="S1")
  res = BIC_sgMCP_mdpd(y, X, W, m_gamma, b1=b1, family=family, Matla=Matla_mdpd)
  beta_mdpd <- res[[1]]
  mdpd_index[i,] <- index_over(beta_mdpd,k_main_nonzero,k_Interaction_nonzero,q,p)
  mdpd_b[,i] <- beta_mdpd
  mdpd_beta_MSE[i] <- norm(beta_mdpd - beta0)
  beta_mdpd_all <- res[[2]]
  
  res = BIC_sgMCP_gamma(y, X, W, gamma, b1=b1, family=family, Matla=Matla_gamma)
  beta_gamma <- res[[1]]
  gamma_index[i,] <- index_over(beta_gamma,k_main_nonzero,k_Interaction_nonzero,q,p)
  gamma_b[,i] <- beta_gamma
  gamma_beta_MSE[i] <- norm(beta_gamma - beta0)
  beta_gamma_all <- res[[2]]
  
  save.image(filename)
  print(paste("!!!replication",i,sep = "-"))
  
}
t1 <- proc.time() - ptm
index_mdpd = c( apply(mdpd_index[,5:8],2,mean),apply(mdpd_index[,5:8],2,sd), mean(mdpd_beta_MSE),sd(mdpd_beta_MSE) )
index_gamma = c( apply(gamma_index[,5:8],2,mean),apply(gamma_index[,5:8],2,sd), mean(gamma_beta_MSE),sd(gamma_beta_MSE) )
index = as.data.frame(rbind(index_mdpd,index_gamma))
names(index) = c(rep("mean",4),rep("sd",4),"mean","sd")
index0 <- trans_table_joint(index, family, gaussian_rownames, binomial_rownames)
index0


