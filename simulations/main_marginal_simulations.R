##################################################################
# Codes for simulation studies: 
# marginal analysis for logistic regression and linear regression
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



# The row names of the output table
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
m_gamma = 0.1
############### S1 & AR1 & n1 ###############
mis <- mis.seq[2]
Gcorr <- corr.seq[1]
n <- 300


Matla_mdpd <- gen_lambda(0.1,0.005,60)
Matla_gamma <- gen_lambda(0.02,0.005,100)
filename <- paste(paste("marginal_binomial",mis,paste("p",p,sep=""),Gcorr,paste("n",n,sep=""),sep = "_"),"RData",sep = ".")
filename <- paste(path,filename,sep = "")
mdpd_index <- list()
mdpd_AUC_top <- as.data.frame(matrix(0,o,10))
gamma_index <- list()
gamma_AUC_top <- as.data.frame(matrix(0,o,10))
ptm <- proc.time()
for (i in 1:o) {
  ##generate data
  dat <- GenGE(n, p, q, Gcorr=Gcorr, Ecorr=Ecorr,mis=mis,
               alpha0=0, a=a, b=b, family=family, E.cat=E.cat, eta1=eta1)
  y <- dat$Y;X <- dat$X;W <- dat$W
  
  res_mdpd <- mdpd_marginal(y, X, W, gamma = m_gamma, Matla_mdpd, family)
  beta_all_mdpd <- res_mdpd[[1]]
  index_mdpd <- index_over(beta_all_mdpd,k_main_nonzero,k_Interaction_nonzero,q,p)
  index_mdpd <- index_mdpd[apply(!is.na(index_mdpd[,5:8]), 1, sum) == 4,]
  mdpd_At <- index_AUC(TFPR=index_mdpd)
  mdpd_index[[i]] <- index_mdpd
  mdpd_AUC_top[i,] <- mdpd_At


  res_gamma <- gamma_marginal(y, X, W, gamma = gamma, Matla_gamma, family)
  beta_all_gamma <- res_gamma[[1]]
  index_gamma <- index_over(beta_all_gamma,k_main_nonzero,k_Interaction_nonzero,q,p)
  index_gamma <- index_gamma[apply(!is.na(index_gamma[,5:8]), 1, sum) == 4,]
  gamma_At <- index_AUC(TFPR=index_gamma)
  gamma_index[[i]] <- index_gamma
  gamma_AUC_top[i,] <- gamma_At
  
  save.image(filename)
  print(paste("!!!replication",i,sep = "-"))
  
}
t1 <- proc.time() - ptm
AUC_top_mdpd = c( apply(mdpd_AUC_top,2,function(x){mean(x, na.rm = T)}),apply(mdpd_AUC_top,2,function(x){sd(x, na.rm = T)}) )
AUC_top_gamma = c( apply(gamma_AUC_top,2,function(x){mean(x, na.rm = T)}),apply(gamma_AUC_top,2,function(x){sd(x, na.rm = T)}) )
AUC_top = as.data.frame(rbind(AUC_top_mdpd,AUC_top_gamma))
names(AUC_top) = c(c("AUC_main", "AUC_inte", "pAUC0.3_main", "pAUC0.3_inte", 
                     "pAUC0.5_main", "pAUC0.5_inte", "top20_main", "top20_inte","top40_main", "top40_inte"),rep("sd",10))
index0 <- trans_table_marginal(AUC_top, family, gaussian_rownames, binomial_rownames, or=10)
index0

