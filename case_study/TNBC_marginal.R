################################
#Codes for TNBC data marginal analysis
################################
rm(list = ls(all = TRUE))
ls()
###############################

library("Rcpp")
library("grpreg")
library("glmnet")
sourceCpp("Rcpp_fun_mdpdlogistic.cpp")
sourceCpp("Rcpp_fun_gammalogistic.cpp")
source("function.R")
# --------------- read data--------------
load("TNBC.RData")
E <- as.data.frame(data$E)
E[,2] <- as.numeric(E[,2])
G <- as.data.frame(data$G)
y <- data$y
TNBC.gene.names <- as.data.frame(data$TNBC.gene.names)
family <- "binomial"
p <- dim(G)[2]
names(G) <- 1:p


# --------- prescreening ------------
P_valueall <- rep(1,length(G[1,]))
for (j in 1:length(G[1,])) {
  x <- G[,j]
  b <- glm(y ~ x,family = binomial)
  P_value <- summary(b)$coefficients
  if (length(P_value[,1]) == 2){
    P_valueall[j] <- P_value[2,4]
  }
  if (j %% 1000 == 0) {
    print(j)
  }
}
sig_n <- which(P_valueall <= sort(P_valueall)[2000])
G <- G[,sig_n]
n <- length(y)
p <- dim(G)[2]
q <- dim(E)[2]

# Standardization of G-E data format
W <- as.data.frame(matrix(0, nrow=n, ncol=p*(q+1)))
X1 <- as.data.frame(cbind(1,E))
for (j in 1:p){
  W[,(1:(q+1))+(q+1)*(j-1)] <- G[,j]*X1
  if(j%%200==0){print(j)}
}



################################
# Codes for marginal analysis of TNBC: identified main effects and interactions.
################################
gen_coe_t <- function(beta, q, G, E, TNBC.gene.names, stab=F){
  # The function for corresponding estimated regression coefficients
  # to identified gene and interaction effects (used in TNBC data)
  nonzero <- which(beta!=0)
  if(length(nonzero)<q+1){
    return(0)
  }
  n_siggene <- nonzero[which(nonzero%%(q+1)==0)]/(q+1)
  siggene <- names(G[,n_siggene])
  gene_coe <- as.data.frame(matrix(0, nrow=length(n_siggene), ncol=q+1))
  for (i in 1:length(n_siggene)) {
    if(stab==T){
      gene_coe[i,] <- beta[(n_siggene[i]*(q+1)):(n_siggene[i]*(q+1)+q),]
    } else{
      gene_coe[i,] <- beta[(n_siggene[i]*(q+1)):(n_siggene[i]*(q+1)+q)]
    }
    
  }
  names(gene_coe) <- c("Main",names(E))
  row.names(gene_coe) <- siggene
  gene_coe[gene_coe==0] <- "NA"
  row.names(gene_coe) <- as.character(TNBC.gene.names[as.numeric(row.names(gene_coe)),])
  return(gene_coe)
}

#-------------gamma-----------
Matla_gamma <- cbind(c(0.050),c(0.050))
gamma <- 0.5
ptm <- proc.time()
res_gamma <- gamma_marginal(as.matrix(y), as.matrix(E), as.matrix(W), gamma = gamma, Matla_gamma, family)
beta_all_gamma <- res_gamma[[1]]
loss_all_gamma <- res_gamma[[2]]
t1 <- proc.time() - ptm
beta_gamma <- mar_TFPR(beta_all_gamma,loss_all_gamma)
gene_coe_gamma <- gen_coe_t(beta_gamma, q, G, E, TNBC.gene.names)
row.names(gene_coe_gamma)


#-------------mdpd-----------
Matla_mdpd <- cbind(c(0.118),c(0.118))
m_gamma <- 0.1
ptm <- proc.time()
res_mdpd <- mdpd_marginal(y, as.matrix(E), as.matrix(W), gamma = m_gamma, Matla_mdpd, family, Armijo=0)
beta_all_mdpd <- res_mdpd[[1]]
loss_all_mdpd <- res_mdpd[[2]]
t1 <- proc.time() - ptm
beta_mdpd <- mar_TFPR(beta_all_mdpd,loss_all_mdpd)
gene_coe_mdpd <- gen_coe_t(beta_mdpd, q, G, E, TNBC.gene.names)
row.names(gene_coe_mdpd)

# ----------------------output Table: marginal analysis of TNBC: identified main effects and interactions----------
coe_mdpd <- as.data.frame(cbind(row.names(gene_coe_mdpd),gene_coe_mdpd))
coe_gamma <- as.data.frame(cbind(row.names(gene_coe_gamma),gene_coe_gamma))
names(coe_mdpd)[1] <- "gene"
names(coe_gamma)[1] <- "gene"
gene_coe <- merge(coe_mdpd,coe_gamma, by = "gene", all = T)
names(gene_coe) <- c("Gene",rep(c("main","age","race"),2))
gene_coe_mdpd_gamma <- gene_coe
write.csv(gene_coe_mdpd_gamma,"gene_coe_gamma_mdpd_TNBC.csv")






################################
# Codes for Joint analysis of TNBC: 
# stability analysis of identified main effects and interactions.
################################
# ----------gamma
o <- 100
ns <- sample(1:n)
beta_gamma_stabality <- as.data.frame(matrix(0, nrow=p*(q+1)+q, ncol=o))
ptm <- proc.time()
for(ttt in 1:o){
  nnn <- ((ttt-1)*round(n/o)+1):(ttt*round(n/o))
  nn <- ns[nnn]
  y_nn <- y[-nn];E_nn <- E[-nn,];W_nn <- W[-nn,]
  res <- gamma_marginal(as.matrix(y_nn), as.matrix(E_nn), as.matrix(W_nn), gamma = gamma, Matla_gamma, family)
  beta_gamma_stabality[,ttt] <- res[[1]]
  # if(ttt %% 5 == 0){save.image("marginal_TNBC_stab.RData")}
}
t3 <- proc.time() - ptm
stab_gamma <- as.data.frame(apply(beta_gamma_stabality,1,function(a){sum(a!=0)/o}))
stabality_gamma <- gen_coe_t(stab_gamma, q, G, E, TNBC.gene.names, stab=T)



# ----------mdpd
beta_mdpd_stabality <- as.data.frame(matrix(0, nrow=p*(q+1)+q, ncol=o))
ptm <- proc.time()
for(ttt in 1:o){
  nnn <- ((ttt-1)*round(n/o)+1):(ttt*round(n/o))
  nn <- ns[nnn]
  y_nn <- y[-nn];E_nn <- E[-nn,];W_nn <- W[-nn,]
  res <- mdpd_marginal(as.matrix(y_nn), as.matrix(E_nn), as.matrix(W_nn), gamma = m_gamma, Matla_mdpd, family, Armijo=0)
  beta_mdpd_stabality[,ttt] <- res[[1]]
  # if(ttt %% 5 == 0){save.image("marginal_TNBC_stab.RData")}
}
t3 <- proc.time() - ptm
stab_mdpd <- as.data.frame(apply(beta_mdpd_stabality,1,function(a){sum(a!=0)/o}))
stabality_mdpd <- gen_coe_t(stab_mdpd, q, G, E, TNBC.gene.names, stab=T)


# ----------------------output Table: marginal analysis of TNBC: stability analysis of identified main effects and interactions-----------------------
stabality_mdpd <- as.data.frame(cbind(row.names(stabality_mdpd),stabality_mdpd))
stabality_gamma <- as.data.frame(cbind(row.names(stabality_gamma),stabality_gamma))
names(stabality_mdpd)[1] <- "gene"
names(stabality_gamma)[1] <- "gene"
stabality <- merge(stabality_mdpd,stabality_gamma, by = "gene", all = T)
names(stabality) <- c("Gene",rep(c("main","age","race"),2))
stabality[,1] <- as.character(stabality[,1])
commongene <- intersect(gene_coe_mdpd_gamma[,1],stabality[,1])
stabality <- stabality[match(commongene,stabality[,1]),]
write.csv(stabality,"stabality_gamma_mdpd_TNBC.csv")