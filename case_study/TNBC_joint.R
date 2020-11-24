################################
#Codes for TNBC data joint analysis
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
group <- c()
group[1:q] <- 1
for (j in 1:p) {
  group[((q+1)*j):((q+1)*j+q)] <- j+1
}




################################
# Codes for Joint analysis of TNBC: identified main effects and interactions.
################################

# ----------gamma---------
ptm <- proc.time()
b1 <- initial_binomial(y, as.matrix(E), as.matrix(W), family=family, methods='divergence')
length(which(b1!=0))
t1 <- proc.time() - ptm
Matla_gamma <- cbind(c(0.025),c(0.025))
gamma <- 0.5
ptm <- proc.time()
res <- BIC_sgMCP_gamma(as.matrix(y), as.matrix(E), as.matrix(W), gamma, b1=b1, family=family, Matla=Matla_gamma)
res_gamma <- res
beta_gamma <- res[[1]]
beta_gamma_all <- res[[2]]
t2 <- proc.time() - ptm
nonzero_gamma <- which(beta_gamma!=0)
n_siggene <- nonzero_gamma[which(nonzero_gamma%%(q+1)==0)]/(q+1)
siggene_gamma <- names(G[,n_siggene])
gene_coe <- as.data.frame(matrix(0, nrow=length(n_siggene), ncol=q+1))
for (i in 1:length(n_siggene)) {
  gene_coe[i,] <- beta_gamma[(n_siggene[i]*(q+1)):(n_siggene[i]*(q+1)+q)]
}
names(gene_coe) <- c("Main",names(E))
row.names(gene_coe) <- siggene_gamma
gene_coe[gene_coe==0] <- "NA"
as.character(TNBC.gene.names[as.numeric(row.names(gene_coe)),])
gene_coe_gamma <- gene_coe
row.names(gene_coe_gamma) <- as.character(TNBC.gene.names[as.numeric(row.names(gene_coe)),])


# ----------mdpd---------
ptm <- proc.time()
b1 <- initial_binomial(y, as.matrix(E), as.matrix(W), family=family, methods='divergence')
length(which(b1!=0))
t1 <- proc.time() - ptm
Matla_mdpd <- cbind(c(0.05),c(0.05))
m_gamma <- 0.1
ptm <- proc.time()
res <- BIC_sgMCP_mdpd(as.matrix(y), as.matrix(E), as.matrix(W), m_gamma, b1=b1, family=family, Matla=Matla_mdpd)
res_mdpd <- res
beta_mdpd <- res[[1]]
beta_mdpd_all <- res[[2]]
t2 <- proc.time() - ptm
nonzero_mdpd <- which(beta_mdpd!=0)
n_siggene <- nonzero_mdpd[which(nonzero_mdpd%%(q+1)==0)]/(q+1)
siggene_mdpd <- names(G[,n_siggene])
gene_coe <- as.data.frame(matrix(0, nrow=length(n_siggene), ncol=q+1))
for (i in 1:length(n_siggene)) {
  gene_coe[i,] <- beta_mdpd[(n_siggene[i]*(q+1)):(n_siggene[i]*(q+1)+q)]
}
names(gene_coe) <- c("Main",names(E))
row.names(gene_coe) <- siggene_mdpd
gene_coe[gene_coe==0] <- "NA"
as.character(TNBC.gene.names[as.numeric(row.names(gene_coe)),])
gene_coe_mdpd <- gene_coe
row.names(gene_coe_mdpd) <- as.character(TNBC.gene.names[as.numeric(row.names(gene_coe)),])

# ----------------------output Table: Joint analysis of TNBC: identified main effects and interactions-----------------------
coe_mdpd <- as.data.frame(cbind(row.names(gene_coe_mdpd),gene_coe_mdpd))
coe_gamma <- as.data.frame(cbind(row.names(gene_coe_gamma),gene_coe_gamma))
names(coe_mdpd)[1] <- "gene"
names(coe_gamma)[1] <- "gene"
gene_coe <- merge(coe_mdpd,coe_gamma, by = "gene", all = T)
names(gene_coe) <- c("Gene",rep(c("main","age","race"),2))
gene_coe_mdpd_gamma <- gene_coe
write.csv(gene_coe_mdpd_gamma,"joint_gene_coe_gamma_mdpd_TNBC.csv")






################################
# Codes for Joint analysis of TNBC: Joint analysis of TNBC data: 
# identified main effects and interactions after removing outliers.
################################

# ----------------gamma-----------------
res <- res_gamma
beta0_gamma <- res[[7]]
data <- cbind(E,W)
gamma_weight <- gammabinomial_weight(y, as.matrix(data), 0.5, as.numeric(beta_gamma), as.numeric(beta0_gamma))
downweight <- which(gamma_weight < 0.5)
delname_gamma <- row.names(E)[downweight]
y_del <- y[-downweight]
E_del <- E[-downweight,]
W_del <- W[-downweight,]




ptm <- proc.time()
b1 <- initial_binomial(y_del, as.matrix(E_del), as.matrix(W_del), family=family, methods='divergence')
t1 <- proc.time() - ptm




ptm <- proc.time()
res_del <- BIC_sgMCP_gamma(as.matrix(y_del), as.matrix(E_del), as.matrix(W_del), gamma, b1=b1, family=family, Matla=Matla_gamma)
beta_gamma_del <- res_del[[1]]
beta_gamma_all_del <- res_del[[2]]
t2 <- proc.time() - ptm




nonzero_gamma_del <- which(beta_gamma_del!=0)
n_siggene_del <- nonzero_gamma_del[which(nonzero_gamma_del%%(q+1)==0)]/(q+1)
siggene_gamma_del <- names(G[,n_siggene_del])
gene_coe_del <- as.data.frame(matrix(0, nrow=length(n_siggene_del), ncol=q+1))
for (i in 1:length(n_siggene_del)) {
  gene_coe_del[i,] <- beta_gamma_del[(n_siggene_del[i]*(q+1)):(n_siggene_del[i]*(q+1)+q)]
}
names(gene_coe_del) <- c("Main",names(E_del))
row.names(gene_coe_del) <- siggene_gamma_del
gene_coe_del[gene_coe_del==0] <- "NA"
as.character(TNBC.gene.names[as.numeric(row.names(gene_coe_del)),])
gene_coe_gamma_del <- gene_coe_del
row.names(gene_coe_gamma_del) <- as.character(TNBC.gene.names[as.numeric(row.names(gene_coe_del)),])
intersect(row.names(gene_coe_gamma_del),row.names(gene_coe_gamma))
# save.image("joint_TNBC_gamma0_delsmallweight.RData")


# ----------------mdpd-----------------
res <- res_mdpd
beta0_mdpd <- res[[7]]
data <- cbind(E,W)
mdpd_weight <- mdpdbinomial_weight(y, as.matrix(data), 0.1, as.numeric(beta_mdpd), as.numeric(beta0_mdpd))
downweight <- which(mdpd_weight < 0.5)
delname_mdpd <- row.names(E)[downweight]
y_del <- y[-downweight]
E_del <- E[-downweight,]
W_del <- W[-downweight,]




ptm <- proc.time()
b1 <- initial_binomial(y_del, as.matrix(E_del), as.matrix(W_del), family=family, methods='divergence')
t1 <- proc.time() - ptm




ptm <- proc.time()
res_del <- BIC_sgMCP_mdpd(as.matrix(y_del), as.matrix(E_del), as.matrix(W_del), m_gamma, b1=b1, family=family, Matla=Matla_mdpd)
beta_mdpd_del <- res_del[[1]]
beta_mdpd_all_del <- res_del[[2]]
t2 <- proc.time() - ptm




nonzero_mdpd_del <- which(beta_mdpd_del!=0)
n_siggene_del <- nonzero_mdpd_del[which(nonzero_mdpd_del%%(q+1)==0)]/(q+1)
siggene_mdpd_del <- names(G[,n_siggene_del])
gene_coe_del <- as.data.frame(matrix(0, nrow=length(n_siggene_del), ncol=q+1))
for (i in 1:length(n_siggene_del)) {
  gene_coe_del[i,] <- beta_mdpd_del[(n_siggene_del[i]*(q+1)):(n_siggene_del[i]*(q+1)+q)]
}
names(gene_coe_del) <- c("Main",names(E_del))
row.names(gene_coe_del) <- siggene_mdpd_del
gene_coe_del[gene_coe_del==0] <- "NA"
as.character(TNBC.gene.names[as.numeric(row.names(gene_coe_del)),])
gene_coe_mdpd_del <- gene_coe_del
row.names(gene_coe_mdpd_del) <- as.character(TNBC.gene.names[as.numeric(row.names(gene_coe_del)),])
intersect(row.names(gene_coe_mdpd_del),row.names(gene_coe_mdpd))
# save.image("joint_TNBC_mdpd0_delsmallweight.RData")

# --output Table: Joint analysis of TNBC data: identified main effects and interactions after removing outliers-----------------------
coe_mdpd_del <- as.data.frame(cbind(row.names(gene_coe_mdpd_del),gene_coe_mdpd_del))
coe_gamma_del <- as.data.frame(cbind(row.names(gene_coe_gamma_del),gene_coe_gamma_del))
a <- as.data.frame(matrix("NA",ncol=4,nrow=length(gene_coe_gamma_del[,1])-length(gene_coe_mdpd_del[,1])))
names(a) <- names(coe_mdpd_del)
row.names(a) <- paste0("ZZZ",1:dim(a)[1])
coe_mdpd_del <- rbind(coe_mdpd_del,a)
mdpd00 <- coe_mdpd_del[order(row.names(coe_mdpd_del)),]
gamma00 <- coe_gamma_del[order(coe_gamma_del[,1]),]
gene_coe_del <- cbind(mdpd00,gamma00)
names(gene_coe_del) <- c(rep(c("Gene","main","age","race"),2))
row.names(gene_coe_del) <- c(1:dim(gene_coe_del)[1])
write.csv(gene_coe_del,"joint_gene_coe_gamma_mdpd_del_TNBC.csv")





################################
# Codes for Joint analysis of TNBC: 
# stability analysis of identified main effects and interactions.
################################

# ----------gamma---------------
o <- 100
beta_gamma_stabality <- as.data.frame(matrix(0, nrow=p*(q+1)+q, ncol=o))
b1 <- initial_binomial(y, as.matrix(E), as.matrix(W), family=family, methods='divergence')
ns <- sample(1:n)
ptm <- proc.time()
for(ttt in 1:o){
  nnn <- ((ttt-1)*floor(n/o)+1):min((ttt*floor(n/o)),n)
  nn <- ns[nnn]
  y_nn <- y[-nn]
  E_nn <- E[-nn,]
  W_nn <- W[-nn,]
  res <- BIC_sgMCP_gamma(as.matrix(y_nn), as.matrix(E_nn), as.matrix(W_nn), gamma, b1=b1, family=family, Matla=Matla_gamma)
  beta_gamma <- res[[1]]
  beta_gamma_stabality[,ttt] <- beta_gamma
  # if(ttt %% 5 == 0){save.image("joint_TNBC_stab.RData")}
}
t3 <- proc.time() - ptm
stab <- as.data.frame(apply(beta_gamma_stabality,1,function(a){sum(a!=0)/o}))
a <- as.data.frame(stab[stab!=0])
stab_gamma <- stab
nonzero <- which(stab!=0)
n_siggene <- nonzero[which(nonzero%%(q+1)==0)]/(q+1)
siggene <- as.character(TNBC.gene.names[as.numeric(names(G[,n_siggene])),])
stabality <- as.data.frame(matrix(0, nrow=length(n_siggene), ncol=q+1))
for (i in 1:length(n_siggene)) {
  stabality[i,] <- stab[(n_siggene[i]*(q+1)):(n_siggene[i]*(q+1)+q),]
}
names(stabality) <- c("Main",names(E))
row.names(stabality) <- siggene
stabality[stabality==0] <- "NA"
stabality_gamma0 <- stabality
# save.image("joint_TNBC_gamma0.RData")


# ----------mdpd-----------------
beta_mdpd_stabality <- as.data.frame(matrix(0, nrow=p*(q+1)+q, ncol=o))
b1 <- initial_binomial(y, as.matrix(E), as.matrix(W), family=family, methods='divergence')
ptm <- proc.time()
for(ttt in 1:o){
  nnn <- ((ttt-1)*floor(n/o)+1):min((ttt*floor(n/o)),n)
  nn <- ns[nnn]
  y_nn <- y[-nn]
  E_nn <- E[-nn,]
  W_nn <- W[-nn,]
  res <- BIC_sgMCP_mdpd(as.matrix(y_nn), as.matrix(E_nn), as.matrix(W_nn), m_gamma, b1=b1, family=family, Matla=Matla_mdpd)
  beta_mdpd <- res[[1]]
  beta_mdpd_stabality[,ttt] <- beta_mdpd
  # if(ttt %% 5 == 0){save.image("joint_TNBC_stab.RData")}
}
t3 <- proc.time() - ptm
stab <- as.data.frame(apply(beta_mdpd_stabality,1,function(a){sum(a!=0)/o}))
a <- as.data.frame(stab[stab!=0])
stab_mdpd <- stab
nonzero <- which(stab!=0)
n_siggene <- nonzero[which(nonzero%%(q+1)==0)]/(q+1)
siggene <- as.character(TNBC.gene.names[as.numeric(names(G[,n_siggene])),])
stabality <- as.data.frame(matrix(0, nrow=length(n_siggene), ncol=q+1))
for (i in 1:length(n_siggene)) {
  stabality[i,] <- stab[(n_siggene[i]*(q+1)):(n_siggene[i]*(q+1)+q),]
}
names(stabality) <- c("Main",names(E))
row.names(stabality) <- siggene
stabality[stabality==0] <- "NA"
stabality_mdpd0 <- stabality
# save.image("joint_TNBC_mdpd0.RData")


# ----------------------output Table: Joint analysis of TNBC: stability analysis of identified main effects and interactions-----------------------
stabality_mdpd <- as.data.frame(cbind(row.names(stabality_mdpd0),stabality_mdpd0))
stabality_gamma <- as.data.frame(cbind(row.names(stabality_gamma0),stabality_gamma0))
names(stabality_mdpd)[1] <- "gene"
names(stabality_gamma)[1] <- "gene"
stabality <- merge(stabality_mdpd,stabality_gamma, by = "gene", all = T)
names(stabality) <- c("Gene",rep(c("main","age","race"),2))
stabality[,1] <- as.character(stabality[,1])
commongene <- intersect(gene_coe_mdpd_gamma[,1],stabality[,1])
stabality <- stabality[match(commongene,stabality[,1]),]
write.csv(stabality,"joint_stabality_gamma_mdpd_commongene_TNBC.csv")