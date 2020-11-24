################################
#Codes for TNBC data joint analysis
################################
rm(list = ls(all = TRUE))
ls()
###############################

#### install brca.data package if needed ####
# library(devtools)
# devtools::install_url('https://github.com/averissimo/brca.data/releases/download/1.0/brca.data_1.0.tar.gz')
library(brca.data)
data("fpkm.per.tissue","fpkm.per.tissue.barcode", 'clinical', package = 'brca.data')
load("TNBC.name.list.RData")
gene.name <- data.name$gene.name
gene.code <- row.names(gene.name)
TNBC.gene.names <- gene.name
sample.name <- data.name$sample.name
sample.name.short <- sample.name
for (i in 1:length(sample.name.short)) {
  a1 <- unlist(strsplit(sample.name.short[i],split='-'))
  sample.name.short[i] <- paste(a1[1],a1[2],a1[3],sep="-")
}

#### Extraction of the dependent variable
tnbc.vars <- c('breast_carcinoma_estrogen_receptor_status',
               'breast_carcinoma_progesterone_receptor_status',
               'her2_immunohistochemistry_level_result',
               'lab_proc_her2_neu_immunohistochemistry_receptor_status',
               'lab_procedure_her2_neu_in_situ_hybrid_outcome_type')
clinical.tnbc <- clinical$primary.solid.tumor[, tnbc.vars]
clinical.tnbc.aux <- clinical.tnbc
names(clinical.tnbc.aux) <- c("ESR1","PGR","HER2_level","HER2_status","HER2_FISH")
clinical.tnbc.aux <- clinical.tnbc.aux[-(which(is.na(clinical.tnbc.aux[,1]))),]
clinical.tnbc.aux$my_HER2_level <- as.character(clinical.tnbc.aux$HER2_level)
clinical.tnbc.aux[which(clinical.tnbc.aux$HER2_level == '0'),6] <- 'Negative'
clinical.tnbc.aux[which(clinical.tnbc.aux$HER2_level == '1+'),6] <- 'Negative'
clinical.tnbc.aux[which(clinical.tnbc.aux$HER2_level == '2+'),6] <- 'Indeterminate'
clinical.tnbc.aux[which(clinical.tnbc.aux$HER2_level == '3+'),6] <- 'Positive'
HER2status.dif.HER2level <- clinical.tnbc.aux[c(rownames(clinical.tnbc.aux[which(clinical.tnbc.aux[,6] == 'Negative' & clinical.tnbc.aux[,4] == 'Positive'),]),rownames(clinical.tnbc.aux[which(clinical.tnbc.aux[,6] == 'Positive' & clinical.tnbc.aux[,4] == 'Negative'),])),]
# getting the individuals short name, as in the clinical data
clinical.tnbc.barcode.pattern <- '(TCGA-[A-Z0-9a-z]{2}-[a-zA-Z0-9]{4})-([0-9]{2}).+'
clinical.tnbc.shortname <- gsub(clinical.tnbc.barcode.pattern,'\\1',rownames(clinical.tnbc))
# individuals with discordant HER2 (IHC) status and HER2 (IHC) level
HER2status.dif.HER2level.names <- clinical.tnbc.shortname[which(clinical.tnbc.shortname %in% rownames(HER2status.dif.HER2level))]
HER2status.dif.FISH <- clinical.tnbc.aux[c(rownames(clinical.tnbc.aux[which(clinical.tnbc.aux[,4] == 'Negative' & clinical.tnbc.aux[,5] == 'Positive'),]),rownames(clinical.tnbc.aux[which(clinical.tnbc.aux[,4] == 'Positive' & clinical.tnbc.aux[,5] == 'Negative'),])),]
HER2status.dif.HER2FISH.names <- clinical.tnbc.shortname[which(clinical.tnbc.shortname %in% rownames( HER2status.dif.FISH))]
suspect <- c(HER2status.dif.HER2level.names,HER2status.dif.HER2FISH.names)
clinical.tnbc.aux$HER2_status_plus_FISH <- clinical.tnbc.aux[,4]
clinical.tnbc.aux[clinical.tnbc.aux[,5] !='',7] <- clinical.tnbc.aux[clinical.tnbc.aux[,5] !='',5]
clinical.tnbc <- clinical.tnbc.aux[,c(1,2,7)]
# replace all '' or NA with 'Indeterminate'
clinical.tnbc[is.na(clinical.tnbc) | clinical.tnbc == '' | clinical.tnbc == 'Equivocal'] <- 'Indeterminate'
# individuals with at least one 'Positive' and keep all are marked as not TNBC
tnbc.status.pos <- clinical.tnbc == 'Positive'
#tnbc.status.pos <- clinical.tnbc[,c(1,2,3)] == 'Positive'
not.tnbc.id <- rownames(tnbc.status.pos[which(sapply(seq(nrow(tnbc.status.pos)), function(ix) { any(tnbc.status.pos[ix,]) })),])
# individuals with at least one 'Indeterminate' and excluded from dataset
tnbc.status.ind <- clinical.tnbc[!(rownames(clinical.tnbc) %in% not.tnbc.id),] == 'Indeterminate'
not.tnbc.id     <- rownames(tnbc.status.ind[which(sapply(seq(nrow(tnbc.status.ind)), function(ix) { any(tnbc.status.ind[ix,]) })),])
clinical.tnbc <- clinical.tnbc[!(rownames(clinical.tnbc) %in% not.tnbc.id), ]
tnbc.status.neg <- clinical.tnbc == 'Negative'
tnbc.ix <- sapply(seq(nrow(clinical.tnbc)), function(ix) { all(tnbc.status.neg[ix,]) })
# set cases with TNBC
clinical.tnbc$tnbc <- 'NO_TNBC'
clinical.tnbc[tnbc.ix, 'tnbc'] <- 'TNBC'
ydata <- data.frame(tnbc = clinical.tnbc$tnbc, row.names = rownames(clinical.tnbc))
ydata$tnbc <- factor(ydata$tnbc)
ydata <- ydata[match(sample.name.short,row.names(ydata)),]
data.Y <- as.matrix(ydata)
data.Y[data.Y =="NO_TNBC"] <- 0
data.Y[data.Y =="TNBC"] <- 1
data.Y <- as.numeric(data.Y)
length(data.Y)

#### Extraction of genetic variables for analysis
xdata <- t(fpkm.per.tissue$primary.solid.tumor)
xdata <- xdata[match(sample.name,row.names(xdata)),match(gene.code,colnames(xdata))]
dim(xdata)
# log-transform
data <- log2(xdata+1)
# normalizing
data <- (data-matrix(apply(data,2,mean),nrow(data),ncol(data),byrow=TRUE))/matrix(apply(data,2,sd),nrow(data),ncol(data),byrow=TRUE)
data <- as.data.frame(data)
row.names(data) <- sample.name
dim(data)

#### Extraction of environmental variables for analysis
clinical.interest <- c('age_at_initial_pathologic_diagnosis','race_list')
TNBC.glm.clinical <- clinical$primary.solid.tumor[,clinical.interest]
TNBC.glm.clinical <- TNBC.glm.clinical[match(sample.name.short,row.names(TNBC.glm.clinical)),]
dim(TNBC.glm.clinical)
#### merge y, G, E
data_GE <- as.data.frame(cbind(data.Y,TNBC.glm.clinical,data))


# coding race variable
race.vect <- as.character(data_GE[,"race_list"])
race.vect[race.vect==""] <- ""
race.vect[race.vect=="AMERICAN INDIAN OR ALASKA NATIVE"] <- "0"
race.vect[race.vect=="WHITE"] <- "0"
race.vect[race.vect=="ASIAN"] <- "0"
race.vect[race.vect=="BLACK OR AFRICAN AMERICAN"] <- "1"
data_GE[,"race_list"] <- as.character(race.vect)
data_GE <- data_GE[data_GE[,"race_list"]!="",]
data_GE <- data_GE[-which(is.na(data_GE[,2])),]
dim(data_GE)
# normalizing age
data_GE[,2] <- (data_GE[,2] - mean(data_GE[,2]))/ sd(data_GE[,2])

#### save data ####
y <- data_GE[,1]
E <- data_GE[,c(2,3)]
G <- data_GE[,-c(1:3)]
data <- list(G=G,E=E,y=y,TNBC.gene.names=TNBC.gene.names)
save(data, file = "TNBC.RData")


