https://www.r-bloggers.com/summer-of-data-science-1-genomic-prediction-machines-sods17/
  
https://github.com/mrtnj/sods17/tree/master/genomic_prediction


https://bookdown.org/rdpeng/RProgDA/working-with-large-datasets.html

#######################################################################################
#                                                                                     #
#                               DATA + FILE MANANGEMENT                               #
#                                                                                     #
#######################################################################################
setwd("c://Biometrics//Palmerston North//David Change//2018-04 Spatial effect analysis")

##### library mgmt #####
library(xlsx)
library(ggplot2)
library(dplyr)
#library(asreml)
#library(asremlPlus)
library(nadiv)
library(gridExtra)
library(knitr)
library(grid)
library(reshape2)
library(gplots)
library(sommer)

##### Setup option #####
my_theme <- theme_bw()+theme(panel.grid = element_blank())
theme_set(my_theme)

options( stringsAsFactors=F,
         max.print = 1000,
         scipen=8,
         digits=4)

##### read data #####
data_D <- read.xlsx("Height phenotype 20 April 2018.xlsx",sheetName = "Phenotypes April 2018",startRow = 1)

snp_D <- read.csv("EC103xEC201_GBS_40clusters_abxaa.csv",header=T)

corr_D <- read.xlsx("Similarity of individuals Joinmap calculated.xlsx",sheetIndex = 2,startRow = 1)

##### data mgmt #####
data_D <- data_D[,1:5]
data_D <- data_D[!is.na(data_D$NA.),]
data_D$id <- paste0("X",data_D$NA.) 
data_D$Rowf <- as.factor(data_D$Row)
data_D$Treef <- as.factor(data_D$Tree)

snp_D[snp_D=="--"] <- NA
snp_D[snp_D=="aa"] <- 0
snp_D[snp_D=="ab"] <- 1
snp_D[snp_D=="bb"] <- 2
snp_D <- snp_D[,-1]
rownames_C <- colnames(snp_D)
snp_M <- data.matrix(snp_D)
snp_M <- t(snp_M)
snp_M <- snp_M[row.names(snp_M) %in% data_D$id,]
A <-  A.mat(snp_M[])
A[1:10,1:10]
colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
heatmap.2(A[1:100,1:100], col = colfunc(100),Colv = "Rowv",revC=TRUE)

# 
# ####=========================================####
# #### random population of 200 lines with 1000 markers
# ####=========================================####
# X <- matrix(rep(0,200*1000),200,1000)
# for (i in 1:200) {
#   X[i,] <- ifelse(runif(1000)<0.5,-1,1)
# }
# 
# A <- A.mat(X)
# 
# ####=========================================####
# #### take a look at the Genomic relationship matrix 
# #### (just a small part)
# ####=========================================####
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap.2(A[1:15,1:15], col = colfunc(100),Colv = "Rowv")
# str(hv)
# (0*0)/1
# 
# var(c(5,6))
# corr_D$Individual1 <- paste0("X",corr_D$Individual1)
# corr_D$Individual2 <- paste0("X",corr_D$Individual2)
# corr_D <- corr_D[corr_D$Individual1 %in% data_D$id,]
# corr_D <- corr_D[corr_D$Individual2 %in% data_D$id,]
# 
# corr_M <- acast(corr_D, Individual1~Individual2, value.var="Similarity")
# corr_M <- data.matrix(corr_M)
# 
# corr_M[1:10,1:10]
# makeSymm <- function(m) {
#   m[lower.tri(m)] <- t(m)[lower.tri(m)]
#   return(m)
# }
# corr_M <- makeSymm(corr_M)
# 
# library(gdata)
# lowerTriangle(corr_M) <- upperTriangle(corr_M) 
# upperTriangle(corr_M) <- lowerTriangle(corr_M)
# A <- corr_M
#######################################################################################
#                                                                                     #
#                                     Analyses                                        #
#                                                                                     #
#######################################################################################
ggplot(data=data_D,aes(x=Row,y=Tree,fill=Height))+
  geom_tile()+  
  scale_fill_gradientn(colours = terrain.colors(10))
citation("sommer")
ggplot(data=data_D,aes(x=Height,group=pruning,fill=pruning))+
  geom_histogram(alpha=0.5)

ggplot(data=data_D,aes(x=Height,group=pruning,fill=pruning))+
  geom_density(alpha=0.5)
Covarrubias-Pazaran G. 2016. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6):1-15.

mix1 <- sommer::mmer(fixed=Height~pruning,
                      
                      random=~vs(id, Gu=A)+ Rowf + Treef + spl2D(Row,Tree),
                      
                      rcov=~units,
                      
                      data=data_D)
summary(mix1)
plot(mix1)
fittedvals <- spatPlots(mix1,row = "Row", range = "Tree")
vm0 <- variogram(mix1, xcoor = "Row", ycoor = "Tree")
plot(vm0$F1)

plot(data_D$Height,mix1$fitted.y)
str(mix1)
data_D$Height_pred <- mix1$u.hat$`g(id)`
  fitted.y

ggplot(data_D,aes(x=Height,y=Height_pred,color=pruning))+
  geom_point()+labs(x="Observed tree height",y="Predicted Relative tree height")

write.csv(data_D,"predicted_Height.csv")

