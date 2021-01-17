library(qtl)
library(tidyr)

D1 = read.csv("/powerplant/workspace/cfpmjt/leafy-green-master/Tidy/Adjusted Data.csv")
M1 = read.csv("/powerplant/workspace/cfpmjt/leafy-green-master/Code Peter Jaksons/EC103xEC201_GBS_40clusters_abxaa.csv")

### Creating a Table 
#
#M2 = as.data.frame(t(M1))
#M2[1:10,1:10]
#colnames(M2) <- unlist(M2[1,])
#M2[1:10,1:10]

#
#M2$ID = row.names(M2)
#D2 = D1
#D2$ID = paste0("X",sep = "",D2$Tree.ID)
#
#D2 = subset(D2, Allocation == "F1")
#Pheno = D2[c("ID","AnthAdj")]
#
#MD = merge(Anth,M2, by = "ID", all = TRUE)
#MD <- MD[-1,]
#
#MD[1:10,1:10]
#
#MD <- MD[-2,]
#MD[1,] <- NA
#
#MD <- MD[-1]
#write.csv(MD,"MD.csv", na = "", row.names = FALSE )
#
####Read in cross
#
#DATA = read.cross(format = "csv", file = "MD.csv",na.strings = "--", genotypes = c("aa","ab","bb"), alleles = c("a","b"),estimate.map=TRUE)
#
#DATA[1:10,1:10]
#summary(DATA)

###############Try seperate files and genotypes

#####Tidy genotype file

M3 = M1
M3[1:10,1:10]
M3 =separate(M3,
               col = "Marker",
               into = c("Chr","Marker"),
               sep = "_",
               remove = TRUE)

M3$M <- paste0(M3$Chr,"m", 1:nrow(M3))
M3$M

M3 = as.data.frame(t(M3))
M3[1:10,1:10]

M3$id <- row.names(M3)

M3 <- M3[c(19128,1:19127)]
M3 <- M3[c(184,1:183),]
M3[1:10,1:10]
M3[1,1] <- ""
M3[2,1] <- ""
M3[1:10,1:10]



write.csv(M3, "Gtype.csv",row.names = FALSE)


#####Tidy phenotype file

D2$id = paste0("X",sep = "",D2$Tree.ID)
D2 = subset(D2, Allocation == "F1")
Anth = D2[c("id","AnthAdj")]
names(Anth)[names(Anth) == "AnthAdj"] <- "pheno"

#write.csv(Anth,"Pheno.csv", row.names = FALSE)

# Read in cross 

Data2 <- read.cross("csvs", genfile ="Gtype.csv", phefile = "Pheno.csv",na.strings = "--", genotypes = c("aa","ab","bb"), alleles = c("a","b"),crosstype = F1,estimate.map=TRUE)
summary(Data2)

nind(Data2)
nphe(Data2)
totmar(Data2)
nchr(Data2)
nmar(Data2)
plot(Data2)

Mappy <- pull.map(Data2)
summary(Mappy)

geno.table(Data2)

Dchek <- est.rf(Data2)
checkAlleles(mapthis)
plotRF(Dchek, alternate.chrid = TRUE)

pull.map(Dchek,"S0")

geno.crosstab(Dchek,"S0m10952","S0m11002")



plotRF(Data2)

nm <- est.map(Dchek,error.prob = 0.001)
plot.map(Dchek,nm)

mapthis <- formLinkageGroups(Data2, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
