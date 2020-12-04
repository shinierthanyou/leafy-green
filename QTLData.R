D1 = read.csv("Adjusted Data.csv")
M1 = read.csv("/powerplant/workspace/cfpmjt/leafy-green-master/Code Peter Jaksons/EC103xEC201_GBS_40clusters_abxaa.csv")

# Creating a Table 

M2 = as.data.frame(t(M1))
colnames(M2) <- M2[1,]
M2$ID = row.names(M2)
D2 = D1
D2$ID = paste0("X",sep = "",D2$Tree.ID)

D2 = subset(D2, Allocation == "F1")
Anth = D2[c("ID","AnthAdj")]


MD = merge(Anth,M2, by = "ID", all = TRUE)
MD <- MD[-1]

MD[1:10,1:10]

MD <- MD[-2,]

write.csv(MD,"MD.csv")

#Read in cross

DATA = read.cross(format = "csv", file = "MD.csv",na.strings = "--", genotypes = c("aa","ab","bb"), alleles = c("a","b"),estimate.map=TRUE,)

DATA[1:10,1:10]
summary(DATA)

#Try seperate files and genotypes

#Tidy genotype file

M3 = as.data.frame(t(M1))
M3$id <- row.names(M3)
M3 <- M3[c(19128,1:19127)]

write.csv(M3, "Gtype.csv",row.names = FALSE)

#Tidy phenotype file

D2$id = paste0("X",sep = "",D2$Tree.ID)
D2 = subset(D2, Allocation == "F1")
Anth = D2[c("id","AnthAdj")]
names(Anth)[names(Anth) == "AnthAdj"] <- "pheno"

write.csv(Anth,"Pheno.csv")

# Read in cross 

Data2 <- read.cross("csvs", genfile ="Gtype.csv", phefile = "Pheno.csv",na.strings = "--", genotypes = c("aa","ab","bb"), alleles = c("a","b"),estimate.map=TRUE)
summary(Data2)
