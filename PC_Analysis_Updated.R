##Lipidomics## PC Analysis

#HC1 Datasets:
#HC1.iPSC = subset(Lipiddata, select = c(1:5)) = 5
#HC1.NI = subset(Lipiddata, select = c(6:8)) = 3
#HC1.DA = subset(Lipiddata, select = c(9:12)) = 4

#HC2 Datasets:
#HC2.iPSC = subset(Lipiddata, select = c(13:16)) = 4
#HC2.NI = subset(Lipiddata, select = c(17:21)) = 5
#HC2.DA = subset(Lipiddata, select = c(22:25)) = 4

#SAD Datasets:
#SAD.iPSC = subset(Lipiddata, select = c(26:29)) = 4
#SAD.NI = subset(Lipiddata, select = c(30:34)) = 5
#SAD.DA = subset(Lipiddata, select = c(35:38)) = 4

#iPSCs = 13
#NI = 13
#DA = 12


#1.Opening files

##Data used for analysis was the PC_iPSC_NI_DA_.csv. 
##Setting working directory (wd) as /Users/monbax/Documents/Science/6._UOW_Science/Lipids_Paper/1._Data/5._Lipidomics/R_Analysis_Lipids

#------#------ : Refers to where you are up to in adjusting this code (you're welcome future Mon)

#Using Working directory on laptop for ease:
setwd("/")
library(ggplot2)
library(gplots)
library(VennDiagram)
library(reshape2)
library(plyr)
install.packages('farver')
install.packages('colorspace')


#Import file:

Lipiddata = PC_iPSC_NI_DA_
Lipidspecies = PC_iPSC_NI_DA

# import the 29 Oct 2020 PC_iPSC_NI_DA
PC.lipid.species = subset(PC_iPSC_NI_DA, select = PC.Species)

#Separating datasets
##Separating out data so that filter works
##Note data sets have N=4,or N=5, a minimum of n=2 should be found (may up this to n=3)


#HC1 Datasets:
#iPSC 5 Ni 3 DA 4
HC1.iPSC = subset(Lipiddata, select = c(1:5))
HC1.NI = subset(Lipiddata, select = c(6:8))
HC1.DA = subset(Lipiddata, select = c(10:13))

#HC2 Datasets:
#iPSC 4 Ni 5 DA 4
HC2.iPSC = subset(Lipiddata, select = c(14:17))
HC2.NI = subset(Lipiddata, select = c(18:22))
HC2.DA = subset(Lipiddata, select = c(23:26))

#SAD Datasets:
#iPSC 4 Ni 5 DA 4
SAD.iPSC = subset(Lipiddata, select = c(27:30))
SAD.NI = subset(Lipiddata, select = c(31:35))
SAD.DA = subset(Lipiddata, select = c(36:39))

#Cell type:
All_iPSC = subset(Lipiddata, select = c(1:5, 14:17, 27:30))
write.csv(All_iPSC, file = "All_iPSC.csv")
All_NI = subset(Lipiddata, select = c(6:8, 18:22, 31:35))
write.csv(All_NI, file = "All_NI.csv")
All_DA = subset(Lipiddata, select = c(10:13, 23:26, 36:39))
write.csv(All_DA, file = "All_DA.csv")


#Cell line:
All_HC1 = subset(Lipiddata, select = c(1:8, 10:13))
write.csv(All_HC1, file = "All_HC1.csv")

All_HC2 = subset(Lipiddata, select = c(14:26))
write.csv(All_HC2, file = "All_HC2.csv")

All_SAD = subset(Lipiddata, select = c(27:39))
write.csv(All_SAD, file = "All_SAD.csv")



#Filtering for species replicates
#HC1 Datasets:
Filtered.HC1.iPSC.Lipids = as.data.frame((HC1.iPSC * ifelse(!(rowSums(HC1.iPSC [1:5]==0)>=2), 1, 0)))
Filtered.HC1.NI.Lipids = as.data.frame((HC1.NI * ifelse(!(rowSums(HC1.NI [1:3]==0)>=2), 1, 0)))
Filtered.HC1.DA.Lipids = as.data.frame((HC1.DA * ifelse(!(rowSums(HC1.DA [1:4]==0)>=2), 1, 0)))

#HC2 Datasets:
Filtered.HC2.iPSC.Lipids = as.data.frame((HC2.iPSC * ifelse(!(rowSums(HC2.iPSC [1:4]==0)>=2), 1, 0)))
Filtered.HC2.NI.Lipids = as.data.frame((HC2.NI * ifelse(!(rowSums(HC2.NI [1:5]==0)>=2), 1, 0)))
Filtered.HC2.DA.Lipids = as.data.frame((HC2.DA * ifelse(!(rowSums(HC2.DA [1:4]==0)>=2), 1, 0)))

#SAD Datasets:
Filtered.SAD.iPSC.Lipids = as.data.frame((SAD.iPSC * ifelse(!(rowSums(SAD.iPSC [1:4]==0)>=2), 1, 0)))
Filtered.SAD.NI.Lipids = as.data.frame((SAD.NI * ifelse(!(rowSums(SAD.NI [1:5]==0)>=2), 1, 0)))
Filtered.SAD.DA.Lipids = as.data.frame((SAD.DA * ifelse(!(rowSums(SAD.DA [1:4]==0)>=2), 1, 0)))

#Putting data back together
#HC1
HC1 = (subset = c(PC.lipid.species[1], Filtered.HC1.iPSC.Lipids, Filtered.HC1.NI.Lipids, Filtered.HC1.DA.Lipids))
write.csv(HC1, file = "HC1.csv") 
#HC2
HC2 = (subset = c(PC.lipid.species[1], Filtered.HC2.iPSC.Lipids, Filtered.HC2.NI.Lipids, Filtered.HC2.DA.Lipids))
write.csv(HC2, file = "HC2.csv") 
#SAD
SAD = (subset = c(PC.lipid.species[1], Filtered.SAD.iPSC.Lipids, Filtered.SAD.NI.Lipids, Filtered.SAD.DA.Lipids))
write.csv(SAD, file = "SAD.csv") 

All_ = (subset = c(PC.lipid.species[1], Filtered.HC1.iPSC.Lipids, Filtered.HC1.NI.Lipids, Filtered.HC1.DA.Lipids, Filtered.HC2.iPSC.Lipids, Filtered.HC2.NI.Lipids, Filtered.HC2.DA.Lipids, Filtered.SAD.iPSC.Lipids, Filtered.SAD.NI.Lipids, Filtered.SAD.DA.Lipids))
All_df = as.data.frame(All_)
All_tf = log(All_df[2:39],2)
write.csv(All_tf, file = "Transformed_Lipids.csv")
All_tf_matrix = as.matrix(All_tf)
All_tf_matrix[!is.finite(All_tf_matrix)] = NA
write.csv(All_tf_matrix, file = "Transformed_Lipids_matrix.csv")


#-------------------------------------------------------------------------------------------------------------------------------------------------
#Came back to the analysis here, so reloading the Transformed_Lipids_matrix.csv", in the file for 4. transformed data, and reloading PC.lipid.species, so I can pop those back on

#trying this:

# import the 29 Oct 2020 PC_iPSC_NI_DA
PC.lipid.species = subset(PC_iPSC_NI_DA, select = PC.Species)
All_tf_matrix = Transformed_Lipids_matrix
View(All_tf_matrix)

Lipids_natransformed = cbind(PC.lipid.species, All_tf_matrix)
View(Lipids_natransformed)


#Now impute - this will make the lower quartile of each replicate NA
Lipids_natransformed[is.na(Lipids_natransformed[,1]),1]=rnorm(n=sum(is.na(Lipids_natransformed[,1])),mean=min(Lipids_natransformed[,1],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,1] <= quantile(Lipids_natransformed[,1], na.rm=T)[1], 1], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,2]),2]=rnorm(n=sum(is.na(Lipids_natransformed[,2])),mean=min(Lipids_natransformed[,2],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,2] <= quantile(Lipids_natransformed[,2], na.rm=T)[2], 2], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,3]),3]=rnorm(n=sum(is.na(Lipids_natransformed[,3])),mean=min(Lipids_natransformed[,3],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,3] <= quantile(Lipids_natransformed[,3], na.rm=T)[2], 3], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,4]),4]=rnorm(n=sum(is.na(Lipids_natransformed[,4])),mean=min(Lipids_natransformed[,4],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,4] <= quantile(Lipids_natransformed[,4], na.rm=T)[2], 4], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,5]),5]=rnorm(n=sum(is.na(Lipids_natransformed[,5])),mean=min(Lipids_natransformed[,5],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,5] <= quantile(Lipids_natransformed[,5], na.rm=T)[2], 5], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,6]),6]=rnorm(n=sum(is.na(Lipids_natransformed[,6])),mean=min(Lipids_natransformed[,6],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,6] <= quantile(Lipids_natransformed[,6], na.rm=T)[2], 6], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,7]),7]=rnorm(n=sum(is.na(Lipids_natransformed[,7])),mean=min(Lipids_natransformed[,7],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,7] <= quantile(Lipids_natransformed[,7], na.rm=T)[2], 7], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,8]),8]=rnorm(n=sum(is.na(Lipids_natransformed[,8])),mean=min(Lipids_natransformed[,8],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,8] <= quantile(Lipids_natransformed[,8], na.rm=T)[2], 8], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,9]),9]=rnorm(n=sum(is.na(Lipids_natransformed[,9])),mean=min(Lipids_natransformed[,9],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,9] <= quantile(Lipids_natransformed[,9], na.rm=T)[2], 9], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,10]),10]=rnorm(n=sum(is.na(Lipids_natransformed[,10])),mean=min(Lipids_natransformed[,10],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,10] <= quantile(Lipids_natransformed[,10], na.rm=T)[2], 10], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,11]),11]=rnorm(n=sum(is.na(Lipids_natransformed[,11])),mean=min(Lipids_natransformed[,11],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,11] <= quantile(Lipids_natransformed[,11], na.rm=T)[2], 11], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,12]),12]=rnorm(n=sum(is.na(Lipids_natransformed[,12])),mean=min(Lipids_natransformed[,12],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,12] <= quantile(Lipids_natransformed[,12], na.rm=T)[2], 12], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,13]),13]=rnorm(n=sum(is.na(Lipids_natransformed[,13])),mean=min(Lipids_natransformed[,13],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,13] <= quantile(Lipids_natransformed[,13], na.rm=T)[2], 13], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,14]),14]=rnorm(n=sum(is.na(Lipids_natransformed[,14])),mean=min(Lipids_natransformed[,14],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,14] <= quantile(Lipids_natransformed[,14], na.rm=T)[2], 14], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,15]),15]=rnorm(n=sum(is.na(Lipids_natransformed[,15])),mean=min(Lipids_natransformed[,15],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,15] <= quantile(Lipids_natransformed[,15], na.rm=T)[2], 15], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,16]),16]=rnorm(n=sum(is.na(Lipids_natransformed[,16])),mean=min(Lipids_natransformed[,16],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,16] <= quantile(Lipids_natransformed[,16], na.rm=T)[2], 16], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,17]),17]=rnorm(n=sum(is.na(Lipids_natransformed[,17])),mean=min(Lipids_natransformed[,17],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,17] <= quantile(Lipids_natransformed[,17], na.rm=T)[2], 17], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,18]),18]=rnorm(n=sum(is.na(Lipids_natransformed[,18])),mean=min(Lipids_natransformed[,18],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,18] <= quantile(Lipids_natransformed[,18], na.rm=T)[2], 18], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,19]),19]=rnorm(n=sum(is.na(Lipids_natransformed[,19])),mean=min(Lipids_natransformed[,19],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,19] <= quantile(Lipids_natransformed[,19], na.rm=T)[2], 19], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,20]),20]=rnorm(n=sum(is.na(Lipids_natransformed[,20])),mean=min(Lipids_natransformed[,20],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,20] <= quantile(Lipids_natransformed[,20], na.rm=T)[2], 20], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,21]),21]=rnorm(n=sum(is.na(Lipids_natransformed[,21])),mean=min(Lipids_natransformed[,21],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,21] <= quantile(Lipids_natransformed[,21], na.rm=T)[2], 21], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,22]),22]=rnorm(n=sum(is.na(Lipids_natransformed[,22])),mean=min(Lipids_natransformed[,22],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,22] <= quantile(Lipids_natransformed[,22], na.rm=T)[2], 22], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,23]),23]=rnorm(n=sum(is.na(Lipids_natransformed[,23])),mean=min(Lipids_natransformed[,23],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,23] <= quantile(Lipids_natransformed[,23], na.rm=T)[2], 23], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,24]),24]=rnorm(n=sum(is.na(Lipids_natransformed[,24])),mean=min(Lipids_natransformed[,24],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,24] <= quantile(Lipids_natransformed[,24], na.rm=T)[2], 24], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,25]),25]=rnorm(n=sum(is.na(Lipids_natransformed[,25])),mean=min(Lipids_natransformed[,25],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,25] <= quantile(Lipids_natransformed[,25], na.rm=T)[2], 25], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,26]),26]=rnorm(n=sum(is.na(Lipids_natransformed[,26])),mean=min(Lipids_natransformed[,26],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,26] <= quantile(Lipids_natransformed[,26], na.rm=T)[2], 26], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,27]),27]=rnorm(n=sum(is.na(Lipids_natransformed[,27])),mean=min(Lipids_natransformed[,27],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,27] <= quantile(Lipids_natransformed[,27], na.rm=T)[2], 27], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,28]),28]=rnorm(n=sum(is.na(Lipids_natransformed[,28])),mean=min(Lipids_natransformed[,28],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,28] <= quantile(Lipids_natransformed[,28], na.rm=T)[2], 28], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,29]),29]=rnorm(n=sum(is.na(Lipids_natransformed[,29])),mean=min(Lipids_natransformed[,29],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,29] <= quantile(Lipids_natransformed[,29], na.rm=T)[2], 29], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,30]),30]=rnorm(n=sum(is.na(Lipids_natransformed[,30])),mean=min(Lipids_natransformed[,30],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,30] <= quantile(Lipids_natransformed[,30], na.rm=T)[2], 30], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,31]),31]=rnorm(n=sum(is.na(Lipids_natransformed[,31])),mean=min(Lipids_natransformed[,31],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,31] <= quantile(Lipids_natransformed[,31], na.rm=T)[2], 31], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,32]),32]=rnorm(n=sum(is.na(Lipids_natransformed[,32])),mean=min(Lipids_natransformed[,32],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,32] <= quantile(Lipids_natransformed[,32], na.rm=T)[2], 32], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,33]),33]=rnorm(n=sum(is.na(Lipids_natransformed[,33])),mean=min(Lipids_natransformed[,33],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,33] <= quantile(Lipids_natransformed[,33], na.rm=T)[2], 33], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,34]),34]=rnorm(n=sum(is.na(Lipids_natransformed[,34])),mean=min(Lipids_natransformed[,34],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,34] <= quantile(Lipids_natransformed[,34], na.rm=T)[2], 34], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,35]),35]=rnorm(n=sum(is.na(Lipids_natransformed[,35])),mean=min(Lipids_natransformed[,35],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,35] <= quantile(Lipids_natransformed[,35], na.rm=T)[2], 35], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,36]),36]=rnorm(n=sum(is.na(Lipids_natransformed[,36])),mean=min(Lipids_natransformed[,36],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,36] <= quantile(Lipids_natransformed[,36], na.rm=T)[2], 36], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,37]),37]=rnorm(n=sum(is.na(Lipids_natransformed[,37])),mean=min(Lipids_natransformed[,37],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,37] <= quantile(Lipids_natransformed[,37], na.rm=T)[2], 37], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,38]),38]=rnorm(n=sum(is.na(Lipids_natransformed[,38])),mean=min(Lipids_natransformed[,38],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,38] <= quantile(Lipids_natransformed[,38], na.rm=T)[2], 38], na.rm=T))
Lipids_natransformed[is.na(Lipids_natransformed[,39]),39]=rnorm(n=sum(is.na(Lipids_natransformed[,39])),mean=min(Lipids_natransformed[,39],na.rm=T),sd=sd(Lipids_natransformed[Lipids_natransformed[,39] <= quantile(Lipids_natransformed[,39], na.rm=T)[2], 39], na.rm=T))

View(Lipids_natransformed)
write.csv(Lipids_natransformed, "Lipids_imputed.csv")
Lipids_Imputed = Lipids_natransformed
#Note that this assigns a value to missing data for stats; if you redo this, you will get different newly assigned values

#Remove the Lipid Species Column:
Lipids_imputed_data_only = subset(Lipids_Imputed, select = -c(PC.Species))
View(Lipids_imputed_data_only)
write.csv(Lipids_imputed_data_only, "Lipids_imputed_data_only.csv")

Lipids_Imputed_dataframe = as.data.frame(Lipids_imputed_data_only)
View(Lipids_Imputed_dataframe)
write.csv(Lipids_imputed_data_only, "Lipids_imputed_data_only.csv")

#Scaling the imputed data **I'm not entirely convinced this needs to be done
Lipids_imputed_data_only_matrix = as.matrix(Lipids_imputed_data_only)
View(Lipids_imputed_data_only_matrix)


#these won't seem to work - saying that the centre is an unused argument:
#Lipids_imputed_scaled  = scale(Lipids_imputed_data_only, centre=TRUE,scale=TRUE)
#Lipids_imputed_scaled = scale(Lipids_imputed_data_only_matrix, centre=FALSE, scale=TRUE)


#going ahead without scaling:

#Splitting up the imputed data for future calcs:

Lipids_Imputed = Lipids_imputed
View(Lipids_Imputed)
Lipids_Imputed = as.matrix(Lipids_Imputed)
Lipids_Imputed = subset(Lipids_Imputed, select =  c(2:39))

#Trying to Labelling the columns before separating
Lipids_Imputed = Lipids_imputed
Lipids_Imputed = subset(Lipids_Imputed, select =  c(2:39))
colnames(Lipids_Imputed)  = c("HC1 iPSC 1", "HC1 iPSC 2", "HC1 iPSC 3", "HC1 iPSC 4", "HC1 iPSC 5", "HC1 NI 1", "HC1 NI 2", "HC1 NI 3", "HC1 DA 1", "HC1 DA 2", "HC1 DA 3", "HC1 DA 4", "HC2 iPSC 1", "HC2 iPSC 2", "HC2 iPSC 3", "HC2 iPSC 4", "HC2 NI 1", "HC2 NI 2", "HC2 NI 3", "HC2 NI 4", "HC2 NI 5", "HC2 DA 1", "HC2 DA 2", "HC2 DA 3", "HC2 DA 4", "SAD iPSC 1", "SAD iPSC 2", "SAD iPSC 3", "SAD iPSC 4", "SAD NI 1", "SAD NI 2", "SAD NI 3", "SAD NI 4", "SAD NI 5" ,"SAD DA 1", "SAD DA 2", "SAD DA 3", "SAD DA 4")
Lipids_Imputed = as.matrix(Lipids_Imputed)
View(Lipids_Imputed)
write.csv(Lipids_Imputed, file = "Lipids_Imputed_with_headings.csv")

#wrote  .csvs for with these headings and without (just needed to run the code again without the colnames )

#Checking the first and last to make sure it's writing properly (had issues with the X1_1 in the previous import)
#first
HC1_iPSC_imputed = subset(Lipids_Imputed, select = c(1:5))
View(HC1_iPSC_imputed)
write.csv(HC1_iPSC_imputed, file = "HC1_iPSC_imputed.csv")
#last
SAD_DA_imputed = subset(Lipids_Imputed, select = c(35:38))
View(SAD_DA_imputed)
write.csv(SAD_DA_imputed, file = "SAD_DA_imputed.csv")

#HC1 Datasetsc#iPSC 5 Ni 3 DA 4
HC1_iPSC_imputed = subset(Lipids_Imputed, select = c(1:5))
write.csv(HC1_iPSC_imputed, file = "HC1_iPSC_imputed.csv")
HC1_NI_imputed = subset(Lipids_Imputed, select = c(6:8))
write.csv(HC1_NI_imputed , file = "HC1_NI_imputed.csv")
HC1_DA_imputed = subset(Lipids_Imputed, select = c(9:12))
write.csv(HC1_DA_imputed, file = "HC1_DA_imputed.csv")

#HC2 Datasets: iPSC 4 Ni 5 DA 4
HC2_iPSC_imputed = subset(Lipids_Imputed, select = c(13:16))
write.csv(HC2_iPSC_imputed, file = "HC2_iPSC_imputed.csv")
HC2_NI_imputed = subset(Lipids_Imputed, select = c(17:21))
write.csv(HC2_NI_imputed , file = "HC2_NI_imputed.csv")
HC2_DA_imputed = subset(Lipids_Imputed, select = c(22:25))
write.csv(HC2_DA_imputed, file = "HC2_DA_imputed.csv")

#SAD Datasets: #iPSC 4 Ni 5 DA 4
SAD_iPSC_imputed = subset(Lipids_Imputed, select = c(26:29))
write.csv(SAD_iPSC_imputed, file = "SAD_iPSC_imputed.csv")
SAD_NI_imputed = subset(Lipids_Imputed, select = c(30:34))
write.csv(SAD_NI_imputed , file = "SAD_NI_imputed.csv")
SAD_DA_imputed = subset(Lipids_Imputed, select = c(35:38))
write.csv(SAD_DA_imputed, file = "SAD_DA_imputed.csv")


#Cell type:
All_iPSC_imputed = subset(Lipids_Imputed, select = c(1:5, 13:16, 26:29))
write.csv(All_iPSC_imputed, file = "All_iPSC_imputed.csv")
All_NI_imputed = subset(Lipids_Imputed, select = c(6:8, 17:21, 30:34))
write.csv(All_NI_imputed , file = "All_NI_imputed.csv")
All_DA_imputed = subset(Lipids_Imputed, select = c(9:12, 22:25, 35:38))
write.csv(All_DA_imputed, file = "All_DA_imputed.csv")


##Calculating Standard deviations ***May not need this**
HC1_iPSC_sd = subset(HC1_iPSC_imputed)
HC1_iPSC_sd$HC1_iPSC_sd = apply(HC1_iPSC_sd, 1, sd)
write.csv(HC1_iPSC_sd, file="HC1_iPSC_sd.csv")
HC1_NI_sd = subset(HC1_NI_imputed)
HC1_NI_sd$HC1_NI_sd = apply(HC1_NI_sd, 1, sd)
write.csv(HC1_NI_sd, file="HC1_NI_sd.csv")
HC1_DA_sd = subset(HC1_DA_imputed)
HC1_DA_sd$HC1_DA_sd = apply(HC1_DA_sd, 1, sd)
write.csv(HC1_DA_sd, file="HC1_DA_sd.csv")

HC2_iPSC_sd = subset(HC2_iPSC_imputed)
HC2_iPSC_sd$HC2_iPSC_sd = apply(HC2_iPSC_sd, 1, sd)
write.csv(HC2_iPSC_sd, file="HC2_iPSC_sd.csv")
HC2_NI_sd = subset(HC2_NI_imputed)
HC2_NI_sd$HC2_NI_sd = apply(HC2_NI_sd, 1, sd)
write.csv(HC2_NI_sd, file="HC2_NI_sd.csv")
HC2_DA_sd = subset(HC2_DA_imputed)
HC2_DA_sd$HC2_DA_sd = apply(HC2_DA_sd, 1, sd)
write.csv(HC2_DA_sd, file="HC2_DA_sd.csv")

SAD_iPSC_sd = subset(SAD_iPSC_imputed)
SAD_iPSC_sd$SAD_iPSC_sd = apply(SAD_iPSC_sd, 1, sd)
write.csv(SAD_iPSC_sd, file="SAD_iPSC_sd.csv")
SAD_NI_sd = subset(SAD_NI_imputed)
SAD_NI_sd$SAD_NI_sd = apply(SAD_NI_sd, 1, sd)
write.csv(SAD_NI_sd, file="SAD_NI_sd.csv")
SAD_DA_sd = subset(SAD_DA_imputed)
SAD_DA_sd$SAD_DA_sd = apply(SAD_DA_sd, 1, sd)
write.csv(SAD_DA_sd, file="SAD_DA_sd.csv")

#Sd values of the cell types
All_iPSC_sd = subset(All_iPSC_imputed)
All_iPSC_sd $All_iPSC_sd = apply(All_iPSC_sd, 1, sd)
write.csv(All_iPSC_sd, file=" All_iPSC_sd.csv")
All_NI_sd = subset(All_NI_imputed)
All_NI_sd $All_NI_sd = apply(All_NI_sd, 1, sd)
write.csv(All_NI_sd, file=" All_NI_sd.csv")
All_DA_sd = subset(All_DA_imputed)
All_DA_sd $All_DA_sd = apply(All_DA_sd, 1, sd)
write.csv(All_DA_sd, file=" All_DA_sd.csv")

## Calculating the mean values
HC1_iPSC_mean = subset(HC1_iPSC_imputed)
HC1_iPSC_mean = data.frame(HC1_iPSC_imputed)
HC1_iPSC_mean$HC1_iPSC_mean = apply(HC1_iPSC_mean, 1, mean)
write.csv(HC1_iPSC_mean, file="HC1_iPSC_mean.csv")
HC1_NI_mean = subset(HC1_NI_imputed)
HC1_NI_mean = data.frame(HC1_NI_imputed)
HC1_NI_mean$HC1_NI_mean = apply(HC1_NI_mean, 1, mean)
write.csv(HC1_NI_mean, file="HC1_NI_mean.csv")
HC1_DA_mean = subset(HC1_DA_imputed)
HC1_DA_mean = data.frame(HC1_DA_imputed)
HC1_DA_mean$HC1_DA_mean = apply(HC1_DA_mean, 1, mean)
write.csv(HC1_DA_mean, file="HC1_DA_mean.csv")


HC2_iPSC_mean = subset(HC2_iPSC_imputed)
HC2_iPSC_mean = data.frame(HC2_iPSC_imputed)
HC2_iPSC_mean$HC2_iPSC_mean = apply(HC2_iPSC_mean, 1, mean)
write.csv(HC2_iPSC_mean, file="HC2_iPSC_mean.csv")
HC2_NI_mean = subset(HC2_NI_imputed)
HC2_NI_mean = data.frame(HC2_NI_imputed)
HC2_NI_mean$HC2_NI_mean = apply(HC2_NI_mean, 1, mean)
write.csv(HC2_NI_mean, file="HC2_NI_mean.csv")
HC2_DA_mean = subset(HC2_DA_imputed)
HC2_DA_mean = data.frame(HC2_DA_imputed)
HC2_DA_mean$HC2_DA_mean = apply(HC2_DA_mean, 1, mean)
write.csv(HC2_DA_mean, file="HC2_DA_mean.csv")

SAD_iPSC_mean = subset(SAD_iPSC_imputed)
SAD_iPSC_mean = data.frame(SAD_iPSC_imputed)
SAD_iPSC_mean$SAD_iPSC_mean = apply(SAD_iPSC_mean, 1, mean)
write.csv(SAD_iPSC_mean, file="SAD_iPSC_mean.csv")
SAD_NI_mean = subset(SAD_NI_imputed)
SAD_NI_mean = data.frame(SAD_NI_imputed)
SAD_NI_mean$SAD_NI_mean = apply(SAD_NI_mean, 1, mean)
write.csv(SAD_NI_mean, file="SAD_NI_mean.csv")
SAD_DA_mean = subset(SAD_DA_imputed)
SAD_DA_mean = data.frame(SAD_DA_imputed)
SAD_DA_mean$SAD_DA_mean = apply(SAD_DA_mean, 1, mean)
write.csv(SAD_DA_mean, file="SAD_DA_mean.csv")



#Mean values of the cell types
All_iPSC_mean = subset(All_iPSC_imputed)
All_iPSC_mean = data.frame(All_iPSC_imputed)
All_iPSC_mean $All_iPSC_mean = apply(All_iPSC_mean, 1, mean)
write.csv(All_iPSC_mean, file=" All_iPSC_mean.csv")
All_NI_mean = subset(All_NI_imputed)
All_NI_mean = data.frame(All_NI_imputed)
All_NI_mean $All_NI_mean = apply(All_NI_mean, 1, mean)
write.csv(All_NI_mean, file=" All_NI_mean.csv")
All_DA_mean = subset(All_DA_imputed)
All_DA_mean = data.frame(All_DA_imputed)
All_DA_mean $All_DA_mean = apply(All_DA_mean, 1, mean)
write.csv(All_DA_mean, file=" All_DA_mean.csv")



## Normalising the data to the mean mean of the other two differentiation values (x rep1-(y+z/2))
##This is the comparison of the cell lines
Normalised.HC1.iPSC = as.data.frame(HC1_iPSC_mean)
Normalised.HC1.iPSC$Scaled.iPSC.1 = (Normalised.HC1.iPSC[1] - ((HC1_NI_mean$HC1_NI_mean + HC1_DA_mean$HC1_DA_mean)/2))
Normalised.HC1.iPSC$Scaled.iPSC.2 = (Normalised.HC1.iPSC[2] - ((HC1_NI_mean$HC1_NI_mean + HC1_DA_mean$HC1_DA_mean)/2))
Normalised.HC1.iPSC$Scaled.iPSC.3 = (Normalised.HC1.iPSC[3] - ((HC1_NI_mean$HC1_NI_mean + HC1_DA_mean$HC1_DA_mean)/2))
Normalised.HC1.iPSC$Scaled.iPSC.4 = (Normalised.HC1.iPSC[4] - ((HC1_NI_mean$HC1_NI_mean + HC1_DA_mean$HC1_DA_mean)/2))
Normalised.HC1.iPSC$Scaled.iPSC.5 = (Normalised.HC1.iPSC[5] - ((HC1_NI_mean$HC1_NI_mean + HC1_DA_mean$HC1_DA_mean)/2))

Normalised.HC1.NI = as.data.frame(HC1_NI_mean)
Normalised.HC1.NI$Scaled.NI.1 = (Normalised.HC1.NI[1] - ((HC1_iPSC_mean$HC1_iPSC_mean + HC1_DA_mean$HC1_DA_mean)/2))
Normalised.HC1.NI$Scaled.NI.2 = (Normalised.HC1.NI[2] - ((HC1_iPSC_mean$HC1_iPSC_mean + HC1_DA_mean$HC1_DA_mean)/2))
Normalised.HC1.NI$Scaled.NI.3 = (Normalised.HC1.NI[3] - ((HC1_iPSC_mean$HC1_iPSC_mean + HC1_DA_mean$HC1_DA_mean)/2))

Normalised.HC1.DA = as.data.frame(HC1_DA_mean)
Normalised.HC1.DA$Scaled.DA.1 = (Normalised.HC1.DA[1] - ((HC1_iPSC_mean$HC1_iPSC_mean + HC1_NI_mean$HC1_NI_mean)/2))
Normalised.HC1.DA$Scaled.DA.2 = (Normalised.HC1.DA[2] - ((HC1_iPSC_mean$HC1_iPSC_mean + HC1_NI_mean$HC1_NI_mean)/2))
Normalised.HC1.DA$Scaled.DA.3 = (Normalised.HC1.DA[3] - ((HC1_iPSC_mean$HC1_iPSC_mean + HC1_NI_mean$HC1_NI_mean)/2))
Normalised.HC1.DA$Scaled.DA.4 = (Normalised.HC1.DA[4] - ((HC1_iPSC_mean$HC1_iPSC_mean + HC1_NI_mean$HC1_NI_mean)/2))


Normalised.HC2.iPSC = as.data.frame(HC2_iPSC_mean)
Normalised.HC2.iPSC$Scaled.iPSC.1 = (Normalised.HC2.iPSC[1] - ((HC2_NI_mean$HC2_NI_mean + HC2_DA_mean$HC2_DA_mean)/2))
Normalised.HC2.iPSC$Scaled.iPSC.2 = (Normalised.HC2.iPSC[2] - ((HC2_NI_mean$HC2_NI_mean + HC2_DA_mean$HC2_DA_mean)/2))
Normalised.HC2.iPSC$Scaled.iPSC.3 = (Normalised.HC2.iPSC[3] - ((HC2_NI_mean$HC2_NI_mean + HC2_DA_mean$HC2_DA_mean)/2))
Normalised.HC2.iPSC$Scaled.iPSC.4 = (Normalised.HC2.iPSC[4] - ((HC2_NI_mean$HC2_NI_mean + HC2_DA_mean$HC2_DA_mean)/2))

Normalised.HC2.NI = as.data.frame(HC2_NI_mean)
Normalised.HC2.NI$Scaled.NI.1 = (Normalised.HC2.NI[1] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_DA_mean$HC2_DA_mean)/2))
Normalised.HC2.NI$Scaled.NI.2 = (Normalised.HC2.NI[2] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_DA_mean$HC2_DA_mean)/2))
Normalised.HC2.NI$Scaled.NI.3 = (Normalised.HC2.NI[3] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_DA_mean$HC2_DA_mean)/2))
Normalised.HC2.NI$Scaled.NI.4 = (Normalised.HC2.NI[4] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_DA_mean$HC2_DA_mean)/2))
Normalised.HC2.NI$Scaled.NI.5 = (Normalised.HC2.NI[5] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_DA_mean$HC2_DA_mean)/2))

Normalised.HC2.DA = as.data.frame(HC2_DA_mean)
Normalised.HC2.DA$Scaled.DA.1 = (Normalised.HC2.DA[1] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_NI_mean$HC2_NI_mean)/2))
Normalised.HC2.DA$Scaled.DA.2 = (Normalised.HC2.DA[2] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_NI_mean$HC2_NI_mean)/2))
Normalised.HC2.DA$Scaled.DA.3 = (Normalised.HC2.DA[3] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_NI_mean$HC2_NI_mean)/2))
Normalised.HC2.DA$Scaled.DA.4 = (Normalised.HC2.DA[4] - ((HC2_iPSC_mean$HC2_iPSC_mean + HC2_NI_mean$HC2_NI_mean)/2))


Normalised.SAD.iPSC = as.data.frame(SAD_iPSC_mean)
Normalised.SAD.iPSC$Scaled.iPSC.1 = (Normalised.SAD.iPSC[1] - ((SAD_NI_mean$SAD_NI_mean + SAD_DA_mean$SAD_DA_mean)/2))
Normalised.SAD.iPSC$Scaled.iPSC.2 = (Normalised.SAD.iPSC[2] - ((SAD_NI_mean$SAD_NI_mean + SAD_DA_mean$SAD_DA_mean)/2))
Normalised.SAD.iPSC$Scaled.iPSC.3 = (Normalised.SAD.iPSC[3] - ((SAD_NI_mean$SAD_NI_mean + SAD_DA_mean$SAD_DA_mean)/2))
Normalised.SAD.iPSC$Scaled.iPSC.4 = (Normalised.SAD.iPSC[4] - ((SAD_NI_mean$SAD_NI_mean + SAD_DA_mean$SAD_DA_mean)/2))

Normalised.SAD.NI = as.data.frame(SAD_NI_mean)
Normalised.SAD.NI$Scaled.NI.1 = (Normalised.SAD.NI[1] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_DA_mean$SAD_DA_mean)/2))
Normalised.SAD.NI$Scaled.NI.2 = (Normalised.SAD.NI[2] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_DA_mean$SAD_DA_mean)/2))
Normalised.SAD.NI$Scaled.NI.3 = (Normalised.SAD.NI[3] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_DA_mean$SAD_DA_mean)/2))
Normalised.SAD.NI$Scaled.NI.4 = (Normalised.SAD.NI[4] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_DA_mean$SAD_DA_mean)/2))
Normalised.SAD.NI$Scaled.NI.5 = (Normalised.SAD.NI[5] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_DA_mean$SAD_DA_mean)/2))

Normalised.SAD.DA = as.data.frame(SAD_DA_mean)
Normalised.SAD.DA$Scaled.DA.1 = (Normalised.SAD.DA[1] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_NI_mean$SAD_NI_mean)/2))
Normalised.SAD.DA$Scaled.DA.2 = (Normalised.SAD.DA[2] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_NI_mean$SAD_NI_mean)/2))
Normalised.SAD.DA$Scaled.DA.3 = (Normalised.SAD.DA[3] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_NI_mean$SAD_NI_mean)/2))
Normalised.SAD.DA$Scaled.DA.4 = (Normalised.SAD.DA[4] - ((SAD_iPSC_mean$SAD_iPSC_mean + SAD_NI_mean$SAD_NI_mean)/2))


##This is the comparison of the cell types

Normalised.All.iPSC = as.data.frame(All_iPSC_mean)
Normalised.All.iPSC$Scaled.iPSC.1 = (Normalised.All.iPSC[1] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.2 = (Normalised.All.iPSC[2] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.3 = (Normalised.All.iPSC[3] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.4 = (Normalised.All.iPSC[4] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.5 = (Normalised.All.iPSC[5] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.6 = (Normalised.All.iPSC[6] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.7 = (Normalised.All.iPSC[7] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.8 = (Normalised.All.iPSC[8] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.9 = (Normalised.All.iPSC[9] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.10 = (Normalised.All.iPSC[10] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.11 = (Normalised.All.iPSC[11] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.12 = (Normalised.All.iPSC[12] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.iPSC$Scaled.iPSC.13 = (Normalised.All.iPSC[13] - ((Normalised.All.NI$All_NI_mean + All_DA_mean$All_DA_mean)/2))

Normalised.All.NI = as.data.frame(All_NI_mean)
Normalised.All.NI$Scaled.NI.1 = (Normalised.All.NI[1] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.2 = (Normalised.All.NI[2] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.3 = (Normalised.All.NI[3] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.4 = (Normalised.All.NI[4] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.5 = (Normalised.All.NI[5] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.6 = (Normalised.All.NI[6] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.7 = (Normalised.All.NI[7] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.8 = (Normalised.All.NI[8] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.9 = (Normalised.All.NI[9] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.10 = (Normalised.All.NI[10] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.11 = (Normalised.All.NI[11] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.12 = (Normalised.All.NI[12] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))
Normalised.All.NI$Scaled.NI.13 = (Normalised.All.NI[13] - ((Normalised.All.iPSC$All_iPSC_mean + All_DA_mean$All_DA_mean)/2))

Normalised.All.DA = as.data.frame(All_DA_mean)
Normalised.All.DA$Scaled.DA.1 = (Normalised.All.DA[1] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.2 = (Normalised.All.DA[2] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.3 = (Normalised.All.DA[3] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.4 = (Normalised.All.DA[4] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.5 = (Normalised.All.DA[5] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.6 = (Normalised.All.DA[6] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.7 = (Normalised.All.DA[7] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.8 = (Normalised.All.DA[8] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.9 = (Normalised.All.DA[9] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.10 = (Normalised.All.DA[10] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.11 = (Normalised.All.DA[11] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))
Normalised.All.DA$Scaled.DA.12 = (Normalised.All.DA[12] - ((Normalised.All.NI$All_NI_mean + All_iPSC_mean$All_iPSC_mean)/2))

####----####----####----####----####----####----####----####----####----####----####----####----####----####----####----

#T-tests
#binding together the normalised values to compare
#T- test to get p value
#adjusted p value
#mean of the normalised data
#reverse transorm
#fold change
#log of the fold change
#-log10 of the adjusted p value
#thresholds of 2 fold change and 0.01 adjusted P value

View(T.Data.HC1.iPSC.NI)
T.Data.HC1.iPSC.NI = cbind(Normalised.HC1.iPSC[7:11], Normalised.HC1.NI[5:7])
T.Data.HC1.iPSC.NI$P.value = apply(T.Data.HC1.iPSC.NI, 1, function(T.Data.HC1.iPSC.NI){t.test (x = T.Data.HC1.iPSC.NI[1:5], y = T.Data.HC1.iPSC.NI[6:8], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.HC1.iPSC.NI$adj.pvalues = p.adjust(T.Data.HC1.iPSC.NI$P.value, method = "BH", n=length(T.Data.HC1.iPSC.NI$P.value))
T.Data.HC1.iPSC.NI$normalised.mean.iPSC = apply(T.Data.HC1.iPSC.NI[1:5], 1, mean)
T.Data.HC1.iPSC.NI$Abundance.iPSC = (2^T.Data.HC1.iPSC.NI["normalised.mean.iPSC"])
T.Data.HC1.iPSC.NI$normalised.mean.NI = apply(T.Data.HC1.iPSC.NI[6:8], 1, mean)
T.Data.HC1.iPSC.NI$Abundance.NI = (2^T.Data.HC1.iPSC.NI["normalised.mean.NI"])
T.Data.HC1.iPSC.NI$Fold.change = (T.Data.HC1.iPSC.NI$Abundance.iPSC/T.Data.HC1.iPSC.NI$Abundance.NI)
T.Data.HC1.iPSC.NI$Log.Fold.change = log2(T.Data.HC1.iPSC.NI[,"Fold.change"])
T.Data.HC1.iPSC.NI$neg.log10.adj.pvalue = -log10(T.Data.HC1.iPSC.NI [,"adj.pvalues"])
T.Data.HC1.iPSC.NI$Threshold = as.factor(T.Data.HC1.iPSC.NI[,"Log.Fold.change"]>=1.5 & T.Data.HC1.iPSC.NI[,"neg.log10.adj.pvalue"] >1.3)
T.Data.HC1.iPSC.NI$Threshold.neg = as.factor(T.Data.HC1.iPSC.NI[,"Log.Fold.change"]<=-1.5 & T.Data.HC1.iPSC.NI[,"neg.log10.adj.pvalue"] >1.3)
View(T.Data.HC1.iPSC.NI)
#need to make it a matrix so it exports friendly
T.Data.HC1.iPSC.NI.matrix = as.matrix(T.Data.HC1.iPSC.NI)
write.csv(T.Data.HC1.iPSC.NI.matrix, file = "T.Data.HC1.iPSC.NI.csv")

View(Normalised.HC1.DA)
View(T.Data.HC1.iPSC.DA)

T.Data.HC1.iPSC.DA = cbind(Normalised.HC1.iPSC[7:11], Normalised.HC1.DA[6:9])
T.Data.HC1.iPSC.DA$P.value = apply(T.Data.HC1.iPSC.DA, 1, function(T.Data.HC1.iPSC.DA){t.test (x = T.Data.HC1.iPSC.DA[1:5], y = T.Data.HC1.iPSC.DA[6:9], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.HC1.iPSC.DA$adj.pvalues = p.adjust(T.Data.HC1.iPSC.DA$P.value, method = "BH", n=length(T.Data.HC1.iPSC.DA$P.value))
T.Data.HC1.iPSC.DA$normalised.mean.iPSC = apply(T.Data.HC1.iPSC.DA[1:5], 1, mean)
T.Data.HC1.iPSC.DA$Abundance.iPSC = (2^T.Data.HC1.iPSC.DA["normalised.mean.iPSC"])
T.Data.HC1.iPSC.DA$normalised.mean.DA = apply(T.Data.HC1.iPSC.DA[6:9], 1, mean)
T.Data.HC1.iPSC.DA$Abundance.DA = (2^T.Data.HC1.iPSC.DA["normalised.mean.DA"])
T.Data.HC1.iPSC.DA$Fold.change = (T.Data.HC1.iPSC.DA$Abundance.iPSC/T.Data.HC1.iPSC.DA$Abundance.DA)
T.Data.HC1.iPSC.DA$Log.Fold.change = log2(T.Data.HC1.iPSC.DA[,"Fold.change"])
T.Data.HC1.iPSC.DA$neg.log10.adj.pvalue = -log10(T.Data.HC1.iPSC.DA [,"adj.pvalues"])
T.Data.HC1.iPSC.DA$Threshold = as.factor(T.Data.HC1.iPSC.DA[,"Log.Fold.change"]>=1.5 & T.Data.HC1.iPSC.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.HC1.iPSC.DA$Threshold.neg = as.factor(T.Data.HC1.iPSC.DA[,"Log.Fold.change"]<=-1.5 & T.Data.HC1.iPSC.DA[,"neg.log10.adj.pvalue"] >1.3)

View(T.Data.HC1.iPSC.DA)
T.Data.HC1.iPSC.DA.matrix = as.matrix(T.Data.HC1.iPSC.DA)
write.csv(T.Data.HC1.iPSC.DA.matrix, file = "T.Data.HC1.iPSC.DA.csv")

View(T.Data.HC1.NI.DA)
T.Data.HC1.NI.DA = cbind(Normalised.HC1.NI[5:7], Normalised.HC1.DA[6:9])
T.Data.HC1.NI.DA$P.value = apply(T.Data.HC1.NI.DA, 1, function(T.Data.HC1.NI.DA){t.test (x = T.Data.HC1.NI.DA[1:3], y = T.Data.HC1.NI.DA[4:7], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.HC1.NI.DA$adj.pvalues = p.adjust(T.Data.HC1.NI.DA$P.value, method = "BH", n=length(T.Data.HC1.NI.DA$P.value))
T.Data.HC1.NI.DA$normalised.mean.NI = apply(T.Data.HC1.NI.DA[1:3], 1, mean)
T.Data.HC1.NI.DA$Abundance.NI = (2^T.Data.HC1.NI.DA["normalised.mean.NI"])
T.Data.HC1.NI.DA$normalised.mean.DA = apply(T.Data.HC1.NI.DA[4:7], 1, mean)
T.Data.HC1.NI.DA$Abundance.DA = (2^T.Data.HC1.NI.DA["normalised.mean.DA"])
T.Data.HC1.NI.DA$Fold.change = (T.Data.HC1.NI.DA$Abundance.NI/T.Data.HC1.NI.DA$Abundance.DA)
T.Data.HC1.NI.DA$Log.Fold.change = log2(T.Data.HC1.NI.DA[,"Fold.change"])
T.Data.HC1.NI.DA$neg.log10.adj.pvalue = -log10(T.Data.HC1.NI.DA [,"adj.pvalues"])
T.Data.HC1.NI.DA$Threshold = as.factor(T.Data.HC1.NI.DA[,"Log.Fold.change"]>=1.5 & T.Data.HC1.NI.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.HC1.NI.DA$Threshold.neg = as.factor(T.Data.HC1.NI.DA[,"Log.Fold.change"]<=1.5 & T.Data.HC1.NI.DA[,"neg.log10.adj.pvalue"] >1.3)

View(T.Data.HC1.NI.DA)
T.Data.HC1.NI.DA.matrix = as.matrix(T.Data.HC1.NI.DA)
write.csv(T.Data.HC1.NI.DA.matrix, file = "T.Data.HC1.NI.DA.csv")


##HC2

View(T.Data.HC2.iPSC.NI)
View(Normalised.HC2.NI)
T.Data.HC2.iPSC.NI = cbind(Normalised.HC2.iPSC[6:9], Normalised.HC2.NI[7:11])
T.Data.HC2.iPSC.NI$P.value = apply(T.Data.HC2.iPSC.NI, 1, function(T.Data.HC2.iPSC.NI){t.test (x = T.Data.HC2.iPSC.NI[1:4], y = T.Data.HC2.iPSC.NI[5:9], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.HC2.iPSC.NI$adj.pvalues = p.adjust(T.Data.HC2.iPSC.NI$P.value, method = "BH", n=length(T.Data.HC2.iPSC.NI$P.value))
T.Data.HC2.iPSC.NI$normalised.mean.iPSC = apply(T.Data.HC2.iPSC.NI[1:4], 1, mean)
T.Data.HC2.iPSC.NI$Abundance.iPSC = (2^T.Data.HC2.iPSC.NI["normalised.mean.iPSC"])
T.Data.HC2.iPSC.NI$normalised.mean.NI = apply(T.Data.HC2.iPSC.NI[5:9], 1, mean)
T.Data.HC2.iPSC.NI$Abundance.NI = (2^T.Data.HC2.iPSC.NI["normalised.mean.NI"])
T.Data.HC2.iPSC.NI$Fold.change = (T.Data.HC2.iPSC.NI$Abundance.iPSC/T.Data.HC2.iPSC.NI$Abundance.NI)
T.Data.HC2.iPSC.NI$Log.Fold.change = log2(T.Data.HC2.iPSC.NI[,"Fold.change"])
T.Data.HC2.iPSC.NI$neg.log10.adj.pvalue = -log10(T.Data.HC2.iPSC.NI [,"adj.pvalues"])
T.Data.HC2.iPSC.NI$Threshold = as.factor(T.Data.HC2.iPSC.NI[,"Log.Fold.change"]>=1.5 & T.Data.HC2.iPSC.NI[,"neg.log10.adj.pvalue"] >1.3)
T.Data.HC2.iPSC.NI$Threshold.neg = as.factor(T.Data.HC2.iPSC.NI[,"Log.Fold.change"]<=-1.5 & T.Data.HC2.iPSC.NI[,"neg.log10.adj.pvalue"] >1.3)
T.Data.HC2.iPSC.NI.matrix = as.matrix(T.Data.HC2.iPSC.NI)
write.csv(T.Data.HC2.iPSC.NI.matrix, file = "T.Data.HC2.iPSC.NI.csv")

View(T.Data.HC2.iPSC.DA)
View(Normalised.HC2.DA)
T.Data.HC2.iPSC.DA = cbind(Normalised.HC2.iPSC[6:9], Normalised.HC2.DA[6:9])
T.Data.HC2.iPSC.DA$P.value = apply(T.Data.HC2.iPSC.DA, 1, function(T.Data.HC2.iPSC.DA){t.test (x = T.Data.HC2.iPSC.DA[1:4], y = T.Data.HC2.iPSC.DA[5:8], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.HC2.iPSC.DA$adj.pvalues = p.adjust(T.Data.HC2.iPSC.DA$P.value, method = "BH", n=length(T.Data.HC2.iPSC.DA$P.value))
T.Data.HC2.iPSC.DA$normalised.mean.iPSC = apply(T.Data.HC2.iPSC.DA[1:4], 1, mean)
T.Data.HC2.iPSC.DA$Abundance.iPSC = (2^T.Data.HC2.iPSC.DA["normalised.mean.iPSC"])
T.Data.HC2.iPSC.DA$normalised.mean.DA = apply(T.Data.HC2.iPSC.DA[5:8], 1, mean)
T.Data.HC2.iPSC.DA$Abundance.DA = (2^T.Data.HC2.iPSC.DA["normalised.mean.DA"])
T.Data.HC2.iPSC.DA$Fold.change = (T.Data.HC2.iPSC.DA$Abundance.iPSC/T.Data.HC2.iPSC.DA$Abundance.DA)
T.Data.HC2.iPSC.DA$Log.Fold.change = log2(T.Data.HC2.iPSC.DA[,"Fold.change"])
T.Data.HC2.iPSC.DA$neg.log10.adj.pvalue = -log10(T.Data.HC2.iPSC.DA [,"adj.pvalues"])
T.Data.HC2.iPSC.DA$Threshold = as.factor(T.Data.HC2.iPSC.DA[,"Log.Fold.change"]>=1.5 & T.Data.HC2.iPSC.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.HC2.iPSC.DA$Threshold.neg = as.factor(T.Data.HC2.iPSC.DA[,"Log.Fold.change"]<=-1.5 & T.Data.HC2.iPSC.DA[,"neg.log10.adj.pvalue"] >1.3)

View(T.Data.HC2.iPSC.DA)
T.Data.HC2.iPSC.DA.matrix = as.matrix(T.Data.HC2.iPSC.DA)
write.csv(T.Data.HC2.iPSC.DA.matrix, file = "T.Data.HC2.iPSC.DA.csv")

View(T.Data.HC2.NI.DA)
View(Normalised.HC2.DA)
T.Data.HC2.NI.DA = cbind(Normalised.HC2.NI[7:11], Normalised.HC2.DA[6:9])
T.Data.HC2.NI.DA$P.value = apply(T.Data.HC2.NI.DA, 1, function(T.Data.HC2.NI.DA){t.test (x = T.Data.HC2.NI.DA[1:5], y = T.Data.HC2.NI.DA[6:9], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.HC2.NI.DA$adj.pvalues = p.adjust(T.Data.HC2.NI.DA$P.value, method = "BH", n=length(T.Data.HC2.NI.DA$P.value))
T.Data.HC2.NI.DA$normalised.mean.NI = apply(T.Data.HC2.NI.DA[1:5], 1, mean)
T.Data.HC2.NI.DA$Abundance.NI = (2^T.Data.HC2.NI.DA["normalised.mean.NI"])
T.Data.HC2.NI.DA$normalised.mean.DA = apply(T.Data.HC2.NI.DA[6:9], 1, mean)
T.Data.HC2.NI.DA$Abundance.DA = (2^T.Data.HC2.NI.DA["normalised.mean.DA"])
T.Data.HC2.NI.DA$Fold.change = (T.Data.HC2.NI.DA$Abundance.NI/T.Data.HC2.NI.DA$Abundance.DA)
T.Data.HC2.NI.DA$Log.Fold.change = log2(T.Data.HC2.NI.DA[,"Fold.change"])
T.Data.HC2.NI.DA$neg.log10.adj.pvalue = -log10(T.Data.HC2.NI.DA [,"adj.pvalues"])
T.Data.HC2.NI.DA$Threshold = as.factor(T.Data.HC2.NI.DA[,"Log.Fold.change"]>=1.5 & T.Data.HC2.NI.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.HC2.NI.DA$Threshold.neg = as.factor(T.Data.HC2.NI.DA[,"Log.Fold.change"]<=-1.5 & T.Data.HC2.NI.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.HC2.NI.DA.matrix = as.matrix(T.Data.HC2.NI.DA)
write.csv(T.Data.HC2.NI.DA.matrix, file = "T.Data.HC2.NI.DA.csv")

##SAD

T.Data.SAD.iPSC.NI = cbind(Normalised.SAD.iPSC[6:9], Normalised.SAD.NI[7:11])
T.Data.SAD.iPSC.NI$P.value = apply(T.Data.SAD.iPSC.NI, 1, function(T.Data.SAD.iPSC.NI){t.test (x = T.Data.SAD.iPSC.NI[1:4], y = T.Data.SAD.iPSC.NI[5:9], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.SAD.iPSC.NI$adj.pvalues = p.adjust(T.Data.SAD.iPSC.NI$P.value, method = "BH", n=length(T.Data.SAD.iPSC.NI$P.value))
T.Data.SAD.iPSC.NI$normalised.mean.iPSC = apply(T.Data.SAD.iPSC.NI[1:4], 1, mean)
T.Data.SAD.iPSC.NI$Abundance.iPSC = (2^T.Data.SAD.iPSC.NI["normalised.mean.iPSC"])
T.Data.SAD.iPSC.NI$normalised.mean.NI = apply(T.Data.SAD.iPSC.NI[5:9], 1, mean)
T.Data.SAD.iPSC.NI$Abundance.NI = (2^T.Data.SAD.iPSC.NI["normalised.mean.NI"])
T.Data.SAD.iPSC.NI$Fold.change = (T.Data.SAD.iPSC.NI$Abundance.iPSC/T.Data.SAD.iPSC.NI$Abundance.NI)
T.Data.SAD.iPSC.NI$Log.Fold.change = log2(T.Data.SAD.iPSC.NI[,"Fold.change"])
T.Data.SAD.iPSC.NI$neg.log10.adj.pvalue = -log10(T.Data.SAD.iPSC.NI [,"adj.pvalues"])
T.Data.SAD.iPSC.NI$Threshold = as.factor(T.Data.SAD.iPSC.NI[,"Log.Fold.change"]>=1.5 & T.Data.SAD.iPSC.NI[,"neg.log10.adj.pvalue"] >1.3)
T.Data.SAD.iPSC.NI$Threshold.neg = as.factor(T.Data.SAD.iPSC.NI[,"Log.Fold.change"]<=-1.5 & T.Data.SAD.iPSC.NI[,"neg.log10.adj.pvalue"] >1.3)
T.Data.SAD.iPSC.NI.matrix = as.matrix(T.Data.SAD.iPSC.NI)
write.csv(T.Data.SAD.iPSC.NI.matrix, file = "T.Data.SAD.iPSC.NI.csv")

View(T.Data.SAD.iPSC.DA)
View(Normalised.SAD.DA)
T.Data.SAD.iPSC.DA = cbind(Normalised.SAD.iPSC[6:9], Normalised.SAD.DA[6:9])
T.Data.SAD.iPSC.DA$P.value = apply(T.Data.SAD.iPSC.DA, 1, function(T.Data.SAD.iPSC.DA){t.test (x = T.Data.SAD.iPSC.DA[1:4], y = T.Data.SAD.iPSC.DA[5:8], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.SAD.iPSC.DA$adj.pvalues = p.adjust(T.Data.SAD.iPSC.DA$P.value, method = "BH", n=length(T.Data.SAD.iPSC.DA$P.value))
T.Data.SAD.iPSC.DA$normalised.mean.iPSC = apply(T.Data.SAD.iPSC.DA[1:4], 1, mean)
T.Data.SAD.iPSC.DA$Abundance.iPSC = (2^T.Data.SAD.iPSC.DA["normalised.mean.iPSC"])
T.Data.SAD.iPSC.DA$normalised.mean.DA = apply(T.Data.SAD.iPSC.DA[5:8], 1, mean)
T.Data.SAD.iPSC.DA$Abundance.DA = (2^T.Data.SAD.iPSC.DA["normalised.mean.DA"])
T.Data.SAD.iPSC.DA$Fold.change = (T.Data.SAD.iPSC.DA$Abundance.iPSC/T.Data.SAD.iPSC.DA$Abundance.DA)
T.Data.SAD.iPSC.DA$Log.Fold.change = log2(T.Data.SAD.iPSC.DA[,"Fold.change"])
T.Data.SAD.iPSC.DA$neg.log10.adj.pvalue = -log10(T.Data.SAD.iPSC.DA [,"adj.pvalues"])
T.Data.SAD.iPSC.DA$Threshold = as.factor(T.Data.SAD.iPSC.DA[,"Log.Fold.change"]>=1.5 & T.Data.SAD.iPSC.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.SAD.iPSC.DA$Threshold.neg = as.factor(T.Data.SAD.iPSC.DA[,"Log.Fold.change"]<=-1.5 & T.Data.SAD.iPSC.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.SAD.iPSC.DA.matrix = as.matrix(T.Data.SAD.iPSC.DA)
write.csv(T.Data.SAD.iPSC.DA.matrix, file = "T.Data.SAD.iPSC.DA.csv")

View(T.Data.SAD.NI.DA)
View(Normalised.SAD.DA)
T.Data.SAD.NI.DA = cbind(Normalised.SAD.NI[7:11], Normalised.SAD.DA[6:9])
T.Data.SAD.NI.DA$P.value = apply(T.Data.SAD.NI.DA, 1, function(T.Data.SAD.NI.DA){t.test (x = T.Data.SAD.NI.DA[1:5], y = T.Data.SAD.NI.DA[6:9], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.SAD.NI.DA$adj.pvalues = p.adjust(T.Data.SAD.NI.DA$P.value, method = "BH", n=length(T.Data.SAD.NI.DA$P.value))
T.Data.SAD.NI.DA$normalised.mean.NI = apply(T.Data.SAD.NI.DA[1:5], 1, mean)
T.Data.SAD.NI.DA$Abundance.NI = (2^T.Data.SAD.NI.DA["normalised.mean.NI"])
T.Data.SAD.NI.DA$normalised.mean.DA = apply(T.Data.SAD.NI.DA[6:9], 1, mean)
T.Data.SAD.NI.DA$Abundance.DA = (2^T.Data.SAD.NI.DA["normalised.mean.DA"])
T.Data.SAD.NI.DA$Fold.change = (T.Data.SAD.NI.DA$Abundance.NI/T.Data.SAD.NI.DA$Abundance.DA)
T.Data.SAD.NI.DA$Log.Fold.change = log2(T.Data.SAD.NI.DA[,"Fold.change"])
T.Data.SAD.NI.DA$neg.log10.adj.pvalue = -log10(T.Data.SAD.NI.DA [,"adj.pvalues"])
T.Data.SAD.NI.DA$Threshold = as.factor(T.Data.SAD.NI.DA[,"Log.Fold.change"]>=1.5 & T.Data.SAD.NI.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.SAD.NI.DA$Threshold.neg = as.factor(T.Data.SAD.NI.DA[,"Log.Fold.change"]<=-1.5 & T.Data.SAD.NI.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.SAD.NI.DA.matrix = as.matrix(T.Data.SAD.NI.DA)
write.csv(T.Data.SAD.NI.DA.matrix, file = "T.Data.SAD.NI.DA.csv")

##-------------------------------------------------------------------------------------------------------------------------

#iPSCs v Ni

View(T.Data.iPSC.NI)
View(Normalised.All.iPSC)
View(Normalised.All.NI)

View(Normalised.All.DA)

T.Data.iPSC.NI = cbind(Normalised.All.iPSC[15:27], Normalised.All.NI[15:27])
T.Data.iPSC.NI$P.value = apply(T.Data.iPSC.NI, 1, function(T.Data.iPSC.NI){t.test (x = T.Data.iPSC.NI[1:13], y = T.Data.iPSC.NI[14:26], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.iPSC.NI$adj.pvalues = p.adjust(T.Data.iPSC.NI$P.value, method = "BH", n=length(T.Data.iPSC.NI$P.value))
T.Data.iPSC.NI$normalised.mean.iPSC = apply(T.Data.iPSC.NI[1:13], 1, mean)
T.Data.iPSC.NI$Abundance.iPSC = (2^T.Data.iPSC.NI["normalised.mean.iPSC"])
T.Data.iPSC.NI$normalised.mean.NI = apply(T.Data.iPSC.NI[14:26], 1, mean)
T.Data.iPSC.NI$Abundance.NI = (2^T.Data.iPSC.NI["normalised.mean.NI"])
T.Data.iPSC.NI$Fold.change = (T.Data.iPSC.NI$Abundance.iPSC/T.Data.iPSC.NI$Abundance.NI)
T.Data.iPSC.NI$Log.Fold.change = log2(T.Data.iPSC.NI[,"Fold.change"])
T.Data.iPSC.NI$neg.log10.adj.pvalue = -log10(T.Data.iPSC.NI [,"adj.pvalues"])
T.Data.iPSC.NI$Threshold = as.factor(T.Data.iPSC.NI[,"Log.Fold.change"]>=1.5 & T.Data.iPSC.NI[,"neg.log10.adj.pvalue"] >1.3)
T.Data.iPSC.NI$Threshold.neg = as.factor(T.Data.iPSC.NI[,"Log.Fold.change"]<=-1.5 & T.Data.iPSC.NI[,"neg.log10.adj.pvalue"] >1.3)
View(T.Data.iPSC.NI)
T.Data.iPSC.NI.matrix = as.matrix(T.Data.iPSC.NI)
write.csv(T.Data.iPSC.NI.matrix, file = "T.Data.iPSC.NI.csv")


#iPSCs v DA

View(T.Data.iPSC.DA)


T.Data.iPSC.DA = cbind(Normalised.All.iPSC[15:27], Normalised.All.DA[14:25])
T.Data.iPSC.DA$P.value = apply(T.Data.iPSC.DA, 1, function(T.Data.iPSC.DA){t.test (x = T.Data.iPSC.DA[1:13], y = T.Data.iPSC.DA[14:25], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.iPSC.DA$adj.pvalues = p.adjust(T.Data.iPSC.DA$P.value, method = "BH", n=length(T.Data.iPSC.DA$P.value))
T.Data.iPSC.DA$normalised.mean.iPSC = apply(T.Data.iPSC.DA[1:13], 1, mean)
T.Data.iPSC.DA$Abundance.iPSC = (2^T.Data.iPSC.DA["normalised.mean.DA"])
T.Data.iPSC.DA$normalised.mean.DA = apply(T.Data.iPSC.DA[14:25], 1, mean)
T.Data.iPSC.DA$Abundance.DA = (2^T.Data.iPSC.DA["normalised.mean.DA"])
T.Data.iPSC.DA$Fold.change = (T.Data.iPSC.DA$Abundance.iPSC/T.Data.iPSC.DA$Abundance.DA)
T.Data.iPSC.DA$Log.Fold.change = log2(T.Data.iPSC.DA[,"Fold.change"])
T.Data.iPSC.DA$neg.log10.adj.pvalue = -log10(T.Data.iPSC.DA [,"adj.pvalues"])
T.Data.iPSC.DA$Threshold = as.factor(T.Data.iPSC.DA[,"Log.Fold.change"]>=1.5 & T.Data.iPSC.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.iPSC.DA$Threshold.neg = as.factor(T.Data.iPSC.DA[,"Log.Fold.change"]<=-1.5 & T.Data.iPSC.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.iPSC.DA.matrix = as.matrix(T.Data.iPSC.DA)
write.csv(T.Data.iPSC.DA.matrix, file = "T.Data.iPSC.DA.csv")

#NI v DA

T.Data.NI.DA = cbind(Normalised.All.NI[7:11], Normalised.All.DA[6:9])
T.Data.NI.DA$P.value = apply(T.Data.NI.DA, 1, function(T.Data.NI.DA){t.test (x = T.Data.NI.DA[1:5], y = T.Data.NI.DA[6:9], alternative = "two.sided", var.equal= FALSE)$p.value})
T.Data.NI.DA$adj.pvalues = p.adjust(T.Data.NI.DA$P.value, method = "BH", n=length(T.Data.NI.DA$P.value))
T.Data.NI.DA$normalised.mean.NI = apply(T.Data.NI.DA[1:5], 1, mean)
T.Data.NI.DA$Abundance.NI = (2^T.Data.NI.DA["normalised.mean.DA"])
T.Data.NI.DA$normalised.mean.DA = apply(T.Data.NI.DA[6:9], 1, mean)
T.Data.NI.DA$Abundance.DA = (2^T.Data.NI.DA["normalised.mean.DA"])
T.Data.NI.DA$Fold.change = (T.Data.NI.DA$Abundance.NI/T.Data.NI.DA$Abundance.DA)
T.Data.NI.DA$Log.Fold.change = log2(T.Data.NI.DA[,"Fold.change"])
T.Data.NI.DA$neg.log10.adj.pvalue = -log10(T.Data.NI.DA [,"adj.pvalues"])
T.Data.NI.DA$Threshold = as.factor(T.Data.NI.DA[,"Log.Fold.change"]>=1.5 & T.Data.NI.DA[,"neg.log10.adj.pvalue"] >1.3)
T.Data.NI.DA$Threshold = as.factor(T.Data.NI.DA[,"Log.Fold.change"]<=-1.5 & T.Data.NI.DA[,"neg.log10.adj.pvalue"] >1.3)

View(T.Data.NI.DA)
T.Data.NI.DA.matrix = as.matrix(T.Data.NI.DA)
write.csv(T.Data.NI.DA.matrix, file = "T.Data.NI.DA.csv")


#Drafts below
##-------------------------------------------------------------------------------------------------------------------------

#TRYING graphic with threshold
T.Data.HC1.iPSC.NI.df = as.data.frame(T.Data.HC1.iPSC.NI)
VP.HC1.iPSC.NI = ggplot(data = T.Data.HC1.iPSC.NI.df, aes(x =Fold.change, y=Log.Fold.change, colour=Threshold))
VP.HC1.iPSC.NI + geom_point(alpha = 0.4, size = 1.75) + xlab ("Log2 Fold Change") + ylab("-Log10 P Values") + scale_colour_manual(values = c("black", "red")) 



p = ggplot(data=T.Data.HC1.iPSC.NI, aes(x=log2(Fold.change), y=-log10(p.value))) + geom_point() + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")



colnames(Lipids_Normalised) = c("HC1 iPSC 1", "HC1 iPSC 2", "HC1 iPSC 3", "HC1 iPSC 4", "HC1 iPSC 5", "HC1 NI 1", "HC1 NI 2", "HC1 NI 3", "HC1 DA 1", "HC1 DA 2", "HC1 DA 3", "HC1 DA 4", "HC2 iPSC 1", "HC2 iPSC 2", "HC2 iPSC 3", "HC2 iPSC 4", "HC2 NI 1", "HC2 NI 2", "HC2 NI 3", "HC2 NI 4", "HC2 NI 5", "HC2 DA 1", "HC2 DA 2", "HC2 DA 3", "HC2 DA 4", "SAD iPSC 1", "SAD iPSC 2", "SAD iPSC 3", "SAD iPSC 4", "SAD NI 1", "SAD NI 2", "SAD NI 3", "SAD NI 4", "SAD NI 5" ,"SAD DA 1", "SAD DA 2", "SAD DA 3", "SAD DA 4")

#Titles with no spaces
colnames(Lipids_Normalised) = c("HC1.iPSC.1", "HC1.iPSC.2", "HC1.iPSC.3", "HC1.iPSC.4", "HC1.iPSC.5", "HC1.NI.1", "HC1.NI.2", "HC1.NI.3", "HC1.DA.1", "HC1.DA.2", "HC1.DA.3", "HC1.DA.4", "HC2.iPSC.1", "HC2.iPSC.2", "HC2.iPSC.3", "HC2.iPSC.4", "HC2.NI.1", "HC2.NI.2", "HC2.NI.3", "HC2.NI.4", "HC2.NI.5", "HC2.DA.1", "HC2.DA.2", "HC2.DA.3", "HC2.DA.4", "SAD.iPSC.1", "SAD.iPSC.2", "SAD.iPSC.3", "SAD.iPSC.4", "SAD.NI.1", "SAD.NI.2", "SAD.NI.3", "SAD.NI.4", "SAD.NI.5" ,"SAD.DA.1", "SAD.DA.2", "SAD.DA.3", "SAD.DA.4")

View(Lipids_Normalised)
write.csv(Lipids_Normalised, file = "Lipids_Normalised.csv")


iPSC.Norm = subset(Lipids_Normalised, select = c(1:5, 13:16, 26:29))
NI.Norm = subset(Lipids_Normalised, select = c(6:8, 17:21, 30:34))
DA.Norm = subset(Lipids_Normalised, select = c(9:12, 22:25, 35:38))

All_Cell.line = subset(Lipids_Normalised, select = c(1:5, 13:16, 26:29, 6:8, 17:21, 30:34, 9:12, 22:25, 35:38))
All_Cell.line_dataframe = as.data.frame(All_Cell.line_matrix)

#multiscatter with the changed order:
All_Cell.line_matrix  = as.matrix(All_Cell.line)
colnames(All_Cell.line_matrix) = c("HC1 iPSC 1", "HC1 iPSC 2", "HC1 iPSC 3", "HC1 iPSC 4", "HC1 iPSC 5","HC2 iPSC 1", "HC2 iPSC 2", "HC2 iPSC 3", "HC2 iPSC 4", "SAD iPSC 1", "SAD iPSC 2", "SAD iPSC 3", "SAD iPSC 4", "HC1 NI 1", "HC1 NI 2", "HC1 NI 3", "HC2 NI 1", "HC2 NI 2", "HC2 NI 3", "HC2 NI 4", "HC2 NI 5", "SAD NI 1", "SAD NI 2", "SAD NI 3", "SAD NI 4", "SAD NI 5", "HC1 DA 1", "HC1 DA 2", "HC1 DA 3", "HC1 DA 4", "HC2 DA 1", "HC2 DA 2", "HC2 DA 3", "HC2 DA 4",  "SAD DA 1", "SAD DA 2", "SAD DA 3", "SAD DA 4")
colnames(All_Cell.line_matrix) = c("HC1.iPSC.1", "HC1 iPSC.2", "HC1.iPSC.3", "HC1.iPSC.4", "HC1.iPSC.5","HC2.iPSC.1", "HC2.iPSC.2", "HC2.iPSC.3", "HC2.iPSC.4", "SAD.iPSC.1", "SAD.iPSC.2", "SAD.iPSC.3", "SAD.iPSC.4", "HC1.NI.1", "HC1.NI.2", "HC1.NI.3", "HC2.NI.1", "HC2.NI.2", "HC2.NI.3", "HC2.NI.4", "HC2.NI.5", "SAD.NI.1", "SAD.NI.2", "SAD.NI.3", "SAD.NI.4", "SAD.NI.5", "HC1.DA.1", "HC1.DA.2", "HC1.DA.3", "HC1.DA.4", "HC2.DA.1", "HC2.DA.2", "HC2.DA.3", "HC2.DA.4",  "SAD.DA.1", "SAD.DA.2", "SAD.DA.3", "SAD.DA.4")

View(All_Cell.line_matrix)

plot(All_Cell.line_matrix, cex = 0.01)

#pearsons - 
All_Pearson=cor(All_Cell.line_matrix, method="pearson", use = "complete.obs")
write.csv(All_Pearson, "All_Pearson.csv")

#pearsons correlation of samples - 
All_Pearson=cor(Lipids_Multiscatterdata, method="pearson", use = "complete.obs")
write.csv(All_Pearson, "All_Pearson.csv")

#heatmap dendrogram
Lipids_Normalised_Matrix = data.matrix(Lipids_Normalised)
heatmap.2(Lipids_Normalised_Matrix, main = "", notecol = "black", density.info = "none", trace = "none", col = "bluered", dendogram = "both", srtCol=45, cexCol=0.75, cexRow=0.01, keysize=1)

##---x---x---x---x---x---x---x---x---x---x

######PCA#####ALL######PCA#####ALL######PCA#####ALL######PCA#####ALL######PCA#####ALL######PCA#####ALL######PCA#####ALL
# ---- 1. Load ----

folder <- "/Users/baxx/Library/CloudStorage/OneDrive-Personal/7.Career/7.4 Victor Chang Era/0. Science Management/5. Writing/20. Bax 2026 iPSC Lipids/1. Raw Data/5.Lipidomics/R Analysis Lipids/6. Normalised Values/Lipids_Normalised.csv"
folder            # print it back — confirm it shows your path, not a function
file.exists(folder)

# now row 1 of `dat` is the sample-name row that was misread; promote it to column names
colnames(dat) <- as.character(unlist(dat[1, ]))
dat <- dat[-1, ]                                   # drop that name row
dat[] <- lapply(dat, function(x) as.numeric(as.character(x)))  # force numeric

# ---- 2. Clean to a numeric matrix (lipids = rows, samples = cols) ----
# If column 1 is an index (1,2,3...) and column 2 is lipid names, adjust indices below.
rownames(dat) <- dat[[1]]                 # assumes col 1 = lipid species names
dat <- dat[, sapply(dat, is.numeric)]     # keep only numeric sample columns
mat <- as.matrix(dat)
dim(mat)                                   # should be [n_lipids x 38]
range(mat, na.rm = TRUE)                   # sanity: normalised data straddles 0 (neg AND pos)
mat <- mat[apply(mat, 1, sd) > 0, ]        # drop any zero-variance lipids (scale needs this)


# ---- Fix the matrix: drop index col + reorder to match your palette/labels ----
sample_order <- c("HC1 iPSC 1","HC1 iPSC 2","HC1 iPSC 3","HC1 iPSC 4","HC1 iPSC 5",
                  "HC2 iPSC 1","HC2 iPSC 2","HC2 iPSC 3","HC2 iPSC 4",
                  "SAD iPSC 1","SAD iPSC 2","SAD iPSC 3","SAD iPSC 4",
                  "HC1 NI 1","HC1 NI 2","HC1 NI 3",
                  "HC2 NI 1","HC2 NI 2","HC2 NI 3","HC2 NI 4","HC2 NI 5",
                  "SAD NI 1","SAD NI 2","SAD NI 3","SAD NI 4","SAD NI 5",
                  "HC1 DA 1","HC1 DA 2","HC1 DA 3","HC1 DA 4",
                  "HC2 DA 1","HC2 DA 2","HC2 DA 3","HC2 DA 4",
                  "SAD DA 1","SAD DA 2","SAD DA 3","SAD DA 4")

mat <- as.matrix(dat[, -1])                 # drop the leftover index column -> 59 x 38
setdiff(sample_order, colnames(mat))        # MUST return character(0) — confirms names match
mat <- mat[, sample_order]                  # reorder columns to match palette + labels
mat <- mat[apply(mat, 1, sd) > 0, ]         # drop zero-variance lipids

dim(mat)        # expect 59 x 38
range(mat)      # expect ~ -2.4 .. +2.x  — the 59 should be GONE



# ============================================================
#  PCA — normalised lipidomics, all samples
#  % variance axes | NR/N relabel | filled CI (no border) | L-shaped axes
# ============================================================
library(ggplot2)

# ---- 1. PCA ----
primp_ALL <- prcomp(mat, center = TRUE, scale = TRUE)
pv <- round(100 * primp_ALL$sdev^2 / sum(primp_ALL$sdev^2), 1)

# ---- 2. Plotting frame ----
loadings <- data.frame(primp_ALL$rotation)

labels_new <- rownames(loadings)
labels_new <- gsub(" NI ", " NR ", labels_new)
labels_new <- gsub(" DA ", " N ",  labels_new)

loadings$State <- sub("^[^ ]+ ([^ ]+).*$", "\\1", labels_new)    # iPSC / NR / N
loadings$State <- factor(loadings$State, levels = c("iPSC","NR","N"))

# display labels: rename cell lines LAST (display only)
labels_disp <- labels_new
labels_disp <- sub("^HC1 ", "Line 1 ", labels_disp)
labels_disp <- sub("^HC2 ", "Line 2 ", labels_disp)
labels_disp <- sub("^SAD ", "Line 3 ", labels_disp)

# ---- 3. Colours ----
prpallete <- c("#b69dff","#b69dff","#b69dff","#b69dff","#b69dff",
               "#8f8dff","#8f8dff","#8f8dff","#8f8dff",
               "#5762ff","#5762ff","#5762ff","#5762ff",
               "#9dd5ff","#9dd5ff","#9dd5ff","#9dd5ff",
               "#8dfffd","#8dfffd","#8dfffd","#8dfffd",
               "#0a6699","#0a6699","#0a6699","#0a6699","#0a6699",
               "#b5ffb5","#b5ffb5","#b5ffb5","#b5ffb5",
               "#5fff57","#5fff57","#5fff57","#5fff57",
               "#1f9e1f","#1f9e1f","#1f9e1f","#1f9e1f")
state_fills <- c("iPSC" = "#8f8dff", "NR" = "#57a8ff", "N" = "#5fff57")

# ---- 4. Plot ----
p <- ggplot(loadings, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(fill = State), geom = "polygon",
               alpha = 0.15, colour = NA, type = "t", level = 0.95) +
  geom_point(size = 5, colour = prpallete) +
  geom_text(aes(label = labels_disp), hjust = 0, vjust = 2, size = 3) +
  scale_fill_manual(values = state_fills) +
  labs(title = "",
       x = paste0("PC1 (", pv[1], "%)"),
       y = paste0("PC2 (", pv[2], "%)"),
       fill = "State") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border     = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 0.8),
        plot.title       = element_text(face = "bold", size = rel(2)))

p   # view it

ggsave("PCA_Lipids_Normalised.pdf", plot = p, width = 7, height = 6, units = "in")
ggsave("PCA_Lipids_Normalised.png", plot = p, width = 7, height = 6, units = "in", dpi = 600)





#-------------------------------------------------------



# ============================================================
#  PCA — Lines 1 & 2 only (Line 3 / SAD removed), recomputed
# ============================================================
library(ggplot2)

# ---- 1. Drop Line 3 (SAD) columns by NAME, then clean ----
mat2 <- mat[, !grepl("^SAD", colnames(mat))]      # keep everything except SAD
mat2 <- mat2[apply(mat2, 1, sd) > 0, ]            # re-drop zero-variance lipids on the subset
dim(mat2)                                          # expect ~25 samples (HC1 + HC2)

# ---- 2. PCA recomputed on the subset ----
primp2 <- prcomp(mat2, center = TRUE, scale = TRUE)
pv2 <- round(100 * primp2$sdev^2 / sum(primp2$sdev^2), 1)

# ---- 3. Plotting frame ----
loadings2 <- data.frame(primp2$rotation)

labels_new <- rownames(loadings2)
labels_new <- gsub(" NI ", " NR ", labels_new)
labels_new <- gsub(" DA ", " N ",  labels_new)

loadings2$Line  <- sub(" .*$", "", labels_new)                  # HC1 / HC2
loadings2$State <- sub("^[^ ]+ ([^ ]+).*$", "\\1", labels_new)  # iPSC / NR / N
loadings2$State <- factor(loadings2$State, levels = c("iPSC","NR","N"))

# display labels: rename lines LAST
labels_disp <- labels_new
labels_disp <- sub("^HC1 ", "Line 1 ", labels_disp)
labels_disp <- sub("^HC2 ", "Line 2 ", labels_disp)

# ---- 4. Colours by Line x State (robust to subsetting) ----
pal <- c(
  "HC1 iPSC" = "#b69dff", "HC2 iPSC" = "#8f8dff",
  "HC1 NR"   = "#9dd5ff", "HC2 NR"   = "#8dfffd",
  "HC1 N"    = "#b5ffb5", "HC2 N"    = "#5fff57"
)
point_cols <- pal[paste(loadings2$Line, loadings2$State)]
stopifnot(!any(is.na(point_cols)))                # stops if any sample didn't match

state_fills <- c("iPSC" = "#8f8dff", "NR" = "#57a8ff", "N" = "#5fff57")

# ---- 5. Plot ----
p2 <- ggplot(loadings2, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(fill = State), geom = "polygon",
               alpha = 0.15, colour = NA, type = "t", level = 0.95) +
  geom_point(size = 5, colour = point_cols) +
  geom_text(aes(label = labels_disp), hjust = 0, vjust = 2, size = 3) +
  scale_fill_manual(values = state_fills) +
  labs(title = "",
       x = paste0("PC1 (", pv2[1], "%)"),
       y = paste0("PC2 (", pv2[2], "%)"),
       fill = "State") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border     = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 0.8),
        plot.title       = element_text(face = "bold", size = rel(2)))

p2

ggsave("PCA_Lipids_Lines1-2_noSAD.pdf", plot = p2, width = 7, height = 6, units = "in")
ggsave("PCA_Lipids_Lines1-2_noSAD.png", plot = p2, width = 7, height = 6, units = "in", dpi = 600)

# Metaboanalyst analysis of pre-imputed data (PCA grid and hierachial clustering)

# PID of current job: 4104180
mSet<-InitDataObjects("conc", "stat", FALSE, 150)
mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-PerformSanityClosure (mSet);
mSet<-CheckContainsBlank(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 150, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 150, width=NA)
mSet<-RemoveMissingByPercent(mSet, percent=1.0, F)
mSet<-FilterVariable(mSet, "F", 20, "iqr", 5, "mean", 0, F,10.0)
mSet<-ANOVA.Anal(mSet, F, 0.05, FALSE)
mSet<-PlotANOVA(mSet, "aov_0_", "png", 150, width=NA)
mSet<-PCA.Anal(mSet)
nrow(mSetObj$dataSet$meta.info)
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_0_", "png", 96, width=NA, 5, "NA", "NA")
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 150, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 150, width=NA, 1,2,0.95,0,0, "na")
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 150, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 150, width=NA, 1,2,10)
mSet<-PlotPCA3DScore(mSet, "pca_score3d_0_", "json", 1,2,3)
mSet<-PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1,2,3)
nrow(mSetObj$dataSet$meta.info)
mSet<-GetGroupNames(mSet, "null")
nrow(mSetObj$dataSet$meta.info)
nrow(mSetObj$dataSet$meta.info)
colVec<-c("##8989c3","##6ca2d7","##87c541")
shapeVec<-c(0,0,0)
mSet<-UpdateGraphSettings(mSet, colVec, shapeVec)
nrow(mSetObj$dataSet$meta.info)
nrow(mSetObj$dataSet$meta.info)
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_1_", "png", 96, width=NA, 4, "NA", "NA")
nrow(mSetObj$dataSet$meta.info)
nrow(mSetObj$dataSet$meta.info)
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_1_", "pdf", 72, width=NA, 4, "NA", "NA")
nrow(mSetObj$dataSet$meta.info)
nrow(mSetObj$dataSet$meta.info)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "pdf", 72, width=NA, 1,2,0.95,0,0, "na")
nrow(mSetObj$dataSet$meta.info)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_1_", "png", 150, width=NA, 2,3,0.95,0,0, "na")
nrow(mSetObj$dataSet$meta.info)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_2_", "png", 150, width=NA, 2,4,0.95,0,0, "na")
nrow(mSetObj$dataSet$meta.info)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_3_", "png", 150, width=NA, 2,4,0.95,0,0, "na")
nrow(mSetObj$dataSet$meta.info)
mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 96, width=NA, 5)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 150, width=NA, 1,2,0.95,0,0, "na")
mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 150, width=NA, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "png", 150, width=NA, 1, 2);
mSet<-PlotPLS3DLoading(mSet, "pls_loading3d_0_", "json", 1,2,3)
mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "png", 150, width=NA, "vip", "Comp. 1", 15,FALSE)
mSet<-PlotPLSBiplot(mSet, "pls_biplot_0_", "png", 150, width=NA, 1,2,10)
mSet<-PlotSubHeatMap(mSet, "heatmap_1_", "png", 150, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8,8, 10.0,0.02,10, 10, "tanova", 20, T, T, T, F, T, T, T,T,"overview")
mSet<-PlotSubHeatMap(mSet, "heatmap_1_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8,8, 10.0,0.02,10, 10, "tanova", 20, T, T, T, F, T, T, T,T,"overview" , T)
mSet<-PlotStaticHeatMap(mSet, "heatmap_2_", "png", 150, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8,8, "overview", T, T, NULL, T, F, T, T, T)
mSet<-PlotStaticHeatMap(mSet, "heatmap_2_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8,8, "overview", T, T, NULL, T, F, T, T, T)
mSet<-PlotStaticCorrHeatMap(mSet, "corr_1_", "png", 150, width=NA, "row", "pearson", "bwm", "overview", F, F, 0.0)
mSet<-PlotStaticCorrHeatMap(mSet, "corr_1_", "pdf", 72, width=NA, "row", "pearson", "bwm", "overview", F, F, 0.0)
