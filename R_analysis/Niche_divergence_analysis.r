

#############Niche divergence analyses across taxonomic ranks ############################################################################################################################

###lode R libraries for this session 
library (scales)
library(base)
library(car)
lirary(ggplot2)



###### Calculate the slop from the relationship of Abundance and Occupancy at the OTU level ######

#load OTU table
OTU_RelAbun <-read.csv("InputFiles/RelOTU_table.csv", row.names=1, header=T) #######RelOTU_Table, without taxonomic information

## Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
OTU_RelAbun.M <-OTU_RelAbun[1:18,]
OTU_RelAbun.S <-OTU_RelAbun[19:36,]

## Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
OTU_RelAbun.M.PA <-OTU_RelAbun.M[2:10,]
OTU_RelAbun.M.FL <-OTU_RelAbun.M[c(1, 11:18), ]
OTU_RelAbun.S.PA <-OTU_RelAbun.S[2:10,]
OTU_RelAbun.S.FL <-OTU_RelAbun.S[c(1, 11:18), ]

## Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
OTU_RelAbun.M.PA <-t(OTU_RelAbun.M.PA)
OTU_RelAbun.M.FL <-t(OTU_RelAbun.M.FL)
OTU_RelAbun.S.PA <-t(OTU_RelAbun.S.PA)
OTU_RelAbun.S.FL <-t(OTU_RelAbun.S.FL)

## Calculate the regional mean abundance for every OTU in each subset.
OTU_MeanRelAbun.M.FL <-rowMeans(OTU_RelAbun.M.FL)
OTU_MeanRelAbun.S.FL <-rowMeans(OTU_RelAbun.S.FL)
OTU_MeanRelAbun.M.PA <-rowMeans(OTU_RelAbun.M.PA)
OTU_MeanRelAbun.S.PA <-rowMeans(OTU_RelAbun.S.PA)

## Count Rows that fit the criteria of Relative abundance =0 for each OTU. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the OTUs. .
OTU_OCCU.M.PA <-rowSums(OTU_RelAbun.M.PA==0)
OTU_OCCU.M.FL <-rowSums(OTU_RelAbun.M.FL==0)
OTU_OCCU.S.PA <-rowSums(OTU_RelAbun.S.PA==0)
OTU_OCCU.S.FL <-rowSums(OTU_RelAbun.S.FL==0)

## Calcualte the number of sites that are occupied by the OTUs, by substracting the number of sites unoccupied from the total number of sites 9
OTU_OCCU.M.PA2 <-9-OTU_OCCU.M.PA
OTU_OCCU.M.FL2 <-9-OTU_OCCU.M.FL
OTU_OCCU.S.PA2 <-9-OTU_OCCU.S.PA
OTU_OCCU.S.FL2 <-9-OTU_OCCU.S.FL

## Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance OTUs did not occur by chance for once.
OTU_OA.M.FL <-cbind(OTU_OCCU.M.FL2, OTU_MeanRelAbun.M.FL)
OTU_OA.S.FL <-cbind(OTU_OCCU.S.FL2, OTU_MeanRelAbun.S.FL)
OTU_OA.M.PA <-cbind(OTU_OCCU.M.PA2, OTU_MeanRelAbun.M.PA)
OTU_OA.S.PA <-cbind(OTU_OCCU.S.PA2, OTU_MeanRelAbun.S.PA)

colnames(OTU_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(OTU_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(OTU_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(OTU_OA.S.PA) <-c("Occupancy","Relative abundance")

## Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
OTU_OA.S.FL[OTU_OA.S.FL <= 0.00002]<-NA
OTU_OA.M.FL[OTU_OA.M.FL <= 0.00002]<-NA
OTU_OA.S.PA[OTU_OA.S.PA <= 0.00002]<-NA
OTU_OA.M.PA[OTU_OA.M.PA <= 0.00002]<-NA


## Remove the OTUs that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
OTU_OA.S.FL_exNA <-na.omit(OTU_OA.S.FL, cols=seq_along(OTU_OA.S.FL))
OTU_OA.M.FL_exNA <-na.omit(OTU_OA.M.FL, cols=seq_along(OTU_OA.M.FL))
OTU_OA.S.PA_exNA <-na.omit(OTU_OA.S.PA, cols=seq_along(OTU_OA.S.PA))
OTU_OA.M.PA_exNA <-na.omit(OTU_OA.M.PA, cols=seq_along(OTU_OA.M.PA))

## Make Mean relative abundance logarithmic
OTU_log.M.FL <-log(OTU_OA.M.FL_exNA[,2])
OTU_log.S.FL <-log(OTU_OA.S.FL_exNA[,2])
OTU_log.M.PA <-log(OTU_OA.M.PA_exNA[,2])
OTU_log.S.PA <-log(OTU_OA.S.PA_exNA[,2])

OTU_log.M.FL <-matrix(OTU_log.M.FL)
OTU_log.S.FL <-matrix(OTU_log.S.FL)
OTU_log.M.PA <-matrix(OTU_log.M.PA)
OTU_log.S.PA <-matrix(OTU_log.S.PA)

####### Calculate the slope of AO relationship #######
OTU_Slope_S.FL <-lm(OTU_log.S.FL ~ OTU_OA.S.FL_exNA[,1])$coeff[[2]]
OTU_Slope_M.FL <-lm(OTU_log.M.FL ~ OTU_OA.M.FL_exNA[,1])$coeff[[2]]
OTU_Slope_S.PA <-lm(OTU_log.S.PA ~ OTU_OA.S.PA_exNA[,1])$coeff[[2]]
OTU_Slope_M.PA <-lm(OTU_log.M.PA ~ OTU_OA.M.PA_exNA[,1])$coeff[[2]]

## Record the value of slope
OTU_Slope_S.FL   # 0.6927573
OTU_Slope_M.FL   # 0.6850764
OTU_Slope_S.PA   # 0.8103446
OTU_Slope_M.PA   # 0.6546258






###### Calculate the slop from the relationship of Abundance and Occupancy at the Genus level ######

## load Genus relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Genus <-read.csv("InputFiles/Genus_Table.csv", row.names=1, header=T) # Genus_RelAbun
Genus <-t(Genus)

## Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
Genus_RelAbun.M <-Genus[1:18,]
Genus_RelAbun.S <-Genus[19:36,]

## Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Genus_RelAbun.M.PA <-Genus_RelAbun.M[2:10,]
Genus_RelAbun.M.FL <-Genus_RelAbun.M[c(1, 11:18), ]
Genus_RelAbun.S.PA <-Genus_RelAbun.S[2:10,]
Genus_RelAbun.S.FL <-Genus_RelAbun.S[c(1, 11:18), ]

## Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Genus_RelAbun.M.PA <-t(Genus_RelAbun.M.PA)
Genus_RelAbun.M.FL <-t(Genus_RelAbun.M.FL)
Genus_RelAbun.S.PA <-t(Genus_RelAbun.S.PA)
Genus_RelAbun.S.FL <-t(Genus_RelAbun.S.FL)

## Calculate the regional mean abundance for every genus in each subset. This step assure 
Genus_MeanRelAbun.M.FL <-rowMeans(Genus_RelAbun.M.FL)
Genus_MeanRelAbun.S.FL <-rowMeans(Genus_RelAbun.S.FL)
Genus_MeanRelAbun.M.PA <-rowMeans(Genus_RelAbun.M.PA)
Genus_MeanRelAbun.S.PA <-rowMeans(Genus_RelAbun.S.PA)

## Count Rows that fit the criteria of Relative abundance =0 for each genus. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the Genus.
Genus_OCCU.M.PA <-rowSums(Genus_RelAbun.M.PA==0)
Genus_OCCU.M.FL <-rowSums(Genus_RelAbun.M.FL==0)
Genus_OCCU.S.PA <-rowSums(Genus_RelAbun.S.PA==0)
Genus_OCCU.S.FL <-rowSums(Genus_RelAbun.S.FL==0)

## Calcualte the number of sites that are occupied by the genus, by substracting the number of sites unoccupied from the total number of sites 9
Genus_OCCU.M.PA2 <-9-Genus_OCCU.M.PA
Genus_OCCU.M.FL2 <-9-Genus_OCCU.M.FL
Genus_OCCU.S.PA2 <-9-Genus_OCCU.S.PA
Genus_OCCU.S.FL2 <-9-Genus_OCCU.S.FL

## Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance genus did not occur by chance for once
Genus_OA.M.FL <-cbind(Genus_OCCU.M.FL2, Genus_MeanRelAbun.M.FL)
Genus_OA.S.FL <-cbind(Genus_OCCU.S.FL2, Genus_MeanRelAbun.S.FL)
Genus_OA.M.PA <-cbind(Genus_OCCU.M.PA2, Genus_MeanRelAbun.M.PA)
Genus_OA.S.PA <-cbind(Genus_OCCU.S.PA2, Genus_MeanRelAbun.S.PA)

colnames(Genus_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Genus_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Genus_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Genus_OA.S.PA) <-c("Occupancy","Relative abundance")

## Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Genus_OA.S.FL[Genus_OA.S.FL <= 0.00002]<-NA
Genus_OA.M.FL[Genus_OA.M.FL <= 0.00002]<-NA
Genus_OA.S.PA[Genus_OA.S.PA <= 0.00002]<-NA
Genus_OA.M.PA[Genus_OA.M.PA <= 0.00002]<-NA

## Remove the genus that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Genus_OA.S.FL_exNA <-na.omit(Genus_OA.S.FL, cols=seq_along(Genus_OA.S.FL))
Genus_OA.M.FL_exNA <-na.omit(Genus_OA.M.FL, cols=seq_along(Genus_OA.M.FL))
Genus_OA.S.PA_exNA <-na.omit(Genus_OA.S.PA, cols=seq_along(Genus_OA.S.PA))
Genus_OA.M.PA_exNA <-na.omit(Genus_OA.M.PA, cols=seq_along(Genus_OA.M.PA))

## Make Mean relative abundance logarithmic
Genus_log.M.FL <-log(Genus_OA.M.FL_exNA[,2])
Genus_log.S.FL <-log(Genus_OA.S.FL_exNA[,2])
Genus_log.M.PA <-log(Genus_OA.M.PA_exNA[,2])
Genus_log.S.PA <-log(Genus_OA.S.PA_exNA[,2])

Genus_log.M.FL <-matrix(Genus_log.M.FL)
Genus_log.S.FL <-matrix(Genus_log.S.FL)
Genus_log.M.PA <-matrix(Genus_log.M.PA)
Genus_log.S.PA <-matrix(Genus_log.S.PA)

####### Calculate the slope of AO relationship #######
Genus_Slope_S.FL <-lm(Genus_log.S.FL ~ Genus_OA.S.FL_exNA[,1])$coeff[[2]]
Genus_Slope_M.FL <-lm(Genus_log.M.FL ~ Genus_OA.M.FL_exNA[,1])$coeff[[2]]
Genus_Slope_S.PA <-lm(Genus_log.S.PA ~ Genus_OA.S.PA_exNA[,1])$coeff[[2]]
Genus_Slope_M.PA <-lm(Genus_log.M.PA ~ Genus_OA.M.PA_exNA[,1])$coeff[[2]]

## Record the value of slope
Genu_Slope_S.FL   # 0.8387731
Genus_Slope_M.FL  # 0.7858986
Genus_Slope_S.PA  # 0.7186197
Genus_Slope_M.PA  # 0.6978413







###### Calculate the slop from the relationship of Abundance and Occupancy at the Family level ######

## load family relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Family <-read.csv("InputFiles/Family_Table.csv", row.names=1, header=T) # Family_RelAbun
Family <-t(Family)

## Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
Family_RelAbun.M <-Family[1:18,]
Family_RelAbun.S <-Family[19:36,]

## Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Family_RelAbun.M.PA <-Family_RelAbun.M[2:10,]
Family_RelAbun.M.FL <-Family_RelAbun.M[c(1, 11:18), ]
Family_RelAbun.S.PA <-Family_RelAbun.S[2:10,]
Family_RelAbun.S.FL <-Family_RelAbun.S[c(1, 11:18), ]

## Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Family_RelAbun.M.PA <-t(Family_RelAbun.M.PA)
Family_RelAbun.M.FL <-t(Family_RelAbun.M.FL)
Family_RelAbun.S.PA <-t(Family_RelAbun.S.PA)
Family_RelAbun.S.FL <-t(Family_RelAbun.S.FL)

## Calculate the regional mean abundance for every family in each subset 
Family_MeanRelAbun.M.FL <-rowMeans(Family_RelAbun.M.FL)
Family_MeanRelAbun.S.FL <-rowMeans(Family_RelAbun.S.FL)
Family_MeanRelAbun.M.PA <-rowMeans(Family_RelAbun.M.PA)
Family_MeanRelAbun.S.PA <-rowMeans(Family_RelAbun.S.PA)

## Count Rows that fit the criteria of Relative abundance =0 for each family. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the family.
Family_OCCU.M.PA <-rowSums(Family_RelAbun.M.PA==0)
Family_OCCU.M.FL <-rowSums(Family_RelAbun.M.FL==0)
Family_OCCU.S.PA <-rowSums(Family_RelAbun.S.PA==0)
Family_OCCU.S.FL <-rowSums(Family_RelAbun.S.FL==0)

## Calcualte the number of sites that are occupied by the genus, by substracting the number of sites unoccupied from the total number of sites 9
Family_OCCU.M.PA2 <-9-Family_OCCU.M.PA
Family_OCCU.M.FL2 <-9-Family_OCCU.M.FL
Family_OCCU.S.PA2 <-9-Family_OCCU.S.PA
Family_OCCU.S.FL2 <-9-Family_OCCU.S.FL

## Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance family did not occur by chance for once
Family_OA.M.FL <-cbind(Family_OCCU.M.FL2, Family_MeanRelAbun.M.FL)
Family_OA.S.FL <-cbind(Family_OCCU.S.FL2, Family_MeanRelAbun.S.FL)
Family_OA.M.PA <-cbind(Family_OCCU.M.PA2, Family_MeanRelAbun.M.PA)
Family_OA.S.PA <-cbind(Family_OCCU.S.PA2, Family_MeanRelAbun.S.PA)

colnames(Family_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Family_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Family_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Family_OA.S.PA) <-c("Occupancy","Relative abundance")

## Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Family_OA.S.FL[Family_OA.S.FL <= 0.00002]<-NA
Family_OA.M.FL[Family_OA.M.FL <= 0.00002]<-NA
Family_OA.S.PA[Family_OA.S.PA <= 0.00002]<-NA
Family_OA.M.PA[Family_OA.M.PA <= 0.00002]<-NA

## Remove the genus that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Family_OA.S.FL_exNA <-na.omit(Family_OA.S.FL, cols=seq_along(Family_OA.S.FL))
Family_OA.M.FL_exNA <-na.omit(Family_OA.M.FL, cols=seq_along(Family_OA.M.FL))
Family_OA.S.PA_exNA <-na.omit(Family_OA.S.PA, cols=seq_along(Family_OA.S.PA))
Family_OA.M.PA_exNA <-na.omit(Family_OA.M.PA, cols=seq_along(Family_OA.M.PA))


## Make Mean relative abundance logarithmic
Family_log.M.FL <-log(Family_OA.M.FL_exNA[,2])
Family_log.S.FL <-log(Family_OA.S.FL_exNA[,2])
Family_log.M.PA <-log(Family_OA.M.PA_exNA[,2])
Family_log.S.PA <-log(Family_OA.S.PA_exNA[,2])

Family_log.M.FL <-matrix(Family_log.M.FL)
Family_log.S.FL <-matrix(Family_log.S.FL)
Family_log.M.PA <-matrix(Family_log.M.PA)
Family_log.S.PA <-matrix(Family_log.S.PA)


####### Calculate the slope of AO relationship #######
Family_Slope_S.FL <-lm(Familylog.S.FL ~ Family_OA.S.FL_exNA[,1])$coeff[[2]]
Family_Slope_M.FL <-lm(Family_log.M.FL ~ Family_OA.M.FL_exNA[,1])$coeff[[2]]
Family_Slope_S.PA <-lm(Family_log.S.PA ~ Familly_OA.S.PA_exNA[,1])$coeff[[2]]
Family_Slope_M.PA <-lm(Family_log.M.PA ~ Family_OA.M.PA_exNA[,1])$coeff[[2]]

## Record the value of slope
Family_Slope_S.FL  #0.811506
Family_Slope_M.FL  #0.8039124
Family_Slope_S.PA  #0.7286986
Family_Slope_M.PA  #0.7243899







###### Calculate the slop from the relationship of Abundance and Occupancy at the Order level ######

## load order relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Order <-read.csv("InputFiles/Order_Table.csv", row.names=1, header=T) # Order_RelAbun
Order <-t(Order)

## Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Order_RelAbun.M.PA <-Order_RelAbun.M[2:10,]
Order_RelAbun.M.FL <-Order_RelAbun.M[c(1, 11:18), ]
Order_RelAbun.S.PA <-Order_RelAbun.S[2:10,]
Order_RelAbun.S.FL <-Order_RelAbun.S[c(1, 11:18), ]

## Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Order_RelAbun.M.PA <-t(Order_RelAbun.M.PA)
Order_RelAbun.M.FL <-t(Order_RelAbun.M.FL)
Order_RelAbun.S.PA <-t(Order_RelAbun.S.PA)
Order_RelAbun.S.FL <-t(Order_RelAbun.S.FL)

## Calculate the regional mean abundance for every Order in each subset.
Order_MeanRelAbun.M.FL <-rowMeans(Order_RelAbun.M.FL)
Order_MeanRelAbun.S.FL <-rowMeans(Order_RelAbun.S.FL)
Order_MeanRelAbun.M.PA <-rowMeans(Order_RelAbun.M.PA)
Order_MeanRelAbun.S.PA <-rowMeans(Order_RelAbun.S.PA)

## Count Rows that fit the criteria of Relative abundance =0 for each Order. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the Order.
Order_OCCU.M.PA <-rowSums(Order_RelAbun.M.PA==0)
Order_OCCU.M.FL <-rowSums(Order_RelAbun.M.FL==0)
Order_OCCU.S.PA <-rowSums(Order_RelAbun.S.PA==0)
Order_OCCU.S.FL <-rowSums(Order_RelAbun.S.FL==0)

## Calcualte the number of sites that are occupied by the orders, by substracting the number of sites unoccupied from the total number of sites 9
Order_OCCU.M.PA2 <-9-Order_OCCU.M.PA
Order_OCCU.M.FL2 <-9-Order_OCCU.M.FL
Order_OCCU.S.PA2 <-9-Order_OCCU.S.PA
Order_OCCU.S.FL2 <-9-Order_OCCU.S.FL

## Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance orders did not occur by chance for once.
Order_OA.M.FL <-cbind(Order_OCCU.M.FL2, Order_MeanRelAbun.M.FL)
Order_OA.S.FL <-cbind(Order_OCCU.S.FL2, Order_MeanRelAbun.S.FL)
Order_OA.M.PA <-cbind(Order_OCCU.M.PA2, Order_MeanRelAbun.M.PA)
Order_OA.S.PA <-cbind(Order_OCCU.S.PA2, Order_MeanRelAbun.S.PA)

colnames(Order_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Order_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Order_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Order_OA.S.PA) <-c("Occupancy","Relative abundance")

## Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Order_OA.S.FL[Order_OA.S.FL <= 0.00002]<-NA
Order_OA.M.FL[Order_OA.M.FL <= 0.00002]<-NA
Order_OA.S.PA[Order_OA.S.PA <= 0.00002]<-NA
Order_OA.M.PA[Order_OA.M.PA <= 0.00002]<-NA

## Remove the OTUs that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Order_OA.S.FL_exNA <-na.omit(Order_OA.S.FL, cols=seq_along(Order_OA.S.FL))
Order_OA.M.FL_exNA <-na.omit(Order_OA.M.FL, cols=seq_along(Order_OA.M.FL))
Order_OA.S.PA_exNA <-na.omit(Order_OA.S.PA, cols=seq_along(Order_OA.S.PA))
Order_OA.M.PA_exNA <-na.omit(Order_OA.M.PA, cols=seq_along(Order_OA.M.PA))

## Make Mean relative abundance logarithmic
Order_log.M.FL <-log(Order_OA.M.FL_exNA[,2])
Order_log.S.FL <-log(Order_OA.S.FL_exNA[,2])
Order_log.M.PA <-log(Order_OA.M.PA_exNA[,2])
Order_log.S.PA <-log(Order_OA.S.PA_exNA[,2])

Order_log.M.FL <-matrix(Order_log.M.FL)
Order_log.S.FL <-matrix(Order_log.S.FL)
Order_log.M.PA <-matrix(Order_log.M.PA)
Order_log.S.PA <-matrix(Order_log.S.PA)

####### Calculate the slope of AO relationship #######
Order_Slope_S.FL <-lm(Order_log.S.FL ~ Order_OA.S.FL_exNA[,1])$coeff[[2]]
Order_Slope_M.FL <-lm(Order_log.M.FL ~ Order_OA.M.FL_exNA[,1])$coeff[[2]]
Order_Slope_S.PA <-lm(Order_log.S.PA ~ Order_OA.S.PA_exNA[,1])$coeff[[2]]
Order_Slope_M.PA <-lm(Order_log.M.PA ~ Order_OA.M.PA_exNA[,1])$coeff[[2]]

## Record the value of slope
Order_Slope_S.FL  #0.8308372
Order_Slope_M.FL  #0.8537709
Order_Slope_S.PA  #0.7801242
Order_Slope_M.PA  #0.7874378





###### Calculate the slop from the relationship of Abundance and Occupancy at the Class level #######

## load order relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Class <-read.csv("InputFiles/Class_Table.csv", row.names=1, header=T) # Class_RelAbun
Class <-t(Class)

## Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
Class_RelAbun.M <-Class_RelAbun[1:18,]
Class_RelAbun.S <-Class_RelAbun[19:36,]

## Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Class_RelAbun.M.PA <-Class_RelAbun.M[2:10,]
Class_RelAbun.M.FL <-Class_RelAbun.M[c(1, 11:18), ]
Class_RelAbun.S.PA <-Class_RelAbun.S[2:10,]
Class_RelAbun.S.FL <-Class_RelAbun.S[c(1, 11:18), ]

## Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Class_RelAbun.M.PA <-t(Class_RelAbun.M.PA)
Class_RelAbun.M.FL <-t(Class_RelAbun.M.FL)
Class_RelAbun.S.PA <-t(Class_RelAbun.S.PA)
Class_RelAbun.S.FL <-t(Class_RelAbun.S.FL)

## Calculate the regional mean abundance for every class in each subset.
Class_MeanRelAbun.M.FL <-rowMeans(Class_RelAbun.M.FL)
Class_MeanRelAbun.S.FL <-rowMeans(Class_RelAbun.S.FL)
Class_MeanRelAbun.M.PA <-rowMeans(Class_RelAbun.M.PA)
Class_MeanRelAbun.S.PA <-rowMeans(Class_RelAbun.S.PA)

## Count Rows that fit the criteria of Relative abundance =0 for each OTU. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the OTUs. .
Class_OCCU.M.PA <-rowSums(Class_RelAbun.M.PA==0)
Class_OCCU.M.FL <-rowSums(Class_RelAbun.M.FL==0)
Class_OCCU.S.PA <-rowSums(Class_RelAbun.S.PA==0)
Class_OCCU.S.FL <-rowSums(Class_RelAbun.S.FL==0)

## Calcualte the number of sites that are occupied by the OTUs, by substracting the number of sites unoccupied from the total number of sites 9
Class_OCCU.M.PA2 <-9-Class_OCCU.M.PA
Class_OCCU.M.FL2 <-9-Class_OCCU.M.FL
Class_OCCU.S.PA2 <-9-Class_OCCU.S.PA
Class_OCCU.S.FL2 <-9-Class_OCCU.S.FL

## Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance OTUs did not occur by chance for once.
Class_OA.M.FL <-cbind(Class_OCCU.M.FL2, Class_MeanRelAbun.M.FL)
Class_OA.S.FL <-cbind(Class_OCCU.S.FL2, Class_MeanRelAbun.S.FL)
Class_OA.M.PA <-cbind(Class_OCCU.M.PA2, Class_MeanRelAbun.M.PA)
Class_OA.S.PA <-cbind(Class_OCCU.S.PA2, Class_MeanRelAbun.S.PA)

colnames(Class_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Class_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Class_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Class_OA.S.PA) <-c("Occupancy","Relative abundance")

## Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Class_OA.S.FL[OTU_OA.S.FL <= 0.00002]<-NA
Class_OA.M.FL[OTU_OA.M.FL <= 0.00002]<-NA
Class_OA.S.PA[OTU_OA.S.PA <= 0.00002]<-NA
Class_OA.M.PA[OTU_OA.M.PA <= 0.00002]<-NA

## Remove the Class that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Class_OA.S.FL_exNA <-na.omit(Class_OA.S.FL, cols=seq_along(Class_OA.S.FL))
Class_OA.M.FL_exNA <-na.omit(Class_OA.M.FL, cols=seq_along(Class_OA.M.FL))
Class_OA.S.PA_exNA <-na.omit(Class_OA.S.PA, cols=seq_along(Class_OA.S.PA))
Class_OA.M.PA_exNA <-na.omit(Class_OA.M.PA, cols=seq_along(Class_OA.M.PA))

## Make Mean relative abundance logarithmic
Class_log.M.FL <-log(Class_OA.M.FL_exNA[,2])
Class_log.S.FL <-log(Class_OA.S.FL_exNA[,2])
Class_log.M.PA <-log(Class_OA.M.PA_exNA[,2])
Class_log.S.PA <-log(Class_OA.S.PA_exNA[,2])

Class_log.M.FL <-matrix(Class_log.M.FL)
Class_log.S.FL <-matrix(Class_log.S.FL)
Class_log.M.PA <-matrix(Class_log.M.PA)
Class_log.S.PA <-matrix(Class_log.S.PA)


####### Calculate the slope of AO relationship #######
Class_Slope_S.FL <-lm(Class_log.S.FL ~ Class_OA.S.FL_exNA[,1])$coeff[[2]]
Class_Slope_M.FL <-lm(Class_log.M.FL ~ Class_OA.M.FL_exNA[,1])$coeff[[2]]
Class_Slope_S.PA <-lm(Class_log.S.PA ~ Class_OA.S.PA_exNA[,1])$coeff[[2]]
Class_Slope_M.PA <-lm(Class_log.M.PA ~ Class_OA.M.PA_exNA[,1])$coeff[[2]]

## Record the value of slope
Class_Slope_S.FL  #0.9653827
Class_Slope_M.FL  #0.8946083
Class_Slope_S.PA  #0.895838
Class_Slope_M.PA  #0.9075367






###### Calculate the slop from the relationship of Abundance and Occupancy at the Phylum level ######

## load order relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Phylum <-read.csv("InputFiles/Phylum_Table.csv", row.names=1, header=T) # Phylum_RelAbun
Phylum <-t(Phylum)

## Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
Phylum_RelAbun.M <-Phylum_RelAbun[1:18,]
Phylum_RelAbun.S <-Phylum_RelAbun[19:36,]

## Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Phylum_RelAbun.M.PA <-Phylum_RelAbun.M[2:10,]
Phylum_RelAbun.M.FL <-Phylum_RelAbun.M[c(1, 11:18), ]
Phylum_RelAbun.S.PA <-Phylum_RelAbun.S[2:10,]
Phylum_RelAbun.S.FL <-Phylum_RelAbun.S[c(1, 11:18), ]

## Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Phylum_RelAbun.M.PA <-t(Phylum_RelAbun.M.PA)
Phylum_RelAbun.M.FL <-t(Phylum_RelAbun.M.FL)
Phylum_RelAbun.S.PA <-t(Phylum_RelAbun.S.PA)
Phylum_RelAbun.S.FL <-t(Phylum_RelAbun.S.FL)

## Calculate the regional mean abundance for every class in each subset.
Phylum_MeanRelAbun.M.FL <-rowMeans(Phylum_RelAbun.M.FL)
Phylum_MeanRelAbun.S.FL <-rowMeans(Phylum_RelAbun.S.FL)
Phylum_MeanRelAbun.M.PA <-rowMeans(Phylum_RelAbun.M.PA)
Phylum_MeanRelAbun.S.PA <-rowMeans(Phylum_RelAbun.S.PA)

## Count Rows that fit the criteria of Relative abundance =0 for each phylum. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the phyla.
Phylum_OCCU.M.PA <-rowSums(Phylum_RelAbun.M.PA==0)
Phylum_OCCU.M.FL <-rowSums(Phylum_RelAbun.M.FL==0)
Phylum_OCCU.S.PA <-rowSums(Phylum_RelAbun.S.PA==0)
Phylum_OCCU.S.FL <-rowSums(Phylum_RelAbun.S.FL==0)

## Calcualte the number of sites that are occupied by the OTUs, by substracting the number of sites unoccupied from the total number of sites 9
Phylum_OCCU.M.PA2 <-9-Phylum_OCCU.M.PA
Phylum_OCCU.M.FL2 <-9-Phylum_OCCU.M.FL
Phylum_OCCU.S.PA2 <-9-Phylum_OCCU.S.PA
Phylum_OCCU.S.FL2 <-9-Phylum_OCCU.S.FL

## Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance phyla did not occur by chance for once.
Phylum_OA.M.FL <-cbind(Phylum_OCCU.M.FL2, Phylum_MeanRelAbun.M.FL)
Phylum_OA.S.FL <-cbind(Phylum_OCCU.S.FL2, Phylum_MeanRelAbun.S.FL)
Phylum_OA.M.PA <-cbind(Phylum_OCCU.M.PA2, Phylum_MeanRelAbun.M.PA)
Phylum_OA.S.PA <-cbind(Phylum_OCCU.S.PA2, Phylum_MeanRelAbun.S.PA)

colnames(Phylum_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Phylum_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Phylum_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Phylum_OA.S.PA) <-c("Occupancy","Relative abundance")

## Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Phylum_OA.S.FL[Phylum_OA.S.FL <= 0.00002]<-NA
Phylum_OA.M.FL[Phylum_OA.M.FL <= 0.00002]<-NA
Phylum_OA.S.PA[Phylum_OA.S.PA <= 0.00002]<-NA
Phylum_OA.M.PA[Phylum_OA.M.PA <= 0.00002]<-NA


## Remove the phylum that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Phylum_OA.S.FL_exNA <-na.omit(Phylum_OA.S.FL, cols=seq_along(Phylum_OA.S.FL))
Phylum_OA.M.FL_exNA <-na.omit(Phylum_OA.M.FL, cols=seq_along(Phylum_OA.M.FL))
Phylum_OA.S.PA_exNA <-na.omit(Phylum_OA.S.PA, cols=seq_along(Phylum_OA.S.PA))
Phylum_OA.M.PA_exNA <-na.omit(Phylum_OA.M.PA, cols=seq_along(Phylum_OA.M.PA))

## Make Mean relative abundance logarithmic
Phylum_log.M.FL <-log(Phylum_OA.M.FL_exNA[,2])
Phylum_log.S.FL <-log(Phylum_OA.S.FL_exNA[,2])
Phylum_log.M.PA <-log(Phylum_OA.M.PA_exNA[,2])
Phylum_log.S.PA <-log(Phylum_OA.S.PA_exNA[,2])

Phylum_log.M.FL <-matrix(Phylum_log.M.FL)
Phylum_log.S.FL <-matrix(Phylum_log.S.FL)
Phylum_log.M.PA <-matrix(Phylum_log.M.PA)
Phylum_log.S.PA <-matrix(Phylum_log.S.PA)


####### Calculate the slope of AO relationship #######
Phylum_Slope_S.FL <-lm(Phylum_log.S.FL ~ Phylum_OA.S.FL_exNA[,1])$coeff[[2]]
Phylum_Slope_M.FL <-lm(Phylum_log.M.FL ~ Phylum_OA.M.FL_exNA[,1])$coeff[[2]]
Phylum_Slope_S.PA <-lm(Phylum_log.S.PA ~ Phylum_OA.S.PA_exNA[,1])$coeff[[2]]
Phylum_Slope_M.PA <-lm(Phylum_log.M.PA ~ Phylum_OA.M.PA_exNA[,1])$coeff[[2]]

## Record the value of slope
Phylum_Slope_S.FL  #1.153625
Phylum_Slope_M.FL  #1.046678
Phylum_Slope_S.PA  #0.9495059
Phylum_Slope_M.PA  #1.064442




####### Plot slope of Abundance-occupancy (AO) relationships across taxonomic ranks, which was used as the estimate of niche divergence rates######################

## load a table which includes all values of the slop of AO derived from subcommunities#######
SLOPE <-read.csv("InputFiles/SLOP.csv")
SLOPE$Ranks <-factor(SLOPE$Ranks, levels=c("Species","Genus", "Family","Order","Class","Phylum"))  ######## setting the levels of grouping according to taxonoic ranks
SLOPE$Lifestyles <-factor(SLOPE$Lifestyles, levels=c("FL","PA"))    ######## setting the levels of grouping according to lifestyles
SLOPE$Subcommunities <-factor(SLOPE$Subcommunities, levels=c("Surface_FL", "Surface_PA","CIL_FL", "CIL_PA"))      ######## setting the levels of grouping according to the community types
SLOPE$Depth <-factor(SLOPE$Depth, levels=c("Surface" ,"CIL"))     ######## setting the levels of grouping according to depths

library (ggplot2)

## Figure 6 

## different shape of symbols with optional for filled colours
p1 <-ggplot(SLOPE, aes(x=Ranks, y=Rates, group=Subcommunities, shape=Depth, fill=Lifestyles), colour="black") + geom_point(size=3.5) + scale_shape_manual(values=c(21,24)) + geom_line(size=0.2) +scale_fill_manual(values=c("orange", "darkolivegreen1")) + theme_bw() + labs(y="Niche divergence rates", x="Taxonomic ranks", title="", size=12) + theme(axis.title.y=element_text(colour = 'black',size = 11, face='bold')) + theme(axis.title.x=element_text(colour = 'black',size = 11, face='bold')) + theme(axis.text.y=element_text(colour = 'black', size = 11)) + theme(axis.text.x=element_text(colour = 'black', size = 11)) + theme(legend.title=element_blank(), legend.position="right") + theme(legend.text=element_text(size=9)) + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())\


## Figuer 6 inner panel

## draw boxplots to compared at which taxonomic ranks the difference in rates between the communities is larger 
Slope_box <-read.csv("InputFiles/boxplot.csv")
Slope_box$Ranks <-factor(Slope_box$Ranks, levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum"))

p2 <-ggplot(Slope_box, aes(x=Ranks, y=Values)) + geom_boxplot() + theme_bw() + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) +  labs(y="", x="", title="", size=11) + theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1), axis.title=element_blank())
p2




## Analysis of variance (ANOVA)-test for the significant differences in niche divergence rates between taxonomic ranks, in which case, the subcommunities function as replicates for the test###############

Slope_box$Ranks <-factor(Slope_box$Ranks, levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum"))
log_values <-log(Slope_box$Values)

## the normal distribution of the residuals of the linear models was tested using the Shapiro-Wilk normality test to fulfill the ANOVA's assumption
library(car)   ####### load this library for Shapiro-Wilk normality test

shapiro.test(log(Slope_box$Values))   
### Shapiro-Wilk normality test
##### output showed:
###### data:  log(Slope_box$Values)
###### W = 0.965, p-value = 0.5467    #######the residual of the model did not significantly differ from the normal distribution as P-value = 0.5467 indicated ###


## next step moves on to the the ANOVA test ##############

fit <-aov(formula=log(Slope_box$Values) ~ Slope_box$Ranks)

#### display the output 
fit
#Analysis of Variance Table

#Response: log(Slope_box$Values)
               # Df Sum Sq  Mean Sq F value    Pr(>F)    
# Slope_box$Ranks  5 0.4219 0.084380  17.437 2.416e-06 ***
# Residuals       18 0.0871 0.004839                      
# ---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## since there was a significant difference detected,move on to post doc test to detect which level the significance occured ######
TukeyHSD(fit)

##################### End #############
