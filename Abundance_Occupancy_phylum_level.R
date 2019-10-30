
#####################################Calculation of AO at the Phylum level###################################
####################################################################################################################################################

## load R libraries for this session 
library(scales)
library(base)
library(graphics)
library('RVAideMemoire')

############################# AO-relationsihp, habitat specialists, habitat generalists, based on Phylum level########################

######################## load order relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Phylum <-read.csv("InputFiles/Phylum_Table.csv", row.names=1, header=T) # Phylum_RelAbun
Phylum <-t(Phylum)


#################Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
Phylum_RelAbun.M <-Phylum_RelAbun[1:18,]
Phylum_RelAbun.S <-Phylum_RelAbun[19:36,]

#################Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Phylum_RelAbun.M.PA <-Phylum_RelAbun.M[2:10,]
Phylum_RelAbun.M.FL <-Phylum_RelAbun.M[c(1, 11:18), ]
Phylum_RelAbun.S.PA <-Phylum_RelAbun.S[2:10,]
Phylum_RelAbun.S.FL <-Phylum_RelAbun.S[c(1, 11:18), ]

#################Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Phylum_RelAbun.M.PA <-t(Phylum_RelAbun.M.PA)
Phylum_RelAbun.M.FL <-t(Phylum_RelAbun.M.FL)
Phylum_RelAbun.S.PA <-t(Phylum_RelAbun.S.PA)
Phylum_RelAbun.S.FL <-t(Phylum_RelAbun.S.FL)

#################Calculate the regional mean abundance for every class in each subset.
Phylum_MeanRelAbun.M.FL <-rowMeans(Phylum_RelAbun.M.FL)
Phylum_MeanRelAbun.S.FL <-rowMeans(Phylum_RelAbun.S.FL)
Phylum_MeanRelAbun.M.PA <-rowMeans(Phylum_RelAbun.M.PA)
Phylum_MeanRelAbun.S.PA <-rowMeans(Phylum_RelAbun.S.PA)

#################Count Rows that fit the criteria of Relative abundance =0 for each phylum. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the phyla.
Phylum_OCCU.M.PA <-rowSums(Phylum_RelAbun.M.PA==0)
Phylum_OCCU.M.FL <-rowSums(Phylum_RelAbun.M.FL==0)
Phylum_OCCU.S.PA <-rowSums(Phylum_RelAbun.S.PA==0)
Phylum_OCCU.S.FL <-rowSums(Phylum_RelAbun.S.FL==0)

#################Calcualte the number of sites that are occupied by the OTUs, by substracting the number of sites unoccupied from the total number of sites 9
Phylum_OCCU.M.PA2 <-9-Phylum_OCCU.M.PA
Phylum_OCCU.M.FL2 <-9-Phylum_OCCU.M.FL
Phylum_OCCU.S.PA2 <-9-Phylum_OCCU.S.PA
Phylum_OCCU.S.FL2 <-9-Phylum_OCCU.S.FL

#################Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance phyla did not occur by chance for once.
Phylum_OA.M.FL <-cbind(Phylum_OCCU.M.FL2, Phylum_MeanRelAbun.M.FL)
Phylum_OA.S.FL <-cbind(Phylum_OCCU.S.FL2, Phylum_MeanRelAbun.S.FL)
Phylum_OA.M.PA <-cbind(Phylum_OCCU.M.PA2, Phylum_MeanRelAbun.M.PA)
Phylum_OA.S.PA <-cbind(Phylum_OCCU.S.PA2, Phylum_MeanRelAbun.S.PA)

colnames(Phylum_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Phylum_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Phylum_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Phylum_OA.S.PA) <-c("Occupancy","Relative abundance")

#################Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Phylum_OA.S.FL[Phylum_OA.S.FL <= 0.00002]<-NA
Phylum_OA.M.FL[Phylum_OA.M.FL <= 0.00002]<-NA
Phylum_OA.S.PA[Phylum_OA.S.PA <= 0.00002]<-NA
Phylum_OA.M.PA[Phylum_OA.M.PA <= 0.00002]<-NA


#################Remove the phylum that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Phylum_OA.S.FL_exNA <-na.omit(Phylum_OA.S.FL, cols=seq_along(Phylum_OA.S.FL))
Phylum_OA.M.FL_exNA <-na.omit(Phylum_OA.M.FL, cols=seq_along(Phylum_OA.M.FL))
Phylum_OA.S.PA_exNA <-na.omit(Phylum_OA.S.PA, cols=seq_along(Phylum_OA.S.PA))
Phylum_OA.M.PA_exNA <-na.omit(Phylum_OA.M.PA, cols=seq_along(Phylum_OA.M.PA))

################Make Mean relative abundance logarithmic
Phylum_log.M.FL <-log(Phylum_OA.M.FL_exNA[,2])
Phylum_log.S.FL <-log(Phylum_OA.S.FL_exNA[,2])
Phylum_log.M.PA <-log(Phylum_OA.M.PA_exNA[,2])
Phylum_log.S.PA <-log(Phylum_OA.S.PA_exNA[,2])

Phylum_log.M.FL <-matrix(Phylum_log.M.FL)
Phylum_log.S.FL <-matrix(Phylum_log.S.FL)
Phylum_log.M.PA <-matrix(Phylum_log.M.PA)
Phylum_log.S.PA <-matrix(Phylum_log.S.PA)





###############################Calculation of Spearman's rank correlation at the Phylum level ########################################
####################################################################################################################################

Phylum_Spearman.S.FL <-cor.test(x=Phylum_OA.S.FL_exNA[,1], y=Phylum_log.S.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Phylum_Spearman.M.FL <-cor.test(x=Phylum_OA.M.FL_exNA[,1], y=Phylum_log.M.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Phylum_Spearman.S.PA <-cor.test(x=Phylum_OA.S.PA_exNA[,1], y=Phylum_log.S.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Phylum_Spearman.M.PA <-cor.test(x=Phylum_OA.M.PA_exNA[,1], y=Phylum_log.M.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)


######################## Vasulize the rho value and P value
Phylum_Spearman.S.FL #  rho=0.9041778 , P value=9.467e-09
Phylum_Spearman.M.FL #  rho=0.8635437 , P value=1.136e-07
Phylum_Spearman.S.PA # rho=0.9360732, P-value=1.106e-12
Phylum_Spearman.M.PA # rho=0.8705646, P-value=1.791e-09





############################# Calculate confidence level of the interval， this is used to estimate the standard error of upper and lower limits of the 95% confidence interval. ########################################################
####################

library('RVAideMemoire')
Phylum_CON.S.FL <-spearman.ci(Phylum_OA.S.FL_exNA[,1], Phylum_log.S.FL, nrep=1000, conf.level=0.95)
Phylum_CON.M.FL <-spearman.ci(Phylum_OA.M.FL_exNA[,1], Phylum_log.M.FL, nrep=1000, conf.level=0.95)
Phylum_CON.S.PA <-spearman.ci(Phylum_OA.S.PA_exNA[,1], Phylum_log.S.PA, nrep=1000, conf.level=0.95)
Phylum_CON.M.PA <-spearman.ci(Phylum_OA.M.PA_exNA[,1], Phylum_log.M.PA, nrep=1000, conf.level=0.95)


## output displaying the lower and upper limits of the 95% confidence interval, which was later used to calculate the standard errors for the coefficient plot ####################################
Phylum_CON.S.FL # 0.7435973 0.9552152
Phylum_CON.M.FL # 0.6752815 0.9534084
Phylum_CON.S.PA # 0.8440948 0.9702454
Phylum_CON.M.PA # 0.6863909 0.9478254





############################### Classification of habitat generalists and habitat specialists###################
############################################################################################################

####################Sorting Generalists， that is phyla occupying more than or equal to 6 sites 
Phylum_Gener.S.FL <-subset(Phylum_OA.S.FL_exNA, Phylum_OA.S.FL_exNA[,1]>=6)
Phylum_Gener.M.FL <-subset(Phylum_OA.M.FL_exNA, Phylum_OA.M.FL_exNA[,1]>=6)
Phylum_Gener.S.PA <-subset(Phylum_OA.S.PA_exNA, Phylum_OA.S.PA_exNA[,1]>=6)
Phylum_Gener.M.PA <-subset(Phylum_OA.M.PA_exNA, Phylum_OA.M.PA_exNA[,1]>=6)

Phylum_count.G.S.FL <-length(which(Phylum_OA.S.FL_exNA[,1]>=6)) # count the number of phyla in that category 
Phylum_count.G.M.FL <-length(which(Phylum_OA.M.FL_exNA[,1]>=6)) # 
Phylum_count.G.S.PA <-length(which(Phylum_OA.S.PA_exNA[,1]>=6)) # 
Phylum_count.G.M.PA <-length(which(Phylum_OA.M.PA_exNA[,1]>=6))  

Phylum_count.G.S.FL # 9
Phylum_count.G.M.FL # 16
Phylum_count.G.S.PA # 13
Phylum_count.G.M.PA # 21

########################Sorting Specialists that is phyla occupying less than or equal to 2 sites 
Phylum_Spec.S.FL <-subset(Phylum_OA.S.FL_exNA, Phylum_OA.S.FL_exNA[,1]<=2)
Phylum_Spec.M.FL <-subset(Phylum_OA.M.FL_exNA, Phylum_OA.M.FL_exNA[,1]<=2)
Phylum_Spec.S.PA <-subset(Phylum_OA.S.PA_exNA, Phylum_OA.S.PA_exNA[,1]<=2)
Phylum_Spec.M.PA <-subset(Phylum_OA.M.PA_exNA, Phylum_OA.M.PA_exNA[,1]<=2)

######################count the number of phyla in that category 
Phylum_count.S.S.FL <-length(which(Phylum_OA.S.FL_exNA[,1]<=2))
Phylum_count.S.M.FL <-length(which(Phylum_OA.M.FL_exNA[,1]<=2))
Phylum_count.S.S.PA <-length(which(Phylum_OA.S.PA_exNA[,1]<=2))
Phylum_count.S.M.PA <-length(which(Phylum_OA.M.PA_exNA[,1]<=2))

Phylum_count.S.S.FL # 5
Phylum_count.S.M.FL # 2
Phylum_count.S.S.PA # 6
Phylum_count.S.M.PA # 1


##################################Proportion of Generalists and Specialists, using the number of phyla in Phylum_OA.S.FL_exNA for example #####################

Phylum_Proportion.G.S.FL <-(Phylum_count.G.S.FL/29)*100
Phylum_Proportion.G.M.FL <-(Phylum_count.G.M.FL/29)*100
Phylum_Proportion.G.S.PA <-(Phylum_count.G.S.PA/29)*100
Phylum_Proportion.G.M.PA <-(Phylum_count.G.M.PA/29)*100

Phylum_Proportion.S.S.FL <-(Phylum_count.S.S.FL/29)*100
Phylum_Proportion.S.M.FL <-(Phylum_count.S.M.FL/29)*100
Phylum_Proportion.S.S.PA <-(Phylum_count.S.S.PA/29)*100
Phylum_Proportion.S.M.PA <-(Phylum_count.S.M.PA/29)*100

Phylum_Proportion.S.S.FL # 17.24138
Phylum_Proportion.S.M.FL # 6.896552
Phylum_Proportion.S.S.PA # 20.68966
Phylum_Proportion.S.M.PA # 3.448276

Phylum_Proportion.G.S.FL # 31.03448
Phylum_Proportion.G.M.FL # 55.17241
Phylum_Proportion.G.S.PA # 44.82759
Phylum_Proportion.G.M.PA # 72.41379






################################################ Graphing for Abundance and Occupancy at the Phylum level ################################################################################################
   
##### Plot mean relative abundance againt Occupancy for the different subcommunities #######################

part(mfrow=c(2,2))  ############## setting the 4 panels of figures in one page

library(scales)


## Supplementary Figure S3 F
#######################Surface -FL community########################################

myColor= c(ifelse(Phylum_OA.S.FL_exNA[,1]<3, "red", ifelse(Phylum_OA.S.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Phylum_OA.S.FL_exNA[,1]), Phylum_OA.S.FL_exNA[,2], log="y", title(main="Surface_FL at phylum level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "17.24%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "31.03%", cex=0.8)
text(7.5, 0.000095, "rho=0.904, P value=9.467e-09", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.FL)
dev.off()


########################CIL -FL community#######################################

myColor= c(ifelse(Phylum_OA.M.FL_exNA[,1]<3, "red", ifelse(Phylum_OA.M.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Phylum_OA.M.FL_exNA[,1]), Phylum_OA.M.FL_exNA[,2], log="y", title(main="CIL_FL at phylum level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "6.90%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "55.17%", cex=0.8)
text(7.5, 0.0003, "rho=0.864, P-value=1.136e-07", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.FL)
dev.off()

########################Surface -PA community#######################################

myColor= c(ifelse(Phylum_OA.S.PA_exNA[,1]<3, "red", ifelse(Phylum_OA.S.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Phylum_OA.S.PA_exNA[,1]), Phylum_OA.S.PA_exNA[,2], log="y", title(main="Surface_PA at phylum level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "20.69%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "44.83%", cex=0.8)
text(7.5, 0.0011, "rho=0.936, P-value=1.106e-12", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.PA)
dev.off()

########################CIL -PA community#######################################

myColor= c(ifelse(Phylum_OA.M.PA_exNA[,1]<3, "red", ifelse(Phylum_OA.M.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Phylum_OA.M.PA_exNA[,1]), Phylum_OA.M.PA_exNA[,2], log="y", title(main="CIL_PA at phylum level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "3.45%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "72.41%", cex=0.8)
text(7.5, 0.00055, "rho=0.871, P-value=1.791e-09", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.PA)
dev.off()


########## End #########