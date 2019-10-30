
#####################################Calculation of AO at the Class level###################################
####################################################################################################################################

## load R libraries for this session 
library(scales)
library(base)
library(graphics)
library('RVAideMemoire')

## AO-relationsihp, habitat specialists, habitat generalists, based on Class level########################

## load order relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Class <-read.csv("InputFiles/Class_Table.csv", row.names=1, header=T) # Class_RelAbun
Class <-t(Class)


#################Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
Class_RelAbun.M <-Class_RelAbun[1:18,]
Class_RelAbun.S <-Class_RelAbun[19:36,]

#################Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Class_RelAbun.M.PA <-Class_RelAbun.M[2:10,]
Class_RelAbun.M.FL <-Class_RelAbun.M[c(1, 11:18), ]
Class_RelAbun.S.PA <-Class_RelAbun.S[2:10,]
Class_RelAbun.S.FL <-Class_RelAbun.S[c(1, 11:18), ]

#################Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Class_RelAbun.M.PA <-t(Class_RelAbun.M.PA)
Class_RelAbun.M.FL <-t(Class_RelAbun.M.FL)
Class_RelAbun.S.PA <-t(Class_RelAbun.S.PA)
Class_RelAbun.S.FL <-t(Class_RelAbun.S.FL)

#################Calculate the regional mean abundance for every class in each subset.
Class_MeanRelAbun.M.FL <-rowMeans(Class_RelAbun.M.FL)
Class_MeanRelAbun.S.FL <-rowMeans(Class_RelAbun.S.FL)
Class_MeanRelAbun.M.PA <-rowMeans(Class_RelAbun.M.PA)
Class_MeanRelAbun.S.PA <-rowMeans(Class_RelAbun.S.PA)

#################Count Rows that fit the criteria of Relative abundance =0 for each OTU. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the OTUs. .
Class_OCCU.M.PA <-rowSums(Class_RelAbun.M.PA==0)
Class_OCCU.M.FL <-rowSums(Class_RelAbun.M.FL==0)
Class_OCCU.S.PA <-rowSums(Class_RelAbun.S.PA==0)
Class_OCCU.S.FL <-rowSums(Class_RelAbun.S.FL==0)

#################Calcualte the number of sites that are occupied by the OTUs, by substracting the number of sites unoccupied from the total number of sites 9
Class_OCCU.M.PA2 <-9-Class_OCCU.M.PA
Class_OCCU.M.FL2 <-9-Class_OCCU.M.FL
Class_OCCU.S.PA2 <-9-Class_OCCU.S.PA
Class_OCCU.S.FL2 <-9-Class_OCCU.S.FL

#################Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance OTUs did not occur by chance for once.
Class_OA.M.FL <-cbind(Class_OCCU.M.FL2, Class_MeanRelAbun.M.FL)
Class_OA.S.FL <-cbind(Class_OCCU.S.FL2, Class_MeanRelAbun.S.FL)
Class_OA.M.PA <-cbind(Class_OCCU.M.PA2, Class_MeanRelAbun.M.PA)
Class_OA.S.PA <-cbind(Class_OCCU.S.PA2, Class_MeanRelAbun.S.PA)

colnames(Class_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Class_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Class_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Class_OA.S.PA) <-c("Occupancy","Relative abundance")

#################Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Class_OA.S.FL[OTU_OA.S.FL <= 0.00002]<-NA
Class_OA.M.FL[OTU_OA.M.FL <= 0.00002]<-NA
Class_OA.S.PA[OTU_OA.S.PA <= 0.00002]<-NA
Class_OA.M.PA[OTU_OA.M.PA <= 0.00002]<-NA


#################Remove the Class that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Class_OA.S.FL_exNA <-na.omit(Class_OA.S.FL, cols=seq_along(Class_OA.S.FL))
Class_OA.M.FL_exNA <-na.omit(Class_OA.M.FL, cols=seq_along(Class_OA.M.FL))
Class_OA.S.PA_exNA <-na.omit(Class_OA.S.PA, cols=seq_along(Class_OA.S.PA))
Class_OA.M.PA_exNA <-na.omit(Class_OA.M.PA, cols=seq_along(Class_OA.M.PA))

################Make Mean relative abundance logarithmic
Class_log.M.FL <-log(Class_OA.M.FL_exNA[,2])
Class_log.S.FL <-log(Class_OA.S.FL_exNA[,2])
Class_log.M.PA <-log(Class_OA.M.PA_exNA[,2])
Class_log.S.PA <-log(Class_OA.S.PA_exNA[,2])

Class_log.M.FL <-matrix(Class_log.M.FL)
Class_log.S.FL <-matrix(Class_log.S.FL)
Class_log.M.PA <-matrix(Class_log.M.PA)
Class_log.S.PA <-matrix(Class_log.S.PA)






###############################Calculation of Spearman's rank correlation at the Class level ########################################
####################################################################################################################################

Class_Spearman.S.FL <-cor.test(x=Class_OA.S.FL_exNA[,1], y=Class_log.S.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Class_Spearman.M.FL <-cor.test(x=Class_OA.M.FL_exNA[,1], y=Class_log.M.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Class_Spearman.S.PA <-cor.test(x=Class_OA.S.PA_exNA[,1], y=Class_log.S.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Class_Spearman.M.PA <-cor.test(x=Class_OA.M.PA_exNA[,1], y=Class_log.M.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)

######################## Vasulize the rho value and P value
Class_Spearman.S.FL #  rho=0.9311273, P value=2.2e-16
Class_Spearman.M.FL #  rho=0.9305386, P value=2.2e-16
Class_Spearman.S.PA # rho=0.945363, P-value=2.2e-16
Class_Spearman.M.PA # rho=0.9392132, P-value=2.2e-16


############################# Calculate confidence level of the interval， this is used to estimate the standard error of upper and lower limits of the 95% confidence interval. ########################################################
####################

library('RVAideMemoire')
Class_CON.S.FL <-spearman.ci(Class_OA.S.FL_exNA[,1], Class_log.S.FL, nrep=1000, conf.level=0.95)
Class_CON.M.FL <-spearman.ci(Class_OA.M.FL_exNA[,1], Class_log.M.FL, nrep=1000, conf.level=0.95)
Class_CON.S.PA <-spearman.ci(Class_OA.S.PA_exNA[,1], Class_log.S.PA, nrep=1000, conf.level=0.95)
Class_CON.M.PA <-spearman.ci(Class_OA.M.PA_exNA[,1], Class_log.M.PA, nrep=1000, conf.level=0.95)

######################################## output displaying the lower and upper limits of the 95% confidence interval, which was later used to calculate the standard errors for the coefficient plot ####################################
Class_CON.S.FL #0.8869516 0.9538740
Class_CON.M.FL # 0.8949453 0.9479686
Class_CON.S.PA #0.9154369 0.9588148
Class_CON.M.PA #0.9063822 0.9546337





############################### Classification of habitat generalists and habitat specialists###################
############################################################################################################

####################Sorting Generalists， that is classes occupying more than or equal to 6 sites 
Class_Gener.S.FL <-subset(Class_OA.S.FL_exNA, Class_OA.S.FL_exNA[,1]>=6)  ###### Gener refers to generalists
Class_Gener.M.FL <-subset(Class_OA.M.FL_exNA, Class_OA.M.FL_exNA[,1]>=6)
Class_Gener.S.PA <-subset(Class_OA.S.PA_exNA, Class_OA.S.PA_exNA[,1]>=6)
Class_Gener.M.PA <-subset(Class_OA.M.PA_exNA, Class_OA.M.PA_exNA[,1]>=6)

######################count the number of classes in that category 
Class_count.G.S.FL <-length(which(Class_OA.S.FL_exNA[,1]>=6))  
Class_count.G.M.FL <-length(which(Class_OA.M.FL_exNA[,1]>=6)) 
Class_count.G.S.PA <-length(which(Class_OA.S.PA_exNA[,1]>=6)) 
Class_count.G.M.PA <-length(which(Class_OA.M.PA_exNA[,1]>=6))  


Class_count.G.S.FL # 23
Class_count.G.M.FL # 39
Class_count.G.S.PA # 28
Class_count.G.M.PA # 48


########################Sorting Specialists that is classes occupying less than or equal to 2 sites 
Class_Spec.S.FL <-subset(Class_OA.S.FL_exNA, Class_OA.S.FL_exNA[,1]<=2)
Class_Spec.M.FL <-subset(Class_OA.M.FL_exNA, Class_OA.M.FL_exNA[,1]<=2)
Class_Spec.S.PA <-subset(Class_OA.S.PA_exNA, Class_OA.S.PA_exNA[,1]<=2)
Class_Spec.M.PA <-subset(Class_OA.M.PA_exNA, Class_OA.M.PA_exNA[,1]<=2)

######################count the number of classes in that category 
Class_count.S.S.FL <-length(which(Class_OA.S.FL_exNA[,1]<=2))
Class_count.S.M.FL <-length(which(Class_OA.M.FL_exNA[,1]<=2))
Class_count.S.S.PA <-length(which(Class_OA.S.PA_exNA[,1]<=2))
Class_count.S.M.PA <-length(which(Class_OA.M.PA_exNA[,1]<=2))

Class_count.S.S.FL # 20
Class_count.S.M.FL # 17
Class_count.S.S.PA # 33
Class_count.S.M.PA # 17


##################################Proportion of Generalists and Specialists, using the number of classes in Class_OA.S.FL_exNA for example #####################
Class_Proportion.G.S.FL <-(Class_count.G.S.FL/94)*100
Class_Proportion.G.M.FL <-(Class_count.G.M.FL/94)*100
Class_Proportion.G.S.PA <-(Class_count.G.S.PA/94)*100
Class_Proportion.G.M.PA <-(Class_count.G.M.PA/94)*100

Class_Proportion.S.S.FL <-(Class_count.S.S.FL/94)*100
Class_Proportion.S.M.FL <-(Class_count.S.M.FL/94)*100
Class_Proportion.S.S.PA <-(Class_count.S.S.PA/94)*100
Class_Proportion.S.M.PA <-(Class_count.S.M.PA/94)*100

Class_Proportion.S.S.FL # 21.2766
Class_Proportion.S.M.FL # 18.08511
Class_Proportion.S.S.PA # 35.10638
Class_Proportion.S.M.PA # 18.08511

Class_Proportion.G.S.FL # 24.46809
Class_Proportion.G.M.FL # 41.48936
Class_Proportion.G.S.PA # 29.78723
Class_Proportion.G.M.PA # 51.06383





################################################ Graphing for Abundance and Occupancy at the Class level ################################################################################################
   
###################### Plot mean relative abundance againt Occupancy for the different subcommunities #######################

part(mfrow=c(2,2))  ############## setting the 4 panels of figures in one page

library(scales)

## Supplementary Figure S3 E
#######################Surface -FL community########################################

myColor= c(ifelse(Class_OA.S.FL_exNA[,1]<3, "red", ifelse(Class_OA.S.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Class_OA.S.FL_exNA[,1]), Class_OA.S.FL_exNA[,2], log="y", title(main="Surface_FL at class level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "21.28%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "24.47%", cex=0.8)
text(7.5, 0.00009, "rho=0.931, P value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.FL)
dev.off()


########################CIL -FL community#######################################

myColor= c(ifelse(Class_OA.M.FL_exNA[,1]<3, "red", ifelse(Class_OA.M.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Class_OA.M.FL_exNA[,1]), Class_OA.M.FL_exNA[,2], log="y", title(main="CIL_FL at class level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "18.09%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "41.49%", cex=0.8)
text(7.5, 0.00009, "rho=0.931, P value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.FL)
dev.off()

########################Surface -PA community#######################################

myColor= c(ifelse(Class_OA.S.PA_exNA[,1]<3, "red", ifelse(Class_OA.S.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Class_OA.S.PA_exNA[,1]), Class_OA.S.PA_exNA[,2], log="y", title(main="Surface_PA at class level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "35.11%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "29.79%", cex=0.8)
text(7.5, 0.0004, "rho=0.945, P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.PA)
dev.off()


########################CIL -PA community#######################################

myColor= c(ifelse(Class_OA.M.PA_exNA[,1]<3, "red", ifelse(Class_OA.M.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Class_OA.M.PA_exNA[,1]), Class_OA.M.PA_exNA[,2], log="y", title(main="CIL_PA at class level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "18.09%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "51.06%", cex=0.8)
text(7.5, 0.0002, "rho=0.939, P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.PA)
dev.off()


############ End ############
