
#####################################Calculation of AO at the OTU level###################################
####################################################################################################################################################
########################## Again: OTU level and species level are interchangable in our case##########################

## load R libraries for this session 
library(scales)
library(base)
library(graphics)
library('RVAideMemoire')


############################# AO-relationsihp, habitat specialists, habitat generalists, based on OTU_level########################

OTU_RelAbun <-read.csv("InputFiles/RelOTU_table.csv", row.names=1, header=T) #######RelOTU_Table, without taxonomic information

#################Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
OTU_RelAbun.M <-OTU_RelAbun[1:18,]
OTU_RelAbun.S <-OTU_RelAbun[19:36,]

#################Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
OTU_RelAbun.M.PA <-OTU_RelAbun.M[2:10,]
OTU_RelAbun.M.FL <-OTU_RelAbun.M[c(1, 11:18), ]
OTU_RelAbun.S.PA <-OTU_RelAbun.S[2:10,]
OTU_RelAbun.S.FL <-OTU_RelAbun.S[c(1, 11:18), ]

#################Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
OTU_RelAbun.M.PA <-t(OTU_RelAbun.M.PA)
OTU_RelAbun.M.FL <-t(OTU_RelAbun.M.FL)
OTU_RelAbun.S.PA <-t(OTU_RelAbun.S.PA)
OTU_RelAbun.S.FL <-t(OTU_RelAbun.S.FL)

#################Calculate the regional mean abundance for every OTU in each subset.
OTU_MeanRelAbun.M.FL <-rowMeans(OTU_RelAbun.M.FL)
OTU_MeanRelAbun.S.FL <-rowMeans(OTU_RelAbun.S.FL)
OTU_MeanRelAbun.M.PA <-rowMeans(OTU_RelAbun.M.PA)
OTU_MeanRelAbun.S.PA <-rowMeans(OTU_RelAbun.S.PA)

#################Count Rows that fit the criteria of Relative abundance =0 for each OTU. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the OTUs. .
OTU_OCCU.M.PA <-rowSums(OTU_RelAbun.M.PA==0)
OTU_OCCU.M.FL <-rowSums(OTU_RelAbun.M.FL==0)
OTU_OCCU.S.PA <-rowSums(OTU_RelAbun.S.PA==0)
OTU_OCCU.S.FL <-rowSums(OTU_RelAbun.S.FL==0)

#################Calcualte the number of sites that are occupied by the OTUs, by substracting the number of sites unoccupied from the total number of sites 9
OTU_OCCU.M.PA2 <-9-OTU_OCCU.M.PA
OTU_OCCU.M.FL2 <-9-OTU_OCCU.M.FL
OTU_OCCU.S.PA2 <-9-OTU_OCCU.S.PA
OTU_OCCU.S.FL2 <-9-OTU_OCCU.S.FL

#################Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance OTUs did not occur by chance for once.
OTU_OA.M.FL <-cbind(OTU_OCCU.M.FL2, OTU_MeanRelAbun.M.FL)
OTU_OA.S.FL <-cbind(OTU_OCCU.S.FL2, OTU_MeanRelAbun.S.FL)
OTU_OA.M.PA <-cbind(OTU_OCCU.M.PA2, OTU_MeanRelAbun.M.PA)
OTU_OA.S.PA <-cbind(OTU_OCCU.S.PA2, OTU_MeanRelAbun.S.PA)

colnames(OTU_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(OTU_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(OTU_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(OTU_OA.S.PA) <-c("Occupancy","Relative abundance")

#################Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
OTU_OA.S.FL[OTU_OA.S.FL <= 0.00002]<-NA
OTU_OA.M.FL[OTU_OA.M.FL <= 0.00002]<-NA
OTU_OA.S.PA[OTU_OA.S.PA <= 0.00002]<-NA
OTU_OA.M.PA[OTU_OA.M.PA <= 0.00002]<-NA


#################Remove the OTUs that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
OTU_OA.S.FL_exNA <-na.omit(OTU_OA.S.FL, cols=seq_along(OTU_OA.S.FL))
OTU_OA.M.FL_exNA <-na.omit(OTU_OA.M.FL, cols=seq_along(OTU_OA.M.FL))
OTU_OA.S.PA_exNA <-na.omit(OTU_OA.S.PA, cols=seq_along(OTU_OA.S.PA))
OTU_OA.M.PA_exNA <-na.omit(OTU_OA.M.PA, cols=seq_along(OTU_OA.M.PA))

################Make Mean relative abundance logarithmic
OTU_log.M.FL <-log(OTU_OA.M.FL_exNA[,2])
OTU_log.S.FL <-log(OTU_OA.S.FL_exNA[,2])
OTU_log.M.PA <-log(OTU_OA.M.PA_exNA[,2])
OTU_log.S.PA <-log(OTU_OA.S.PA_exNA[,2])

OTU_log.M.FL <-matrix(OTU_log.M.FL)
OTU_log.S.FL <-matrix(OTU_log.S.FL)
OTU_log.M.PA <-matrix(OTU_log.M.PA)
OTU_log.S.PA <-matrix(OTU_log.S.PA)





 
###############################Calculation of Spearman's rank correlation at the OTU level ########################################
####################################################################################################################################

OTU_Spearman.S.FL <-cor.test(x=OTU_OA.S.FL_exNA[,1], y=OTU_log.S.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
OTU_Spearman.M.FL <-cor.test(x=OTU_OA.M.FL_exNA[,1], y=OTU_log.M.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
OTU_Spearman.S.PA <-cor.test(x=OTU_OA.S.PA_exNA[,1], y=OTU_log.S.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
OTU_Spearman.M.PA <-cor.test(x=OTU_OA.M.PA_exNA[,1], y=OTU_log.M.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)

######################## Vasulize the rho value and P value
OTU_Spearman.S.FL # rho 0.6473082 
OTU_Spearman.M.FL # rho 0.7162025
OTU_Spearman.S.PA # rho 0.6260733
OTU_Spearman.M.PA # rho 0.7194695




#############################Install.packages "RVAideMemoire" testing and plotting Procedures for clacuating confidence level of the interval， this is used to estimate the standard error of upper and lower limits of the 95% confidence interval. ########################################################
####################

Install.packages('RVAideMemoire')
library('RVAideMemoire')
OTU_CON.S.FL <-spearman.ci(OTU_OA.S.FL_exNA[,1], OTU_log.S.FL, nrep=1000, conf.level=0.95)
OTU_CON.M.FL <-spearman.ci(OTU_OA.M.FL_exNA[,1], OTU_log.M.FL, nrep=1000, conf.level=0.95)
OTU_CON.S.PA <-spearman.ci(OTU_OA.S.PA_exNA[,1], OTU_log.S.PA, nrep=1000, conf.level=0.95)
OTU_CON.M.PA <-spearman.ci(OTU_OA.M.PA_exNA[,1], OTU_log.M.PA, nrep=1000, conf.level=0.95)


######################################## output displaying the lower and upper limits of the 95% confidence interval, which was later used to calculate the standard errors for the coefficient plot ####################################
OTU_CON.S.FL   # 95 percent confidence interval:0.6370073 0.6573729
OTU_CON.M.FL   # 95 percent confidence interval:0.7071199 0.7250515
OTU_CON.S.PA   # 95 percent confidence interval:0.6112686 0.6405946
OTU_CON.M.PA   # 95 percent confidence interval:0.7084961 0.7302933





############################### Classification of habitat generalists and habitat specialists###################
############################################################################################################

####################Sorting Generalists， that is OTUs occupying more than or equal to 6 sites 
OTU_Gener.S.FL <-subset(OTU_OA.S.FL_exNA, OTU_OA.S.FL_exNA[,1]>=6)  ###### Gener refers to generalists
OTU_Gener.M.FL <-subset(OTU_OA.M.FL_exNA, OTU_OA.M.FL_exNA[,1]>=6)
OTU_Gener.S.PA <-subset(OTU_OA.S.PA_exNA, OTU_OA.S.PA_exNA[,1]>=6)
OTU_Gener.M.PA <-subset(OTU_OA.M.PA_exNA, OTU_OA.M.PA_exNA[,1]>=6)

######################count the number of OTUs in that category 
OTU_count.G.S.FL <-length(which(OTU_OA.S.FL_exNA[,1]>=6))  
OTU_count.G.M.FL <-length(which(OTU_OA.M.FL_exNA[,1]>=6)) 
OTU_count.G.S.PA <-length(which(OTU_OA.S.PA_exNA[,1]>=6)) 
OTU_count.G.M.PA <-length(which(OTU_OA.M.PA_exNA[,1]>=6))  

OTU_count.G.S.FL   # 443
OTU_count.G.M.FL   # 640
OTU_count.G.S.PA   # 280
OTU_count.G.M.PA   # 741


########################Sorting Specialists that is OTUs occupying less than or equal to 2 sites 
OTU_Spec.S.FL <-subset(OTU_OA.S.FL_exNA, OTU_OA.S.FL_exNA[,1]<=2)
OTU_Spec.M.FL <-subset(OTU_OA.M.FL_exNA, OTU_OA.M.FL_exNA[,1]<=2)
OTU_Spec.S.PA <-subset(OTU_OA.S.PA_exNA, OTU_OA.S.PA_exNA[,1]<=2)
OTU_Spec.M.PA <-subset(OTU_OA.M.PA_exNA, OTU_OA.M.PA_exNA[,1]<=2)

######################count the number of OTUs in that category 
OTU_count.S.S.FL <-length(which(OTU_OA.S.FL_exNA[,1]<=2))
OTU_count.S.M.FL <-length(which(OTU_OA.M.FL_exNA[,1]<=2))
OTU_count.S.S.PA <-length(which(OTU_OA.S.PA_exNA[,1]<=2))
OTU_count.S.M.PA <-length(which(OTU_OA.M.PA_exNA[,1]<=2))

OTU_count.S.S.FL  #12126
OTU_count.S.M.FL  #11480
OTU_count.S.S.PA  #6773
OTU_count.S.M.PA  #7934



##################################Proportion of Generalists and Specialists, using the number of OTU in OTU_OA.S.FL_exNA for example #####################
OTU_Proportion.G.S.FL <-(OTU_count.G.S.FL/OTU_OA.S.FL_exNA)*100   ####### calcualte for generalists
OTU_Proportion.G.M.FL <-(OTU_count.G.M.FL/OTU_OA.M.FL_exNA)*100
OTU_Proportion.G.S.PA <-(OTU_count.G.S.PA/OTU_OA.S.PA_exNA)*100
OTU_Proportion.G.M.PA <-(OTU_count.G.M.PA/OTU_OA.M.PA_exNA)*100

OTU_Proportion.S.S.FL <-(OTU_count.S.S.FL/OTU_OA.S.FL_exNA)*100    ####### calcualte for specialists   
OTU_Proportion.S.M.FL <-(OTU_count.S.M.FL/OTU_OA.M.FL_exNA)*100 
OTU_Proportion.S.S.PA <-(OTU_count.S.S.PA/OTU_OA.S.PA_exNA)*100
OTU_Proportion.S.M.PA <-(OTU_count.S.M.PA/OTU_OA.M.PA_exNA)*100

OTU_Proportion.S.S.FL  # 83.96344
OTU_Proportion.S.M.FL  # 79.89422
OTU_Proportion.S.S.PA  # 87.74453
OTU_Proportion.S.M.PA  # 75.0402

OTU_Proportion.G.S.FL  # 3.067442
OTU_Proportion.G.M.FL  # 4.454033
OTU_Proportion.G.S.PA  # 3.627413
OTU_Proportion.G.M.PA  # 7.008418





################################################ Graphing for Abundance and Occupancy at the OTU level ################################################################################################
   
######################### Plot mean relative abundance againt Occupancy for the different subcommunities #######################
part(mfrow=c(2,2))  ############## setting the 4 panels of figures in one page

library(scales)

## Supplementary Figure S3 A
#######################Surface -FL community ########################################

myColor= c(ifelse(OTU_OA.S.FL_exNA[,1]<3, "red", ifelse(OTU_OA.S.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(OTU_OA.S.FL_exNA[,1]), OTU_OA.S.FL_exNA[,2], log="y", title(main="Surface_FL at OTU level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "83.96%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "3.06%", cex=0.8)
text(7.5, 0.00012, "rho=0.647,  P=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
dev.off()


########################CIL -FL community#######################################

myColor= c(ifelse(OTU_OA.M.FL_exNA[,1]<3, "red", ifelse(OTU_OA.M.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(OTU_OA.M.FL_exNA[,1]), OTU_OA.M.FL_exNA[,2], log="y", title(main="CIL_FL at OTU level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "79.89%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "4.45%", cex=0.8)
text(7.5, 0.00012, "rho=0.716,  P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
dev.off()


########################Surface -PA community#######################################

myColor= c(ifelse(OTU_OA.S.PA_exNA[,1]<3, "red", ifelse(OTU_OA.S.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(OTU_OA.S.PA_exNA[,1]), OTU_OA.S.PA_exNA[,2], log="y", title(main="Surface_PA at OTU level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "87.74%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "3.63%", cex=0.8)
text(7.5, 0.00012, "rho=0.626,  P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
dev.off()

########################CIL -PA community#######################################

myColor= c(ifelse(OTU_OA.M.PA_exNA[,1]<3, "red", ifelse(OTU_OA.M.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(OTU_OA.M.PA_exNA[,1]), OTU_OA.M.PA_exNA[,2], log="y", title(main="CIL_PA at OTU level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "75.04%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "7.01%", cex=0.8)
text(7.5, 0.00012, "rho=0.719,  P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
dev.off()


############## End ############
