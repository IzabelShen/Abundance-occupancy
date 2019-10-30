
#####################################Calculation of AO at the Family level###################################
####################################################################################################################################################

## load R libraries for this session 
library(scales)
library(base)
library(graphics)
library('RVAideMemoire')

############################# AO-relationsihp, habitat specialists, habitat generalists, based on Family level########################

######################## load family relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Family <-read.csv("InputFiles/Family_Table.csv", row.names=1, header=T) # Family_RelAbun
Family <-t(Family)


#################Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
Family_RelAbun.M <-Family[1:18,]
Family_RelAbun.S <-Family[19:36,]

#################Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Family_RelAbun.M.PA <-Family_RelAbun.M[2:10,]
Family_RelAbun.M.FL <-Family_RelAbun.M[c(1, 11:18), ]
Family_RelAbun.S.PA <-Family_RelAbun.S[2:10,]
Family_RelAbun.S.FL <-Family_RelAbun.S[c(1, 11:18), ]

#################Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Family_RelAbun.M.PA <-t(Family_RelAbun.M.PA)
Family_RelAbun.M.FL <-t(Family_RelAbun.M.FL)
Family_RelAbun.S.PA <-t(Family_RelAbun.S.PA)
Family_RelAbun.S.FL <-t(Family_RelAbun.S.FL)

#################Calculate the regional mean abundance for every family in each subset 
Family_MeanRelAbun.M.FL <-rowMeans(Family_RelAbun.M.FL)
Family_MeanRelAbun.S.FL <-rowMeans(Family_RelAbun.S.FL)
Family_MeanRelAbun.M.PA <-rowMeans(Family_RelAbun.M.PA)
Family_MeanRelAbun.S.PA <-rowMeans(Family_RelAbun.S.PA)

#################Count Rows that fit the criteria of Relative abundance =0 for each family. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the family.
Family_OCCU.M.PA <-rowSums(Family_RelAbun.M.PA==0)
Family_OCCU.M.FL <-rowSums(Family_RelAbun.M.FL==0)
Family_OCCU.S.PA <-rowSums(Family_RelAbun.S.PA==0)
Family_OCCU.S.FL <-rowSums(Family_RelAbun.S.FL==0)

#################Calcualte the number of sites that are occupied by the genus, by substracting the number of sites unoccupied from the total number of sites 9
Family_OCCU.M.PA2 <-9-Family_OCCU.M.PA
Family_OCCU.M.FL2 <-9-Family_OCCU.M.FL
Family_OCCU.S.PA2 <-9-Family_OCCU.S.PA
Family_OCCU.S.FL2 <-9-Family_OCCU.S.FL

#################Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance family did not occur by chance for once
Family_OA.M.FL <-cbind(Family_OCCU.M.FL2, Family_MeanRelAbun.M.FL)
Family_OA.S.FL <-cbind(Family_OCCU.S.FL2, Family_MeanRelAbun.S.FL)
Family_OA.M.PA <-cbind(Family_OCCU.M.PA2, Family_MeanRelAbun.M.PA)
Family_OA.S.PA <-cbind(Family_OCCU.S.PA2, Family_MeanRelAbun.S.PA)

colnames(Family_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Family_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Family_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Family_OA.S.PA) <-c("Occupancy","Relative abundance")

#################Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Family_OA.S.FL[Family_OA.S.FL <= 0.00002]<-NA
Family_OA.M.FL[Family_OA.M.FL <= 0.00002]<-NA
Family_OA.S.PA[Family_OA.S.PA <= 0.00002]<-NA
Family_OA.M.PA[Family_OA.M.PA <= 0.00002]<-NA


#################Remove the genus that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Family_OA.S.FL_exNA <-na.omit(Family_OA.S.FL, cols=seq_along(Family_OA.S.FL))
Family_OA.M.FL_exNA <-na.omit(Family_OA.M.FL, cols=seq_along(Family_OA.M.FL))
Family_OA.S.PA_exNA <-na.omit(Family_OA.S.PA, cols=seq_along(Family_OA.S.PA))
Family_OA.M.PA_exNA <-na.omit(Family_OA.M.PA, cols=seq_along(Family_OA.M.PA))


################Make Mean relative abundance logarithmic
Family_log.M.FL <-log(Family_OA.M.FL_exNA[,2])
Family_log.S.FL <-log(Family_OA.S.FL_exNA[,2])
Family_log.M.PA <-log(Family_OA.M.PA_exNA[,2])
Family_log.S.PA <-log(Family_OA.S.PA_exNA[,2])

Family_log.M.FL <-matrix(Family_log.M.FL)
Family_log.S.FL <-matrix(Family_log.S.FL)
Family_log.M.PA <-matrix(Family_log.M.PA)
Family_log.S.PA <-matrix(Family_log.S.PA)





###############################Calculation of Spearman's rank correlation at the Family level ########################################
####################################################################################################################################

Family_Spearman.S.FL <-cor.test(x=Family_OA.S.FL_exNA[,1], y=Family_log.S.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Family_Spearman.M.FL <-cor.test(x=Family_OA.M.FL_exNA[,1], y=Family_log.M.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Family_Spearman.S.PA <-cor.test(x=Family_OA.S.PA_exNA[,1], y=Family_log.S.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Family_Spearman.M.PA <-cor.test(x=Family_OA.M.PA_exNA[,1], y=Family_log.M.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)

######################## Vasulize the rho value and P value
Family_Spearman.S.FL #  rho=0.9074308, P-value=2.2e-16
Family_Spearman.M.FL #  rho=0.8974094, P-value=2.2e-16
Family_Spearman.S.PA # rho=0.9216762, P-value=2.2e-16
Family_Spearman.M.PA # rho=0.9351546, P-value= 2.2e-16





############################# Calculate confidence level of the interval， this is used to estimate the standard error of upper and lower limits of the 95% confidence interval. ########################################################
####################
library('RVAideMemoire')

Family_CON.S.FL <-spearman.ci(Family_OA.S.FL_exNA[,1], Family_log.S.FL, nrep=1000, conf.level=0.95)
Family_CON.M.FL <-spearman.ci(Family_OA.M.FL_exNA[,1], Family_log.M.FL, nrep=1000, conf.level=0.95)
Family_CON.S.PA <-spearman.ci(Family_OA.S.PA_exNA[,1], Family_log.S.PA, nrep=1000, conf.level=0.95)
Family_CON.M.PA <-spearman.ci(Family_OA.M.PA_exNA[,1], Family_log.M.PA, nrep=1000, conf.level=0.95)

######################################## output displaying the lower and upper limits of the 95% confidence interval, which was later used to calculate the standard errors for the coefficient plot ####################################
Family_CON.S.FL # 0.8813844 0.9271879
Family_CON.M.FL # 0.8640309 0.9209236
Family_CON.S.PA # 0.9044862 0.9344704
Family_CON.M.PA # 0.9192241 0.9459054




############################### Classification of habitat generalists and habitat specialists###################
############################################################################################################

####################Sorting Generalists， that is Family occupying more than or equal to 6 sites 
Family_Gener.S.FL <-subset(Family_OA.S.FL_exNA, Family_OA.S.FL_exNA[,1]>=6)
Family_Gener.M.FL <-subset(Family_OA.M.FL_exNA, Family_OA.M.FL_exNA[,1]>=6)
Family_Gener.S.PA <-subset(Family_OA.S.PA_exNA, Family_OA.S.PA_exNA[,1]>=6)
Family_Gener.M.PA <-subset(Family_OA.M.PA_exNA, Family_OA.M.PA_exNA[,1]>=6)

######################count the number of genus in that category 
Family_count.G.S.FL <-length(which(Family_OA.S.FL_exNA[,1]>=6)) # count the number of genus in that category 
Family_count.G.M.FL <-length(which(Family_OA.M.FL_exNA[,1]>=6)) # 
Family_count.G.S.PA <-length(which(Family_OA.S.PA_exNA[,1]>=6)) # 
Family_count.G.M.PA <-length(which(Family_OA.M.PA_exNA[,1]>=6))  

Family_count.G.S.FL # 96
Famiy_count.G.M.FL # 134
Family_count.G.S.PA # 102
Family_count.G.M.PA # 167



########################Sorting Specialists that is Family occupying less than or equal to 2 sites 
Family_Spec.S.FL <-subset(Family_OA.S.FL_exNA, Family_OA.S.FL_exNA[,1]<=2)
Family_Spec.M.FL <-subset(Family_OA.M.FL_exNA, Family_OA.M.FL_exNA[,1]<=2)
Family_Spec.S.PA <-subset(Family_OA.S.PA_exNA, Family_OA.S.PA_exNA[,1]<=2)
Family_Spec.M.PA <-subset(Family_OA.M.PA_exNA, Family_OA.M.PA_exNA[,1]<=2)

######################count the number of genus in that category genus_count.S.S.FL <-length(which(genus_OA.S.FL_exNA[,1]<=2))
Family_count.S.S.FL <-length(which(Family_OA.S.FL_exNA[,1]<=2))
Family_count.S.M.FL <-length(which(Family_OA.M.FL_exNA[,1]<=2))
Family_count.S.S.PA <-length(which(Family_OA.S.PA_exNA[,1]<=2))
Family_count.S.M.PA <-length(which(Family_OA.M.PA_exNA[,1]<=2))

Family_count.G.S.FL # 96
Family_count.G.M.FL # 134
Family_count.G.S.PA # 102
Family_count.G.M.PA # 167



##################################Proportion of Generalists and Specialists, using the number of family in family_OA.S.FL_exNA for example #####################
Family_Proportion.G.S.FL <-(Family_count.G.S.FL/354)*100     ####### calcualte for generalists
Family_Proportion.G.M.FL <-(Family_count.G.M.FL/354)*100
Family_Proportion.G.S.PA <-(Family_count.G.S.PA/354)*100
Family_Proportion.G.M.PA <-(Family_count.G.M.PA/354)*100

Family_Proportion.S.S.FL <-(Family_count.S.S.FL/354)*100    ####### calcualte for specialists 
Family_Proportion.S.M.FL <-(Family_count.S.M.FL/354)*100
Family_Proportion.S.S.PA <-(Family_count.S.S.PA/354)*100
Family_Proportion.S.M.PA <-(Family_count.S.M.PA/354)*100

Family_Proportion.S.S.FL # 20.0565
Family_Proportion.S.M.FL # 18.36158
Family_Proportion.S.S.PA # 30.79096
Family_Proportion.S.M.PA # 21.46893

Family_Proportion.G.S.FL # 27.11864
Family_Proportion.G.M.FL # 37.85311
Family_Proportion.G.S.PA # 28.81356
Family_Proportion.G.M.PA # 47.17514








################################################ Graphing for Abundance and Occupancy at the family level ################################################################################################
   
######################### Plot mean relative abundance againt Occupancy for the different subcommunities #######################

part(mfrow=c(2,2))  ############## setting the 4 panels of figures in one page


#Plot mean relative abundance againt Occupancy for the different communities 
library(scales)


## Supplementary Figure S3 C
#######################Surface -FL community########################################

myColor= c(ifelse(Family_OA.S.FL_exNA[,1]<3, "red", ifelse(Family_OA.S.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Family_OA.S.FL_exNA[,1]), Family_OA.S.FL_exNA[,2], log="y", title(main="Surface_FL at family level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "20.06%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "27.12%", cex=0.8)
text(7.5, 0.00012, "rho=0.907, P-value=2.2e-16 ", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.FL)
dev.off()



########################CIL -FL community#######################################

myColor= c(ifelse(Family_OA.M.FL_exNA[,1]<3, "red", ifelse(Family_OA.M.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Family_OA.M.FL_exNA[,1]), Family_OA.M.FL_exNA[,2], log="y", title(main="CIL_FL at family level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "18.36%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "37.85%", cex=0.8)
text(7.5, 0.00012, "rho=0.897, P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.FL)
dev.off()



########################Surface -PA community#######################################

myColor= c(ifelse(Family_OA.S.PA_exNA[,1]<3, "red", ifelse(Family_OA.S.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)
plot(jitter(Family_OA.S.PA_exNA[,1]), Family_OA.S.PA_exNA[,2], log="y", title(main="Surface_PA at family level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "30.79%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "28.81%", cex=0.8)
text(7.5, 0.00012, "rho=0.922, P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.PA)
dev.off()


########################CIL -PA community#######################################

myColor= c(ifelse(Family_OA.M.PA_exNA[,1]<3, "red", ifelse(Family_OA.M.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Family_OA.M.PA_exNA[,1]), Family_OA.M.PA_exNA[,2], log="y", title(main="CIL_PA at family level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "21.47%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "47.18%", cex=0.8)
text(7.5, 0.00012, "rho=0.935, P-value= 2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.PA)
dev.off()


########### End ##########