
#####################################Calculation of AO at the Order level###################################
####################################################################################################################################################

## load R libraries for this session 
library(scales)
library(base)
library(graphics)
library('RVAideMemoire')

############################# AO-relationsihp, habitat specialists, habitat generalists, based on Order level########################

######################## load order relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Order <-read.csv("InputFiles/Order_Table.csv", row.names=1, header=T) # Order_RelAbun
Order <-t(Order)

#################Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Order_RelAbun.M.PA <-Order_RelAbun.M[2:10,]
Order_RelAbun.M.FL <-Order_RelAbun.M[c(1, 11:18), ]
Order_RelAbun.S.PA <-Order_RelAbun.S[2:10,]
Order_RelAbun.S.FL <-Order_RelAbun.S[c(1, 11:18), ]

#################Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Order_RelAbun.M.PA <-t(Order_RelAbun.M.PA)
Order_RelAbun.M.FL <-t(Order_RelAbun.M.FL)
Order_RelAbun.S.PA <-t(Order_RelAbun.S.PA)
Order_RelAbun.S.FL <-t(Order_RelAbun.S.FL)

#################Calculate the regional mean abundance for every Order in each subset.
Order_MeanRelAbun.M.FL <-rowMeans(Order_RelAbun.M.FL)
Order_MeanRelAbun.S.FL <-rowMeans(Order_RelAbun.S.FL)
Order_MeanRelAbun.M.PA <-rowMeans(Order_RelAbun.M.PA)
Order_MeanRelAbun.S.PA <-rowMeans(Order_RelAbun.S.PA)

#################Count Rows that fit the criteria of Relative abundance =0 for each Order. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the Order.
Order_OCCU.M.PA <-rowSums(Order_RelAbun.M.PA==0)
Order_OCCU.M.FL <-rowSums(Order_RelAbun.M.FL==0)
Order_OCCU.S.PA <-rowSums(Order_RelAbun.S.PA==0)
Order_OCCU.S.FL <-rowSums(Order_RelAbun.S.FL==0)

#################Calcualte the number of sites that are occupied by the orders, by substracting the number of sites unoccupied from the total number of sites 9
Order_OCCU.M.PA2 <-9-Order_OCCU.M.PA
Order_OCCU.M.FL2 <-9-Order_OCCU.M.FL
Order_OCCU.S.PA2 <-9-Order_OCCU.S.PA
Order_OCCU.S.FL2 <-9-Order_OCCU.S.FL

#################Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance orders did not occur by chance for once.
Order_OA.M.FL <-cbind(Order_OCCU.M.FL2, Order_MeanRelAbun.M.FL)
Order_OA.S.FL <-cbind(Order_OCCU.S.FL2, Order_MeanRelAbun.S.FL)
Order_OA.M.PA <-cbind(Order_OCCU.M.PA2, Order_MeanRelAbun.M.PA)
Order_OA.S.PA <-cbind(Order_OCCU.S.PA2, Order_MeanRelAbun.S.PA)

colnames(Order_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Order_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Order_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Order_OA.S.PA) <-c("Occupancy","Relative abundance")

#################Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Order_OA.S.FL[Order_OA.S.FL <= 0.00002]<-NA
Order_OA.M.FL[Order_OA.M.FL <= 0.00002]<-NA
Order_OA.S.PA[Order_OA.S.PA <= 0.00002]<-NA
Order_OA.M.PA[Order_OA.M.PA <= 0.00002]<-NA


#################Remove the OTUs that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Order_OA.S.FL_exNA <-na.omit(Order_OA.S.FL, cols=seq_along(Order_OA.S.FL))
Order_OA.M.FL_exNA <-na.omit(Order_OA.M.FL, cols=seq_along(Order_OA.M.FL))
Order_OA.S.PA_exNA <-na.omit(Order_OA.S.PA, cols=seq_along(Order_OA.S.PA))
Order_OA.M.PA_exNA <-na.omit(Order_OA.M.PA, cols=seq_along(Order_OA.M.PA))

################Make Mean relative abundance logarithmic
Order_log.M.FL <-log(Order_OA.M.FL_exNA[,2])
Order_log.S.FL <-log(Order_OA.S.FL_exNA[,2])
Order_log.M.PA <-log(Order_OA.M.PA_exNA[,2])
Order_log.S.PA <-log(Order_OA.S.PA_exNA[,2])

Order_log.M.FL <-matrix(Order_log.M.FL)
Order_log.S.FL <-matrix(Order_log.S.FL)
Order_log.M.PA <-matrix(Order_log.M.PA)
Order_log.S.PA <-matrix(Order_log.S.PA)





###############################Calculation of Spearman's rank correlation at the Order level ########################################
####################################################################################################################################

Order_Spearman.S.FL <-cor.test(x=Order_OA.S.FL_exNA[,1], y=Order_log.S.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Order_Spearman.M.FL <-cor.test(x=Order_OA.M.FL_exNA[,1], y=Order_log.M.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Order_Spearman.S.PA <-cor.test(x=Order_OA.S.PA_exNA[,1], y=Order_log.S.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Order_Spearman.M.PA <-cor.test(x=Order_OA.M.PA_exNA[,1], y=Order_log.M.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)

######################## Vasulize the rho value and P value
Order_Spearman.S.FL  #  rho=0.8932855, P value=2.2e-16 
Order_Spearman.M.FL  #  rho=0.9015657, P value=2.2e-16
Order_Spearman.S.PA  # rho=0.9338354, P-value=2.2e-16
Order_Spearman.M.PA  # rho=0.9269942, P-value=2.2e-16





############################# Calculate confidence level of the interval， this is used to estimate the standard error of upper and lower limits of the 95% confidence interval. ########################################################
####################

library('RVAideMemoire')
Order_CON.S.FL <-spearman.ci(Order_OA.S.FL_exNA[,1], Order_log.S.FL, nrep=1000, conf.level=0.95)
Order_CON.M.FL <-spearman.ci(Order_OA.M.FL_exNA[,1], Order_log.M.FL, nrep=1000, conf.level=0.95)
Order_CON.S.PA <-spearman.ci(Order_OA.S.PA_exNA[,1], Order_log.S.PA, nrep=1000, conf.level=0.95)
Order_CON.M.PA <-spearman.ci(Order_OA.M.PA_exNA[,1], Order_log.M.PA, nrep=1000, conf.level=0.95)

######################################## output displaying the lower and upper limits of the 95% confidence interval, which was later used to calculate the standard errors for the coefficient plot ####################################
Order_CON.S.FL #0.8418080 0.9247332
Order_CON.M.FL #0.8522198 0.9292167
Order_CON.S.PA # 0.9110638 0.9473803
Order_CON.M.PA # 0.8959659 0.9445520





############################### Classification of habitat generalists and habitat specialists###################
############################################################################################################

####################Sorting Generalists， that is orders occupying more than or equal to 6 sites 
Order_Gener.S.FL <-subset(Order_OA.S.FL_exNA, Order_OA.S.FL_exNA[,1]>=6)  ###### Gener refers to generalists
Order_Gener.M.FL <-subset(Order_OA.M.FL_exNA, Order_OA.M.FL_exNA[,1]>=6)
Order_Gener.S.PA <-subset(Order_OA.S.PA_exNA, Order_OA.S.PA_exNA[,1]>=6)
Order_Gener.M.PA <-subset(Order_OA.M.PA_exNA, Order_OA.M.PA_exNA[,1]>=6)

######################count the number of Orders in that category 
Order_count.G.S.FL <-length(which(Order_OA.S.FL_exNA[,1]>=6))  
Order_count.G.M.FL <-length(which(Order_OA.M.FL_exNA[,1]>=6)) 
Order_count.G.S.PA <-length(which(Order_OA.S.PA_exNA[,1]>=6)) 
Order_count.G.M.PA <-length(which(Order_OA.M.PA_exNA[,1]>=6))  


Order_count.G.S.FL # 54 
Order_count.G.M.FL # 78
Order_count.G.S.PA # 62
Order_count.G.M.PA # 92



########################Sorting Specialists that is OTUs occupying less than or equal to 2 sites 
Order_Spec.S.FL <-subset(Order_OA.S.FL_exNA, Order_OA.S.FL_exNA[,1]<=2)
Order_Spec.M.FL <-subset(Order_OA.M.FL_exNA, Order_OA.M.FL_exNA[,1]<=2)
Order_Spec.S.PA <-subset(Order_OA.S.PA_exNA, Order_OA.S.PA_exNA[,1]<=2)
Order_Spec.M.PA <-subset(Order_OA.M.PA_exNA, Order_OA.M.PA_exNA[,1]<=2)

######################count the number of OTUs in that category 
Order_count.S.S.FL <-length(which(Order_OA.S.FL_exNA[,1]<=2))
Order_count.S.M.FL <-length(which(Order_OA.M.FL_exNA[,1]<=2))
Order_count.S.S.PA <-length(which(Order_OA.S.PA_exNA[,1]<=2))
Order_count.S.M.PA <-length(which(Order_OA.M.PA_exNA[,1]<=2))

Order_count.S.S.FL # 32
Order_count.S.M.FL # 29
Order_count.S.S.PA # 52
Order_count.S.M.PA # 34



##################################Proportion of Generalists and Specialists, using the number of OTU in OTU_OA.S.FL_exNA for example #####################
Order_Proportion.G.S.FL <-(Order_count.G.S.FL/Order_OA.S.FL_exNA)*100   ####### calcualte for generalists
Order_Proportion.G.M.FL <-(Order_count.G.M.FL/Order_OA.M.FL_exNA)*100
Order_Proportion.G.S.PA <-(Order_count.G.S.PA/Order_OA.S.PA_exNA)*100
Order_Proportion.G.M.PA <-(Order_count.G.M.PA/Order_OA.M.PA_exNA)*100

Order_Proportion.S.S.FL <-(Order_count.S.S.FL/Order_OA.S.FL_exNA)*100    ####### calcualte for specialists   
Order_Proportion.S.M.FL <-(Order_count.S.M.FL/Order_OA.M.FL_exNA)*100 
Order_Proportion.S.S.PA <-(Order_count.S.S.PA/Order_OA.S.PA_exNA)*100
Order_Proportion.S.M.PA <-(Order_count.S.M.PA/Order_OA.M.PA_exNA)*100


Order_Proportion.S.S.FL # 18.0791
Order_Proportion.S.M.FL # 16.38418
Order_Proportion.S.S.PA # 29.37853
Order_Proportion.S.M.PA # 19.20904

Order_Proportion.G.S.FL # 30.50847
Order_Proportion.G.M.FL # 44.0678
Order_Proportion.G.S.PA # 35.02825
Order_Proportion.G.M.PA # 51.9774






################################################ Graphing for Abundance and Occupancy at the order level ################################################################################################
   
######################### Plot mean relative abundance againt Occupancy for the different subcommunities #######################

part(mfrow=c(2,2))  ############## setting the 4 panels of figures in one page

library(scales)


## Supplementary Figure S3 D
#######################Surface -FL community########################################

myColor= c(ifelse(Order_OA.S.FL_exNA[,1]<3, "red", ifelse(Order_OA.S.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Order_OA.S.FL_exNA[,1]), Order_OA.S.FL_exNA[,2], log="y", title(main="Surface_FL at order level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "18.08%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "30.51%", cex=0.8)
text(7.5, 0.00011, "rho=0.893, P value=2.2e-16 ", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.FL)
dev.off()


########################CIL -FL community#######################################

myColor= c(ifelse(Order_OA.M.FL_exNA[,1]<3, "red", ifelse(Order_OA.M.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Order_OA.M.FL_exNA[,1]), Order_OA.M.FL_exNA[,2], log="y", title(main="CIL_FL at order level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "16.38%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "44.07%", cex=0.8)
text(7.5, 0.00012, "rho=0.902, P value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.FL)
dev.off()


########################Surface -PA community#######################################

myColor= c(ifelse(Order_OA.S.PA_exNA[,1]<3, "red", ifelse(Order_OA.S.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)
plot(jitter(Order_OA.S.PA_exNA[,1]), Order_OA.S.PA_exNA[,2], log="y", title(main="Surface_PA at order level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "29.38%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "35.03%", cex=0.8)
text(7.5, 0.0005, "rho=0.934, P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.PA)
dev.off()


########################CIL -PA community#######################################

myColor= c(ifelse(Order_OA.M.PA_exNA[,1]<3, "red", ifelse(Order_OA.M.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Order_OA.M.PA_exNA[,1]), Order_OA.M.PA_exNA[,2], log="y", title(main="CIL_PA at order level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "19.21%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "51.98%", cex=0.8)
text(7.5, 0.00035, "rho=0.927, P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.PA)
dev.off()


########### End ###########