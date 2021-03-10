
#####################################Calculation of AO at the Genus level###################################
####################################################################################################################################################

## load R libraries for this session 
library(scales)
library(base)
library(graphics)
library('RVAideMemoire')

############################# AO-relationsihp, habitat specialists, habitat generalists, based on Genus level########################

######################## load Genus relative abundance table, pleaes see the first few lines of scripts in the 'Abundance-occupancy at OTU level'   
Genus <-read.csv("InputFiles/Genus_Table.csv", row.names=1, header=T) # Genus_RelAbun
Genus <-t(Genus)


#################Divide dataset into Cold Intermidate layer (however, for simplicity, we used 'M' instead of 'CIL' in all scripts)  and Surface sub-datasets
Genus_RelAbun.M <-Genus[1:18,]
Genus_RelAbun.S <-Genus[19:36,]

#################Divide dataset into PA (Particle-attached) and FL (Free-living) sub-datasets
Genus_RelAbun.M.PA <-Genus_RelAbun.M[2:10,]
Genus_RelAbun.M.FL <-Genus_RelAbun.M[c(1, 11:18), ]
Genus_RelAbun.S.PA <-Genus_RelAbun.S[2:10,]
Genus_RelAbun.S.FL <-Genus_RelAbun.S[c(1, 11:18), ]

#################Transpose Ordination. Depth: M refers to cold intermidate layer samples; S to surface samples. Size fractionated assemblages: PA refers to particle-attached communities; FL to free-living communities.
Genus_RelAbun.M.PA <-t(Genus_RelAbun.M.PA)
Genus_RelAbun.M.FL <-t(Genus_RelAbun.M.FL)
Genus_RelAbun.S.PA <-t(Genus_RelAbun.S.PA)
Genus_RelAbun.S.FL <-t(Genus_RelAbun.S.FL)

#################Calculate the regional mean abundance for every genus in each subset. This step assure 
Genus_MeanRelAbun.M.FL <-rowMeans(Genus_RelAbun.M.FL)
Genus_MeanRelAbun.S.FL <-rowMeans(Genus_RelAbun.S.FL)
Genus_MeanRelAbun.M.PA <-rowMeans(Genus_RelAbun.M.PA)
Genus_MeanRelAbun.S.PA <-rowMeans(Genus_RelAbun.S.PA)

#################Count Rows that fit the criteria of Relative abundance =0 for each genus. Here, it means that if relative abundance is 0, considered as absent; therefore, The output returns the a set of number indicating the number of site unoccupied by the Genus.
Genus_OCCU.M.PA <-rowSums(Genus_RelAbun.M.PA==0)
Genus_OCCU.M.FL <-rowSums(Genus_RelAbun.M.FL==0)
Genus_OCCU.S.PA <-rowSums(Genus_RelAbun.S.PA==0)
Genus_OCCU.S.FL <-rowSums(Genus_RelAbun.S.FL==0)

#################Calcualte the number of sites that are occupied by the genus, by substracting the number of sites unoccupied from the total number of sites 9
Genus_OCCU.M.PA2 <-9-Genus_OCCU.M.PA
Genus_OCCU.M.FL2 <-9-Genus_OCCU.M.FL
Genus_OCCU.S.PA2 <-9-Genus_OCCU.S.PA
Genus_OCCU.S.FL2 <-9-Genus_OCCU.S.FL

#################Create Combined Matrix with Occupancy and Mean relative abundance. This step ensures that the low-abundance genus did not occur by chance for once
Genus_OA.M.FL <-cbind(Genus_OCCU.M.FL2, Genus_MeanRelAbun.M.FL)
Genus_OA.S.FL <-cbind(Genus_OCCU.S.FL2, Genus_MeanRelAbun.S.FL)
Genus_OA.M.PA <-cbind(Genus_OCCU.M.PA2, Genus_MeanRelAbun.M.PA)
Genus_OA.S.PA <-cbind(Genus_OCCU.S.PA2, Genus_MeanRelAbun.S.PA)

colnames(Genus_OA.M.FL) <-c("Occupancy","Relative abundance")
colnames(Genus_OA.S.FL) <-c("Occupancy","Relative abundance")
colnames(Genus_OA.M.PA) <-c("Occupancy","Relative abundance")
colnames(Genus_OA.S.PA) <-c("Occupancy","Relative abundance")

#################Set all values <0.00002 to NA, because the regional mean abundance <0.002% might have been due to the sequence errors according to Pandit et al., 2009
Genus_OA.S.FL[Genus_OA.S.FL <= 0.00002]<-NA
Genus_OA.M.FL[Genus_OA.M.FL <= 0.00002]<-NA
Genus_OA.S.PA[Genus_OA.S.PA <= 0.00002]<-NA
Genus_OA.M.PA[Genus_OA.M.PA <= 0.00002]<-NA


#################Remove the genus that have mean relative abundance lower than 0.00002, this step further reduced the noise for the downstream analyses
Genus_OA.S.FL_exNA <-na.omit(Genus_OA.S.FL, cols=seq_along(Genus_OA.S.FL))
Genus_OA.M.FL_exNA <-na.omit(Genus_OA.M.FL, cols=seq_along(Genus_OA.M.FL))
Genus_OA.S.PA_exNA <-na.omit(Genus_OA.S.PA, cols=seq_along(Genus_OA.S.PA))
Genus_OA.M.PA_exNA <-na.omit(Genus_OA.M.PA, cols=seq_along(Genus_OA.M.PA))


################Make Mean relative abundance logarithmic
Genus_log.M.FL <-log(Genus_OA.M.FL_exNA[,2])
Genus_log.S.FL <-log(Genus_OA.S.FL_exNA[,2])
Genus_log.M.PA <-log(Genus_OA.M.PA_exNA[,2])
Genus_log.S.PA <-log(Genus_OA.S.PA_exNA[,2])

Genus_log.M.FL <-matrix(Genus_log.M.FL)
Genus_log.S.FL <-matrix(Genus_log.S.FL)
Genus_log.M.PA <-matrix(Genus_log.M.PA)
Genus_log.S.PA <-matrix(Genus_log.S.PA)





###############################Calculation of Spearman's rank correlation at the Genus level ########################################
####################################################################################################################################

Genus_Spearman.S.FL <-cor.test(x=Genus_OA.S.FL_exNA[,1], y=Genus_log.S.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Genus_Spearman.M.FL <-cor.test(x=Genus_OA.M.FL_exNA[,1], y=Genus_log.M.FL, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Genus_Spearman.S.PA <-cor.test(x=Genus_OA.S.PA_exNA[,1], y=Genus_log.S.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
Genus_Spearman.M.PA <-cor.test(x=Genus_OA.M.PA_exNA[,1], y=Genus_log.M.PA, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)

######################## Vasulize the rho value and P value
Genus_Spearman.S.FL #  rho=0.9021975, P value=2.2e-16 
Genus_Spearman.M.FL #  rho=0.9024015, P value=2.2e-16 
Genus_Spearman.S.PA # rho=0.8862197,  P-value=2.2e-16
Genus_Spearman.M.PA # rho=0.9377624, P-value=2.2e-16





############################# Calculate confidence level of the interval， this is used to estimate the standard error of upper and lower limits of the 95% confidence interval. ########################################################
####################

library('RVAideMemoire')

Genus_CON.S.FL <-spearman.ci(Genus_OA.S.FL_exNA[,1], Genus_log.S.FL, nrep=1000, conf.level=0.95)
Genus_CON.M.FL <-spearman.ci(Genus_OA.M.FL_exNA[,1], Genus_log.M.FL, nrep=1000, conf.level=0.95)
Genus_CON.S.PA <-spearman.ci(Genus_OA.S.PA_exNA[,1], Genus_log.S.PA, nrep=1000, conf.level=0.95)
Genus_CON.M.PA <-spearman.ci(Genus_OA.M.PA_exNA[,1], Genus_log.M.PA, nrep=1000, conf.level=0.95)


######################################## output displaying the lower and upper limits of the 95% confidence interval, which was later used to calculate the standard errors for the coefficient plot ####################################

Genus_CON.S.FL #95 percent confidence interval: 0.8795701 0.9197097
Genus_CON.M.FL #95 percent confidence interval:0.8828237 0.9166612
Genus_CON.S.PA #0.8647913 0.9048455
Genus_CON.M.PA # 0.9288185 0.9446666





############################### Classification of habitat generalists and habitat specialists###################
############################################################################################################

####################Sorting Generalists， that is Genus occupying more than or equal to 6 sites 
Genus_Gener.S.FL <-subset(Genus_OA.S.FL_exNA, Genus_OA.S.FL_exNA[,1]>=6)
Genus_Gener.M.FL <-subset(Genus_OA.M.FL_exNA, Genus_OA.M.FL_exNA[,1]>=6)
Genus_Gener.S.PA <-subset(Genus_OA.S.PA_exNA, Genus_OA.S.PA_exNA[,1]>=6)
Genus_Gener.M.PA <-subset(Genus_OA.M.PA_exNA, Genus_OA.M.PA_exNA[,1]>=6)

######################count the number of genus in that category 
Genus_count.G.S.FL <-length(which(Genus_OA.S.FL_exNA[,1]>=6)) # count the number of genus in that category 
Genus_count.G.M.FL <-length(which(Genus_OA.M.FL_exNA[,1]>=6)) # 
Genus_count.G.S.PA <-length(which(Genus_OA.S.PA_exNA[,1]>=6)) # 
Genus_count.G.M.PA <-length(which(Genus_OA.M.PA_exNA[,1]>=6))  

Genus_count.G.S.FL #  139
Genus_count.G.M.FL # 198
Genus_count.G.S.PA # 153
Genus_count.G.M.PA # 272




########################Sorting Specialists that is genus occupying less than or equal to 2 sites 
Genus_Spec.S.FL <-subset(Genus_OA.S.FL_exNA, Genus_OA.S.FL_exNA[,1]<=2)
Genus_Spec.M.FL <-subset(Genus_OA.M.FL_exNA, Genus_OA.M.FL_exNA[,1]<=2)
Genus_Spec.S.PA <-subset(Genus_OA.S.PA_exNA, Genus_OA.S.PA_exNA[,1]<=2)
Genus_Spec.M.PA <-subset(Genus_OA.M.PA_exNA, Genus_OA.M.PA_exNA[,1]<=2)

######################count the number of genus in that category genus_count.S.S.FL <-length(which(genus_OA.S.FL_exNA[,1]<=2))
Genus_count.S.S.FL <-length(which(Genus_OA.S.FL_exNA[,1]<=2))
Genus_count.S.M.FL <-length(which(Genus_OA.M.FL_exNA[,1]<=2))
Genus_count.S.S.PA <-length(which(Genus_OA.S.PA_exNA[,1]<=2))
Genus_count.S.M.PA <-length(which(Genus_OA.M.PA_exNA[,1]<=2))

Genus_count.S.S.FL # 156
Genus_count.S.M.FL # 165
Genus_count.S.S.PA # 272
Genus_count.S.M.PA # 196



##################################Proportion of Generalists and Specialists, using the number of genus in genus_OA.S.FL_exNA for example #####################
Genus_Proportion.G.S.FL <-(Genus_count.G.S.FL/733)*100     ####### calcualte for generalists
Genus_Proportion.G.M.FL <-(Genus_count.G.M.FL/733)*100
Genus_Proportion.G.S.PA <-(Genus_count.G.S.PA/733)*100
Genus_Proportion.G.M.PA <-(Genus_count.G.M.PA/733)*100

Genus_Proportion.S.S.FL <-(Genus_count.S.S.FL/733)*100    ####### calcualte for specialists 
Genus_Proportion.S.M.FL <-(Genus_count.S.M.FL/733)*100
Genus_Proportion.S.S.PA <-(Genus_count.S.S.PA/733)*100
Genus_Proportion.S.M.PA <-(Genus_count.S.M.PA/733)*100

Genus_Proportion.S.S.FL # 21.2824
Genus_Proportion.S.M.FL # 22.51023
Genus_Proportion.S.S.PA # 37.10778
Genus_Proportion.S.M.PA # 26.73943

Genus_Proportion.G.S.FL # 18.96317
Genus_Proportion.G.M.FL # 27.01228
Genus_Proportion.G.S.PA # 20.87312
Genus_Proportion.G.M.PA # 37.10778





################################################ Graphing for Abundance and Occupancy at the genus level ################################################################################################
   
######################### Plot mean relative abundance againt Occupancy for the different subcommunities #######################

part(mfrow=c(2,2))  ############## setting the 4 panels of figures in one page


#Plot mean relative abundance againt Occupancy for the different communities 
library(scales)

## Supplementary Figure S3 B
#######################Surface -FL community########################################

myColor= c(ifelse(Genus_OA.S.FL_exNA[,1]<3, "red", ifelse(Genus_OA.S.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Genus_OA.S.FL_exNA[,1]), Genus_OA.S.FL_exNA[,2], log="y", title(main="Surface_FL at genus level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "21.28%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "18.96%", cex=0.8)
text(7.5, 0.00011, "rho=0.902 , P value=2.2e-16 ", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.FL)
dev.off()


########################CIL -FL community#######################################

myColor= c(ifelse(Genus_OA.M.FL_exNA[,1]<3, "red", ifelse(Genus_OA.M.FL_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Genus_OA.M.FL_exNA[,1]), Genus_OA.M.FL_exNA[,2], log="y", title(main="CIL_FL at genus level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "22.51%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "27.01%", cex=0.8)
text(7.5, 0.00012, "rho=0.902, P value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.FL)
dev.off()


########################Surface -PA community#######################################

myColor= c(ifelse(Genus_OA.S.PA_exNA[,1]<3, "red", ifelse(Genus_OA.S.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)
plot(jitter(Genus_OA.S.PA_exNA[,1]), Genus_OA.S.PA_exNA[,2], log="y", title(main="Surface_PA at genus level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "37.11%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "20.87%", cex=0.8)
text(7.5, 0.00012, "rho=0.886, P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.S.PA)
dev.off()


########################CIL -PA community#######################################

myColor= c(ifelse(Genus_OA.M.PA_exNA[,1]<3, "red", ifelse(Genus_OA.M.PA_exNA[,1]>5, "blue", "grey")))
myColorAlpha <-alpha(myColor, 0.4)

plot(jitter(Genus_OA.M.PA_exNA[,1]), Genus_OA.M.PA_exNA[,2], log="y", title(main="CIL_PA at genus level"), xlab="Occupancy", ylab="Mean relative abundance", col= myColorAlpha, cex.lab=0.8, cex.axis=0.8, pch=16)
text(1.5, 6,"Specialists", cex=0.8)
text(1.5,2.5, "26.74%", cex=0.8)
text(7.5, 6,"Generalists", cex=0.8)
text(7.5,2.5, "37.11%", cex=0.8)
text(7.5, 0.00012, "rho=0.938, P-value=2.2e-16", cex=0.8)
abline(v=2.5, col="red", lty=3)
abline(v=5.5, col="blue", lty=3)
#abline(Linear_Reg.M.PA)
dev.off()


########### End ##########
