
########################Graphical display for Spearman rank's coefficient analyses###################

### load R library for this session
library(base)
library(graphics)
library(ggplot2)
library(gridExtra)

#########################################

coe <-read.csv("InputFiles/CorePlot.csv", header=T) 
############## including Variable that is taxonomic ranks, coefficient value, lifestyles, Depth, and Standard Errors (SE) of the uppler and lower of the 95% confidence interval #################### the values of upper and lower of the 95% confidence interval were calculated in the Spearman's rank analysis for each taxonomic rank ######## 

################ substracting surface relevant data ################
coe1 <-coe[1:12,]       

################ substracting cold intermediate layer relevant data ################
coe2 <-coe[13:24,]       

interval1 <-qnorm((1-0.9)/2)  # inner coefficient intervals 
interval2 <-qnorm((1-0.95)/2) # outer coefficient intervals


################## setting the factors to order the sequence of taxonomic ranks ###########################
coe1$Variable <-factor(coe1$Variable, levels=c("OTU","Genus", "Family","Order","Class","Phylum"))
coe2$Variable <-factor(coe1$Variable, levels=c("OTU","Genus", "Family","Order","Class","Phylum"))
coe1$Lifestyles <-factor(coe1$Lifestyles, levels=c("PA","FL"))
coe2$Lifestyles <-factor(coe2$Lifestyles, levels=c("PA","FL"))


## Figure 1 Surface panel
############################################ plot the surface samples #######################################
zp1 <-ggplot(coe1, aes(Variable, Coefficient,colour=Lifestyles)) + geom_hline(yintercept=0, colour=gray (1/2), lty=2)  + scale_color_manual(values=c("darkolivegreen1","orange"))+ geom_linerange(aes(x=Variable, ymin=Coefficient - SE*interval1, ymax=Coefficient+SE*interval1), lwd=1, position= position_dodge(width=1/2)) +geom_pointrange(aes(x=Variable, ymin=Coefficient-SE*interval2, ymax=Coefficient+SE*interval2), lwd=1/2, position= position_dodge(width=1/2), shape=21, fill="WHITE") +  theme_bw() +ggtitle("Surface") + coord_flip() 


## Figure 1 CIL panel
############################ plot the cold intermediate layer samles ####################################################
zp2 <-ggplot(coe2, aes(Variable, Coefficient,colour=Lifestyles)) + geom_hline(yintercept=0, colour=gray (1/2), lty=2) + scale_color_manual(values=c("darkolivegreen1","orange")) + geom_linerange(aes(x=Variable, ymin=Coefficient - SE*interval1, ymax=Coefficient+SE*interval1), lwd=1, position= position_dodge(width=1/2)) +geom_pointrange(aes(x=Variable, ymin=Coefficient -SE*interval2, ymax=Coefficient+SE*interval2), lwd=1/2, position= position_dodge(width=1/2), shape=21, fill="WHITE") + theme_bw() + ggtitle("CIL") + coord_flip() 
pdf("/Users/Dandan/Desktop/coefficentPlot.pdf", width=4.5, height=5.5)
grid.arrange(zp1, zp2, nrow=2)
dev.off()


########################### end #############################