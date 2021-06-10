#################################################
#Alpha and beta diversity plotting and statistics 
#################################################


####Community alpha diversity

####Phylogenetic diversity 
library("picante")
Otutable<-read.csv("/Users/Dandan/Desktop/IOW/MicroFun/PA_FL_manuscript/OTU_subampled.csv",row.names = 1, header = T)
Environ <-read.csv("/Users/Dandan/Desktop/IOW/MicroFun/PA_FL_manuscript/PcoA/EnvironData.csv", row.names=1, header=T)
rownames(Otutable) <-rownames(Environ) ##sample ID in row, OTU in column

## load the community data matrix that used for UniFrac, because the tips ID on the tree must be the same names as the in the community data 
UniFrac<-read.csv("InputFiles/UniFrac_new.csv",row.names = 1, header = T)
colnames(Otutable) <-colnames(UniFrac)

#load tree file
tree <-read.tree("InputFiles/rep.Tree.exdouble.tre")
tree.p <-as.phylo(tree)

###Perform phylogenetic diversity analyses 
PD <-pd(Otutable, tree.p, include.root=TRUE)
write.csv(PD, file="./Community/PD.csv")  #to you directory



library(vegan)
library(ggplot2)

#values for species richness, evennes, and phylogenetic diversity are already combined into one data set for plotting
Comm_Diversity <-read.csv("InputFiles/CommunDiver.csv", header=T, row.names = 1)


setwd("/Users/Dandan/Desktop/IOW/MicroFun/PA_FL_manuscript/")
##drop off the shannon i
Comm_Diversity$Diversity <-factor(Comm_Diversity$Diversity, levels=c("Phylogenetic diversity", "Species richness","Evenness"))
Comm_Diversity$Depth <-factor(Comm_Diversity$Depth, levels =c("Surface","CIL"))

Comm_Diversity$Color <-as.character(Comm_Diversity$Color)
 
Comm_Diversity$Lifestyle<-factor(Comm_Diversity$Lifestyle, levels=c("FL", "PA"))


###plot Figure 2A
p1 <-ggplot(Comm_Diversity, aes(x=Depth,y=Values, fill=Lifestyle)) + 
geom_boxplot(fatten=0.65,alpha=0.8) + 
scale_fill_manual(values=c("FL"="orange","PA"="darkolivegreen1")) + 
theme_bw()+ 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + 
theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + 
labs(title="Alpha diversity") + theme(axis.title.y=element_text(colour = 'black',size = 11)) + 
theme(axis.title.x=element_text(colour = 'black',size = 11)) + 
theme(axis.text.y=element_text(colour = 'black', size = 10)) + 
theme(plot.title=element_text(colour= 'black', size=11, face='bold')) + 
theme(axis.text.x=element_text(colour = 'black', size = 10)) + 
facet_grid(Comm_Diversity$Diversity ~., scales="free") + 
theme(strip.text=element_text(size=11,face='bold')) + 
theme(strip.background=element_rect(size=1))+geom_point(position=jitter,alpha=0.45)

pdf("./AlphaDiversity.pdf", width =5, height=7) #save to your working directory
plot(p1)
dev.off()


################################
########Student's t-test #####
#################################

######Divide the dataset into different diversity matrix for statistics 
PD_Surface <-Comm_Diversity[1:18,]
PD_CIL <-Comm_Diversity[19:36,]
Richness_Surface <-Comm_Diversity[37:54,]
Richness_CIL<-Comm_Diversity[55:72,]
Evenness_Surface <-Comm_Diversity[73:90,]
Evenness_CIL<-Comm_Diversity[91:108,]

#########       
Richness_Surface$Lifestyle <-factor(Richness_Surface$Lifestyle)
t3 <-t.test(Richness_Surface$log ~Richness_Surface$Lifestyle, alternative="less")
t3
#	Welch Two Sample t-test

#data:  Richness_CIL$log by Richness_CIL$Type
#t = -2.3256, df = 15.901, p-value = 0.0168
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
    #    -Inf -0.01671068
#sample estimates:
#mean in group FL mean in group PA 
 #       3.004055         3.071170 
        
Richness_CIL$Lifestyle<-factor(Richness_CIL$Lifestyle)
t4 <-t.test(Richness_CIL$log ~Richness_CIL$Lifestyle, alternative="less")
t4
#	Welch Two Sample t-test

#data:  Richness_CIL$log by Richness_CIL$Lifestyle
#t = 1.3059, df = 10.985, p-value = 0.1091
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
# -0.0217428        Inf
#sample estimates:
#mean in group FL mean in group PA 
 #       2.928877         2.870950 

Evenness_Surface$Lifestyle <-factor(Evenness_Surface$Lifestyle)
t5<-t.test(Evenness_Surface$log ~Evenness_Surface$Lifestyle, alternative="less")
t5

#	Welch Two Sample t-test
#data:  Evenness_CIL$log by Evenness_CIL$Type
#t = -0.58599, df = 15.717, p-value = 0.2831
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
 #      -Inf 0.00996807
#sample estimates:
#mean in group FL mean in group PA 
 #     -0.1743562       -0.1693286 
         
      
Evenness_CIL$Lifestyle <-factor(Evenness_CIL$Lifestyle)
t6<-t.test(Evenness_CIL$log ~Evenness_CIL$Lifestyle, alternative="less")
t6
#Welch Two Sample t-test

#data:  Evenness_Surface$log by Evenness_Surface$Type
#t = -3.3173, df = 15.522, p-value = 0.002254
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
 #       -Inf -0.01585151
#sample estimates:
#mean in group FL mean in group PA 
 #     -0.2147442       -0.1812112 	
	

PD_Surface$Lifestyle <-factor(PD_Surface$Lifestyle)
t7<-t.test(PD_Surface$log ~PD_Surface$Lifestyle, alternative="less")
t7
#
 Welch Two Sample t-test

#data:  PD_Surface$log by PD_Surface$Lifestyle
#t = 6.369, df = 11.228, p-value = 2.415e-05
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
# 0.1392815       Inf
#sample estimates:
#mean in group FL mean in group PA 
 #       2.068567         1.874730 
        
PD_CIL$Lifestyle <-factor(PD_CIL$Lifestyle)
t8<-t.test(PD_CIL$log ~PD_CIL$Lifestyle, alternative="greater")
t8
##
#	Welch Two Sample t-test

#data:  PD_CIL$log by PD_CIL$Lifestyle
#t = 2.1644, df = 15.915, p-value = 0.02299
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
# 0.009018127         Inf
#sample estimates:
#mean in group FL mean in group PA 
 #       2.093459         2.046760          
      


############################
#### NMDS plot Betadiversity ##############
########################
Otutable<-read.csv("./OTU_subampled.csv",row.names = 1, header = T)
#Environ <-read.csv("./EnvironData.csv", row.names=1, header=T)
rownames(Otutable) <-rownames(Environ) ##sample ID in row, OTU in column
rownames(Otutable) <-c("DM1","DM10", "DM11","DM12", "DM13", "DM14","DM15","DM16","DM17", "DM18","DM2",	"DM3",	"DM4","DM5",	"DM6",	"DS7",	"DM8",	"DM9",	"DS1",	"DS10",	"DS11",	"DS12",	"DS13",	"DS14",	"DS15",	"DS16",	"DS17",	"DS18",	"DS2",	"DS3",	"DS4",	"DS5",	"DS6",	"DM7",	"DS8",	"DS9")

library("vegan")
library("gridExtra")
library("ggplot2")


##metadata for plotting 
Symbol <-read.csv("./Symbol.csv",header=T, row.names=1)
Symbol$Depth <-factor(Symbol$Depth, levels=c("Surface_FL", "Surface_PA", "CIL_FL","CIL_PA"))

#make grouping for anosim test later
Symbol$Layer <-factor(Symbol$Layer, levels=c("Surface", "CIL"))
Symbol$Lifestyle <-factor(Symbol$Lifestyle, levels=c("FL","PA"))

#Calculate the dissimilarity matrix 
nmds.otu<-metaMDS(Otutable, distance="bray", k=2, trymax=200, autotransform=T)
names(nmds.otu)
nmds.otu$stress
#0.1265208

#make NMDS table for plot
###Plot Figure 2B
NMDS = data.frame(NMDS1 =nmds.otu$points[,1], NMDS2 =nmds.otu$points[,2], Depth=Symbol$Layer, Lifestyle=Symbol$Lifestyle)

####look at the output
NMDS
p1 <-ggplot(data=NMDS, aes(NMDS1, NMDS2)) + geom_point(aes(fill=Lifestyle, shape=Depth), size=2.8) + scale_fill_manual(values=c("orange","darkolivegreen1")) + scale_shape_manual(values=c(21,24)) + theme_bw()+ theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())+ theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + labs(title="Bray-Curtis NMDS") + theme(axis.title.y=element_text(colour = 'black',size = 9)) + theme(axis.title.x=element_text(colour = 'black',size = 9)) + theme(axis.text.y=element_text(colour = 'black', size = 9)) + theme(plot.title=element_text(colour= 'black', size=10, face='bold')) + theme(axis.text.x=element_text(colour = 'black', size = 9)) + theme(legend.text=element_text(size=9)) + theme(legend.key=element_rect(fill="transparent", colour=NA)) + theme(legend.background=element_rect(fill=NA)) + theme(legend.title=element_text(color="black", size=9, face='bold')) + guides(color=guide_legend(override.aes=list(shape=15, size=3))) + guides(shape=guide_legend(override.aes=list(size=2.5))) + theme(legend.key.size=unit("0.32","cm"))  + theme(legend.title=element_blank())+ theme(plot.margin=margin(0.08,0.05,0.08,0.08,"cm"))
p1
#### combined Figure 2A and B in Illustrator



############END#############

