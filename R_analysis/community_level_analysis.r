community level analysis

####PCA principle Component Analysis analysis to visualize the similarity among the sampling stations and their relation with the environmental variables
####Figure 1B###

EnvironPCA <-read.csv("./EnvironPCA.csv", row.names=1, header=T) # samples in raw, and environmental parameters in column

#look at the Environ, drop other parameters not used for calculation of dist of environmental variables

Environ1 <-EnvironPCA[,-1:-2]
Environ1 <-Environ1[,-12:-13]

#perform pca
pca<-prcomp(Environ1, center=TRUE, scale=TRUE)

#look at the output
summary(pca)


prop1<-summary(pca)$importance[2]*100
prop2<-summary(pca)$importance[5]*100

class(pca)

#look at the parameters return
names(pca)

pca.scores<-pca$x[,1:2]
pca.loadings<-pca$rotation[,1:2]

#The first thing is to scale our eigenvectors
 #to unit length...this is what is done in the packaged
#biplot call - this also has the same effect as normlaizing
 #the PC axes: normalize PC axes to unit length (McCune and Grace)
 
pca.scores[,1]<-pca.scores[,1]*(1/sqrt(sum(pca.scores[,1]^2)))
x=pca.scores[,1]

pca.scores[,2]<-pca.scores[,2]*(1/sqrt(sum(pca.scores[,2]^2)))
y=pca.scores[,2]

sc<-1
unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)),abs(max(x, na.rm = TRUE)))
x.scores = unsigned.range(pca.scores[,1])
y.scores = unsigned.range(pca.scores[,2])
x.loadings = unsigned.range(pca.loadings[,1])
y.loadings = unsigned.range(pca.loadings[,2])
xlim <- ylim <- x.scores <- x.loadings <- range(x.scores, x.loadings)
ratio <- max(y.scores/x.scores, y.loadings/x.loadings)/sc


par(mar=c(5, 6, 4, 4))
plot(pca.scores[,1],pca.scores[,2], main="PCA", ylab=paste("PC2 - ",round(prop2,0),"%"), xlab=paste("PC1 - ",round(prop1,0),"%"),cex=1, cex.axis=0.8, cex.lab=0.76, pch=EnvironPCA[,15])

library(calibrate)
par(new=TRUE)
p1<-plot(pca.loadings,axes=FALSE,type="n",xlim=xlim*ratio,ylim=ylim*ratio,xlab = "", ylab = "")
axis(3, col = "gray40", cex.axis=0.8)
axis(4, col = "gray40", cex.axis=0.8)
#we can also remove the axis
#text(pca.loadings * 0.9, labels = rownames(pca.loadings), cex =0.7, col = "gray40")
textxy(X=pca.scores[,1], Y=pca.scores[,2], labs=EnvironPCA[,1], cex=0.6) # as the labs on points are too small to be seen
arrows(0, 0, pca.loadings[, 1] * 0.6, pca.loadings[, 2] * 0.6, col = "gray40",length = 0.09)
abline(v=0,h=0,lty="dotted")
legend("topright", inset=c(0.006,0), y.intersp=1.20, c("Bb_ND","Bb_DT","Bm_ND","Bm_DT","Mb_ND","Mb_DT","Mm_ND","Mm_DT"), cex=0.65, col=c(rep("orange",2), rep("darkolivegreen1", 2), rep ("darkorange4",2),rep("darkgreen",2)), pch=rep (c(4,17)), pt.cex=1, xpd = TRUE)
pdf("netapp1/homes/shen/Profile/Desktop/PCA/PCAnew.pdf", width = 5.7, height=6.2)
plot()
dev.off()



#######ANOSIM analysis ########
#anosim analysis to test the grouping factor 
Otutable<-read.csv("./OTU_subampled.csv",row.names = 1, header = T)

####
##Test the effct of depth on community similarity
T1<-anosim(Otutable, Symbol$Layer, permutations = 999, distance = "bray", strat=NULL)
# grouping by depth
#Call:
#anosim(x = Otutable, grouping = Symbol$Layer, permutations = 999,      distance = "bray", strata = NULL) 
#Dissimilarity: bray 

ANOSIM statistic R: 0.5882 
      Significance: 0.001 

Permutation: free
Number of permutations: 999

####
###Test the effect of size fraction on community similarity
T2 <-anosim(Otutable, Symbol$Lifestyle, permutations = 999, distance = "bray", strat=NULL)
#Call:
#anosim(x = Otutable, grouping = Symbol$Lifestyle, permutations = 999,      distance = "bray", strata = NULL) 
#Dissimilarity: bray 

#ANOSIM statistic R: 0.672 
 #     Significance: 0.001 

#Permutation: free
#Number of permutations: 999






########Principle Coordinates Analysis CPCoA analysis
####PcoA analysis ####
#########################


######PcoA analyses
library("ape")
library("vegan")

#### PCA on the surface_Fl, surface_PA, CIL_FL, and CIL_PA

Otu_FL <-read.csv("./Otu_FL_subsampled.csv", header=T, row.names=1)
Otu_PA <-read.csv("./Otu_PA_subsampled.csv", header=T, row.names=1)

##### separate the FL communities by depth ####
Otu_FL_CIL <-Otu_FL[1:9,]
Otu_FL_Sur <-Otu_FL[10:18,]

##### separate the FL communities by depth 
Otu_PA_CIL <-Otu_PA[1:9,]
Otu_PA_Sur <-Otu_PA[10:18,]

#### making fitting environmental variables
Environ_Sur <-read.csv("./Environ_Sur.csv", row.names=1, header=T)
Environ_CIL <-read.csv("./Environ_CIL.csv", row.names=1, header=T)


###remove the first column including the depth as it is not there is lifestrategy bacteria from each depth###
Environ_SurNew <-Environ_Sur[,-1:-4]
Environ_SurNew <-Environ_SurNew[,-12:-13]


Environ_CILNew <-Environ_CIL[,-1:-4]
Environ_CILNew <-Environ_CILNew[,-12:-13]

###############
########surface_FL#########
###############
bc.d<-vegdist(Otu_FL_Sur, method="bray")
bc.pcoa=cmdscale(bc.d, eig=TRUE)

#Look at the output
bc.pcoa

###Second, calculate percent variance explained by each axes1 and 2, then add to plot later 
ax1.v.bc=bc.pcoa$eig[1]/sum(bc.pcoa$eig)
ax2.v.bc=bc.pcoa$eig[2]/sum(bc.pcoa$eig)
# look at the output
ax1.v.bc
ax2.v.bc

##Third, correlate measured environmental variables (env) to the axes
envEF.bc=envfit(bc.pcoa, Environ_SurNew, permu=999, na.rm=TRUE)
#look at the output
envEF.bc


################
########surface_PA#######
##################
bc.d<-vegdist(Otu_PA_Sur, method="bray")
bc.pcoa=cmdscale(bc.d, eig=TRUE)

#Look at the output
bc.pcoa

###Second, calculate percent variance explained by each axes1 and 2, then add to plot later 
ax1.v.bc=bc.pcoa$eig[1]/sum(bc.pcoa$eig)
ax2.v.bc=bc.pcoa$eig[2]/sum(bc.pcoa$eig)
# look at the output
ax1.v.bc
ax2.v.bc

##Third, correlate measured environmental variables (env) to the axes
envEF.bc=envfit(bc.pcoa, Environ_SurNew, permu=999, na.rm=TRUE)
#look at the output
envEF.bc


####################
########CIL_FL##########
####################
bc.d<-vegdist(Otu_FL_CIL, method="bray")
bc.pcoa=cmdscale(bc.d, eig=TRUE)

#Look at the output
bc.pcoa

###Second, calculate percent variance explained by each axes1 and 2, then add to plot later 
ax1.v.bc=bc.pcoa$eig[1]/sum(bc.pcoa$eig)
ax2.v.bc=bc.pcoa$eig[2]/sum(bc.pcoa$eig)
# look at the output
ax1.v.bc
ax2.v.bc

##Third, correlate measured environmental variables (env) to the axes
envEF.bc=envfit(bc.pcoa, Environ_CILNew, permu=999, na.rm=TRUE)
#look at the output
envEF.bc


##################
######## CIL_PA#########
###################
bc.d<-vegdist(Otu_PA_CIL, method="bray")
bc.pcoa=cmdscale(bc.d, eig=TRUE)

#Look at the output
bc.pcoa

###Second, calculate percent variance explained by each axes1 and 2, then add to plot later 
ax1.v.bc=bc.pcoa$eig[1]/sum(bc.pcoa$eig)
ax2.v.bc=bc.pcoa$eig[2]/sum(bc.pcoa$eig)
# look at the output
ax1.v.bc
ax2.v.bc

##Third, correlate measured environmental variables (env) to the axes
envEF.bc=envfit(bc.pcoa, Environ_CILNew, permu=999, na.rm=TRUE)
#look at the output
envEF.bc

