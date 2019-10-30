
###### This script is to calcualte community dissimilarity metric based on Horn-distance and UniFrac distance which is visualized in the non-metric multidimensional scaling

## load R libraries for this session
library(vegan)
library(base)
library(ape)
library(GUniFrac)


###################### NMDS plotS for Horn distance##############################

## load community data
community <-read.csv ("InputFiles/OTU_table_exDoubleton.csv", row.names=1, header=T)
nmds.otu <-metaMDS(community, distance="horn", trymax=200)
nmds.otu$stress # 0.1539898

## load symbol table functions as the color code and symbols for NMDS plot. 
symbols <-read.csv("InputFiles/symbols_NMDS.csv", row.names=1, header=T) 



################ UniFrac distance##################################


## load all data for UniFrac distance 
UniFrac <-read.csv("InputFiles/UniFrac.csv", header=T, row.names=2) 
UniFrac <-UniFrac[,-1] ####### this OTU table should contain the same sequence name for each OTU as in the Treefile 
UniFrac <-t(UniFrac) ########## OTU count table, row -n sample, column -1 OTU 
write.csv(UniFrac, file="InputFiles/UniFrac_new.csv")

## load tree file
#the tree was unrooted for UniFrac, and it was generated in QIIME using the usage of 'make_phylogeny.py'
tree <-read.tree("InputFiles/rep.Tree.exdouble.tre") 

## tranform to phylo class file
tree.p <-as.phylo(tree)  

## load a table containing all information for setting the color and symbols for the figures
symbols <-read.csv("InputFiles/symbols_NMDS.csv", row.names=1, header=T)


## Calculate the UniFrac distance######################
UniFrac_distance <-read.csv("InputFiles/UniFrac_new.csv", row.names=1, header=T)
UniFrac_distance <-GUniFrac(UniFrac, tree.p, alpha=c(0,0.5,1))$UniFrac_distance
d5 <-UniFrac_distance[,,"d_0.5"] ######## combined unweight and weighted UniFrac


## NMDS of UniFrac matrix##################
nmds_wUnifrac <-metaMDS(as.dist(d5))
nmds_wUnifrac$stress ###### 0.108631

## Switch axis for a easier comparison between analysis#######################
nmds_wUnifrac$points<-nmds_wUnifrac$points/(-1)



####################Plot Horn and Unifrac NMDS in one page######################

## Display plots in one row but two columns
par(mfrow=c(1,2)) 

## margins of plot
par(mar=c(5,5,5,3)) 

## Graphical features for legend 
Communities<-c("Surface_FL / Rep 1", "Surface_FL / Rep 2", "Surface_PA / Rep 1", "Surface_PA / Rep 2", "CIL_FL / Rep 1", "CIL_FL / Rep 2", "CIL_PA / Rep 1","CIL_PA / Rep 2")

Colors <-c("orange","orange","darkolivegreen1","darkolivegreen1","orange","orange","darkolivegreen1","darkolivegreen1")

## Supplementary Figure S2 NDMS (Horn-Morisita)
plot(nmds.otu$points, pch=symbols$Symbols, col=as.character(symbols$Colors), xlab="NMDS1", ylab="NMDS2", cex.axis=0.8, cex.lab=0.8, cex=1.2) #fontsize
title(main="NMDS (Horn-Morisita)", adj=0, cex.main=1.1) ##############add title to the plot
text(1.4,2, "stress-level=0.154", cex=0.8) ###################add stress level as text 

# margins of plot
par(mar=c(5,1,5,7))


##Supplementary Figure S2 NDMS (UniFrac)
plot(nmds_wUnifrac$points, pch=symbols$Symbols, col=as.character(symbols$Colors), xlab="NMDS1", ylab="NMDS2", cex.axis=0.8, cex.lab=0.8, cex=1.2)
title(main="NMDS (UniFrac)", adj=0, cex.main=1.1)
text(0.12,0.2, "stress-level=0.109",cex=0.8)
Legend("bottom", inset=c(-0.27, 0), y.intersp=2.4, c(Communities), col=c(Colors), pch=c(1,16,1,16,2,17,2,17), pt.cex=15, title="Subcommunities", cex=1, 
xpd=TRUE) ## Plotting outside of margins 


########### End ########