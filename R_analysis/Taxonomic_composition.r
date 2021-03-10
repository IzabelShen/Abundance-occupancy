##Taxonomic composition plot#####
###Overall taxonomic composition across stations 
####Specialist and generalist taxonomic composition


####
###Make plot of taxonomic affiliation of All phyla across stations and detph and lifestrategy >1%

Phyl <-read.csv("./Phylum.csv", header=T)
Phyl$Phylum <-factor(Phyl$Phylum, levels=c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria","Actinobacteria","Bacteroidetes",   "Chlamydiae",   "Chloroflexi",  "Cyanobacteria","Deferribacteres",  "Lentisphaerae", "Plantomycetes", "Verrucomicrobia", "Unclassified_Bacteria"))
Phyl$Type <-factor(Phyl$Type, levels=c("Surface_FL","Surface_PA", "CIL_FL", "CIL_PA"))
Phyl$Station <-factor(Phyl$Station)
 
 
###stacked bar plot
#####Figure 3 in the manuscript
#####
p6<-ggplot(Phyl, aes(factor(Station), Abundance, fill=Phylum, order=as.numeric(Phylum))) + geom_bar(stat="identity", width=0.64) +
scale_fill_manual(values=c("Alphaproteobacteria"="darkcyan", "Betaproteobacteria"="tan4", "Gammaproteobacteria"="dodgerblue3", "Deltaproteobacteria"="darkmagenta","Actinobacteria"="gold2","Bacteroidetes"="darkolivegreen2","Chlamydiae"="cyan3",   "Chloroflexi"="brown2",   "Cyanobacteria"="pink","Deferribacteres"="khaki2",  "Lentisphaerae"="lightskyblue2", "Plantomycetes"="darkslategray","Verrucomicrobia"="chocolate4",  "Unclassified_Bacteria"="gray50",  )) + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + theme(panel.background =element_rect(fil="white", color="black")) + theme(axis.text.x =element_text(colour = 'black', size = 10.5, family='serif')) + theme(axis.text.y=element_text(colour = 'black', size = 10.5,family='serif')) + theme(axis.title.x=element_text(colour = 'black',size = 12, family='serif')) + theme(axis.ticks.x = element_blank()) + theme(axis.title.y=element_text(colour = 'black',size = 12, family='serif')) +
theme(legend.title=element_text(size=9, face='bold', family='serif'))  + theme(legend.text=element_text(size=9, family='serif')) + theme(legend.key.size=unit("0.38","cm"))+ facet_grid(. ~ Type, scales="free", space="free_x") +  labs(title="", y="Relative abundance %", x="Stations")
 
pdf("./PhylumRel_1.pdf", width =8.5, height=5)
plot(p6)
dev.off()



####################################################################
####Statistic analyses for phylum level abundance###
###test the abundance of phyla across PA and FL communities at CIL
####################################################################

Phyla <-read.csv("./Taxo_Phylum/PhylaStat.csv", header=1)

####Test indiviudal taxonomic groups
##The results are reported in SUpplementary Table S2######
###########

######
Acidobacteria_CIL <-subset(Phyla, Phyla$Phylum=='Acidobacteria' & Phyla$Depth=='CIL')
Acidobacteria_Sur <-subset(Phyla, Phyla$Phylum=='Acidobacteria' & Phyla$Depth=='Surface')
#droplevels(Acidobacteria_CIL$Type)
#droplevels(Acidobacteria_Sur$Type)

T1<-wilcox.test(Acidobacteria_CIL$Relative ~ Acidobacteria_CIL$Type, alternative = c("less"), paired=F)
T2<-wilcox.test(Acidobacteria_Sur$Relative ~ Acidobacteria_Sur$Type, paired=F, alternative = c("less"))
T1
T2

#####
Actinobacteria_CIL <-subset(Phyla, Phyla$Phylum=='Actinobacteria' & Phyla$Depth=='CIL')
Actinobacteria_Sur <-subset(Phyla, Phyla$Phylum=='Actinobacteria' & Phyla$Depth=='Surface')

T3 <-wilcox.test(Actinobacteria_CIL$Relative ~ Actinobacteria_CIL$Type, paired=F,alternative = c("greater"))
T4 <-wilcox.test(Actinobacteria_Sur$Relative ~ Actinobacteria_Sur$Type, paired=F,alternative = c("greater") )
T3
T4

####BD1-5
p1 <-subset(Phyla, Phyla$Phylum=='BD1-5' & Phyla$Depth=='CIL')
p2 <-subset(Phyla, Phyla$Phylum=='BD1-5' & Phyla$Depth=='Surface')

T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
T1
T2


##Unclassified_Bacteria
p1 <-subset(Phyla, Phyla$Phylum=='Unclassified_Bacteria' & Phyla$Depth=='CIL')
p2 <-subset(Phyla, Phyla$Phylum=='Unclassified_Bacteria' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less") )
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F,alternative = c("less") )

##look at results
T1
T2

#Bacteroidetes
p1 <-subset(Phyla, Phyla$Phylum=='Bacteroidetes' & Phyla$Depth=='CIL')
p2 <-subset(Phyla, Phyla$Phylum=='Bacteroidetes' & Phyla$Depth=='Surface')

T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F, alternative = c("less") )
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less") )

##look at results
T1
T2

#Candidate_BRC1
p1 <-subset(Phyla, Phyla$Phylum=='Candidate_division_BRC1' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Candidate_division_BRC1' & Phyla$Depth=='Surface')

T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F, alternative = c("less") )
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less") )

##look at results
T1
T2


###Candidate_division_OD1
p1 <-subset(Phyla, Phyla$Phylum=='Candidate_division_OD1' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Candidate_division_OD1' & Phyla$Depth=='Surface')

#Candidate_OD1
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F, alternative = c("less") )
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F,  alternative = c("less") )

##look at results
T1
T2

#### TM7
p1 <-subset(Phyla, Phyla$Phylum=='Candidate_division_TM7' & Phyla$Depth=='CIL')
p2 <-subset(Phyla, Phyla$Phylum=='Candidate_division_TM7' & Phyla$Depth=='Surface')

T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F, alternative = c("less") )
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less") )

#Candidate_TM7
##look at results
T1
T2

##Chlamydiae
p1 <-subset(Phyla, Phyla$Phylum=='Chlamydiae' & Phyla$Depth=='CIL')
p2 <-subset(Phyla, Phyla$Phylum=='Chlamydiae' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F, alternative = c("less") )
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less") )

##look at results
#Chlamydiae
T1
T2

##Chloroflexi
p1 <-subset(Phyla, Phyla$Phylum=='Chloroflexi' & Phyla$Depth=='CIL')
p2 <-subset(Phyla, Phyla$Phylum=='Chloroflexi' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("greater"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
#Chloroflexi
T1
T2

###Cyanobacteria
p1 <-subset(Phyla, Phyla$Phylum=='Cyanobacteria' & Phyla$Depth=='CIL')
p2 <-subset(Phyla, Phyla$Phylum=='Cyanobacteria' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
#Cyanobacteria
T1
T2

##Deferribacteres
p1<-subset(Phyla, Phyla$Phylum=='Deferribacteres' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Deferribacteres' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("greater"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F,alternative = c("greater"))

##look at results
##Deferribacteres
T1
T2

##
p1<-subset(Phyla, Phyla$Phylum=='Firmicutes' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Firmicutes' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
#Firmicutes
T1
T2

###
p1<-subset(Phyla, Phyla$Phylum=='Fusobacteria' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Fusobacteria' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
#Fusobacteria
T1
T2

##
p1<-subset(Phyla, Phyla$Phylum=='Gemmatimonadetes' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Gemmatimonadetes' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
#Gemmatimonadetes
T1
T2

##
p1<-subset(Phyla, Phyla$Phylum=='Lentisphaerae' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Lentisphaerae' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F,alternative = c("less"))

##look at results
#Lentisphaerae
T1
T2

##
p1<-subset(Phyla, Phyla$Phylum=='Planctomycetes' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Planctomycetes' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=T,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=T, alternative = c("less"))

##look at results
#Planctomycetes
T1
T2

###
p1<-subset(Phyla, Phyla$Phylum=='Alphaproteobacteria' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Alphaproteobacteria' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("greater"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("greater"))

##look at results
#Alphaproteobacteria
T1
T2


###
p1<-subset(Phyla, Phyla$Phylum=='Betaproteobacteria' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Betaproteobacteria' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("greater"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("greater"))

##look at results
#Betaproteobacteria
T1
T2



###
p1<-subset(Phyla, Phyla$Phylum=='Deltaproteobacteria' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Deltaproteobacteria' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
#Deltaproteobacteria
T1
T2


###
p1<-subset(Phyla, Phyla$Phylum=='Epsilonproteobacteria' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Epsilonproteobacteria' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("greater"))

##look at results
#Epsilonproteobacteria 
T1
T2

###
p1<-subset(Phyla, Phyla$Phylum=='Gammaproteobacteria' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Gammaproteobacteria' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
#Gammaproteobacteria 
T1
T2



###
p1<-subset(Phyla, Phyla$Phylum=='TM6' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='TM6' & Phyla$Depth=='Surface')
p1$Type<-factor(p1$Type)
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F,alternative = c("less"))

##look at results
#TM6
T1
T2


###
p1<-subset(Phyla, Phyla$Phylum=='Tenericutes' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Tenericutes' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
##Tenericutes
T1
T2

##
p1<-subset(Phyla, Phyla$Phylum=='Verrucomicrobia' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Verrucomicrobia' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F,alternative = c("less"))

##look at results
Verrucomicrobia
T1
T2

##
p1<-subset(Phyla, Phyla$Phylum=='Other' & Phyla$Depth=='CIL')
p2<-subset(Phyla, Phyla$Phylum=='Other' & Phyla$Depth=='Surface')
T1 <-wilcox.test(p1$Relative ~ p1$Type, paired=F,alternative = c("less"))
T2 <-wilcox.test(p2$Relative ~ p2$Type, paired=F, alternative = c("less"))

##look at results
#Other
T1
T2






##################################################################
####Taxonomic composition for specialists and generalists
########Make plot of taxonomic affiliation of generalists (top 20 classes)
##################################################################


########Make plot of taxonomic affiliation of specialists (top 20 classes)
RelSpe <-read.csv("./RelSpeplot.csv", header=T)
RelSpe$Class <-factor(RelSpe$Class, levels=c("SAR202_clade", "SL56_marine_group","JG30-KF-CM66","Ardenticatenia", "Anaerolineae","SPG12-401-411-B72", "SGST604", "Planctomycetes_unclassified", "BD7-11",
"Spartobacteria","SS1-B-03-39","Verrucomicrobia_Incertae_Sedis","Arctic97B.4_marine_group","DEV055","Verrucomicrobia_unclassified","Acidobacteria","Lentisphaeria","Oligosphaeria","Armatimonadetes_unclassified","Candidate_division_BRC1_unclassified"))
RelSpe$Type <-factor(RelSpe$Type, levels=c("Surface_FL","Surface_PA", "CIL_FL", "CIL_PA"))

###stacked bar plot
p5<-ggplot(RelSpe, aes(factor(Type), Abundance, fill=Class, order=as.numeric(Class))) + 
geom_bar(stat="identity", width=0.64) +
scale_fill_manual(values=c("SAR202_clade"="darkorange", "SL56_marine_group"="darkorange1","JG30-KF-CM66"="darkorange2","Ardenticatenia"="darkorange3", "Anaerolineae"="darkorange4","SPG12-401-411-B72"="lightskyblue", "SGST604"="lightskyblue1", "Planctomycetes_unclassified"="lightskyblue3", "BD7-11"="lightskyblue4",
"Spartobacteria"="darkolivegreen","SS1-B-03-39"="darkolivegreen1","Verrucomicrobia_Incertae_Sedis"="darkolivegreen2","Arctic97B.4_marine_group"="darkolivegreen3","DEV055"="darkolivegreen4","Verrucomicrobia_unclassified"="darkgreen","Acidobacteria"="darkmagenta","Lentisphaeria"="gold","Oligosphaeria"="gold4","Armatimonadetes_unclassified"="firebrick4","Candidate_division_BRC1_unclassified"="darksalmon")) + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + theme(panel.background =element_rect(fil="white", color="black")) + theme(axis.text.x =element_text(colour = 'black', size = 10.5, angle=45, hjust = 1, family='serif')) + theme(axis.text.y=element_text(colour = 'black', size = 10.5,family='serif')) + theme(axis.title.x=element_text(colour = 'black',size = 12, family='serif')) + theme(axis.ticks.x = element_blank()) + theme(axis.title.y=element_text(colour = 'black',size = 12, family='serif')) + 
theme(legend.title=element_text(size=10, face='bold', family='serif'))  + theme(legend.text=element_text(size=10, family='serif')) + theme(legend.key.size=unit("0.38","cm"))  + labs(title="", y="Relative proportion %", x="Subcommunities") 

pdf("./Class_spe.pdf", width =6, height=6.5)
plot(p5)
dev.off()


#####Generalists #####
RelGene <-read.csv("./RelGeneplot.csv", header=T)
RelGene$Class <-factor(RelGene$Class, levels=c("Acidimicrobiia",	"Actinobacteria","Alphaproteobacteria",	"Betaproteobacteria",	"Gammaproteobacteria",	"Deltaproteobacteria",	"Cytophagia",	"Flavobacteriia",	"Sphingobacteriia",	"Bacteroidetes_unclassified",	"Opitutae",	"Verrucomicrobiae",	"OM190","Phycisphaerae",	"Planctomycetacia",	"Cyanobacteria","Deferribacteres",	"Chlamydiae",	"Bacteria_unclassified"))
RelGene$Type <-factor(RelGene$Type, levels=c("Surface_FL","Surface_PA", "CIL_FL", "CIL_PA"))

###stacked bar plot
p7<-ggplot(RelGene, aes(factor(Type), Abundance, fill=Class, order=as.numeric(Class))) + 
geom_bar(stat="identity", width=0.64) + 
scale_fill_manual(values=c("Acidimicrobiia"="lightpink",	"Actinobacteria"="mediumorchid1","Alphaproteobacteria"="mediumorchid4",	"Betaproteobacteria"="mediumpurple1",	"Gammaproteobacteria"="mediumpurple4",	"Deltaproteobacteria"="midnightblue",	"Cytophagia"="tan",	"Flavobacteriia"="tan2",	"Sphingobacteriia"="indianred4",	"Bacteroidetes_unclassified"="wheat4",	"Opitutae"="slategray3",	"Verrucomicrobiae"="slategray4","OM190"="snow3","Phycisphaerae"="peachpuff1",	"Planctomycetacia"="palevioletred3",	"Cyanobacteria"="seagreen4","Deferribacteres"="salmon3",	"Chlamydiae"="paleturquoise4",	"Bacteria_unclassified"="black"))+ theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + theme(panel.background =element_rect(fil="white", color="black")) + theme(axis.text.x =element_text(colour = 'black', size = 10.5, angle=45, hjust = 1, family='serif')) + theme(axis.text.y=element_text(colour = 'black', size = 10.5,family='serif')) + theme(axis.title.x=element_text(colour = 'black',size = 12, family='serif')) + theme(axis.ticks.x = element_blank()) + theme(axis.title.y=element_text(colour = 'black',size = 12, family='serif')) + 
theme(legend.title=element_text(size=10, face='bold', family='serif'))  + theme(legend.text=element_text(size=10, family='serif')) + theme(legend.key.size=unit("0.38","cm"))  + labs(title="", y="Relative proportion %", x="Subcommunities") 

#For Supplementary Figure S4 B
pdf("./Class_Gene.pdf", width =6, height=6.5)
plot(p7)
dev.off()


#############END #########
