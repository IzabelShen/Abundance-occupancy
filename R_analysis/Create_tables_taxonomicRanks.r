
###### This script is to create subtables which were later used as inputfiles, by spliting the OTU table according to respective taxonomic ranks, including species, genus, family, order, class, phylum.

## load R libraries for this session 
library(scales)
library(base)

## load OTU table containg relative abundance and taxonomic affiliation of each OTU
Tax_OTU_Rel <-read.csv("InputFiles/Tax_OTU_Rel.csv", row.names=1, header=T) ################in our study, The threshold of  ≥ 99% similarity in 16S rRNA gene sequences was used in our study to cluster OTUs, and these finely resolved microbial taxa correspond to bacterial “species”. Therefore, this OTU table is 'equivalent' to Species_Table##########################
########################## OTU level and species level are interchangable in our case##########################


## Aggregate Abundances at different taxonomic Levels by spliting the whole OTUs according to the taxonomic ranks ###############
Phylum_Table <-aggregate(Tax_OTU_Rel[, 6:41], by=list(Category=Tax_OTU_Rel$Phylum), FUN=sum)
Class_Table <-aggregate(Tax_OTU_Rel[, 6:41], by=list(Category=Tax_OTU_Rel$Class), FUN=sum)
Order_Table <-aggregate(Tax_OTU_Rel[, 6:41], by=list(Category=Tax_OTU_Rel$Order), FUN=sum)
Family_Table <-aggregate(Tax_OTU_Rel[, 6:41], by=list(Category=Tax_OTU_Rel$Family), FUN=sum)
Genus_Table <-aggregate(Tax_OTU_Rel[, 6:41], by=list(Category=Tax_OTU_Rel$Genus), FUN=sum)

## Save the relative abudnance count table based at different taxonomic ranks, these tables were later used as inputfiles #############
write.csv(Phylum_Table, file='InputFiles/Phylum_Table.csv')
write.csv(Class_Table, file='InputFiles/Class_Table.csv')
write.csv(Order_Table, file='InputFiles/Order_Table.csv')
write.csv(Family_Table, file='InputFiles/Family_Table.csv')
write.csv(Genus_Table, file='InputFiles/Genus_Table.csv')

## Check the number of taxononic classification ################
dim (Phylum_Table)    ##[1] 29 37
dim (Class_Table)     ##[1] 94 37
dim (Order_Table)     ##[1] 177 37
dim (Family_Table)    ##[1] 354 37
dim (Genus_Table)     ##[1] 733 37


