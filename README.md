# Abundance-Occupancy-along-taxonomicRanks

Abundance-occupancy relationships (AO) is a relation describing the abundance and the number of sites that each speices attain on a regional scale. Quantifying the strength of AO relationships in microbial communities in natural systems incorporated with phylogenetic framework and life-history styles is currently not known. In this study, we first assesed the strength of AO relationships along taxonomic ranks within particle-attached and free-living assemblages across depths in St. Lawrence Estuary. Furthermore, the slope of AO relationships along taxonomic ranks was calculated and served as the estimate of niche divergence rates of bacterial communities surveyed.

_This work has been published in Frontiers in Microbiology http://dx.doi.org/10.3389/fmicb.2021.690712

The codes and computing notes in this Github repository are also published in Zenodo: https://doi.org/10.5281/zenodo.4743168. 



# To cite this work or code:

Izabel-Shen, D., Höger, A., Jürgens, K. (2021, In Press). Abundance-occupancy relationships along taxonomic ranks reveal a consistency of niche differentiation in marine bacterioplankton with distinct lifestyles. Frontiers in Microbiology. doi: 10.3389/fmicb.2021.690712 .


Corresponding author: Dandan Izabel-Shen, Department of Ecology, Environment and Plant Sciences, Stockholm University; email: dand.shen@gmail.com


## Installation

These scripts have been tested with R 4.0.0, running under: macOS 10.14.2.

The following packages (and their dependencies) are required to run the whole analyses

R version 4.0.0
```
ggplot2 (3.3.3)
vegan (2.5.7)
ape (5.4.1)
GUniFrac (1.1)
gridExtra (2.3)
RVAideMemoire (0.9.79)
scales (1.1.1)
graphics (3.4.0)
base (4.0.4)
picante (1.8.2)
dada2 (1.18.0)
readxl (1.3.1)
Biostrings (2.58.0)
ShortRead (1.48.0)


```

## How to run
The R scripts in the base folder should be run interactively (for instance with RStudio).

'InputFiles' contains all data files required to run the analyses. In the repository, the input files are compressed with gzip. Before you proceed with the analyses, the files should be uncompressed.
