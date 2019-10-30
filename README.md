# Abundance-Occupancy-along-taxonomicRanks

A repository for R code, data, figures, scripts from the manuscript "Abundance-occupancy relationships along taxonomic ranks reveals a high degree of niche divergence at broad taxonomic levels" (in preparation). Abundance-occupancy relationships (AO) is a relation describing the abundance and the number of sites that each speices attain on a regional scale. Quantifying the strength of AO relationships in microbial communities in natural systems incorporated with phylogenetic framework and life-history styles is currently not known. In this study, we first assesed the strength of AO relationships along taxonomic ranks within particle-attached and free-living assemblages across depths in St. Lawrence Estuary. Furthermore, the slope of AO relationships along taxonomic ranks was calculated and served as the estimate of niche divergence rates of bacterial communities surveyed.

by D Shen, A Höger, K Jürgens

Prepared 30 January 2019

Author: Dandan Shen, Leibniz Institute for Baltic Sea Research Warnemünde (IOW), Germany; dand.shen@gmail.com

## Installation

These scripts have been tested with R 3.4.0, running under: macOS 10.14.2.

The following packages (and their dependencies) are required to run the whole analyses

```
ggplot2 2.2.1
vegan 2.4-3
ape 4.1
GUniFrac 1.0
gridExtra 2.2.1
RVAideMemoire 0.9-69-3
scales 0.4.1
graphics 3.4.0
base 3.4.0

```

## How to run
The R scripts in the base folder should be run interactively (for instance with RStudio).

The folder 'InputFiles' contains all data files required to run the analyses. In the repository, the input files are compressed with gzip. Before you proceed with the analyses, the files should be uncompressed.
