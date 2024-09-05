# Landscape Olea europea

This page is created to track progresses on my postdoctoral research in modelling genomic offset in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.

I started by entering the vcf into R using the vcfR package as follow
```
library(vcfR)
library(adegenet)

setwd("/lustre/rocchettil")
genoLAND.VCF <- read.vcfR("Oe9_admixed_and_genuine_wild_clean.mac1.maf.vcf")#import vcf file
gl.genoLAND <- vcfR2genind(genoLAND.VCF)#transfrom file in genind object
genotype<-as.data.frame(gl.genoLAND)
```
Considering that downstream analysis like RDA do not work with NA values I found the following R for cycle for genetic data imputation. This code can be found in https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd

```
for (i in 1:ncol(genotype))
{
  genotype[which(is.na(genotype[,i])),i] <- median(genotype[-which(is.na(genotype[,i])),i], na.rm=TRUE)
}

write.table(genotype, "geno_data_west_olive_imputated.txt")

#once the dataset was created I can enter the data table with read.table
genotype<- read.table("geno_data_west_olive_imputated.txt", header=TRUE)
```
The following excel file data359 which includes uncorrelated (r<0.7) bioclimatic variables as well as latitude and longitude was uploaded

```
data359<- read.csv("dataset_359_olive.csv", header = TRUE)
```
In the next chuck of codes we are going to define the Euclidean distqnce for geogrqphy, enviroment and genetic with the final aim to conduct a Mantel test
```
#bioclim PCdata frame
PCbio = data359[,13:22]
Env <- scale(PCbio, center=TRUE, scale=TRUE)
dist.PCbio = dist(Env, method = "euclidean")

#geographic data
geo = data.frame(data359$long, data359$lat)
dist.geo = dist(geo, method = "euclidean")

#genetic distance
distgenEUCL <- dist(gl.genoLAND, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
```



