## Landscape Olea europea

This page is created to track progresses on my postdoctoral research in modelling genomic offset in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.
## data input
I started by entering the vcf into R using the vcfR package as follow
```
library(vcfR)
library(adegenet)

setwd("/lustre/rocchettil")
genoLAND.VCF <- read.vcfR("359_Olive_west_MAF005.vcf.recode.vcf")#import vcf file
gl.genoLAND <- vcfR2genind(genoLAND.VCF)#transfrom file in genind object
genotype<-as.data.frame(gl.genoLAND)
genotype<-tibble::rownames_to_column(genotype, "geno") #transform raw name in column

```
Considering that downstream analysis like RDA do not work with NA values I found the following R for cycle for genetic data imputation. This code can be found in https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd

```
for (i in 1:ncol(genotype))
{
  genotype[which(is.na(genotype[,i])),i] <- median(genotype[-which(is.na(genotype[,i])),i], na.rm=TRUE)
}

write.table(genotype, "geno_359_west_olive_MAF005__imputated.txt")

#once the dataset was created I can enter the data table with read.table
genotype<- read.table("geno_359_west_olive_MAF005__imputated.txt", header=TRUE)
```
The following excel file data359 which includes uncorrelated (r<0.7) bioclimatic variables as well as latitude and longitude was uploaded

```
data359<- read.csv("dataset_359_olive.csv", header = TRUE)
```
# Mantel test
The Mantel test allows to conduct a linear regression analysis between the genetic distance and environmental distance. The significance of this regression suggests a Isolation by Environment IBE effect, where individual in ecologically similar locations they are more genetically similar compare to individuals in ecologically diverse locations.

As first step I'm going to estimate the genetic distance as pairwise FST among the 27 populations using the R package Hierfstat.
```
library("hierfstat")
pops<-read.table("27_Pops.txt", header=T) #one column table wih pop info for each individual
#convert genInd and pop in hierfstat
hierfstat<-genind2hierfstat(gl.genoLAND ,pop=pops)
genet.dist(hierfstat,diploid=TRUE,method="WC84")
```


In the next chuck of codes we are going to define the Euclidean distqnce for geogrqphy, enviroment and genetic with the final aim to conduct a Mantel test
```
#bioclim PCdata frame
PCbio = data359[,15:24]
Env <- scale(PCbio, center=TRUE, scale=TRUE)
dist.PCbio = dist(Env, method = "euclidean")

#geographic data
geo = data.frame(data359$long, data359$lat)
dist.geo = dist(geo, method = "euclidean")

#genetic distance
distgenEUCL <- dist(gl.genoLAND, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
```
#Mantel test
With the aim if there to underline Isolation by Environment (IBE) I used a Mantel test to see if there is a linear correlation between ecological distance and genetic distance matrices
```
# mantel test Genetic distance-ecological distance
geno_eco = mantel(distgenEUCL, dist.PCbio, method="spearman", permutations=1000,  na.rm = TRUE)
geno_eco
summary(lm(dist.PCbio~distgenEUCL))
graph = mantel.correlog(distgenEUCL, dist.PCbio, XY=NULL, n.class=0, break.pts=NULL, 
                        cutoff=TRUE, r.type="pearson", nperm=999, mult="holm", progressive=TRUE)


xx = as.vector(distgenEUCL) #convert distance matrix into a vector
zz = as.vector(dist.PCbio)
manatelmatrix = data.frame(zz,xx)
mm = ggplot(manatelmatrix, aes(y = xx, x = zz))+
  geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21,fill = "grey") + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2)+
  labs(y = "Euclidean genetic distance", x = "Euclidean ecological distance")+
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18, colour = "black"), 
         axis.title= element_text(face = "bold", size = 18, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =12, face = "bold", colour = "black"),
         legend.text = element_text(size = 10, face = "bold", colour = "black"), 
         legend.position = "top", strip.background = element_rect(fill = "grey90", colour = "black"),
         strip.text = element_text(size = 9, face = "bold"))
jpeg(file = "/lustre/rocchettil/mantel_olive_geno_eco.jpeg", width = 350, height = "350")
plot(mm)
dev.off()

```

![mantel_olive_geno_eco](https://github.com/user-attachments/assets/9f58ecf7-48fb-4bfa-8181-8006f5a6b857)

Results suggests a moderate (r: 0.21) though significant (P<0.01) correlation between genetic euclidean distance and ecological euclidean distances. 
Considering the potential effect of geography in the ecological distance I used partial Mantel test which test for correlation between genetic distance and environmental distance considering geographic distance as covariate


```
#partial Mantel test 
partial_mantel = mantel.partial(distgenEUCL, dist.PCbio, dist.geo, method = "spearman", permutations = 1000,
                                na.rm = TRUE)
partial_mantel
summary(lm(distgenEUCL~dist.PCbio|dist.geo))
#plotting partial Mantel test
xx = as.vector(distgenEUCL) #convert distance matrix into a vector
yy= as.vector(dist.geo)
zz = as.vector(dist.PCbio)
partial_mantel_matrix = data.frame(xx,zz,yy)#new data frame with vectorize distance matrix

mp = ggplot(partial_mantel_matrix, aes(y = xx, x = zz)) + 
  geom_point(size = 2.5, alpha = 0.75, colour = "black",shape = 21, aes(fill = yy)) + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) + 
  scale_fill_continuous(high = "navy", low = "lightblue")+
  labs(y = "Eucledian genetic distance", x = "Eucledian ecological distance", fill= "geographic distance")+
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18, colour = "black"), 
         axis.title= element_text(face = "bold", size = 18, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =12, face = "bold", colour = "black"),
         legend.text = element_text(size = 10, face = "bold", colour = "black"), 
         legend.position = "top", strip.background = element_rect(fill = "grey90", colour = "black"),
         strip.text = element_text(size = 9, face = "bold"))


jpeg(file = "/lustre/rocchettil/partial_mantel_olive.jpeg")
plot(mp)
dev.off()
```
The result confirmes the significant (P<0.01) and moderate correlation (r:0.16) between environmental distance and genetic distance. Overall this result suggest a potential Isolation-by-envirnoment effect IBE on the sampled population.

![partial_mantel_olive](https://github.com/user-attachments/assets/dd34f232-d6d1-4e13-81d4-a06992c4e8f9)


# Redundancy analysis

Within the landscape genomic framework, Redundancy analysis (RDA) represent a useful tool that allows to dissect the the total genetic variance among the environment, geographic and demographic components. 
In this first analysis I used RDA on the following linear model to see if we can detect specif environmental variables diverging Wild vs Admixed genotypes or geographic regions.
 $` Gen \sim Environment `$

The follwing chuck of code illustrates the step undertaken for assembly the dataset for RDA. The main step is the standardization of environmental variables
```
#standardize bioclim variable
PCbio = data359[ ,16:29]
Env <- scale(PCbio, center=TRUE, scale=TRUE)
Env <- as.data.frame(Env)

#combining geographic, Popstructure, environmental (scaled) variables
Variables <- data.frame(data359$IDSample, data359$long, data359$lat, data359$group,data359$latitude_range, data359$region, data359$PC1, Env)
names(Variables)[1]<- paste("geno")
names(Variables)[2]<- paste("long")
names(Variables)[3]<- paste("lat")
names(Variables)[4]<- paste("group")
names(Variables)[5]<- paste("latitude_range")
names(Variables)[6]<- paste("region")
names(Variables)[7]<- paste("PC1")
 ```
To reduce collinearity, I want to check if the selected environmental variance have low VIF variance inflation factor

```
RDAgeo_env <- rda(genotype ~ bio1+bio2+bio4+ bio5 + bio6+bio8	+ bio9 + bio10+bio11+ bio12 + bio14+	bio15	+ bio18 + bio19, Variables)

sqrt(vif.cca(RDAgeo_env))
```
| bio1    |  bio2    |  bio4   |   bio5   |   bio6    |  bio8   |   bio9 | bio10 |   bio11      |   bio12  |   bio14 | bio15 | bio18 | bio19 |
|---------|----------|---------|----------|-----------|----------|---------|-----------|----------|---------|--------|-------|-----|---------|
|29.004607 | 24.449079  |49.46922| 15.210764 | 29.510669| 3.320524 |30.741673 |22.032603|39.983301| 6.3260372 |5.625297|4.254793|3.822378 | 5.766060|


I attempted to select significant ecological variables with a Variance Inflation Factor (VIF) lower than 10. Considering the termopluviometric graph for the study area, which includes France, Spain, and Morocco, we observe a common bell-shaped temperature distribution. Temperatures are higher in the spring and summer, while lower during the autumn and winter. However, precipitation patterns differ across regions. For example, the internal regions of France experience a more uniform distribution of precipitation throughout the year.

To address specific temperature and precipitation adaptations, I carefully selected bioclimatic variables exluding those reffering temperature and precipitation to the wettest or driest quarters. I focused instead on variables that reflect temperature and precipitation variations during the warmest and coldest quarters, allowing to precisly asses temperature and precipitation during winter and summer quarters. Following this selection, I checked the Variance Inflation Factor (VIF) for the chosen variables to ensure they were appropriate.


I selected bio2, bio6, bio8, bio12, bio14, bio15 and bio 19 and run again the VIF analysis

```
RDAgeo_env <- rda(genotype ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19, Variables)

sqrt(vif.cca(RDAgeo_env))
```

|  bio2    |   bio10   |   bio11  |  bio15  |   bio18  |   bio19  | 
|---------|----------|---------|----------|-----------|----------|
1.655992| 1.287833 |2.315390 |2.776476 |1.942473| 1.371975|

The variable selected showed low inflation factor due to correlation


In the next part I'm preparing the data for plotting
```
#write.table(score, "Genotypevalue_RDAgeo_env")
summary(eigenvals(RDAgeo_env, model = "constrained"))
#data_RDAgeo_env <- read.table(file = "clipboard", sep = "\t", header=TRUE)
score<-scores(RDAgeo_env , display = "sites")
TAB_gen <- data.frame(geno = row.names(score), score)
dataRDA<-merge(Variables, TAB_gen, by = "geno")
#install.packages("ggrepel")
#library(ggrepel)
TAB_var <- as.data.frame(scores(RDAgeo_env, choices=c(1,2), display="bp"))
```
I would like to prepare three plot highlighting WLDvsADM, latitude ranges and geographic regions

> Wild vs Admixed
```
##color ADM vs WLD
loading_RDAgeo_env<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = dataRDA, aes(x=RDA1, y=RDA2, color=group), size = 4.5) +
  scale_color_manual(values = c("blue", "darkorange")) + 
  geom_segment(data = TAB_var, aes(xend=RDA1*10, yend=RDA2*10, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*10, y=RDA2*11, label = row.names(TAB_var)), size = 4.5, family = "Times") +
  xlab("RDA 1: 74 %") + ylab("RDA 2: 8 %") +
  guides(color=guide_legend(title="Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_env.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```

![RDA_env](https://github.com/user-attachments/assets/8e2b2e1a-b648-48d4-ac8f-bf072b731df5)



>Latitude ranges
```
#color latitude range
loading_RDAgeo_env<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = dataRDA, aes(x=RDA1, y=RDA2, color=latitude_range), size = 4.5) +
  scale_color_manual(values = c("darkgreen","red", "darkorange")) + 
  geom_segment(data = TAB_var, aes(xend=RDA1*10, yend=RDA2*10, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*10, y=RDA2*11, label = row.names(TAB_var)), size = 4.5, family = "Times") +
  xlab("RDA 1: 74%") + ylab("RDA 2: 8 %") +
  guides(color=guide_legend(title="latitude range")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_geo_env_lat_range.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```


![RDA_geo_env_lat_range](https://github.com/user-attachments/assets/77470479-2423-4acd-8781-eb9ec1f32a3d)



>Geographic regions: France, Corse, Spain, Morocco

```
#color regions
loading_RDAgeo_env<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = dataRDA, aes(x=RDA1, y=RDA2, color=region), size = 4.5) +
  scale_color_manual(values = c("darkgreen","purple", "darkorange", "blue")) + 
  geom_segment(data = TAB_var, aes(xend=RDA1*10, yend=RDA2*10, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*10, y=RDA2*11, label = row.names(TAB_var)), size = 4.5, family = "Times") +
  xlab("RDA 1: 74%") + ylab("RDA 2: 8 %") +
  guides(color=guide_legend(title="latitude range")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_geo_env_region.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```


![RDA_geo_env_region](https://github.com/user-attachments/assets/31019223-6a4c-4687-bbc6-d4965e3315a8)



The result show a clear differentiation between Wild and Admixed populations. The two groups are mainly divided along the RDA1 component which is positively correlated with bio6 (Min Temperature of Coldest Month) and bio 15 (Precipitation Seasonality) . The result suggest that the wild populations can trive in warmer winters, and drier summers, compared to the admixed group. I would speculate from this outcome that the introgression of cultivated genepool can decrease the potential adaptation in future environmental scenarios were temperature levels are forecast to increase.

RDA can be used for variance partitioning




#RDA for GEA discovery
Redundancy analysis can be used to identify GEA based on the Mhallanoise distance of SNPs in the RDA-biplot. Within the RDA model we can effectively correct for population structure (PC1) and Isolation by distanc (lqtitude and longitude) using them as covariates in the RDA model
As first attempt I decided to run the anlysis seperate for temperature and precipitation variables.

>Temperature

```
RDA_temp <- rda(genotype ~ bio2+bio10+bio11 +  Condition(PC1 + lat + long), Variables)
summary(eigenvals(RDA_temp, model = "constrained"))
library(robust)
remotes::install_github("koohyun-kwon/rdadapt")
source("./src/rdadapt.R")
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

rdadapt_env<- rdadapt(RDA_temp, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genotype)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))
write.table(top_outliers, "Bonferroni_temp")
qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values, q.value = rdadapt_env$q.value)
outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$q.values<0.05)], p.value = rdadapt_env$p.values[which(rdadapt_env$q.values<0.05)])

locus_scores <- scores(RDA_temp, choices=c(1:2), display="species", scaling="none")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
TAB_var <- as.data.frame(scores(RDA_temp, choices=c(1,2), display="bp"))
loading_temp<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*40, y=RDA2*40, colour = type), size = 2.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1: 46%") + ylab("RDA 2: 29%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_temp
jpeg(file = "/lustre/rocchettil/RDA_temp_biplot.jpeg")
plot(loading_temp)
dev.off()

write.table(qvalue, "Temp_GEA_Olive")

#plotting Mhanattan plot using the library qqman

library(qqman)
Manhattan_temp <- read.csv(file = "tempGEA.csv", header=TRUE) #import the p value result for temperature
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.000909433), genomewideline = -log10(3.797084e-07))
jpeg(file = "/lustre/rocchettil/Manh_RDA_temp.jpeg")
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.000909433), genomewideline = -log10(3.797084e-07))
dev.off()
```
![RDA_temp_biplot](https://github.com/user-attachments/assets/aa4cc238-7d40-418d-bb22-924c609748fa)
![Manh_RDA_temp](https://github.com/user-attachments/assets/383f0850-efb7-4754-92de-2eca06ccc75b)
![Phist_Manh_RDA_temp](https://github.com/user-attachments/assets/1c07b997-b249-4e1a-b43b-f0881ce38239)


> Precipitation
```
RDA_prec <- rda(genotype ~ 	bio15	+ bio18 + bio19 +  Condition(PC1 + lat + long), Variables)
summary(eigenvals(RDA_prec, model = "constrained"))
library(robust)
remotes::install_github("koohyun-kwon/rdadapt")
source("./src/rdadapt.R")
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

rdadapt_env<- rdadapt(RDA_prec, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genotype) [which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))
write.table(top_outliers, "Bonferroni_prec")
qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values, q.value = rdadapt_env$q.value)
outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$q.values<0.05)], p.value = rdadapt_env$p.values[which(rdadapt_env$q.values<0.05)])

locus_scores <- scores(RDA_prec, choices=c(1:2), display="species", scaling="none")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
TAB_var <- as.data.frame(scores(RDA_prec, choices=c(1,2), display="bp"))
loading_prec<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*40, y=RDA2*40, colour = type), size = 2.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1: 44%") + ylab("RDA 2: 32%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_prec
jpeg(file = "/lustre/rocchettil/RDA_prec_biplot.jpeg")
plot(loading_prec)
dev.off()


write.table(qvalue, "Prec_GEA_Olive")

#plotting Mhanattan plot using the library qqman

library(qqman)
Manhattan_prec <- read.csv(file = "precGEA.csv", header=TRUE) #import the p value result for precipitation
jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
manhattan(Manhattan_prec, col = c("blue", "gray60"),suggestiveline = -log10(0.000100972), genomewideline = -log10(3.797084e-07))
dev.off()
```
![RDA_prec_biplot](https://github.com/user-attachments/assets/2f6a3ecc-9928-4c98-95f7-7781aea76802)
![Manh_RDA_prec](https://github.com/user-attachments/assets/8c493917-6826-4495-8ece-25f1076237f0)
![Phist_Manh_RDA_prec](https://github.com/user-attachments/assets/b539de2e-1114-4c09-8ca4-16e627b53b27)

>All together

In this attempt I am going to run the GEA analysis considering all the bioclimatic variable selected. The derived GEA will be used for the adaptive index projection and Genomic offset estimation

```
RDA_all <- rda(genotype ~ 	bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 +  Condition(PC1 + lat + long), Variables)
summary(eigenvals(RDA_all, model = "constrained"))
library(robust)
remotes::install_github("koohyun-kwon/rdadapt")
source("./src/rdadapt.R")
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

rdadapt_env<- rdadapt(RDA_all, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genotype) [which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values, q.value = rdadapt_env$q.value)
outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$q.values<0.05)], p.value = rdadapt_env$p.values[which(rdadapt_env$q.values<0.05)])

locus_scores <- scores(RDA_all, choices=c(1:2), display="species", scaling="none")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
TAB_var <- as.data.frame(scores(RDA_all, choices=c(1,2), display="bp"))
loading_all<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*40, y=RDA2*40, colour = type), size = 2.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1: 44%") + ylab("RDA 2: 32%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_all
jpeg(file = "/lustre/rocchettil/RDA_all_biplot.jpeg")
plot(loading_all)
dev.off()
```
![RDA_all_biplot](https://github.com/user-attachments/assets/8f76c89e-6453-477a-9273-b3ca15eaf00e)



To visualize the adaptive differentiation among genotypes, I conducted an additional Redundancy Analysis (RDA) using only the previously identified GEA SNPs. In this analysis, I did not include geography and population structure as covariates for two reasons: First, I aimed to observe the differentiation between wild and admixed genotypes. Second, the GEA SNPs used have already been identified with corrections for population structure and geography.

```
#partial redundancy analysis (RDA only with GEA QTL)
geno_all_enrich<-genotype[which(rdadapt_env$q.values<0.05)]
RDA_all_enriched<-rda(geno_all_enrich ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19, Variables)
summary(eigenvals(RDA_all_enriched, model = "constrained"))


#plot genotypes

TAB_gen <- data.frame(geno = row.names(score), scores(RDA_all_enriched , display = "sites"))
Geno <- merge(TAB_gen, Variables[, 1:7] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, colour = group), size = 2.5) +
  scale_color_manual(values = c("blue", "darkorange")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*10, yend=RDA2*10, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=10*RDA1, y=10*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1: 46%") + ylab("RDA 2: 26%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_all_enriched
jpeg(file = "/lustre/rocchettil/RDA_all_geno_biplot.jpeg")
plot(loading_geno_all_enriched)
dev.off()
```
![RDA_all_geno_biplot](https://github.com/user-attachments/assets/5c7bdf6e-a3fc-4900-acf9-69ecfd965ac2)



# Adaptive index projection
Adaptive indeix function Capblach.
Using the raster environmental data the function allows to predict the adaptie value of each pixel

```
adaptive_index <- function(RDA, K, env_pres, range = NULL, method = "loadings", scale_env, center_env){
  
  # Formatting environmental rasters for projection
  var_env_proj_pres <- as.data.frame(rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]]))
  
  # Standardization of the environmental variables
  var_env_proj_RDA <- as.data.frame(scale(var_env_proj_pres[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))
  
  # Predicting pixels genetic component based on RDA axes
  Proj_pres <- list()
  if(method == "loadings"){
    for(i in 1:K){
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(apply(var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Prediction with RDA model and linear combinations
  if(method == "predict"){ 
    pred <- predict(RDA, var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], type = "lc")
    for(i in 1:K){
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Returning projections for current climates for each RDA axis
  return(Proj_pres = Proj_pres)
}
```
Run a RDA analysis using only the GEA marker previously identify using a FDR (q<0.05) threshold

```
geno_enrich<-genotype[which(rdadapt_env$q.values<0.05)]
RDA_temp_enriched<-rda(geno_enrich ~ bio2+bio6+bio8 +  Condition(PC1 + lat + long), Variables)
summary(eigenvals(RDA_temp_enriched, model = "constrained"))
jpeg(file = "/lustre/rocchettil/RDA_tempenriched_biplot.jpeg")
plot(RDA_temp_enriched)
dev.off()
```
![RDA_tempenriched_biplot](https://github.com/user-attachments/assets/c595862b-8329-43bc-b30c-0db4915bc0ca)

> Recovering scaling coefficients
Create a table of scaled temperature variable
```

library('dplyr')
PCbio = data359[ ,16:25]
tempvar = PCbio %>% select('bio2','bio6', 'bio8')
Temp <- scale(tempvar, center=TRUE, scale=TRUE)
scale_Temp <- attr(Temp, 'scaled:scale')
center_Temp <- attr(Temp, 'scaled:center')


Env <- scale(PCbio, center=TRUE, scale=TRUE)

scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')
```
Enter the raster file for the specific bioclimatic variable for current climatic situation
```
# ras temperature
library(raster)
library("readxl")


bio2<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio2_studyarea_ext.tif"))
bio6<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio6_studyarea_ext.tif"))
bio8<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio8_studyarea_ext.tif"))
names(bio2) = 'bio2'
names(bio6) = 'bio6'
names(bio8) = 'bio8'
#stack the different raster file
ras_current_temp<-stack(c(bio2,bio6, bio8))
```
Predict tha adaptive index for each pixel grid

```
## Function to predict the adaptive index across the landscape
source("./src/adaptive_index.R")

res_RDA_temp_proj_current <- adaptive_index(RDA = RDA_temp_enriched, K = 2, env_pres = ras_current_temp, range = range, method = "loadings", scale_env = scale_Temp, center_env = center_Temp)
projection<- stack(c(res_RDA_temp_proj_current$RDA1, res_RDA_temp_proj_current$RDA2))
plot(projection)

## Vectorization of the climatic rasters for ggplot
RDA_proj <- list(res_RDA_temp_proj_current$RDA1, res_RDA_temp_proj_current$RDA2)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}

## Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
distrib_temp<- ggplot(data = TAB_RDA) + 
  geom_tile(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  #coord_sf(xlim = c(-148, -98), ylim = c(35, 64), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw(base_size = 7, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=7))
jpeg(file = "/lustre/rocchettil/Adaptive_temp_proj.jpeg", height=1500, width=3000, res=600)
plot(distrib_temp)
dev.off()
```

![Adaptive_temp_proj](https://github.com/user-attachments/assets/97f90439-c8c8-419e-b158-620d281a63a1)


>Genomic offset RDA based

```
#future temp scenario
# ras temperature
library(raster)
library("readxl")


bio2<- raster(paste("/lustre/rocchettil/bio2_ext_studyarea.tif"))
bio6<- raster(paste("/lustre/rocchettil/bio6_ext_studyarea.tif"))
bio8<- raster(paste("/lustre/rocchettil/bio8_ext_studyarea.tif"))
names(bio2) = 'bio2'
names(bio6) = 'bio6'
names(bio8) = 'bio8'
#stack the different raster file
ras_2050_temp<-stack(c(bio2,bio6, bio8))
```
function genomic offset
```
#### Function to predict genomic offset from a RDA model
genomic_offset <- function(RDA, K, env_pres, env_fut, range = NULL, method = "loadings", scale_env, center_env){
  
  
  # Formatting and scaling environmental rasters for projection
  var_env_proj_pres <- as.data.frame(scale(rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]])[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))
  var_env_proj_fut <- as.data.frame(scale(rasterToPoints(env_fut[[row.names(RDA$CCA$biplot)]])[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))

  # Predicting pixels genetic component based on the loadings of the variables
  if(method == "loadings"){
    # Projection for each RDA axis
    Proj_pres <- list()
    Proj_fut <- list()
    Proj_offset <- list()
    for(i in 1:K){
      # Current climates
      ras_pres <- env_pres[[1]]
      ras_pres[!is.na(ras_pres)] <- as.vector(apply(var_env_proj_pres[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      # Future climates
      ras_fut <- env_fut[[1]]
      ras_fut[!is.na(ras_fut)] <- as.vector(apply(var_env_proj_fut[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))
      Proj_fut[[i]] <- ras_fut
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Predicting pixels genetic component based on predict.RDA
  if(method == "predict"){ 
    # Prediction with the RDA model and both set of envionments 
    pred_pres <- predict(RDA, var_env_proj_pres[,-c(1,2)], type = "lc")
    pred_fut <- predict(RDA, var_env_proj_fut[,-c(1,2)], type = "lc")
    # List format
    Proj_offset <- list()    
    Proj_pres <- list()
    Proj_fut <- list()
    for(i in 1:K){
      # Current climates
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_pres[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      # Future climates
      ras_fut <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_fut[,i])), crs = crs(env_pres))
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      Proj_fut[[i]] <- ras_fut
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Weights based on axis eigen values
  weights <- RDA$CCA$eig/sum(RDA$CCA$eig)
  
  # Weighing the current and future adaptive indices based on the eigen values of the associated axes
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_pres[[x]])[,-c(1,2)]))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x]*weights[x])))
  Proj_offset_fut <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut[[x]])[,-c(1,2)]))
  Proj_offset_fut <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut[,x]*weights[x])))
  
  # Predict a global genetic offset, incorporating the K first axes weighted by their eigen values
  ras <- Proj_offset[[1]]
  ras[!is.na(ras)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut[x,]), method = "euclidean")))
  names(ras) <- "Global_offset"
  Proj_offset_global <- ras
  
  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset 
  return(list(Proj_pres = Proj_pres, Proj_fut = Proj_fut, Proj_offset = Proj_offset, Proj_offset_global = Proj_offset_global, weights = weights[1:K]))
}
```

Genomic offset projection
```
res_RDA_temp_proj_2050 <- genomic_offset(RDA = RDA_temp_enriched, K = 2, env_pres = ras_current_temp, env_fut = ras_2050_temp, range = range, method = "loadings", scale_env = scale_Temp, center_env = center_Temp)

```



# Gradient Forest
Gradient Forest is an alternative approach widely use in landscape genomics studies, where the relation between genetic component and environmental component is constructed using the random forest machine learning approach.
This code is still under construction. In this part I'm keeping track of the progresses achived.

In this example I used a genotipic data file only from the Wild group
```
#Genomic offset runGF only on wild (file from Lison)
install.packages("gradientForest", repos="http://R-Forge.R-project.org")

library(vcfR)
library(adegenet)

setwd("/lustre/rocchettil")

genoLAND.VCF <- read.vcfR("Oe9_genuine_clean_outliers.annotated.vcf")#import vcf file
geno_genuine <- vcfR2genind(genoLAND.VCF)#transfrom file in genind object
geno_wild<-as.data.frame(geno_genuine)

for (i in 1:ncol(geno_wild))
{
  geno_wild[which(is.na(geno_wild[,i])),i] <- median(geno_wild[-which(is.na(geno_wild[,i])),i], na.rm=TRUE)
}
```
I created with excell a table with the indivdual name and the different bioclimatic variable 

```
Env_tab<- read.table("Env_tab_WLD.txt")
```
To run the Gradient Forest function I used the gradient forest package.
In this link there is a guide https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf

Once we are sure that the geno_wild dataset and Env_tab dataset are in the same order for sample IDs we can apply the following code

```
library(gradientForest)


gf <- gradientForest(cbind(geno_wild, Env_tab), 
                     predictor.vars=colnames(Env_tab),
                     response.vars=colnames(geno_wild), 
                     ntree=500, #set the number of individual decision tree
                     trace=TRUE)
```
From this function we can print:

>the Environmental variable importance.

>split (split node of decison tree; their order in the decision tree reflects the variable importance) density plot.

```
most_important <- names(importance(gf))[1:25]
par(mgp = c(2, 0.75, 0))
plot(gf)
plot(gf, plot.type = "S", imp.vars = most_important,leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6, cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
```
Once the gradient forest model is created we can use it to estimate the adaptive component of each environmental pixel data. This function allows to map at the geographic level the biological adaptive space.
To do so we use raster file from CHELSA dataset previosuly clipped for our study area using QGIS. Each raster represent a biovariable value for each pixel. 
The first step is to stack all the raster information
```
library(raster)
library("readxl")

bio1<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio1_studyarea_ext.tif"))
bio2<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio2_studyarea_ext.tif"))
bio4<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio4_studyarea_ext.tif"))
bio6<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio6_studyarea_ext.tif"))
bio8<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio8_studyarea_ext.tif"))
bio9<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio9_studyarea_ext.tif"))
bio12<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio12masked_studyarea_ext.tif"))
bio14<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio14_studyarea_ext.tif"))
bio15<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio15_studyarea_ext.tif"))
bio19<- raster(paste("/lustre/rocchettil/biovar_studyarea/bio19_studyarea_ext.tif"))
names(bio1) = 'bio1'
names(bio2) = 'bio2'
names(bio4) = 'bio4'
names(bio6) = 'bio6'
names(bio8) = 'bio8'
names(bio9) = 'bio9'
names(bio12) = 'bio12'
names(bio14) = 'bio14'
names(bio15) = 'bio15'
names(bio19) = 'bio19'
#stack the different raster file
ras_current<-stack(c(bio1, bio2, bio4, bio6, bio8, bio9, bio12, bio14, bio15, bio19))
```
The next step will be to trasforme the raster file into a centered spatial grid (x,y) where each environmental variable is given by the average of each pixel.

```
#spatial grid

coord_r<-rasterToPoints(ras_current, spatial = TRUE)
map_pts<-data.frame(x = coordinates(coord_r)[,1], y=coordinates(coord_r)[,2], coord_r@data)
```
Subsequently we are going to use the GF function to estimate the adaptive value for each environmental data point.
```
library(tidyr)
newmap_pts <- map_pts %>% drop_na()
imp.var<- names(importance(gf))
Trns_grid<-cbind(newmap_pts[, c("x", "y")], predict(gf, newmap_pts[, imp.var]))
```
The resulted data table is a list of grid cells with specific lat/long values and the estimate genetic adaptive value for each environmental variable.
The multi-dimentinal adaptive space can be efficiently plotted using a PCA. The resuts provides a biplot for each cell grid colored follwing a palette scale. The same color palette will be used to plot the cell grid in lat/long space.
```
#color PCs
PCs<- prcomp(Trns_grid[, imp.var])
a1<- PCs$x[,1]
a2<- PCs$x[,2]
a3<- PCs$x[,3]
r<- a1+a2
g<- -a2
b<- a3+a2-a1
r<- (r-min(r))/(max(r)-min(r))*255
g<- (g-min(g))/(max(g)-min(g))*255
b<- (b-min(b))/(max(b)-min(b))*255

#color biplot

nvs <- dim(PCs$rotation)[1]
vec<-c("bio15", "bio12", "bio14", "bio19", "bio1",  "bio2", "bio8",  "bio6",  "bio4", "bio9")
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 5
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) 
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal)
plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1)
points(PCs$rotation[!vind, 1:2]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec,1]), PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
                                                         
  

#plot grid
plot(Trns_grid[, c("x", "y")], pch = ".", cex = 3, col= rgb(r,g,b,max=255))
```

![adaptive_PCA_space](https://github.com/user-attachments/assets/bc91fde3-600a-4eee-97a1-d11ab677eb2f)

![adaptive_geo_space](https://github.com/user-attachments/assets/fe373f9c-8c3a-4b74-af5c-3592949d9ffa)
