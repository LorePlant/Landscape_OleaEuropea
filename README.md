## Landscape Olea europea

This page is created to track progresses on my postdoctoral research in modelling genomic offset in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.
# data input
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
PCbio = data359[ ,15:24]
Env <- scale(PCbio, center=TRUE, scale=TRUE)
Env <- as.data.frame(Env)

#combining geographic, Popstructure, environmental (scaled) variables
Variables <- data.frame(data359$IDSample, data359$long, data359$lat, data359$group,data359$latitude_range, data359$region, Env)
names(Variables)[1]<- paste("geno")
names(Variables)[2]<- paste("long")
names(Variables)[3]<- paste("lat")
names(Variables)[4]<- paste("group")
names(Variables)[5]<- paste("latitude_range")
names(Variables)[6]<- paste("region")
 ```
To reduce collinearity, I want to check if the selected environmental variance have low VIF variance inflation factor

```
RDAgeo_env <- rda(genotype ~ bio1+bio2+bio4+bio6+bio8	+ bio9 + bio12 + bio14+	bio15	+ bio19, Variables)

sqrt(vif.cca(RDAgeo_env))
```
| bio1    |  bio2    |  bio4   |   bio6   |   bio8    |  bio9   |  bio12  |   bio14      |   bio15  |   bio19 |
|---------|----------|---------|----------|-----------|----------|---------|-----------|----------|-----------|
|22.380360 | 8.374068 |25.856036| 15.688084 | 2.594527| 22.626527 | 5.328737 | 4.204566|3.296078 | 5.147522 |


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
  xlab("RDA 1: 71 %") + ylab("RDA 2: 11 %") +
  guides(color=guide_legend(title="Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_env.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```

![RDA_env](https://github.com/user-attachments/assets/af247ac7-2998-48a8-8f72-eed424309d55)


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
  xlab("RDA 1: 71%") + ylab("RDA 2: 11 %") +
  guides(color=guide_legend(title="latitude range")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_geo_env_lat_range.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```

![RDA_env_region](https://github.com/user-attachments/assets/cfead7dc-e084-4460-a171-9fdc1cde5323)


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
  xlab("RDA 1: 71%") + ylab("RDA 2: 11 %") +
  guides(color=guide_legend(title="latitude range")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_geo_env_region.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```

![RDA_geo_env_lat_range](https://github.com/user-attachments/assets/06eb813c-7f89-44e9-af05-3a03d3890e15)


The result show a clear differentiation between Wild and Admixed populations. The two groups are mainly divided along the RDA1 component which is positively correlated with bio6 (Min Temperature of Coldest Month) and bio 15 (Precipitation Seasonality) . The result suggest that the wild populations can trive in warmer winters, and drier summers, compared to the admixed group. I would speculate from this outcome that the introgression of cultivated genepool can decrease the potential adaptation in future environmental scenarios were temperature levels are forecast to increase.





