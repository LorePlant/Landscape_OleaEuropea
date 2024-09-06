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
PCbio = data359[,13:22]
Env <- scale(PCbio, center=TRUE, scale=TRUE)
dist.PCbio = dist(Env, method = "euclidean")

#geographic data
geo = data.frame(data359$long, data359$lat)
dist.geo = dist(geo, method = "euclidean")

#genetic distance
distgenEUCL <- dist(gl.genoLAND, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
```
#Mantel test
With the aim if there to underline Isolation by Environment (IBE) I used a Mantel test to see if thre is a linear correlation betaeen ecological distqnce and genetic distance matrices
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
Considering the potential effect of geogrqphy in the ecological distance I used partial Mantel test which consider the correlation betaeen genetic distance and environmental distance considering geographic distance as covariate
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
The result show a correlation between environmental distance and geogrqphic distance, where geographic distqnce increases ( darker dots color) environmental distqnce also increase. Considering the correlation with genetic distance there is no correlation suggesting a non significant IBE

![mantel_olive](https://github.com/user-attachments/assets/314ccdd4-0ada-4c21-a1ae-f48330044726)

# Redundancy analysis

Within the landscape genomic framework, Redundancy analysis (RDA) represent a useful tool that allo to dissect the the total genetic variance among the environment, geographic and demographic components. 
In this first analysis I used the following RDA model to see if we can detect specif environmental variable diverging Wild vs Admixed genotypes or geographic regions.
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

```
##color ADM vs WLD
loading_RDAgeo_env<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = data_RDAgeo_env, aes(x=RDA1, y=RDA2, color=group), size = 4.5) +
  scale_color_manual(values = c("blue", "darkorange")) + 
  geom_segment(data = TAB_var, aes(xend=RDA1*10, yend=RDA2*10, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*10, y=RDA2*11, label = row.names(TAB_var)), size = 4.5, family = "Times") +
  xlab("RDA 1: 67 %") + ylab("RDA 2: 7 %") +
  guides(color=guide_legend(title="Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_geo_env.jpeg")
loading_RDAgeo_env
dev.off()
```
![RDA_env](https://github.com/user-attachments/assets/b2576bb0-411f-4d0b-919e-481dd8153557)


The result show a clear differentiation between Wild and Admixed populations. The two groups are mainly divided along the RDA1 component which is positively correlated with bio2 (Mean Diurnal Range (Mean of monthly (max temp - min temp)) and negatively with bio6 (Min Temperature of Coldest Month). It see,s that the wild population require cold temperature especially in winter, while the admixed group seems it reduces this need. From this early result I might suspect that the introgression of cultivated material can increase adaptation in the future climatic scenario where temperature level will rise.


![RDA_geo_env_lat_range](https://github.com/user-attachments/assets/5296d3ed-e6dc-4880-b246-6a0a3cd16fc0)

