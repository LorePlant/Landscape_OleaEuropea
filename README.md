
## Landscape Olea europea

This page is created to track progresses on my postdoctoral research in modelling genomic offset in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.

## Data input
I started by entering the vcf into R using the vcfR package as follow
```
library(vcfR)
library(adegenet)

setwd("/lustre/rocchettil")
genoLAND.VCF <- read.vcfR("202_Olive_west_MAF005.vcf.recode.vcf")#import vcf file
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

As first step I'm going to estimate the genetic distance as pairwise FST (WC 1984) among the 27 populations using the R package Hierfstat. To run this, I have used a thinned 250Kb version of the vcf file.
```
library("hierfstat")
pops<-read.table("27_Pops.txt", header=T) #one column table wih pop info for each individual
#convert genInd and pop in hierfstat
hierfstat<-genind2hierfstat(gl.genoLAND ,pop=pops)
m<-genet.dist(hierfstat,diploid=TRUE,method="WC84")
```
We obtained a pairwide distance matrix of 325 pairwise comparisons.
For each FST value we calcolated the transformed value as FST/(1-FST). To do so we converted the dist object in matrix, saved the text to ultimetly calcolate the FST/(1-FST) for each data point. Ther resulted matrix was then entered in R as matrix and trasnformed in dist object to run the subseauent mantel test.
```
p<-as.matrix(m)
write.table(p, 'FSTpop27.txt')
# calcolated FST/(1-FST) for each data point in excel

f<-read.table('FSTtransf.txt')
FSTpp<- as.matrix(f)
distFST<- as.dist(FSTpp)
```



In the next chuck of codes we are going to define the Euclidean distance for geography and enviromental variables
```
#population environmental and geographic variables
pop_env<- read.csv("pop_data.csv", header = TRUE)

#bioclim PCdata frame
Env <- scale(pop_env[,4:17], center=TRUE, scale=TRUE)
dist.PCbio = dist(Env, method = "euclidean")

#geographic data
geo = data.frame(pop_env$long, pop_env$lat)
dist.geo = dist(geo, method = "euclidean")


```
#Mantel test
With the aim if there to underline Isolation by Environment (IBE) I used a Mantel test to see if there is a linear correlation between ecological distance and genetic distance matrices
```
# mantel test Genetic distance-ecological distance
geno_eco = mantel(distFST, dist.PCbio, method="spearman", permutations=1000,  na.rm = TRUE)
geno_eco
summary(lm(dist.PCbio~m))
graph = mantel.correlog(m, dist.PCbio, XY=NULL, n.class=0, break.pts=NULL, 
                        cutoff=TRUE, r.type="pearson", nperm=999, mult="holm", progressive=TRUE)


xx = as.vector(distFST) #convert distance matrix into a vector
zz = as.vector(dist.PCbio)
manatelmatrix = data.frame(zz,xx)
mm = ggplot(manatelmatrix, aes(y = xx, x = zz))+
  geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21,fill = "grey") + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2)+
  labs(y = "FST/(1-FST)", x = "Euclidean ecological distance")+
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18, colour = "black"), 
         axis.title= element_text(face = "bold", size = 18, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =12, face = "bold", colour = "black"),
         legend.text = element_text(size = 10, face = "bold", colour = "black"), 
         legend.position = "top", strip.background = element_rect(fill = "grey90", colour = "black"),
         strip.text = element_text(size = 9, face = "bold"))
jpeg(file = "/lustre/rocchettil/mantel_olive_geno_eco.jpeg", width = 350, height = 350)
plot(mm)
dev.off()

```

![mantel_olive_geno_eco](https://github.com/user-attachments/assets/90200675-9ae9-4ee3-be18-e21248c72ec7)



Results suggests a moderate (r: 0.32) though significant (P<0.01) correlation between genetic euclidean distance and ecological euclidean distances. 
Considering the potential effect of geography in the ecological distance I used partial Mantel test which test for correlation between genetic distance and environmental distance considering geographic distance as covariate


```
#partial Mantel test 
partial_mantel = mantel.partial(distFST, dist.PCbio, dist.geo, method = "spearman", permutations = 1000,
                                na.rm = TRUE)
partial_mantel
summary(lm(m~dist.PCbio|dist.geo))
#plotting partial Mantel test
xx = as.vector(distFST) #convert distance matrix into a vector
yy= as.vector(dist.geo)
zz = as.vector(dist.PCbio)
partial_mantel_matrix = data.frame(xx,zz,yy)#new data frame with vectorize distance matrix

mp = ggplot(partial_mantel_matrix, aes(y = xx, x = zz)) + 
  geom_point(size = 2.5, alpha = 0.75, colour = "black",shape = 21, aes(fill = yy)) + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) + 
  scale_fill_continuous(high = "navy", low = "lightblue")+
  labs(y = "FST/(1-FST)", x = "Eucledian ecological distance", fill= "geographic distance")+
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
The result confirmes the significant (P<0.01) and moderate correlation (r:0.24) between environmental distance and genetic distance. Overall this result suggest a potential Isolation-by-envirnoment effect IBE on the sampled population.

![partial_mantel_olive](https://github.com/user-attachments/assets/32430757-a42f-4b22-b496-ec5b2480ab37)

# indentify recent hybrids vs past introgression
The significant admixture present in our collection between wild and cultivated material can be derived from recent crossing forming F1 hybrids or from past generations of crossing where natural selection had the possibilty to act. The GEA identification presume that the associated QTL derived from processes of local adaptation where environmental selection had the generational time to act. In our case the cultivated genome from cultivars vegetatevly progated mainly evolved in the eastern part of the mediterrenan basin, so it is paramount to identify the presence of recent F1 hybrids where the cultivated genome have been recently introduced and where selections did not have the generational time to act.

To investigate the presence of F1 hybrids I identified a recent devoped Rpackage that allow to identified ancestry-informative markers and estimate their hybrids index with the relative presence of F1, BC1, BC2 or past introgression.  https://omys-omics.github.io/triangulaR/index.html

From the PCA analysis I've diveded the germplasm in Wild and Cultivated using PC1 values <-20 and >40 respectively. The intermidiate genotypes represent admixed genotypes for which hybrid index will be checked

![image](https://github.com/user-attachments/assets/3e2f6668-8e47-4df7-9652-07d6fc575484)

```
library(triangulaR)
# make a pop map
popmap<-read.table("popmap.txt")

# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
vcfR.diff <- alleleFreqDiff(vcfR = genoLAND.VCF, pm = popmap, p1 = " P1", p2 = "P2", difference = 0.7)
#"3042 sites passed allele frequency difference threshold"

# Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
hi.het <- hybridIndex(vcfR = vcfR.diff, pm = popmap, p1 = "P1", p2 = "P2")

# Generate colors (or leave blank to use default)
cols <- c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837")
# View triangle plot
jpeg(file = "/lustre/rocchettil/triangular_plot.jpeg")
triangle.plot(hi.het, colors = cols)
dev.off()
```
![triangular_plot](https://github.com/user-attachments/assets/a65ae610-7812-4c6f-9181-244587150fe2)

In the following figure we can see the  theoretical expectations for combinations of hybrid index and interclass heterozygosity under Hardy-Weinberg Equilibrium (HWE). In Larson et al 2013 F1 1 hybrids (hybrid index = 0.5, interspecificheterozygosity ≥ 85%), multi-generation hybrids (hybridindex 0.25–0.75, interspecific heterozygosity < 85%). 
Following this I used a more stringent selection selecting individuals with ****_interclass heterozygosity_**  < 0.7** and ****_hybrid index_** <  **0.5**** selecting individuals with at least one generation of segregation and selection, obtaining a total of 202 genotypes 62 admixed and 140 Wild.


![image](https://github.com/user-attachments/assets/82d441d1-70e4-432c-a08a-c5dd92ea617d)


# Principal component analysis
Principal Component Analysis (PCA) is the amongst the most common multivariate analyses
used in genetics.
Let's first upload the new genotypic datafile of 202 individuals and applying the MAF 0.05.

```
setwd("/lustre/rocchettil")
genoLAND.VCF <- read.vcfR("202_Olive_west_MAF005.vcf.recode.vcf")#import vcf file
gl.genoLAND <- vcfR2genind(genoLAND.VCF)#transfrom file in genind object
genotype<-as.data.frame(gl.genoLAND)
#genotype<-tibble::rownames_to_column(genotype, "geno") #transform raw name in column

```
Enter the population informations. For simplicity I al going to use a popolation differentiations based on geographic origins: (Morocco, Spain, France, Corse)
```
pop_region<-read.table("pops_geo_regions.txt") #one column table wih pop info for each individual
geneIndpop <- vcfR2genind(genoLAND.VCF, pop=pop_region)#transfrom file in genind object

x.olive <- tab(geneIndpop, freq=TRUE, NA.method="mean")
pca.olive <- dudi.pca(x.olive, center=TRUE, scale=FALSE)
popfac<-as.factor(pop_region$V1)
s.class(pca.olive$li, fac=popfac,col=c("darkorange", "darkgreen", "blue", "red"))
jpeg(file = "/lustre/rocchettil/PCA_202.jpeg")
s.class(pca.olive$li, fac=popfac,col=c("darkorange", "darkgreen", "blue", "red"))
dev.off()

eig.perc <- 100*pca.olive$eig/sum(pca.olive$eig)
head(eig.perc)
```
We can compare the result of PCA 202 with the PCA 359 to see if by remouving the F& hybrids and largely cultivated material we can better explain structure of the Wild
```
genoLAND359.VCF <- read.vcfR("359_Olive_west_MAF005.vcf.recode.vcf")#import vcf file
pop_region_359<-read.table("pops_regions_359.txt") #one column table wih pop info for each individual
geneIndpop359 <- vcfR2genind(genoLAND359.VCF, pop=pop_region_359)#transfrom file in genind object

x.olive_359 <- tab(geneIndpop359, freq=TRUE, NA.method="mean")
pca.olive359 <- dudi.pca(x.olive_359, center=TRUE, scale=FALSE)
popfac359<-as.factor(pop_region_359$V1)
s.class(pca.olive359$li, fac=popfac359,col=c("darkorange", "darkgreen", "blue", "red"))
jpeg(file = "/lustre/rocchettil/PCA_359.jpeg")
s.class(pca.olive359$li, fac=popfac359,col=c("darkorange", "darkgreen", "blue", "red"))
dev.off()

eig.perc <- 100*pca.olive$eig/sum(pca.olive$eig)
head(eig.perc)
```
The result show that by remouving the F1 hybrids and genotypes with high membership from cultivated material we can asses in a better way the populations group among the wild germplasm. In details, in the PCA202 we can appreciate the group differentiation between wild of the south (Morocco) and wild of the north (Corse)


![PCA_359](https://github.com/user-attachments/assets/ca186106-4d65-43e0-bcf4-f59a1b3203f8)
>PCA from the whole population of 359 individuals

![PCA_202](https://github.com/user-attachments/assets/048cd8fa-5cd8-431a-b982-1076a52f124b)
>PCA from the filtered population of 202 individuals

Considering that downstream analysis like RDA do not work with NA values I found the following R for cycle for genetic data imputation. This code can be found in https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd


```
for (i in 1:ncol(genotype))
{
  genotype[which(is.na(genotype[,i])),i] <- median(genotype[-which(is.na(genotype[,i])),i], na.rm=TRUE)
}

write.table(genotype, "geno_202_west_olive_MAF005__imputated.txt")
genotype<- read.table("geno_202_west_olive_MAF005__imputated.txt", header=TRUE)
```

# Redundancy analysis

Within the landscape genomic framework, Redundancy analysis (RDA) represent a useful tool that allows to dissect the the total genetic variance among the environment, geographic and demographic components. 
In this first analysis I used RDA on the following linear model to see if we can detect specif environmental variables diverging Wild vs Admixed genotypes or geographic regions.
 $` Gen \sim Environment `$

let's first upload the new genotypic datafile of 202 individuals and applying the MAF 0.05.

```
setwd("/lustre/rocchettil")
genoLAND.VCF <- read.vcfR("202_Olive_west_MAF005.vcf.recode.vcf")#import vcf file
gl.genoLAND <- vcfR2genind(genoLAND.VCF)#transfrom file in genind object
genotype<-as.data.frame(gl.genoLAND)
#genotype<-tibble::rownames_to_column(genotype, "geno") #transform raw name in column

```
Considering that downstream analysis like RDA do not work with NA values I found the following R for cycle for genetic data imputation. This code can be found in https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd

```
for (i in 1:ncol(genotype))
{
  genotype[which(is.na(genotype[,i])),i] <- median(genotype[-which(is.na(genotype[,i])),i], na.rm=TRUE)
}

write.table(genotype, "geno_202_west_olive_MAF005__imputated.txt")
genotype<- read.table("geno_202_west_olive_MAF005__imputated.txt", header=TRUE)
```
The following chuck of code illustrates the step undertaken for assembly the dataset for RDA. The main step is the standardization of environmental variables
```
#standardize bioclim variable
data202<- read.csv("dataset_202_west.csv", header = TRUE)
bio = data202[ ,17:30]
Env <- scale(bio, center=TRUE, scale=TRUE)
Env <- as.data.frame(Env)

#combining geographic, Popstructure, environmental (scaled) variables
Variables <- data.frame(data202$IDSample, data202$long, data202$lat, data202$group,data202$latitude_range, data202$region, data202$PC1, data202$PC2, data202$PC3,  Env)
names(Variables)[1]<- paste("geno")
names(Variables)[2]<- paste("long")
names(Variables)[3]<- paste("lat")
names(Variables)[4]<- paste("group")
names(Variables)[5]<- paste("latitude_range")
names(Variables)[6]<- paste("region")
names(Variables)[7]<- paste("PC1")
names(Variables)[8]<- paste("PC2")
names(Variables)[9]<- paste("PC3")
 ```
To reduce collinearity, I want to check if the selected environmental variance have low VIF variance inflation factor

```
RDAgeo_env <- rda(genotype ~ bio1+bio2+bio4+ bio5 + bio6+bio8	+ bio9 + bio10+bio11+ bio12 + bio14+	bio15	+ bio18 + bio19, Variables)

sqrt(vif.cca(RDAgeo_env))
```
| bio1    |  bio2    |  bio4   |   bio5   |   bio6    |  bio8   |   bio9 | bio10 |   bio11      |   bio12  |   bio14 | bio15 | bio18 | bio19 |
|---------|----------|---------|----------|-----------|----------|---------|-----------|----------|---------|--------|-------|-----|---------|
|29.004607 | 24.449079  |49.46922| 15.210764 | 29.510669| 3.320524 |30.741673 |22.032603|39.983301| 6.3260372 |5.625297|4.254793|3.822378 | 5.766060|


I attempted to select significant ecological variables with a Variance Inflation Factor (VIF) lower than 10. Considering the termopluviometric graph for the study area, which includes France, Spain, and Morocco, we can observe a bell-shaped temperature distribution in common among all locations. Temperatures are higher in the spring and summer, while lower during the autumn and winter. However, precipitation patterns differ across regions. For example, the internal regions of France experience a more uniform distribution of precipitation throughout the year.

To address specific temperature and precipitation adaptations, I carefully selected bioclimatic variables exluding those reffering to temperature and precipitation to the wettest or driest quarters. I focused instead on variables that reflect temperature and precipitation variations during the warmest and coldest quarters, allowing to precisly asses temperature and precipitation during winter and summer quarters. Following this selection, I checked the Variance Inflation Factor (VIF) for the chosen variables to ensure they were appropriate.


I selected bio2, bio6, bio8, bio12, bio14, bio15 and bio 19 and run again the VIF analysis

```
RDAgeo_env <- rda(genotype ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19, Variables)

sqrt(vif.cca(RDAgeo_env))
```

|  bio2    |   bio10   |   bio11  |  bio15  |   bio18  |   bio19  | 
|---------|----------|---------|----------|-----------|----------|
1.939547| 1.368211 |2.865102 |3.484063 |1.964664| 1.438614|


The variable selected showed low inflation factor due to correlation


In the next part I'm preparing the data for plotting the RDA biplot, where genotypes will be presented by dots and loadings will be the environmental variable
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
  xlab("RDA 1: 46 %") + ylab("RDA 2: 17 %") +
  guides(color=guide_legend(title="Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_env.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```
![RDA_env](https://github.com/user-attachments/assets/5e602783-0122-4339-833c-00e45be20fd0)



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
  xlab("RDA 1: 46%") + ylab("RDA 2: 17 %") +
  guides(color=guide_legend(title="latitude range")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_geo_env_lat_range.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```

![RDA_geo_env_lat_range](https://github.com/user-attachments/assets/d040ec52-5a43-4cf6-801c-5a8294f92dc8)



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
  xlab("RDA 1: 46%") + ylab("RDA 2: 17 %") +
  guides(color=guide_legend(title="latitude range")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
jpeg(file = "/lustre/rocchettil/RDA_geo_env_region.jpeg")
plot(loading_RDAgeo_env)
dev.off()
```

![RDA_geo_env_region](https://github.com/user-attachments/assets/cf1b91aa-11d0-4055-821d-f84b20e400c1)




The result show a clear differentiation between Wild and Admixed populations. The two groups are mainly divided along the RDA1 component which is positively correlated with temperature variable bio10 and bio11 and precipitation seasonality bio15. Among the wild group the group of Corse is distinguished by the southern Wild for winter precipitation bio19.
The result suggest that the wild populations can trive in warmer winters, and drier summers, compared to the admixed group. I would speculate from this outcome that the introgression of cultivated genepool can decrease the potential adaptation in future environmental scenarios were temperature levels are forecast to increase.

## RDA for variance partitioning

```
RDAwhole_model <- rda(genotype ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19+  PC1 + PC2 + PC3 + lat + long, Variables)
RsquareAdj(RDAwhole_model)
anova(RDAwhole_model)
## Pure climate model
pRDAclim <- rda(genotype ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19+ Condition(PC1 + PC2 + PC3 + lat + long), Variables)
RsquareAdj(pRDAclim)
anova.cca(pRDAclim)
## Pure neutral population structure model  
pRDAstruct <- rda(genotype ~ PC1 + PC2 + PC3 + Condition(long + lat + bio2+bio10+bio11+	bio15	+ bio18 + bio19), Variables)
RsquareAdj(pRDAstruct)
#anova(pRDAstruct)
##Pure geography model
pRDAgeog <- rda(genotype ~ long + lat + Condition(PC1 + PC2 + PC3 +bio2+bio10+bio11+	bio15	+ bio18 + bio19), Variables)
RsquareAdj(pRDAgeog)
#anova(pRDAgeog)
```

|Partial RDA models |  variance | ADJ R2 | P(<F) | Proportion of explainable variance | Proportion of total variance |
|-------------------------------|--------|----------|--------|--------|----------------------------------------------|
| Full model Y = G+E+Geo+Struct|   |    0.107        |       |         |                                              |
| climate Y = G + E:( Geo + Struct)|  |    0.0139       |       |         |                                              |
| geo Y = G + Geo:(E + Struct)| |        0.006    |      |         |                                              |
| Struct Y = G + Struct:(E + geo)| |  0.041        |       |         |                                              |
| Total unexplained|
| Total variance | |           |       |         |                                              |

## RDA for Genotype Environment Associations (GEA)

Redundancy analysis can be used to identify GEA based on the Mhallanoise distance of SNPs in the RDA-biplot. Within the RDA model we can effectively correct for population structure (PC1 + PC2 + PC3) and geography (latitude and longitude) using them as covariates in the RDA model
As first attempt I decided to run the anlysis seperate for temperature and precipitation variables.

>Temperature

```
RDA_temp <- rda(genotype ~ bio2+bio10+bio11 +  Condition(PC1 + PC2 + PC3 + lat + long), Variables)
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

rdadapt_temp<- rdadapt(RDA_temp, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_temp$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_temp$p.values<thres_env)], p.value = rdadapt_temp$p.values[which(rdadapt_temp$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genotype)[which(rdadapt_temp$p.values<thres_env)], split = "_"), function(x) x[1])))
write.table(top_outliers, "Bonferroni_temp")
qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_temp$p.values, q.value = rdadapt_temp$q.value)
outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_temp$q.values<0.05)], p.value = rdadapt_temp$p.values[which(rdadapt_temp$q.values<0.05)])

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
  xlab("RDA 1: 40%") + ylab("RDA 2: 31%") +
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
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.000414996314981081), genomewideline = -log10(4.45474e-07))
jpeg(file = "/lustre/rocchettil/Manh_RDA_temp.jpeg")
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.000414996314981081), genomewideline = -log10(4.45474e-07))
dev.off()

#P distribution
jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_temp")
hist(Manhattan_temp$P)
dev.off()
```
![RDA_temp_biplot](https://github.com/user-attachments/assets/be1e7698-04e6-4b6e-bac8-4799bb3b6582)
![Manh_RDA_temp](https://github.com/user-attachments/assets/d4baf3a4-187a-4fac-9c46-7e0a6cb34ec5)
![Phist_Manh_RDA_temp](https://github.com/user-attachments/assets/9e4ede77-6459-4482-8a47-a2361a5f7b6a)




> Precipitation
```
RDA_prec <- rda(genotype ~ 	bio15	+ bio18 + bio19 +  Condition(PC1 + PC2 + PC3 + lat + long), Variables)
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

rdadapt_prec<- rdadapt(RDA_prec, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_prec$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_prec$p.values<thres_env)], p.value = rdadapt_prec$p.values[which(rdadapt_prec$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genotype) [which(rdadapt_prec$p.values<thres_env)], split = "_"), function(x) x[1])))
write.table(top_outliers, "Bonferroni_prec")
qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_prec$p.values, q.value = rdadapt_prec$q.value)
outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_prec$q.values<0.05)], p.value = rdadapt_prec$p.values[which(rdadapt_prec$q.values<0.05)])

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
  xlab("RDA 1: 36%") + ylab("RDA 2: 34%") +
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
manhattan(Manhattan_prec, col = c("blue", "gray60"),suggestiveline = -log10(0.000339951095413677), genomewideline = -log10(4.45474e-07))
dev.off()

#P distribution
jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
hist(Manhattan_prec$P)
dev.off()

```
![RDA_prec_biplot](https://github.com/user-attachments/assets/34b802a1-956c-4e34-b410-64ac52d57f28)
![Manh_RDA_prec](https://github.com/user-attachments/assets/c7b9fe24-edc7-42f3-a9b7-8af09640c2d8)
![Phist_Manh_RDA_prec](https://github.com/user-attachments/assets/bed97f85-2103-4e53-aa41-15a3f52f0eab)


>All together

Not Updated!!!!!!

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


## Enriched RDA

To visualize the adaptive differentiation among genotypes, I conducted an additional Redundancy Analysis (RDA) using only the 755 previously identified GEA SNPs for the two seperate analysis for temperature and precipitation (FDR, q<0.05). In this analysis, I did not include geography and population structure as covariates for two reasons: First, I aimed to observe the differentiation between wild and admixed genotypes. Second, the GEA SNPs used have already been identified with corrections for population structure and geography.


```
#partial redundancy analysis (RDA only with GEA QTL)
geno_all_enrich<-genotype[which((rdadapt_temp$q.values<0.05)|(rdadapt_prec$q.values<0.05))]
write.table(geno_all_enrich, "geno_202_GEA_QTLs.txt") #save the new GEA genotype data

#once the dataset was created I can enter the data table with read.table
geno_all_enrich<- read.table("geno_202_GEA_QTLs.txt", header=TRUE)

RDA_all_enriched<-rda(geno_all_enrich ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19, Variables)
summary(eigenvals(RDA_all_enriched, model = "constrained"))
```
As confront I'm going to run the same analysis using the most stringent threshold of Bonferroni

```
geno_all_enrich<-genotype[which((rdadapt_temp$p.values<thres_env)|(rdadapt_prec$p.values<thres_env))]
RDA_all_enriched<-rda(geno_all_enrich ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19, Variables)
summary(eigenvals(RDA_all_enriched, model = "constrained"))
```

>ADM vs WLD using FDR

```
#plot genotypes

TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites"))

Geno <- merge(TAB_gen, Variables[, 1:7] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, colour = group), size = 2.5) +
  scale_color_manual(values = c("blue", "darkorange")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*5, yend=RDA2*5, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=5*RDA1, y=5*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1: 30%") + ylab("RDA 2: 23%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_all_enriched
jpeg(file = "/lustre/rocchettil/RDA_all_geno_biplot_WLD_adm.jpeg")
plot(loading_geno_all_enriched)
dev.off()
write.table(TAB_gen, "geno_all_adaptive_values.txt")
```
> Geographic regions

```
#plot genotypes

TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites"))

Geno <- merge(TAB_gen, Variables[, 1:7] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched_region<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, colour = region), size = 2.5) +
  scale_color_manual(values = c("darkgreen","purple", "darkorange", "blue")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*5, yend=RDA2*5, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=5*RDA1, y=5*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1: 30%") + ylab("RDA 2: 23%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_all_enriched_region
jpeg(file = "/lustre/rocchettil/RDA_all_geno_biplot_region.jpeg")
plot(loading_geno_all_enriched_region)
dev.off()
```

![RDA_all_geno_biplot_WLD_adm](https://github.com/user-attachments/assets/5861860c-c2aa-4d54-b52b-6c43662f91e6)

![RDA_all_geno_biplot_region](https://github.com/user-attachments/assets/ce7d409e-7e60-4180-8560-c08192654459)

>ADM vs WLD using Bonferroni

```
#plot genotypes

TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites"))

Geno <- merge(TAB_gen, Variables[, 1:7] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, colour = group), size = 2.5) +
  scale_color_manual(values = c("blue", "darkorange")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*5, yend=RDA2*5, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=5*RDA1, y=5*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1: 51%") + ylab("RDA 2: 20%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_all_enriched
jpeg(file = "/lustre/rocchettil/RDA_all_geno_biplot_WLD_adm_Bonferroni.jpeg")
plot(loading_geno_all_enriched)
dev.off()
```
![RDA_all_geno_biplot_WLD_adm_Bonferroni](https://github.com/user-attachments/assets/b11db27b-e089-467c-8156-f3017e735cfc)


The enriched RDA based on GEA-QTLs shows different environmental adaptation of the WILD pop. Specfically, we can distringuish between North Morocco (Mediterrnenan) similar to south Spain from south atlantic Morocco and Corse.
By doing the same RDA using only the most significant GEA QTLs (Bonferroni threshold) we can not see any differentiation across genetic group or geographic areas. This shows that by lowering the GEA treshold we inevitably include GEA associations that covary with demography and geography.
The RDA biplot shows that wild and adm have different adaptive landscapes.
The GEA analysis supposes that this two groups are adapted to the condition where they germinated. COnsidering that the introgressions pressure are different among geographic regions depending on the agriculture olive system, we can not say anything about (mal)adaptive introgression. 
I belive the only way to hypothesize this is by doing the GEA analysis only on admixed population and check if the GEA QTLs fall into specific wild genomic windows.


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

> Recovering scaling coefficients
Create a table of scaled temperature variable
```

library('dplyr')
PCbio = data202[ ,17:30]
vars = PCbio %>% select('bio2','bio10', 'bio11', 'bio15', 'bio18', 'bio19')
Var_scale <- scale(vars, center=TRUE, scale=TRUE)
scale_var <- attr(Var_scale, 'scaled:scale')
center_var <- attr(Var_scale, 'scaled:center')

```
Enter the raster file for the specific bioclimatic variable for current climatic situation
```
# ras temperature
library(raster)
library("readxl")


bio2<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_ENM_clipped_biova/bio2_current_masked.tif"))
bio10<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_ENM_clipped_biova/bio10_current_masked.tif"))
bio11<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_ENM_clipped_biova/bio11_current_masked.tif"))
bio15<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_ENM_clipped_biova/bio15_current_masked.tif"))
bio18<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_ENM_clipped_biova/bio18_current_masked.tif"))
bio19<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_ENM_clipped_biova/bio19_current_masked.tif"))
names(bio2) = 'bio2'
names(bio10) = 'bio10'
names(bio11) = 'bio11'
names(bio15) = 'bio15'
names(bio18) = 'bio18'
names(bio19) = 'bio19'
#stack the different raster file
ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19))
```
Predict tha adaptive index for each pixel grid

```
## Function to predict the adaptive index across the landscape
source("./src/adaptive_index.R")

res_RDA_all_proj_current <- adaptive_index(RDA = RDA_all_enriched, K = 2, env_pres = ras_current_var, range = range, method = "loadings", scale_env = scale_var, center_env = center_var)
projection<- stack(c(res_RDA_all_proj_current$RDA1, res_RDA_all_proj_current$RDA2))
plot(projection)
writeRaster(projection,'projection_RDA_current.tif',options=c('TFW=YES'))#save raster for QGIS


## Vectorization of the climatic rasters for ggplot
RDA_proj <- list(res_RDA_all_proj_current$RDA1, res_RDA_all_proj_current$RDA2)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}

## Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
distrib_all<- ggplot(data = TAB_RDA) + 
  geom_tile(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  #coord_sf(xlim = c(-148, -98), ylim = c(35, 64), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw(base_size = 7, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=7))
jpeg(file = "/lustre/rocchettil/adaptive_current.jpeg", height=1500, width=3000, res=600)
plot(distrib_all)
dev.off()
```

![RDA_all_geno_biplot_WLD_adm](https://github.com/user-attachments/assets/5861860c-c2aa-4d54-b52b-6c43662f91e6)
![image](https://github.com/user-attachments/assets/43f8f775-d8c1-4abe-ab80-9e0d7789d8e3)




# Local Genomic offset RDA based

Extract future bioclimatic raster file from CHELSA database. In this first step I used predictions 2071-2100 IPSL ssp585.
```
#future temp scenario

library(raster)
library("readxl")


bio2<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/future_clim_2071_2100/IPSL(france)/IPSLssp585/bio2IPSL_2100_masked_enm.tif"))
bio10<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/future_clim_2071_2100/IPSL(france)/IPSLssp585/bio10IPSL_2100_masked_enm.tif"))
bio11<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/future_clim_2071_2100/IPSL(france)/IPSLssp585/bio11IPSL_2100_masked_enm.tif"))
bio15<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/future_clim_2071_2100/IPSL(france)/IPSLssp585/bio15IPSL_2100_masked_enm.tif"))
bio18<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/future_clim_2071_2100/IPSL(france)/IPSLssp585/bio18IPSL_2100_masked_enm.tif"))
bio19<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/future_clim_2071_2100/IPSL(france)/IPSLssp585/bio19IPSL_2100_masked_enm.tif"))
names(bio2) = 'bio2'
names(bio10) = 'bio10'
names(bio11) = 'bio11'
names(bio15) = 'bio15'
names(bio18) = 'bio18'
names(bio19) = 'bio19'
#stack the different raster file
ras_2100_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19))
ras_2100_var_genotypes<- data.frame(data202$IDSample, data202$group, extract(ras_2100_var, data202[,9:10]))
write.csv(ras_2100_var_genotypes, "ras_2100_var_genotypes.csv")

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
local_offeset_proj_2100_ssp585 <- genomic_offset(RDA = RDA_all_enriched, K = 2, env_pres = ras_current_var, env_fut = ras_2100_var, range = range, method = "loadings", scale_env = scale_var, center_env = center_var)

#local_offeset_proj_2100_ssp585$Proj_offset_global is the raster file for overall GO

plot(local_offeset_proj_2100_ssp585$Proj_offset_global)

writeRaster(local_offeset_proj_2100_ssp585$Proj_offset_global,'Local_GO_2100_ssp585.tif',options=c('TFW=YES'))#save raster for QGIS
```
I prefere to plot the raster using QGIS.

![image](https://github.com/user-attachments/assets/a9c385f8-271b-4820-b6f7-0a5074cbe64c)


The result show a lower GO at low latitude compared to higher latitude level, suggesting that the current adaptive value of souther genotypes will allow them to continue to grow in future climatic condition.
The Local Genomic offsets reflects the adaptive genomic landscape where the southern part of Spain and the costal area of Morocco already presented adaptive GEA for high temperature and low precipitation. In the 2100future scenarios where temperature are going to rise following the latitude gradient will threat the norther part of the olive niche where GEA for higher temperature are not yet present.


>z score

Tentative to normalize the raster file using R ongoingthe zscore

```
library(spatialEco)
raster.Zscore(x, p.value = FALSE, file.name = NULL, ...)
```
In the meanwhile I have done a zscore calcolation using raster calculation in QGIS 

![image](https://github.com/user-attachments/assets/02e6484d-449b-4eef-997f-ab2a6704bfb8)

The results highlight the area of Occitanie and central Spain with potential high genomic offsets in 2100 climatic scenario.








The raster obtained has a GO value for each pixel. In theory I can extract genotype values using latitude and longitude info to ultimetly run an ANOVA between the WLD and ADM group.

```
## Extracting environmental values for each source population
GO_2100_ssp585_genotypes<- data.frame(data359$IDSample, data359$group, extract(res_RDA_all_proj_2100_ssp585$Proj_offset_global, data359[,9:10]))
write.table(GO_2100_ssp585_genotypes, "GO_2100_ssp585_genotypes;txt")
write.csv(GO_2100_ssp585_genotypes, "GO_2100_ssp585_genotypes.csv")

#One-way anova with boxplots. Mainly used to plot the environmental effect for a specific trait

GO_2100_ssp585_geno<- read.csv("GO_2100_ssp585_genotypes.csv", header = TRUE)
names(GO_2100_ssp585_geno)[2]<- paste("geno")
names(GO_2100_ssp585_geno)[3]<- paste("group")
names(GO_2100_ssp585_geno)[4]<- paste("GO_2100_ssp585")

D<- ggplot(GO_2100_ssp585_geno, aes(x=group, y=GO_2100_ssp585, fill=group)) + 
  geom_boxplot()+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(hjust=1,face="bold", size=14, angle=45))+ stat_compare_means(method = "anova")+
scale_fill_manual(values=c( "blue", "darkorange"))+
  labs(y = "GO_2100_ssp585_geno")
jpeg(file = "/lustre/rocchettil/GO_2100_ssp585_WLD_ADM.jpeg")
plot(D)
dev.off()
```
![GO_2100_ssp585_WLD_ADM](https://github.com/user-attachments/assets/98a8d41d-2763-42ef-8934-6266c5ee4910)

tentative of correction for specific genotype climatic distance

```
#wild vs ADM
GO_2100_ssp585_geno<- read.csv("GO_2100_ssp585_genotypes.csv", header = TRUE)
model <- lm(GO_2100_ssp585 ~ group + cov, data = GO_2100_ssp585_geno)
summary(model)
value<-lsmeans(model,~ group|cov)

#latitude range
GO_2100_ssp585_geno<- read.csv("GO_2100_ssp585_genotypes.csv", header = TRUE)
model <- lm(GO_2100_ssp585 ~ latitude_range + cov, data = GO_2100_ssp585_geno)
summary(model)
value<-lsmeans(model,~ latitude_range|cov)
```
The script applied the follwing model 
$Y = group + clim.dist + e$

where _Y_ is the estimated genomic offsets, _group_ the ADM and WLD differentiation and _clim.dist_ the climatic distance between current and future scenario for each genotype.


_clim.dist_ is the value calcolated for each genotype that averages all the coefficients of variation calcolated  for the same clim variable between future and current clim scenario


 |            | Estimate Std |  Error | t value | P| 
 |-------------|--------------|--------|---------|---|
 | Intercept |  2.28728   |  0.15488 |  14.768 |  < 2e-16 ***| 
 | groupWild   | -0.38014  |   0.04019  | -9.459 |  < 2e-16 ***| 
 | cov         | -7.42711  |   1.43705 |  -5.168 | 3.95e-07 ***| 
 



 | group |  lsmean  |   SE | df |lower.CL| upper.CL|
 |--------|---------|-------|---|---------|----------|
 |  Admixed |  1.51 |0.0245| 356 |    1.46 |    1.55|
 |  Wild    |  1.13 |0.0317 |356  |   1.06 |    1.19|

```
df<- as.data.frame(value)
p<- ggplot(df, aes(x=group, y=lsmean, color = group)) + 
  geom_line() +
  geom_point()+
scale_color_manual(values=c( "blue", "darkorange"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                 position=position_dodge(0.05))+
 theme_bw(base_size = 14)
jpeg(file = "/lustre/rocchettil/GO_2100_ssp585_adjusted_WLD_ADM.jpeg")
plot(p)
dev.off()


df<- as.data.frame(value)
l<- ggplot(df, aes(x=latitude_range, y=lsmean, color = latitude_range)) + 
  geom_line() +
  geom_point()+
scale_color_manual(values = c("darkgreen","red", "darkorange"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                 position=position_dodge(0.05))+
 theme_bw(base_size = 14)
jpeg(file = "/lustre/rocchettil/GO_2100_ssp585_adjusted_Lat_range.jpeg")
plot(l)
dev.off()
 
```
![GO_2100_ssp585_adjusted_WLD_ADM](https://github.com/user-attachments/assets/9946f745-7f66-4966-99b7-f5e3164280be)
![GO_2100_ssp585_adjusted_Lat_range](https://github.com/user-attachments/assets/899f67ec-be2c-4eef-b0f8-1d656a2b9287)

Even applying a correction using genotype climatic distance as covariate we can still see a significant difference between Wild and Admixed.



# RDA on candidate introgression zones
I identified populations where both wild and admixed occure together. This subset of population can be used to investigate the potetial presence of GEA related to admixture event.

Let's first upload the new genotypic datafile of 85 individuals and applying the MAF 0.05.

```
setwd("/lustre/rocchettil")
  genoINTRO.VCF <- read.vcfR("introgressed_region_Olive_west_MAF005.vcf.recode.vcf")#import vcf file
gl.genointro <- vcfR2genind(genoINTRO.VCF)#transfrom file in genind object
genotype_intro<-as.data.frame(gl.genointro)
write.table(genotype_intro, "geno_86_west_olive_introgressed_MAF005__imputated.txt")

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
