
# Landscape Olea europea

This page is created to track progresses on my postdoctoral research in landscape genomics in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.

## Data input
I started by entering the vcf into R using the vcfR package as follow
```

#### loading genotypic file and filtering individuals and MAF
setwd("C:/Users/rocchetti/Desktop/running RDA GO")    #### cluster  setwd("/lustre/rocchettil")

geno359.VCF <- read.vcfR("359_Olive_west_MAF005.vcf.recode.vcf")#import vcf file
gl.genoLAND <- vcfR2genind(geno359.VCF)#transfrom file in genind object
geno359<-as.data.frame(gl.genoLAND)
geno359<-geno359%>% select(ends_with(".0"))
for (i in 1:ncol(geno359))
{
  geno359[which(is.na(geno359[,i])),i] <- median(geno359[-which(is.na(geno359[,i])),i], na.rm=TRUE)
}
```
## Diversity analysis and evaluation of introgression

From the work of Lison 2024 135 truly wild genotypes were selected with ancestry _q>0.70_. The rest of 224 genotype were classified as admixed. To distinguish between historical vs recent introgression, we analyzed the hybrid index using ancestry-informative SNPs. To define the two parental hybrid sources, we conducted a PCA. 61 genotypes with the most extreame postion along the PC1 from the truly wild were classified as cultivated. These individulas most likely represent cultivated seeds that migrated in the wild environment from the near agricultural field. 

>PCA on the whole dataset

```
res.pca359<-PCA(geno359, scale.unit = TRUE, ncp = 5, graph = TRUE)
fviz_eig(res.pca359, addlabels = TRUE, ylim = c(0, 50))#  scree plot
# Create a data frame for PCA results
ind359 <- get_pca_ind(res.pca359)
ind359
pca_data359 <- as.data.frame(ind359$coord)
pca_data359$group <- "Admixed"
pca_data359$group[rownames(pca_data359)%in%rownames(geno_Wild_GEA)] <- "wild"
pca_data359$group[rownames(pca_data359)%in%rownames(geno_Cul_GEA)] <- "cultivated"

qq<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = pca_data359, aes(x=Dim.1, y=Dim.2, colour = group), size = 2.5) +
  scale_color_manual(values = c("grey", "darkgreen","purple")) +
  xlab("PC1: 16%") + ylab("PC2: 2%") +
  guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
qq

```
![image](https://github.com/user-attachments/assets/479d6c01-d30e-4073-9c2b-a615d72589c2)


> Indentifing recent hybrids vs past introgression

The significant admixture present in our collection between wild and cultivated material can be derived from recent crossing forming F1 hybrids or from past generations of crossing where natural selection had the possibilty to act. The GEA identification presume that the associated QTL derived from processes of local adaptation where environmental selection had the generational time to act. In our case the cultivated genome from cultivars vegetatevly progated mainly evolved in the eastern part of the mediterrenan basin, so it is paramount to identify the presence of recent F1 hybrids where the cultivated genome have been recently introduced and where selections did not have the generational time to act.

To investigate the presence of F1 hybrids I identified a recent devoped Rpackage that allow to identified ancestry-informative markers and estimate their hybrids index with the relative presence of F1, BC1, BC2 or past introgression.  https://omys-omics.github.io/triangulaR/index.html

In this analysis I used the vcf file that was not filtered for MAF 5%.


```
library(triangulaR)
# make a pop map
popmap<-read.table("popmap.txt")

genoLAND.VCF <- read.vcfR("359_Olive_west.vcf.recode.vcf")#import vcf file

# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
vcfR.diff <- alleleFreqDiff(vcfR = genoLAND.VCF, pm = popmap, p1 = " P1", p2 = "P2", difference = 0.7)
#"3175 sites passed allele frequency difference threshold"

# Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
hi.het <- hybridIndex(vcfR = vcfR.diff, pm = popmap, p1 = "P1", p2 = "P2")

# Generate colors (or leave blank to use default)
cols <- c("darkgrey", "purple", "darkgreen")
# View triangle plot
jpeg(file = "/lustre/rocchettil/triangular_plot.jpeg")
triangle.plot(hi.het, colors = cols)
dev.off()
```
![image](https://github.com/user-attachments/assets/ab6673f7-6bab-4c90-bf3a-fe592ddf2c52)


In the following figure we can see the  theoretical expectations for combinations of hybrid index and interclass heterozygosity under Hardy-Weinberg Equilibrium (HWE). In Larson et al 2013 F1 1 hybrids (hybrid index = 0.5, interspecificheterozygosity ≥ 85%), multi-generation hybrids (hybridindex 0.25–0.75, interspecific heterozygosity < 85%). 
The results highlight the large presence of rencet hybrids like F1 and BC1. From this results we confidentially used the 135 wild previously identified by Lison for the Genotype Environment Association GEA analysis.


![image](https://github.com/user-attachments/assets/82d441d1-70e4-432c-a08a-c5dd92ea617d)



>filtering wild individual and filtering sites for MAF 0.05

```
##### Filter Wild individual 
wild_list<-read.table("Wild_list.txt", header = F)
genoWild <-  geno359[rownames(geno359)%in% wild_list$V1, ]

### Filter wild dataset for MAF 005
freq_mean <- colMeans(genoWild)
genoWild_MAF005 <- genoWild[,-which(freq_mean>=0.95 | freq_mean<=0.05)]

write.table(genoWild_MAF005, "genoWild_MAF005_imputated.txt")
genoWild_MAF005<- read.table("genoWild_MAF005_imputated.txt", header = T)
```
>prepare datafile with standardized environmental variables

```
# Wild Environment datafile

#standardize bioclim variable
data_wild<- read.csv("WILD_135.csv", header = TRUE)
test_env <- data_wild%>% select(long, lat, bio2, bio10, bio11, bio15, bio18, bio19)
Env <- scale(test_env, center=TRUE, scale=TRUE)
# Extract the centering values
env_center <- attr(Env, "scaled:center") #mean of each variable
# Extract the scaling values
env_scale <- attr(Env, "scaled:scale") #standard deviation of each variable
#transform into dataset
Env <- as.data.frame(Env)


#combining geographic, Popstructure, environmental (scaled) variables
Variables <- data.frame(data_wild$IDSample, data_wild$group,data_wild$latitude_range, Env)
names(Variables)[1]<- paste("geno")
names(Variables)[2]<- paste("group")
names(Variables)[3]<- paste("latitude_range")
```
## Population structure analysis of 135 Wild
To asses the population structure of the 135 wild we used a Principal component analysis, searching for potential differentiation by countries or latitude gradient.
```
## population structure for 135 WILD
library(FactoMineR)
library(factoextra)

res.pca<-PCA(genoWild_MAF005, scale.unit = TRUE, ncp = 5, graph = TRUE)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))#  scree plot
# Create a data frame for PCA results
library(tibble)
ind <- get_pca_ind(res.pca)
pca_data <- as.data.frame(ind$coord)
eig.val <- get_eigenvalue(res.pca) #selected 5 PCs
# Create a data frame for PCA results
ind135 <- get_pca_ind(res.pca)
ind135
pca_data135 <- as.data.frame(ind135$coord)
pca_data135<- rownames_to_column(pca_data135, var = "IDSample")
PCA_wild<-left_join(pca_data135, data_wild, by ="IDSample")

#regions
qq<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = PCA_wild, aes(x=Dim.1, y=Dim.2, colour = regions), size = 2.5) +
  scale_color_manual(values = c("blue", "darkgreen","purple", "darkorange")) +
  xlab("PC1: 4.2%") + ylab("PC2: 2.5%") +
  guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
qq

#latitude gradient
ee<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = PCA_wild, aes(x=Dim.1, y=Dim.2, colour = lat), size = 2.5) +
  #scale_color_manual(values = c("darkgreen","red", "darkorange")) +
  scale_color_viridis_c(option = "D", name = "Latitude")+
  xlab("PC1: 4.2%") + ylab("PC2: 2.5%") +
  #guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
ee
library(ggpubr)
ggarrange(qq,ee,nrow = 1 , ncol = 2)

library(metan)
dff<-PCA_wild[,c("Dim.1","bio2", "bio10", "bio11","bio15", "bio18", "bio19", "long", "lat")]
a<-corr_plot(dff)
a
```
![image](https://github.com/user-attachments/assets/32292952-f935-4d53-bfec-abeb11299bde)


## Redundancy analysis

Within the landscape genomic framework, Redundancy analysis (RDA) represent a useful tool that allows to dissect the the total genetic variance among the environment, geographic and demographic components. 


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
2.314209| 1.491136 |3.308390 |3.656316 |2.216436| 1.632229|


The variable selected showed low inflation factor due to correlation



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

Redundancy analysis can be used to identify GEA based on the Mhallanoise distance of SNPs in the RDA-biplot. Within the RDA model we can effectively correct for population structure  and geography (latitude and longitude) using them as covariates in the RDA model. As population structure correction we used latent factor derived from the LEA package.

As first attempt I decided to run the anlysis seperate for temperature and precipitation variables.

>Temperature

```
## Use latent factor for covariable correction
# latent factor temperature variable
Y <- genoWild_MAF005
sel_temp<- data.frame(Env%>% dplyr::select(bio2, bio10, bio11))
write.env(sel_temp, "Temp_variable.env")
X = read.table("Temp_variable.env")

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 4)
str(mod.lfmm2)
mod.lfmm2@U
#Merge latent factor to Variable
latent_temp<-data.frame(rownames(genoWild_MAF005), mod.lfmm2@U)
Temp_Var<-cbind(Variables,latent_temp)


#GEA Temperature
RDA_temp <- rda(genoWild_MAF005 ~ bio2+bio10+bio11 +  Condition(X1 + X2 +X3 +X4 + lat + long), Temp_Var)
summary(eigenvals(RDA_temp, model = "constrained"))
library(robust)
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
top_outliers <- data.frame(Loci = colnames(genoWild_MAF005)[which(rdadapt_temp$p.values<thres_env)], p.value = rdadapt_temp$p.values[which(rdadapt_temp$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genoWild_MAF005)[which(rdadapt_temp$p.values<thres_env)], split = "_"), function(x) x[1])))
qvalue <- data.frame(Loci = colnames(genoWild_MAF005), p.value = rdadapt_temp$p.values, q.value = rdadapt_temp$q.value)
outliers <- data.frame(Loci = colnames(genoWild_MAF005)[which(rdadapt_temp$q.values<0.05)], p.value = rdadapt_temp$p.values[which(rdadapt_temp$q.values<0.05)])


#plot GEA temp

locus_scores <- scores(RDA_temp, choices=c(1:2), display="species", scaling="none")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
TAB_var <- as.data.frame(scores(RDA_temp, choices=c(1,2), display="bp"))
loading_temp<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = TAB_loci, aes(x=20*RDA1, y=20*RDA2, colour = type), size = 2.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=1.1*RDA1, yend=1.1*RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.8, family = "Times") +
  xlab("RDA 1: 37%") + ylab("RDA 2: 33%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_temp
jpeg(file = "/lustre/rocchettil/RDA_temp_biplot.jpeg")
plot(loading_temp)
dev.off()

write.table(qvalue, "Temp_GEA_Olive.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

Manhattan_temp <- read.csv(file = "Temp_GEA_Olive.csv", header=TRUE) #import the p value result for temperature
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.000364239), genomewideline = -log10(3.907776e-06))
jpeg(file = "/lustre/rocchettil/Manh_RDA_temp.jpeg")
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.000364239), genomewideline = -log10(3.907776e-06))
dev.off()

#P distribution
jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_temp")
hist(Manhattan_temp$P)
dev.off()

hist(qvalue$p.value)

```
![image](https://github.com/user-attachments/assets/c17ee38b-e07f-47e5-8526-424b7a1acb93)
![image](https://github.com/user-attachments/assets/7998d98a-10dd-4338-bd5f-686b64c7a380)
![image](https://github.com/user-attachments/assets/c369b9c3-4b8f-4de1-a776-7898e72d7a88)





> Precipitation
```
#latent factor precipitation variable

Y <- genoWild_MAF005
sel_prec<- data.frame(Env%>% select(bio15, bio18, bio19))
write.env(sel_prec, "prec_variable.env")
X = read.table("prec_variable.env")
mod.lfmm2 <- lfmm2(input = Y, env = X, K = 4)
str(mod.lfmm2)
mod.lfmm2@U
#Merge latent factor to Variable
latent_prec<-data.frame(rownames(genoWild_MAF005), mod.lfmm2@U)
Prec_Var<-cbind(Variables,latent_prec)



## GEA Precipitation
RDA_prec <- rda(genoWild_MAF005 ~ 	bio15	+ bio18 + bio19 +  Condition(X1 + X2 + X3 + X4 + lat + long), Prec_Var)
summary(eigenvals(RDA_prec, model = "constrained"))

rdadapt_prec<- rdadapt(RDA_prec, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_prec$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(genoWild_MAF005)[which(rdadapt_prec$p.values<thres_env)], p.value = rdadapt_prec$p.values[which(rdadapt_prec$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genoWild_MAF005)[which(rdadapt_prec$p.values<thres_env)], split = "_"), function(x) x[1])))

qvalue <- data.frame(Loci = colnames(genoWild_MAF005), p.value = rdadapt_prec$p.values, q.value = rdadapt_prec$q.value)
outliers <- data.frame(Loci = colnames(genoWild_MAF005)[which(rdadapt_prec$q.values<0.05)], p.value = rdadapt_prec$p.values[which(rdadapt_prec$q.values<0.05)])

#plot GEA precipitation

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
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 2.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.8, family = "Times") +
  xlab("RDA 1: 40%") + ylab("RDA 2: 34%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_prec
jpeg(file = "/lustre/rocchettil/RDA_prec_biplot.jpeg")
plot(loading_prec)
dev.off()


write.table(qvalue, "Prec_GEA_Olive.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

#plotting Mhanattan plot using the library qqman

library(qqman)
Manhattan_prec <- read.csv(file = "Prec_GEA_Olive.csv", header=TRUE) #import the p value result for precipitation
jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
manhattan(Manhattan_prec, col = c("blue", "gray60"),suggestiveline = -log10(0.00030798), genomewideline = -log10(3.907776e-06))
dev.off()

#P distribution
jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
hist(Manhattan_prec$P)
dev.off()

```
![image](https://github.com/user-attachments/assets/43607b9e-2d80-47d3-b637-bcbf89990e7a)
![image](https://github.com/user-attachments/assets/84320d9d-2814-437a-867c-b976c14b273e)
![image](https://github.com/user-attachments/assets/d7f03d5d-7872-43a9-967d-a0ddbfe90ee9)



## Enriched RDA

To visualize the adaptive differentiation among genotypes, I conducted an additional Redundancy Analysis (RDA) using only the 163 previously identified GEA SNPs for the two seperate analysis for temperature and precipitation (FDR, q<0.05). 

```
#partial redundancy analysis (RDA only with GEA QTL)
geno_Wild_GEA<-genoWild_MAF005[which((rdadapt_temp$q.values<0.05)|(rdadapt_prec$q.values<0.05))]
write.table(geno_Wild_GEA, "geno_Wild_GEA.txt") #save the new GEA genotype data

#scaled variable
data_wild<- read.csv("WILD_135.csv", header = TRUE)
test_env <- data_wild%>% select(bio2, bio10, bio11, bio15, bio18, bio19)
Env <- scale(test_env, center=TRUE, scale=TRUE)
# Extract the centering values
env_center <- attr(Env, "scaled:center")
# Extract the scaling values
env_scale <- attr(Env, "scaled:scale")
#transform into dataset
Env <- as.data.frame(Env)

RDA_all_enriched<-rda(geno_Wild_GEA ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19, Variables)
summary(eigenvals(RDA_all_enriched, model = "constrained"))
```

>geographic region differentiation

```
# plot Geographic regions


TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites"))

Geno <- merge(TAB_gen, Variables[, 1:7] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched_region<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, colour = latitude_range), size = 2.5) +
  scale_color_manual(values = c("darkgreen","darkred", "darkorange")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*5, yend=RDA2*5, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=5*RDA1, y=5*RDA2, label = row.names(TAB_var)), size = 2.8, family = "Times") +
  xlab("RDA 1: 28%") + ylab("RDA 2: 24%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_all_enriched_region
jpeg(file = "/lustre/rocchettil/RDA_all_geno_biplot_region.jpeg")
plot(loading_geno_all_enriched_region)
dev.off()

```
![image](https://github.com/user-attachments/assets/059b8d3b-87ff-442e-8746-c185d38afde6)


# Adaptive index projection
Adaptive indeix function Capblach.
Using the raster environmental data the function allows to predict the adaptie value at each pixel

```
#adaptive index function

adaptive_index <- function(RDA, K, env_pres, range = NULL, method = "loadings", scale_env, center_env){
  
  # Formatting environmental rasters for projection
  var_env_proj_pres <- as.data.frame(rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]]))
  
  # Standardization of the environmental variables
  var_env_proj_RDA <- as.data.frame(scale(var_env_proj_pres[,-c(1,2)], env_center[row.names(RDA$CCA$biplot)], env_scale[row.names(RDA$CCA$biplot)]))
  
  # Prediction with RDA model and linear combinations
  if(method == "predict"){ 
    pred <- predict(RDA, var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], type = "lc", scaling =2)
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

Raster files dowloaded from CHELSA were clipped using Environmental niche modelling masked generated with the package bioclim2. The raster preparation was conducted in QGIS
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

res_RDA_all_proj_current <- adaptive_index(RDA = RDA_all_enriched, K = 2, env_pres = ras_current_var, range = range, method = "predict", scale_env = env_scale, center_env = env_center)
projection<- stack(c(res_RDA_all_proj_current$RDA1, res_RDA_all_proj_current$RDA2))
plot(projection)
  writeRaster(projection,'Wild_adaptive_landscape.tif',options=c('TFW=YES'))#save raster for QGIS

```
![image](https://github.com/user-attachments/assets/b98fd2af-48f6-430e-bfb6-d7e16ca10702)
![image](https://github.com/user-attachments/assets/8b694ec0-5e89-4b83-9d4a-d3e013e1000e)


# estimation of cultivars offsets

The previously identified GEA from the wild dataset represent QTLs involved in adaptation in the western Mediterrenean. Can we use this information to inform cultivar adaptive values? To answer this question we can enter in the RDA model the domesticated genotype info filtered for the wild GEA QTL. The idea is to plot in the RDA space the position of cultivated genotypes based on their GEA information, as weel as, the "best suited" genotypes using the enviromental information from the sample location. The Euclidean distance in the RDA adaptive space between this two point will represent the cultivar offset.

As cultivated genotypes I selected nine different genotype representative of ç different population selecting the ones with higher _hybrid index_ and lower _interclass heterozygosity
_.

>prepare dataset with cultivated genotype and GEA QTLs

```

#filter cultivar individuals
cul_list<-read.table("Dom_9_list.txt", header = F)
genoCul <-  geno359[rownames(geno359)%in% cul_list$V1, ]
geno_Cul_GEA <- genoCul[, colnames(genoCul) %in% colnames(geno_Wild_GEA)]

```
Predict the RDA scores fro the domesticated genotypes
```
fitted_dom_RDAscores <- predict(RDA_all_enriched, newdata=geno_Cul_GEA, type="wa", scaling="sites")



ordiplot(fitted_dom_RDAscores)

TAB_gen135 <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites", choices=c(1:2), scaling = "sites"))
TAB_gen135$type <- "wild"
TAB_gen61 <- data.frame(geno = row.names(fitted_dom_RDAscores), fitted_dom_RDAscores[,1:2])
TAB_gen61$type <- "cultivated"
wild_cult <- rbind(TAB_gen135, TAB_gen61)

TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1:2), display="bp"))

```
Prediction of cultivated best fit
```
data_cul<- read.csv("Dom_9.csv", header = TRUE)
env_c = data_cul[ ,2:15]
enc_c <- env_c%>% select(bio2, bio10, bio11, bio15, bio18, bio19)

#scale environmental variable using the scaling factor from 135 wild
var_env_proj_RDA <- scale(enc_c, env_center[row.names(RDA_all_enriched$CCA$biplot)], env_scale[row.names(RDA_all_enriched$CCA$biplot)])
scaled_env<-as.data.frame(var_env_proj_RDA)

fitted_dom_Envscores <- predict(RDA_all_enriched, newdata=scaled_env, type="lc", scaling = "sites")

plot(fitted_dom_Envscores)

TAB_best_cul<- data.frame(geno = data_cul$IDSample, fitted_dom_Envscores[,1:2])

TAB_best_cul$type <- "best cultivated"

wild_cult_best<-rbind(wild_cult, TAB_best_cul)

```
>plot the result in the RDA space

```

hh<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = wild_cult_best, aes(x=RDA1, y=RDA2, shape=type, color=type, size=type)) +
  scale_shape_manual(values=c(16, 17, 3))+
  scale_color_manual(values=c('#56B4E9', '#E69F00','grey48'))+
  scale_size_manual(values=c(3,3,2))+
  geom_segment(data = TAB_var, aes(xend=RDA1*4, yend=RDA2*4, x=0, y=0), colour="black", size=0.15, linetype=1, arrow = arrow(length=unit(0.20,"cm"),type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x=4*RDA1, y=4*RDA2, label = row.names(TAB_var)), size = 3, family = "Times") +
  xlab("RDA 1: 28%") + ylab("RDA 2: 23%") +
  #guides(legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
hh
```
![image](https://github.com/user-attachments/assets/9f541476-180c-4a28-9764-ba3b311db657)


Calculate the euclidean distance between cultivated and best cultivated

```
dist_data<-merge(TAB_gen61,TAB_best_cul, by = "geno" )
dist_data$offset<-sqrt((dist_data$RDA1.x - dist_data$RDA1.y)^2 + (dist_data$RDA2.x - dist_data$RDA2.y)^2)
hist(dist_data$offset)
colnames(dist_data)[colnames(dist_data) == "geno"] <- "IDSample"
dist_data<-left_join(dist_data, data_359, by ="IDSample")
write.table(dist_data, file = "cultivar_mismatch_scaling1.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

cultivar_offset<- read.csv("cultivar_mismatch_scaling1.csv")

```

Among the cultivated population I selected two contrasting cultivars _F3_ from south France and _E19_ from south Spain. For this two genotypes the Euclidean distance was calcolated between the RDA coordinates from their GEA based predictions and the position of each environmental pixel in the RDA space. The results gives a rapresentation of the adaptive landscape of a specific cultivars base on the GEA QTL from wild genotypes.
The final result was graphically presented using QGIS.

```
#Estimation cultivars adaptive landscape

#load pixel data

pixel<-read.csv("current_env_pixel.csv", sep = " ")
pixel_env<- pixel%>% select(bio2, bio10, bio11, bio15, bio18, bio19)
scaled_pixel <- scale(pixel_env, env_center[row.names(RDA_all_enriched$CCA$biplot)], env_scale[row.names(RDA_all_enriched$CCA$biplot)])
scaled_pixel<-as.data.frame(scaled_pixel)



scaled_pixel_LC <- predict(RDA_all_enriched, newdata=scaled_pixel, type="lc", scaling = "sites")

plot(scaled_pixel_LC)

TAB_pixel_LC<- data.frame(lat = pixel$y, long = pixel$x, scaled_pixel_LC[,1:2])
#distace calculated for the genotype E19 with "wc" prediction of RDA1 = -0.13, RDA2 = 0.15
TAB_pixel_LC$offset<-sqrt((-0.13 - TAB_pixel_LC$RDA1)^2 + (0.15 - TAB_pixel_LC$RDA2)^2)
hist(TAB_pixel_LC$offset)
#The output table was uploaded as vector in QGIS andthen interpolation (IDW method) was used to plot the raster file. The obtained raster was ultilmately standardized using _zscore_.
write.csv(TAB_pixel_LC, "E19_pixel_offset.csv", sep = " ")


#distace calculated for the genotype F3 with "wc" prediction of RDA1 = -0.13, RDA2 = 0.15
TAB_pixel_LC$offset<-sqrt((-0.1676 - TAB_pixel_LC$RDA1)^2 + (-0.1417 - TAB_pixel_LC$RDA2)^2)
hist(TAB_pixel_LC$offset)
#The output table was uploaded as vector in QGIS andthen interpolation (IDW method) was used to plot the raster file. The obtained raster was ultilmately standardized using _zscore_.
write.csv(TAB_pixel_LC, "F3_pixel_offset.csv", sep = " ")
```
![image](https://github.com/user-attachments/assets/54624bb9-dd75-4b9c-aaf5-3a8216585bdd)

