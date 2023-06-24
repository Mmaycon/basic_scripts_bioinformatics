# Find co-correlated DNA methylation probes 

In this tutorial, we gonna find features (DNAmet probes) that might catch some phenotype on Glioblastoma Stem Cells (GSCs).
We will use the WGCNA package to compute the weighted correlation within GSCs.

## PC setup 
asdsadsadasdsadasdasdasdasd

## Preparing packages

```r
library(WGCNA)
library(meffil)
library(impute)

```

## Load GSC DNAmet matrix - I collected them from the GEO database (GSE109330 and GSE119774)

#Obs: You should process those IDAT files, get beta-values, and then transform beta-values into M-values in order to get the DNAmet matrix below

```r
DNAmatrix = readRDS('./loads/67_GSC_Mvalue_DNAmet_matrix_GSE109330_GSE119774.rds')
dim(DNAmatrix) #452567     67
```

## Filter only Most Variable Features (features = DNAmet probes)

```r
require(meffil)
MostVF_GSC <- meffil.most.variable.cpgs(
  as.matrix(DNAmatrix),
  n = 10000,
  sites = NULL,
  samples = NULL,
  autosomal = T,
  winsorize.pct = NA,
  outlier.iqr.factor = NA
)

length(MostVF_GSC) #10000

GSC_matrix_MVF <- DNAmatrix[MostVF_GSC,]
```

## Fix NAs 

Obs: for the sake of practice, we're gonna use only one of these two 'NA removal' strategies below (na.omit)
### Remove probes which got an NA (no M-value) at least for one sample

```r
GSC_matrix_MVF_naomit <- na.omit(GSC_matrix_MVF)
dim(MostVF_GSC_naomit) #9396   67
```

### Replace NA through KNN (K-nearest neighbors)

```r
library(impute)
GSC_matrix_MVF__knn <- impute.knn(as.matrix(GSC_matrix_MVF), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]
```

## Identify outliers 

### Probes outliers
```r
gsg <- goodSamplesGenes(t(GSC_matrix_MVF_naomit))
summary(gsg)
gsg$allOK # no outliers present
```

If gsg$allOK returns FALSE, do the following 

```r
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(GSC_matrix_MVF_naomit)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(GSC_matrix_MVF_naomit)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  GSC_matrix_MVF_naomit <- GSC_matrix_MVF_naomit[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}
```

### Samples outliers

```r
sampleTree <- hclust(dist(t(GSC_matrix_MVF_naomit)), method = "average") 
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
rect.hclust(tree = sampleTree, k = 3, which = 1:3, border = 1:3)
```

[IMAGE] [IMAGE] [IMAGE] [IMAGE] [IMAGE] [IMAGE] [IMAGE] [IMAGE]

If there is any sample outleir 

```r
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10) 
GSC_matrix_MVF_naomit <- GSC_matrix_MVF_naomit[,cut.sampleTree==1] # Remove outlier
```


## Network Construction (Pairwise Gene Co-expression similarity)

### Soft-thresholding
```r
spt <- pickSoftThreshold(t(GSC_matrix_MVF_naomit)) 
```

### R^2 values as a function of the soft thresholds (power)
```r
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")
```

### mean connectivity as a function of soft thresholds (power)
```r
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")
```

#Obs: We should maximize the R2 value and minimize mean connectivity
#Choose a Power higher then a R^2 of 0.8 & the first lowest mean connectivity
#In this case a Power = 5 is within a R^2 of 0.80 & also is the first lowest mean connectivity


### calling the Adjacency Function

```r
softPower <- 5
adjacency <- adjacency(t(GSC_matrix_MVF_naomit), power = softPower)
```


## Module Construction










