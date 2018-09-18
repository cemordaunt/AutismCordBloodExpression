# MARBLES QC preprocessing ####
# Charles Mordaunt
# 6/7/18

setwd()

# Packages ####
# bioconductor
library(oligo)
library(arrayQualityMetrics)
library(pd.hugene.2.0.st)
library(limma)
library(oligoData)
library(pd.hg.u95av2)
library(sva)

# CRAN
library(Rcpp)
library(RSQLite)
library(DBI)
library(VennDiagram)
library(ggplot2)
library(scales)
library(reshape)

# Read in data ####
sampleInfo <- read.delim(file = "tables/COI_Dx_alg_Pheno_Min.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE) # Only includes cord samples with algorithm diagnosis
directory_B1 <- "raw_files/batch_1"
directory_B2 <- "raw_files/batch_2"
directory_B3 <- "raw_files/batch_3"
childSamples_B1 <- paste(directory_B1, sampleInfo$CEL_FileName[sampleInfo$Batch == 1], sep = "/")
childSamples_B2 <- paste(directory_B2, sampleInfo$CEL_FileName[sampleInfo$Batch == 2], sep = "/")
childSamples_B3 <- paste(directory_B3, sampleInfo$CEL_FileName[sampleInfo$Batch == 3], sep = "/")
childSamples_all <- c(childSamples_B1, childSamples_B2, childSamples_B3)
rawData <- read.celfiles(childSamples_all)
saveRDS(object = rawData, file = "R_objects/child_rawData_B123.rds")

# Normalize with RMA ####
sumData <- rma(rawData, target = "core")
saveRDS(object = sumData, file = "R_objects/child_sumData_B123.rds")

# Add Feature Data ####
# Add features
featureData(sumData) <- getNetAffx(sumData, "transcript")

# Remove control probes
tmp <- featureData(sumData)@data$category
keep <- grep("main", tmp)
sumData2 <- sumData[keep,]

# Take out genes without annotation
tmp <- featureData(sumData2)@data$geneassignment
drop <- which(is.na(tmp))
sumData2 <- sumData2[-drop,]

# Output Normalized Expression Data with Annotations
gene <- sumData2
colnames(gene) <- substr(colnames(gene), 1, 17)
genes <- matrix(unlist(lapply(strsplit(featureData(sumData2)@data$geneassignment, split = "//", fixed = TRUE), function(x) paste(x[1:2]))), ncol = 2, byrow = TRUE)
outdat <- cbind(rownames(gene), genes, exprs(gene))
colnames(outdat)[c(1:3)] <- c("Probe", "GeneID", "GeneName")
write.table(outdat, "tables/Normalized_MARBLES_Child_Expression_Data_B123.txt", quote = FALSE, sep = "\t", row.names = FALSE)
rm(gene, genes, keep, tmp, outdat, sumData2, drop)

# Array Quality Metrics ####
# RMA Normalized data
pData(sumData)$Batch <- c(rep("Batch_1", length(childSamples_B1)), rep("Batch_2", length(childSamples_B2)), rep("Batch_3", length(childSamples_B3)))
arrayQualityMetrics(expressionset=sumData, outdir="Child Report RMA-Normalized with Batch B123",force=TRUE, do.logtransform = FALSE, intgroup = "Batch")

# PCA Plots ####
# Setup pData
pData <- pData(sumData)
pData$CEL_FileName <- rownames(pData)
pData <- merge(pData, sampleInfo, by="CEL_FileName")
pData <- pData[,c("index", "IBC", "CEL_FileName", "Batch.x", "Dx_alg", "COI_GENDER")]
colnames(pData) <- c("index", "IBC", "CEL_FileName", "Batch", "Dx_alg", "COI_GENDER")
pData(sumData) <- pData

# Make plots
RMAeset <- exprs(sumData)
RMAeset.data<-as.data.frame(RMAeset)
order.RMA<-RMAeset.data[,as.character(pData$CEL_FileName)]
prin <- princomp(order.RMA)
pervar<-prin$sdev^2 / sum(prin$sdev^2)
# pervar[1:5]
# 1     2       3       4       5
# 0.973 0.002   0.002   0.001   0.001

loadingsPC1 <- prin$loadings[,"Comp.1"]
thresholdPC1 <- quantile(loadingsPC1, 0.95)
outliersPC1 <- loadingsPC1[loadingsPC1 > thresholdPC1]
names(outliersPC1)
# [1] "IHer-P1B4-B1_M114600_541804-HuGene2.0_(HuGene-2_0-st).CEL"
# [2] "IHer-P1B6-B1_M113173_539103-HuGene2.0_(HuGene-2_0-st).CEL"
# [3] "IHer-P1H5-B1_M110289_503904-HuGene2.0_(HuGene-2_0-st).CEL"
# [4] "IHer-P3B4-B3_M109016_527505-HuGene2.0_(HuGene-2_0-st).CEL"
# [5] "IHer-P5A5-B5_M110790_531404-HuGene2.0_(HuGene-2_0-st).CEL"
# [6] "RSCh-B11-P2_M104827_517106-HuGene2_(HuGene-2_0-st).CEL"
# [7] "RSCh-E2-P1_M107053_502404-HuGene2_(HuGene-2_0-st).CEL"
# [8] "RSCh-F4-P3_M113489_536104-HuGene2_(HuGene-2_0-st).CEL"
# [9] "RSCh-F5-P3_M113555_538704-HuGene2_(HuGene-2_0-st).CEL"
# 9 Outliers: 114600, 113173, 110289, 109016, 110790, 104827, 107053, 113489, 113555

pdf("figures/Child RMA Normalized PCA Plots B123.pdf")
screeplot(prin, col="dodgerblue", xlab="Principal Components of RMA Normalized Beta Values", main="", cex.lab=1.3)

plot.new()
myColors<-c("seagreen3","dodgerblue", "darkorchid", "firebrick1", "darkorange", "khaki1", "azure4", "lightcyan1", "brown", "black")
palette(myColors)
par(mar=c(0, 0, 0, 0))
legend("bottom", legend=levels(as.factor(pData$Batch)), fill=myColors, title="Principal Components by Batch")
pairs(prin$loadings[,1:6], col=as.factor(pData$Batch), labels=c("PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch=1, cex=0.6)

plot.new()
legend("bottom", legend=levels(as.factor(pData$Dx_alg)), fill=myColors, title="Principal Components by Diagnosis")
pairs(prin$loadings[,1:6], col=as.factor(pData$Dx_alg), labels=c("PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch=1, cex=0.6)

plot.new()
legend("bottom", legend=levels(as.factor(pData$COI_GENDER)), fill=myColors, title="Principal Components by Sex")
pairs(prin$loadings[,1:6], col=as.factor(pData$COI_GENDER), labels=c("PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch=1, cex=0.6)

dev.off()

# KS Test ####
preparedData <- prepdata(expressionset = sumData, do.logtransform = FALSE, intgroup = c("Batch"))
boxplot <- aqm.boxplot(preparedData, subsample = 20000, outlierMethod = "KS")
box.outliers <- boxplot@outliers
box.list <- box.outliers@which
box.stat <- box.outliers@statistic
threshold <- box.outliers@threshold #0.03243091
write.csv(box.list,"tables/Child KS boxplot outlier RMA B123.csv")
write.csv(box.stat,"tables/Child KS boxplot KS stat RMA B123.csv")

# Outliers: 110289, 109016, 114600

# MA diagnostics, d statistic ####
ma <- aqm.maplot(preparedData, subsample = 20000, Dthresh = 0.15, maxNumArrays = 8, nrColumns = 4)
ma.outliers <- ma@outliers
ma.list <- ma.outliers@which #No outliers
ma.stat <- ma.outliers@statistic
ma.threshold <- ma.outliers@threshold #0.15
write.csv(ma.stat,"tables/Child MA stat RMA B123.csv")

# Upper Quartile Outliers ####
out.upper <- outliers(order.RMA, method = "upperquartile")
stat.upper <- out.upper@statistic
out.threshold <- out.upper@threshold # 6.067925
upper.list <- out.upper@which # No Outliers
write.csv(stat.upper, "tables/Child upper quartile stat RMA B123.csv")

# Prep Data ####
saveRDS(pData, file = "R_objects/child_pData_B123.rds")

# NUSE Values and Boxplots (Laptop) ####
# New R session
setwd()
library(oligo)
library(VennDiagram)
memory.limit(50000)

rawData <- readRDS("R Objects/child_rawData_B123.rds")
fit <- fitProbeLevelModel(rawData)
rm(rawData)
NUSEvalues <- NUSE(fit, type = "values")
write.csv(NUSEvalues, "Tables/Child NUSE Values B123.csv")
rm(fit)
NUSEvalues <- as.data.frame(NUSEvalues)
medians <- sapply(NUSEvalues, median, na.rm = TRUE)
names(medians[medians > 1.05])
# No outliers

summary(medians)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9842  0.9915  1.0005  1.0028  1.0105  1.0464 
rm(medians, NUSEvalues)

# Outlier venn diagram ####
distance <- as.character(c(110289, 109016, 110790, 114600, 113489, 113555)) # from arrayQualityMetrics
ks <- as.character(c(110289, 109016, 114600)) # from KS test (above)
pc1 <- as.character(c(114600, 113173, 110289, 109016, 110790, 104827, 107053, 113489, 113555))
# NUSE, MA plot d stat, and upper quartile measures detected no outliers

venn.diagram(list("Distance" = distance, "KS Stat" = ks, "PC1" = pc1), 
             filename = "Figures/Child Outlier Venn B123.png", cat.cex = 1.5, cex = 2, fill = c("lightblue", "lightpink", "lightgreen"), resolution = 600, 
             imagetype = "png", cat.pos = c(0,0,0), fontfamily = "sans", cat.fontfamily = "sans", rotation.degree = 0, 
             cat.dist = c(-0.015,-0.015,-0.025), margin=0.05)

outliers <- intersect(intersect(distance, ks), pc1)
#"110289" "109016" "114600"
subset(sampleInfo, IBC %in% outliers)
#    IBC Dx_alg COI_GENDER Batch                                              
# 110289     TD          M     2 
# 109016     TD          M     2 
# 114600    ASD          M     2 

sampleInfo_noOutliers <- subset(sampleInfo, !IBC %in% outliers)
table(sampleInfo_noOutliers$Dx_alg, sampleInfo_noOutliers$COI_GENDER)
#        F  M All
# ASD   11 30  41
# NonTD 17 27  44
# TD    37 40  77

# Normalize without outliers ####
setwd()

library(oligo)
library(pd.hugene.2.0.st)
library(arrayQualityMetrics)

# Remove Outliers
rawData <- readRDS("R_objects/child_rawData_B123.rds")
samples <- rownames(rawData@phenoData@data)
outliers <- c(grep("110289", samples), grep("109016", samples), grep("114600", samples))
rawData <- rawData[,-outliers]
samples2 <- rownames(rawData@phenoData@data)
table(grepl("110289", samples2) | grepl("109016", samples2) | grepl("114600", samples2)) # All False, outliers are excluded
saveRDS(rawData, "R_objects/child_rawData_noOutliers_B123.rds")
rm(samples, outliers, samples2)

# Normalize with RMA
sumData2 <- rma(rawData, target = "core")
rm(rawData)

# Add Features
featureData(sumData2) <- getNetAffx(sumData2, "transcript")

# Remove control probes
tmp <- featureData(sumData2)@data$category
keep <- grep("main", tmp)
sumData3 <- sumData2[keep,]

# Take out genes without annotation
tmp <- featureData(sumData3)@data$geneassignment
drop <- which(is.na(tmp))
sumData3 <- sumData3[-drop,]
sampleInfo <- read.delim(file = "tables/COI_Dx_alg_Pheno_Min.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE) # Only includes cord samples with algorithm diagnosis
sampleInfo_noOutliers <- subset(sampleInfo, !IBC %in% c(110289, 109016, 114600))
pData <- pData(sumData3)
pData$CEL_FileName <- row.names(pData)
pData <- merge(pData, sampleInfo_noOutliers, by="CEL_FileName", sort=FALSE, all.x=TRUE, all.y=FALSE)
pData(sumData3) <- pData
saveRDS(sumData3, "R_objects/child_sumData_noOutliers_filtered_B123.rds")

# Output Normalized Expression Data with Annotations
gene <- sumData3
colnames(gene) <- substr(colnames(gene), 1, 17)
genes <- matrix(unlist(lapply(strsplit(featureData(sumData3)@data$geneassignment, split = "//", fixed = TRUE), function(x) paste(x[1:2]))), ncol = 2, byrow = TRUE)
outdat <- cbind(rownames(gene), genes, exprs(gene))
colnames(outdat)[c(1:3)] <- c("Probe", "GeneID", "GeneName")
write.table(outdat, "tables/Normalized_MARBLES_Child_Expression_Data_noOutliers_B123.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Array Quality Metrics
pData <- pData(sumData2)
pData$CEL_FileName <- row.names(pData)
pData <- merge(pData, sampleInfo_noOutliers, by="CEL_FileName", sort=FALSE, all.x=TRUE, all.y=FALSE)
pData(sumData2) <- pData
arrayQualityMetrics(expressionset=sumData2, outdir="Child Report RMA-Normalized with Batch no Outliers B123",force=TRUE, do.logtransform = FALSE, intgroup = "Batch")
