# MARBLES Differential Gene Expression Analysis ####
# Charles Mordaunt
# 6/9/18

setwd()

# Packages ####
library(sva)
library(limma)
library(oligo)
library(ggplot2)
library(scales)
library(reshape2)
library(VennDiagram)

# Functions ####
# Data ####
sumData3 <- readRDS("R Objects/child_sumData_noOutliers_filtered_B123.rds")
pheno <- pData(sumData3)
pheno$Dx_alg <- factor(pheno$Dx_alg, levels=c("TD", "ASD", "NonTD"), ordered=FALSE)
write.table(pheno, "Tables/Pheno Data for MARBLES Dx_alg Samples B123.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
edata <- exprs(sumData3)
write.table(edata, "Tables/MARBLES Child Expression Data for Dx_alg B123.txt", sep="\t", quote=FALSE)

# Differential Expression Analysis with One Model for Diagnosis ####
# Setup model matrices
pheno$Dx_alg <- relevel(as.factor(pheno$Dx_alg), ref = "TD")
mod <- model.matrix(~Dx_alg, data = pheno)
colnames(mod) # "(Intercept)" "Dx_algASD"   "Dx_algNonTD"
colnames(mod) <- c("Intercept", "Dx_alg_ASD", "Dx_alg_NonTD")
mod0 <- model.matrix(~1, data = pheno)
colnames(mod0) <- "Intercept"

# Run SVA
n.sv <- num.sv(edata, mod, method = "be") #21
svobj <- sva(edata, mod, mod0, n.sv = n.sv)
sv <- svobj$sv

# Differential Expression Analysis with limma
modSv <- cbind(mod, svobj$sv)
mod0Sv <- cbind(mod0, svobj$sv)
fit <- lmFit(edata, modSv)
colnames(modSv)[4:length(colnames(modSv))] <- paste("Sv", 1:21, sep = "")
colnames(mod0Sv)[2:length(colnames(mod0Sv))] <- paste("Sv", 1:21, sep = "")
contrast.matrix <- makeContrasts(Dx_alg_ASD, Dx_alg_NonTD, levels = modSv)
fitContrasts <- contrasts.fit(fit,contrast.matrix)
eb <- eBayes(fitContrasts)

# Make results table
ftable <- topTable(eb, adjust="BH", confint=TRUE, sort.by="none", number=Inf)
ftable$Variance <- apply(edata, 1, var)
asd <- topTable(eb, coef="Dx_alg_ASD", adjust="BH", confint=TRUE, sort.by="none", number=Inf)
NonTD <- topTable(eb, coef="Dx_alg_NonTD", adjust="BH", confint=TRUE, sort.by="none", number=Inf)
outTable <- cbind(ftable, asd[,c("CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B")], NonTD[,c("CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B")])
colnames(outTable)[8:13] <- paste("ASD", colnames(outTable)[8:13], sep="_")
colnames(outTable)[14:19] <- paste("NonTD", colnames(outTable)[14:19], sep="_")
write.table(outTable, "Tables/Child_Dx_alg_limma_results_adjSVAonly_B123.txt", sep="\t", quote=FALSE)

# Add GeneNames
genes <- as.data.frame(cbind(as.character(row.names(sumData3)), 
                             matrix(unlist(lapply(strsplit(featureData(sumData3)@data$geneassignment, split = "//", fixed = TRUE), function(x) paste(x[1:2]))), 
                ncol = 2, byrow = TRUE)))
colnames(genes) <- c("Probe", "GeneID", "GeneName")
stats <- outTable
stats$Probe <- row.names(stats)
stats <- merge(stats, genes, by="Probe", all.x=TRUE, all.y=FALSE, sort=FALSE)
write.table(stats, "Tables/MARBLES All Genes Stats Dx_alg SVAonly B123.txt", sep = "\t", quote = FALSE)

# Differential Probes ####
# ASD Differential Probes
ASD_diff <- as.character(stats$Probe[stats$ASD_P.Value < 0.01 & abs(stats$Dx_alg_ASD) > 0.1])
ASD_diff <- sort(unique(ASD_diff))
write.table(x=ASD_diff, file="Differentially Expressed Probes Dx_alg B123/MARBLES_ASD_diff_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

ASD_up <- as.character(stats$Probe[stats$ASD_P.Value < 0.01 & stats$Dx_alg_ASD > 0.1])
ASD_up <- sort(unique(ASD_up))
write.table(x=ASD_up, file="Differentially Expressed Probes Dx_alg B123/MARBLES_ASD_up_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

ASD_down <- as.character(stats$Probe[stats$ASD_P.Value < 0.01 & stats$Dx_alg_ASD < -0.1])
ASD_down <- sort(unique(ASD_down))
write.table(x=ASD_down, file="Differentially Expressed Probes Dx_alg B123/MARBLES_ASD_down_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

# NonTD Differential Probes
NonTD_diff <- as.character(stats$Probe[stats$NonTD_P.Value < 0.01 & abs(stats$Dx_alg_NonTD) > 0.1])
NonTD_diff <- sort(unique(NonTD_diff))
write.table(x=NonTD_diff, file="Differentially Expressed Probes Dx_alg B123/MARBLES_NonTD_diff_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

NonTD_up <- as.character(stats$Probe[stats$NonTD_P.Value < 0.01 & stats$Dx_alg_NonTD > 0.1])
NonTD_up <- sort(unique(NonTD_up))
write.table(x=NonTD_up, file="Differentially Expressed Probes Dx_alg B123/MARBLES_NonTD_up_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

NonTD_down <- as.character(stats$Probe[stats$NonTD_P.Value < 0.01 & stats$Dx_alg_NonTD < -0.1])
NonTD_down <- sort(unique(NonTD_down))
write.table(x=NonTD_down, file="Differentially Expressed Probes Dx_alg B123/MARBLES_NonTD_down_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

# All probes
all <- sort(unique(as.character(stats$Probe)))
write.table(x=all, file="Differentially Expressed Probes Dx_alg B123/MARBLES_all_probes_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
