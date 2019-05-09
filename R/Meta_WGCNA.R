# WGCNA Consensus Module Meta-Analysis ---------------------------------------
# ASD Cord Expression Revisions
# Charles Mordaunt
# 4/8/19

# Packages ####
setwd("~/CM_Microarray_CordBlood")
sapply(c("tidyverse", "oligo", "sva","pd.hugene.2.0.st", "matrixStats", "biomaRt", "WGCNA", "VennDiagram"), require, character.only = TRUE)

# Data ####
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
enableWGCNAThreads()
exp_marbles <- readRDS("R Objects/MARBLES_Exp_BatchAdj.rds")
exp_earli <- readRDS("R Objects/EARLI_Exp_noBatchAdj.rds")
exp <- list(marbles = list(data = as.data.frame(exp_marbles)),
            earli = list(data = as.data.frame(exp_earli)))
checkSets(exp)
# $nSets 2
# $nGenes 36459
# $nSamples 161 108
# $structureOK TRUE
rm(exp_earli, exp_marbles)

# Get Consensus Modules ####
consensusMods <- blockwiseConsensusModules(exp, checkMissingData = FALSE, maxBlockSize = 40000, corType = "bicor",
                                           maxPOutliers = 0.1, power = 10, networkType = "signed", 
                                           checkPower = FALSE, TOMType = "signed", 
                                           networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                                           deepSplit = 4, mergeCutHeight = 0.1, verbose = 5)
table(consensusMods$colors) %>% sort(decreasing = TRUE)

# Plot Merged Gene Dendrogram with Modules
pdf("Figures/Revisions Consensus Dx_alg WGCNA Merged Gene Module Dendrogram.pdf", width = 10, height = 5)
sizeGrWindow(10, 5)
plotDendroAndColors(dendro = consensusMods$dendrograms[[1]], colors = consensusMods$colors, 
                    groupLabels = "Modules", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, marAll = c(1, 5, 1, 0), main = "", cex.colorLabels = 1.3)
dev.off()

# Cluster Modules by MARBLES Eigengenes and Plot Dendrogram
METree <- (1 - bicor(consensusMods$multiMEs$marbles$data, maxPOutliers = 0.1)) %>% as.dist %>% 
        hclust(method = "average")
pdf("Figures/Revisions Consensus Dx_alg WGCNA MARBLES Module Eigengene Dendrogram.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()

# Cluster Modules by EARLI Eigengenes and Plot Dendrogram
METree <- (1 - bicor(consensusMods$multiMEs$earli$data, maxPOutliers = 0.1)) %>% as.dist %>% 
        hclust(method = "average")
pdf("Figures/Revisions Consensus Dx_alg WGCNA EARLI Module Eigengene Dendrogram.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()
rm(METree)

# Compare Eigengene Networks Between MARBLES and EARLI
consensusMEs <- consensusOrderMEs(consensusMods$multiMEs)
pdf(file = "Figures/Revisions Consensus Dx_alg WGCNA Eigengene Networks.pdf", width = 8, height = 7)
sizeGrWindow(width = 8, height = 7)
par(cex = 0.8)
plotEigengeneNetworks(consensusMEs, setLabels = c("MARBLES", "EARLI"), 
                      plotDendrograms = FALSE, marHeatmap = c(3, 3, 2, 1), zlimPreservation = c(0.5, 1), 
                      xLabelsAngle = 90)
dev.off()

# Calculate Module Membership ####
moduleMembership <- mtd.mapply(bicorAndPvalue, exp, consensusMEs, 
                               MoreArgs = list(alternative = "two.sided", use = "pairwise.complete.obs", 
                                               maxPOutliers = 0.1))
MM_marbles <- as.data.frame(moduleMembership$marbles$data$bicor)
colnames(MM_marbles) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_marbles), fixed = TRUE)
MM_marbles$Probe <- rownames(MM_marbles)
MM_marbles$Module <- consensusMods
write.table(MM_marbles, "Tables/Revision Consensus Modules MARBLES Probe Module Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_earli <- as.data.frame(moduleMembership$earli$data$bicor)
colnames(MM_earli) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_earli), fixed = TRUE)
MM_earli$Probe <- rownames(MM_earli)
MM_earli$Module <- consensusMods
write.table(MM_earli, "Tables/Revision Consensus Modules EARLI Probe Module Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Get Module Hub Probes and Genes ####
hubProbes_marbles <- sapply(colnames(MM_marbles)[!colnames(MM_marbles) %in% c("Probe", "Module")], function(x){
        temp <- MM_marbles[MM_marbles$Module == x,]
        temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})
hubProbes_earli <- sapply(colnames(MM_earli)[!colnames(MM_earli) %in% c("Probe", "Module")], function(x){
        temp <- MM_earli[MM_earli$Module == x,]
        temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
hubGenes_marbles <- lapply(hubProbes_marbles, function(x){
        getBM(attributes = "external_gene_name", filters = "affy_hugene_2_0_st_v1", values = x, mart = ensembl, 
              verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
hubGenes_earli <- lapply(hubProbes_earli, function(x){
        getBM(attributes = "external_gene_name", filters = "affy_hugene_2_0_st_v1", values = x, mart = ensembl, 
              verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist

# Combine Covariates ####
cov_marbles <- read.delim("Tables/Revisions MARBLES Covariates for WGCNA.txt", sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE)
cov_earli <- read.delim("Tables/Revisions EARLI Covariates for WGCNA.txt", sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)
cov_marbles <- cov_marbles[,c("IBC", "Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", 
                              "gest_age_deliv_wk", "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth",
                              "MomAgeYr", "DadAgeYr", "Mat_BMI_PrePreg", "DeliveryMethod", 
                              "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", "MomEduBachelors", 
                              "DadEduBachelors", "HomeOwn")]
colnames(cov_marbles)[colnames(cov_marbles) == "IBC"] <- "sampleID"
table(colnames(cov_marbles) == colnames(cov_earli)) # All TRUE

rownames(cov_marbles) <- cov_marbles$sampleID
cov_marbles <- cov_marbles[,c("Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", "gest_age_deliv_wk", 
                              "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth", "MomAgeYr", "DadAgeYr", 
                              "Mat_BMI_PrePreg", "DeliveryMethod", "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", 
                              "MomEduBachelors", "DadEduBachelors", "HomeOwn")]
cov_marbles <- as.matrix(cov_marbles)
table(rownames(exp$marbles$data) == rownames(cov_marbles)) # All TRUE

rownames(cov_earli) <- cov_earli$sampleID
cov_earli <- cov_earli[,c("Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", "gest_age_deliv_wk", 
                          "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth", "MomAgeYr", "DadAgeYr", 
                          "Mat_BMI_PrePreg", "DeliveryMethod", "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", 
                          "MomEduBachelors", "DadEduBachelors", "HomeOwn")]
cov_earli <- as.matrix(cov_earli)
table(rownames(exp$earli$data) == rownames(cov_earli)) # All TRUE
pheno <- list(marbles = list(data = cov_marbles),
              earli = list(data = cov_earli))

# Get Meta-Analysis Correlations ####
moduleMembership <- read.delim("Tables/Revision Consensus Modules MARBLES Probe Module Membership.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE)
consensusMods <- moduleMembership$Module
MEs_marbles <- moduleEigengenes(exp_marbles, colors = consensusMods)$eigengenes
rownames(MEs_marbles) <- rownames(exp_marbles)
MEs_earli <- moduleEigengenes(exp_earli, colors = consensusMods)$eigengenes
rownames(MEs_earli) <- rownames(exp_earli)
consensusMEs <- list(marbles = list(data = MEs_marbles), earli = list(data = MEs_earli))
consensusMEs <- orderMEs(consensusMEs)

MEMAs <- list()
for (t in 1:20){
        MEMAs[[t]] = metaAnalysis(consensusMEs, mtd.subset(pheno, colIndex = t), useRankPvalue = FALSE, 
                                  corFnc = bicor,
                                  corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs"), 
                                  getQvalues = TRUE)
}

zscores <- sapply(MEMAs, function(x) x[["Z.RootDoFWeights"]])
rownames(zscores) <- colnames(consensusMEs$marbles$data)
colnames(zscores) <- colnames(pheno$marbles$data)
pvalues <- sapply(MEMAs, function(x) x[["p.RootDoFWeights"]])
dimnames(pvalues) <- dimnames(zscores)
qvalues <- sapply(MEMAs, function(x) x[["q.RootDoFWeights"]])
dimnames(qvalues) <- dimnames(zscores)

# Plot Correlations ####
# Plot All Correlations for Meta Analysis (Z-scores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Figures/Revisions Consensus Modules Meta Covariate Correlation Plot.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = zscores, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Plot Significant Correlations for Meta Analysis (Z-scores)
sigRows <- rowSums2(qvalues < 0.05) >= 1
sigCols <- colSums2(qvalues < 0.05) >= 1
zscores_sub <- zscores[sigRows, sigCols]
qvalues_sub <- qvalues[sigRows, sigCols]
colnames(zscores_sub) <- c("Diagnosis", "ASD vs TD", "Sex", "Gestational Age", "Birthweight", "Hispanic Ethnicity",
                           "Father Age", "Cesarean Delivery", "Mother Smoking")
star <- apply(qvalues_sub, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Figures/Revisions Consensus Modules Meta Covariate Correlation Plot Significant Only.pdf", width = 8, height = 11)
sizeGrWindow(width = 8, height = 11)
par(mar = c(6, 8, 1, 1))
labeledHeatmap(Matrix = zscores_sub, xLabels = colnames(zscores_sub), yLabels = rownames(zscores_sub), 
               ySymbols = gsub("ME", "", rownames(zscores_sub)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Plot All Correlations for MARBLES (Z-scores)
zscores <- sapply(MEMAs, function(x) x[["Z.marbles"]])
rownames(zscores) <- colnames(consensusMEs$marbles$data)
colnames(zscores) <- colnames(pheno$marbles$data)
pvalues <- sapply(MEMAs, function(x) x[[which(names(x) %in% c("pvalueStudent.marbles", "pvalueStudent.1.vs.2.marbles"))]])
dimnames(pvalues) <- dimnames(zscores)
qvalues <- sapply(MEMAs, function(x) x[[which(names(x) %in% c("qvalueStudent.marbles", "q.Student.marbles"))]])
dimnames(qvalues) <- dimnames(zscores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Figures/Revisions Consensus Modules MARBLES Covariate Correlation Plot zscores.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = zscores, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Plot All Correlations for EARLI (Z-scores)
zscores <- sapply(MEMAs, function(x) x[["Z.earli"]])
rownames(zscores) <- colnames(consensusMEs$earli$data)
colnames(zscores) <- colnames(pheno$earli$data)
pvalues <- sapply(MEMAs, function(x) x[[which(names(x) %in% c("pvalueStudent.earli", "pvalueStudent.1.vs.2.earli"))]])
dimnames(pvalues) <- dimnames(zscores)
qvalues <- sapply(MEMAs, function(x) x[[which(names(x) %in% c("qvalueStudent.earli", "q.Student.earli"))]])
dimnames(qvalues) <- dimnames(zscores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Figures/Revisions Consensus Modules EARLI Covariate Correlation Plot zscores.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = zscores, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Put Together Meta-Analysis Stats Table ####
n_Probes <- sapply(gsub(pattern = "ME", replacement = "", 
                        x = colnames(MM_marbles)[!colnames(MM_marbles) %in% c("Probe", "Module")], fixed = TRUE),
                   function(x) length(MM_marbles$Module[MM_marbles$Module == x]))
moduleStats_meta <- cbind(colnames(consensusMEs$marbles$data),
                          n_Probes,
                          sapply(hubProbes_marbles, paste, collapse = ", "),
                          sapply(hubProbes_earli, paste, collapse = ", "),
                          hubGenes_marbles,
                          hubGenes_earli,
                          sapply(MEMAs, function(x) x[, "Z.marbles"]),
                          sapply(MEMAs, function(x) x[, colnames(x)[colnames(x) %in% c("pvalueStudent.marbles", "pvalueStudent.1.vs.2.marbles")]]),
                          sapply(MEMAs, function(x) x[, colnames(x)[colnames(x) %in% c("qvalueStudent.marbles", "q.Student.marbles")]]),
                          sapply(MEMAs, function(x) x[, "Z.earli"]),
                          sapply(MEMAs, function(x) x[, colnames(x)[colnames(x) %in% c("pvalueStudent.earli", "pvalueStudent.1.vs.2.earli")]]),
                          sapply(MEMAs, function(x) x[, colnames(x)[colnames(x) %in% c("qvalueStudent.earli", "q.Student.earli")]]),
                          sapply(MEMAs, function(x) x[, "Z.RootDoFWeights"]),
                          sapply(MEMAs, function(x) x[, "p.RootDoFWeights"]),
                          sapply(MEMAs, function(x) x[, "q.RootDoFWeights"]))
moduleStats_meta <- as.data.frame(moduleStats_meta)
colnames(moduleStats_meta) <- c("Module", "n_Probes", "hubProbes_MARBLES", "hubProbes_EARLI", "hubGenes_MARBLES", "hubGenes_EARLI",
                                paste(rep(c("Zscore_MARBLES", "pvalue_MARBLES", "qvalue_MARBLES", "Zscore_EARLI", "pvalue_EARLI",
                                            "qvalue_EARLI", "Zscore_Meta", "pvalue_Meta", "qvalue_Meta"), 
                                          each = length(colnames(pheno$marbles$data))), 
                                      rep(colnames(pheno$marbles$data), 5), sep = "_"))
write.table(moduleStats_meta, "Tables/Revision Consensus Modules Meta Covariate Correlation Stats with Hub Genes.txt", sep = "\t", quote = FALSE)

# Skyblue1 Module Expression ####
skyblue1_probes <- colnames(exp$marbles$data)[consensusMods$colors == "skyblue1"]
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
skyblue1_genes <- getBM(attributes = "external_gene_name", filters = "affy_hugene_2_0_st_v1", 
                        values = skyblue1_probes, mart = ensembl) %>% unlist %>% as.character %>% unique %>% sort

# skyblue1 module eigengene plots
skyblue1_ME_marbles <- consensusMEs$marbles$data$MEskyblue1
skyblue1_ME_earli <- consensusMEs$earli$data$MEskyblue1
skyblue1_ME <- as.data.frame(cbind(c(rownames(pheno$marbles$data), rownames(pheno$earli$data)),
                                   c(pheno$marbles$data[, "Dx_alg"], pheno$earli$data[, "Dx_alg"]),
                                   c(pheno$marbles$data[, "Sex"], pheno$earli$data[, "Sex"]),
                                   c(skyblue1_ME_marbles, skyblue1_ME_earli)))
colnames(skyblue1_ME) <- c("sampleID", "Diagnosis", "Sex", "ME")
skyblue1_ME$Diagnosis[skyblue1_ME$Diagnosis == 1] <- "TD"
skyblue1_ME$Diagnosis[skyblue1_ME$Diagnosis == 2] <- "NonTD"
skyblue1_ME$Diagnosis[skyblue1_ME$Diagnosis == 3] <- "ASD"
skyblue1_ME$Sex[skyblue1_ME$Sex == 1] <- "M"
skyblue1_ME$Sex[skyblue1_ME$Sex == 2] <- "F"
skyblue1_ME$Study <- c(rep("MARBLES", length(skyblue1_ME_marbles)), rep("EARLI", length(skyblue1_ME_earli)))
skyblue1_ME$Diagnosis <- factor(skyblue1_ME$Diagnosis, levels = c("TD", "NonTD", "ASD"))
skyblue1_ME$Sex <- factor(skyblue1_ME$Sex, levels = c("M", "F"))
skyblue1_ME$Study <- factor(skyblue1_ME$Study, levels = c("MARBLES", "EARLI"))
skyblue1_ME$ME <- as.numeric(skyblue1_ME$ME)
skyblue1_ME <- skyblue1_ME[, c("sampleID", "Study", "Diagnosis", "Sex", "ME")]
gg <- ggplot(data = skyblue1_ME)
gg <- gg +
        geom_boxplot(aes(x = Sex, y = ME, fill = Diagnosis), size = 0.8, outlier.size = 0.8) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.87, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_text(size = 22),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0, 1, 1, 0.5), "lines"), axis.title.x = element_blank(), 
              axis.text = element_text(size = 22, color = "black")) +
        ylab("skyblue1 Module Eigengene") +
        scale_fill_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        scale_color_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        coord_cartesian(ylim = c(-0.13, 0.13)) +
        facet_wrap(facets = vars(Study), nrow = 2)
ggsave(filename = "Figures/skyblue1 Module Eigengene Boxplot.png", plot = gg, dpi = 600, width = 8, height = 9, units = "in")

# skyblue1 hub probes plots
moduleStats_meta$hubProbes_MARBLES[moduleStats_meta$Module == "MEskyblue1"] #"17117126"
moduleStats_meta$hubProbes_EARLI[moduleStats_meta$Module == "MEskyblue1"] #"17117126"
skyblue1_hub_marbles <- exp$marbles$data[,"17117126"]
skyblue1_hub_earli <- exp$earli$data[,"17117126"]
skyblue1_hub <- as.data.frame(cbind(skyblue1_ME[, c("sampleID", "Study", "Diagnosis", "Sex")],
                                    c(skyblue1_hub_marbles, skyblue1_hub_earli)))
colnames(skyblue1_hub) <- c("sampleID", "Study", "Diagnosis", "Sex", "KDM5D")
gg <- ggplot(data = skyblue1_hub)
gg <- gg +
        geom_boxplot(aes(x = Sex, y = KDM5D, fill = Diagnosis), size = 0.8, outlier.size = 0.8) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.87, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_text(size = 22),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0, 1, 1, 0.5), "lines"), axis.title.x = element_blank(), 
              axis.text = element_text(size = 22, color = "black")) +
        ylab(expression(log[2]*"("*italic(KDM5D)*" Expression)", parse = TRUE)) +
        scale_fill_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        scale_color_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        coord_cartesian(ylim = c(3, 8)) +
        facet_wrap(facets = vars(Study), nrow = 2)
ggsave(filename = "Figures/skyblue1 Hub KDM5D Expression Boxplot.png", plot = gg, dpi = 600, width = 8, height = 9, units = "in")

skyblue1_hub_stats <- list(Diagnosis = metaAnalysis(mtd.subset(exp, colIndex = which(colnames(exp$marbles$data) %in% c("17117126", "16657436"))), 
                                                    mtd.subset(pheno, colIndex = 1), useRankPvalue = FALSE, corFnc = bicor, 
                                                    corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs")),
                           Sex = metaAnalysis(mtd.subset(exp, colIndex = which(colnames(exp$marbles$data) %in% c("17117126", "16657436"))), 
                                              mtd.subset(pheno, colIndex = 5), useRankPvalue = FALSE, corFnc = bicor, 
                                              corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs")))
skyblue1_hub_stats$Diagnosis$Z.RootDoFWeights[2] # Diagnosis Z = 4.851036
skyblue1_hub_stats$Diagnosis$p.RootDoFWeights[2] # Diagnosis p = 1.228185e-06
skyblue1_hub_stats$Sex$Z.RootDoFWeights[2] # Sex Z = -27.28498
skyblue1_hub_stats$Sex$p.RootDoFWeights[2] # Sex p = 6.394725e-164

# Overlap of Diagnosis Modules with Meta-Analysis Differentially-Expressed Genes ####
# Data
moduleStats_meta <- read.delim("Tables/Revision Consensus Modules Meta Covariate Correlation Stats with Hub Genes.txt",
                               sep = "\t", header = TRUE, stringsAsFactors = FALSE)
moduleMembership <- read.delim("Tables/Revision Consensus Modules MARBLES Probe Module Membership.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE)
ASD <- read.delim("Tables/ASD Meta-Analysis with METAL All Genes B123.txt", sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
NonTD <- read.delim("Tables/NonTD Meta-Analysis with METAL All Genes B123.txt", sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
table(ASD$Probe == NonTD$Probe) # All TRUE

# Probes in Diagnosis Associated Modules
DxMods <- moduleStats_meta$Module[moduleStats_meta$pvalue_Meta_Dx_alg < 0.05 | 
                                          moduleStats_meta$pvalue_Meta_ASDvsTD < 0.05 |
                                          moduleStats_meta$pvalue_Meta_NonTDvsTD < 0.05 |
                                          moduleStats_meta$pvalue_Meta_ASDvsNonTD < 0.05]
# "MEgrey60"   "MEskyblue1"
probes_grey60 <- moduleMembership$Probe[moduleMembership$Module == "grey60"]
probes_skyblue1 <- moduleMembership$Probe[moduleMembership$Module == "skyblue1"]

# Probes from Differential Expression Analysis
probes_ASDdiff <- ASD$Probe[abs(ASD$Meta_logFC) > 0.1 & ASD$Meta_pValue < 0.01 & 
                                    ASD$marbles_logFC * ASD$earli_logFC > 0]
probes_NonTDdiff <- NonTD$Probe[abs(NonTD$Meta_logFC) > 0.1 & NonTD$Meta_pValue < 0.01 & 
                                        NonTD$marbles_logFC * NonTD$earli_logFC > 0]
probes_ASDup <- ASD$Probe[ASD$Meta_logFC > 0.1 & ASD$Meta_pValue < 0.01 & ASD$marbles_logFC * ASD$earli_logFC > 0]

# Overlap
allProbes <- ASD$Probe
intersect(probes_grey60, probes_ASDdiff) #  16927756 (IGLV1-40), 16927853 (IGLC2)
fisher.test(allProbes %in% probes_grey60, allProbes %in% probes_ASDdiff) # odds ratio = 1.80474, p-value = 0.3078
intersect(probes_grey60, probes_NonTDdiff) # 16927756 (IGLV1-40)
fisher.test(allProbes %in% probes_grey60, allProbes %in% probes_NonTDdiff) # odds ratio = 2.440187, p-value = 0.3405
intersect(probes_skyblue1, probes_ASDdiff) # 17115953 (ZFY), 17117173 (TTTY10)
fisher.test(allProbes %in% probes_skyblue1, allProbes %in% probes_ASDdiff) # odds ratio = 17.15977, p-value = 0.007132
intersect(probes_skyblue1, probes_ASDup) # 17115953 (ZFY), 17117173 (TTTY10)
fisher.test(allProbes %in% probes_skyblue1, allProbes %in% probes_ASDup) # odds ratio 33.60367, p-value = 0.001968
intersect(probes_skyblue1, probes_NonTDdiff) # None

# Overlapping Genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_overlap <- getBM(attributes = c("affy_hugene_2_0_st_v1", "external_gene_name"), filters = "affy_hugene_2_0_st_v1", 
                                      values = c(16927756, 16927853, 17115953, 17117173), mart = ensembl)
ASD[ASD$Probe %in% c(16927756, 16927853, 17115953, 17117173), c("Probe", "Meta_logFC")]
NonTD[NonTD$Probe %in% 16927756, c("Probe", "Meta_logFC")]

venn.diagram(list("ASD" = probes_ASDdiff, "skyblue1" = probes_skyblue1), 
             rotation.degree = 0, margin = 0.02, cat.cex = 1, cex = 1.4, fill = c("lightblue", "lightpink"), 
             filename = "Figures/Revisions ASD Diff and skyblue1 Module Probe Overlap Venn.png", 
             resolution = 600, imagetype = "png", cat.pos = c(180,180), fontfamily = "Arial", cat.fontfamily = "Arial", 
             height = 2100, width = 2200, cat.dist = c(0.02, 0.02), ext.dist = -0.15, ext.length = 0.8)

venn.diagram(list("ASD Up" = probes_ASDup, "skyblue1" = probes_skyblue1), 
             rotation.degree = 0, margin = 0.02, cat.cex = 1, cex = 1.4, fill = c("lightblue", "lightpink"), 
             filename = "Figures/Revisions ASD Up and skyblue1 Module Probe Overlap Venn.png", 
             resolution = 600, imagetype = "png", cat.pos = c(180,180), fontfamily = "Arial", cat.fontfamily = "Arial", 
             height = 2100, width = 2200, cat.dist = c(0.02, 0.02), ext.dist = -0.15, ext.length = 0.8)
rm(ASD, ensembl, genes_overlap, NonTD, allProbes, DxMods, probes_ASDdiff, probes_ASDup, probes_grey60,
   probes_NonTDdiff, probes_skyblue1)

# Grey60 Module Expression ####
# Expression Data
exp_marbles <- readRDS("R Objects/MARBLES_Exp_BatchAdj.rds")
exp_earli <- readRDS("R Objects/EARLI_Exp_noBatchAdj.rds")
consensusMods <- moduleMembership$Module
MEs_marbles <- moduleEigengenes(exp_marbles, colors = consensusMods)$eigengenes
rownames(MEs_marbles) <- rownames(exp_marbles)
MEs_earli <- moduleEigengenes(exp_earli, colors = consensusMods)$eigengenes
rownames(MEs_earli) <- rownames(exp_earli)

# Plot Module Eigengene
grey60_ME <- as.data.frame(cbind(c(rownames(pheno$marbles$data), rownames(pheno$earli$data)),
                                   c(pheno$marbles$data[, "Dx_alg"], pheno$earli$data[, "Dx_alg"]),
                                   c(pheno$marbles$data[, "DeliveryMethod"], pheno$earli$data[, "DeliveryMethod"]),
                                   c(MEs_marbles$MEgrey60, MEs_earli$MEgrey60)), stringsAsFactors = FALSE)
colnames(grey60_ME) <- c("sampleID", "Diagnosis", "DeliveryMethod", "ME")
grey60_ME$Diagnosis[grey60_ME$Diagnosis == 1] <- "TD"
grey60_ME$Diagnosis[grey60_ME$Diagnosis == 2] <- "NonTD"
grey60_ME$Diagnosis[grey60_ME$Diagnosis == 3] <- "ASD"
grey60_ME$DeliveryMethod[grey60_ME$DeliveryMethod == 1] <- "Vaginal"
grey60_ME$DeliveryMethod[grey60_ME$DeliveryMethod == 2] <- "Cesarean"
grey60_ME$Study <- c(rep("MARBLES", nrow(MEs_marbles)), rep("EARLI", nrow(MEs_earli)))
grey60_ME$Diagnosis <- factor(grey60_ME$Diagnosis, levels = c("TD", "NonTD", "ASD"))
grey60_ME$DeliveryMethod <- factor(grey60_ME$DeliveryMethod, levels = c("Vaginal", "Cesarean"))
grey60_ME$Study <- factor(grey60_ME$Study, levels = c("MARBLES", "EARLI"))
grey60_ME$ME <- as.numeric(grey60_ME$ME)
grey60_ME <- grey60_ME[, c("sampleID", "Study", "Diagnosis", "DeliveryMethod", "ME")]
grey60_ME <- subset(grey60_ME, !is.na(grey60_ME$DeliveryMethod))

gg <- ggplot(data = grey60_ME)
gg <- gg +
        geom_boxplot(aes(x = DeliveryMethod, y = ME, fill = Diagnosis), size = 0.8, outlier.size = 0.8) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.87, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_text(size = 22),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0, 1, 1, 0.5), "lines"), axis.title.x = element_blank(), 
              axis.text = element_text(size = 22, color = "black")) +
        ylab("grey60 Module Eigengene") +
        scale_fill_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        scale_color_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        coord_cartesian(ylim = c(-0.4, 0.4)) +
        facet_wrap(facets = vars(Study), nrow = 2)
ggsave(filename = "Figures/grey60 Module Eigengene by Diagnosis and Delivery Method Boxplot.png", plot = gg, dpi = 600, width = 8, height = 9, units = "in")

# Plot Hub Probe
moduleStats_meta$hubProbes_MARBLES[moduleStats_meta$Module == "MEgrey60"] #"16862604"
moduleStats_meta$hubProbes_EARLI[moduleStats_meta$Module == "MEgrey60"] #"16862604"
grey60_hub <- as.data.frame(cbind(grey60_ME[, c("sampleID", "Study", "Diagnosis", "DeliveryMethod")],
                                    c(exp_marbles[!is.na(pheno$marbles$data[, "DeliveryMethod"]),"16862604"], 
                                      exp_earli[!is.na(pheno$earli$data[, "DeliveryMethod"]),"16862604"])))
colnames(grey60_hub) <- c("sampleID", "Study", "Diagnosis", "DeliveryMethod", "CD79A")

gg <- ggplot(data = grey60_hub)
gg <- gg +
        geom_boxplot(aes(x = DeliveryMethod, y = CD79A, fill = Diagnosis), size = 0.8, outlier.size = 0.8) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.87, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_text(size = 22),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0, 1, 1, 0.5), "lines"), axis.title.x = element_blank(), 
              axis.text = element_text(size = 22, color = "black")) +
        ylab(expression(log[2]*"("*italic(CD79A)*" Expression)", parse = TRUE)) +
        scale_fill_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        scale_color_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        coord_cartesian(ylim = c(8, 13)) +
        facet_wrap(facets = vars(Study), nrow = 2)
ggsave(filename = "Figures/grey60 Hub CD79A Expression by Diagnosis and Delivery Method Boxplot.png", plot = gg, dpi = 600, width = 8, height = 9, units = "in")

exp <- list(marbles = list(data = exp_marbles), earli = list(data = exp_earli))
grey60_hub_stats <- list(ASDvsNonTD = metaAnalysis(mtd.subset(exp, colIndex = which(colnames(exp$marbles$data) %in% c("17117126", "16862604"))), 
                                                    mtd.subset(pheno, colIndex = 4), useRankPvalue = FALSE, corFnc = bicor, 
                                                    corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs")),
                           DeliveryMethod = metaAnalysis(mtd.subset(exp, colIndex = which(colnames(exp$marbles$data) %in% c("17117126", "16862604"))), 
                                              mtd.subset(pheno, colIndex = 15), useRankPvalue = FALSE, corFnc = bicor, 
                                              corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs")))
grey60_hub_stats$ASDvsNonTD$Z.RootDoFWeights[1] # ASDvsNonTD Z = -1.988504
grey60_hub_stats$ASDvsNonTD$p.RootDoFWeights[1] # ASDvsNonTD p = 0.04675598
grey60_hub_stats$DeliveryMethod$Z.RootDoFWeights[1] # DeliveryMethod Z = 3.527298
grey60_hub_stats$DeliveryMethod$p.RootDoFWeights[1] # DeliveryMethod p = 0.0004198241
