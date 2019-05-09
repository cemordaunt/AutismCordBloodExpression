# CIBERSORT Meta-Analysis --------------------------------------
# Cord Expression Revisions
# Charles Mordaunt
# 4/21/19

# Packages ####
sapply(c("tidyverse", "scales", "variancePartition", "sm", "biomaRt", "reshape2", "wesanderson", "WGCNA", "matrixStats"), 
       require, character.only = TRUE)

# Format Mixture Files ####
# Expression Data
exp_marbles <- readRDS("R Objects/MARBLES_Exp_BatchAdj.rds") %>% t %>% as.data.frame
exp_earli <- readRDS("R Objects/EARLI_Exp_noBatchAdj.rds") %>% t %>% as.data.frame
table(rownames(exp_marbles) == rownames(exp_earli)) # All TRUE

# Match Probe IDs with HGNC Gene Symbols
probes <- rownames(exp_marbles)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("affy_hugene_2_0_st_v1", "hgnc_symbol"), filters = "affy_hugene_2_0_st_v1", 
               values = probes, mart = ensembl)
table(probes %in% genes$affy_hugene_2_0_st_v1) # 1671 probes missing gene
table(duplicated(genes$affy_hugene_2_0_st_v1)) # 5628 duplicates
probes <- as.data.frame(probes, stringsAsFactors = FALSE)
matched <- merge(x = probes, y = genes, by.x = "probes", by.y = "affy_hugene_2_0_st_v1", all = FALSE, sort = FALSE)
matched <- subset(matched, !duplicated(matched$probes) & !matched$hgnc_symbol == "")

# MARBLES File
exp_marbles$probes <- rownames(exp_marbles)
exp_marbles <- merge(x = matched, y = exp_marbles, by = "probes", all = FALSE)
exp_marbles <- exp_marbles[, !colnames(exp_marbles) == "probes"]
write.table(exp_marbles, file = "Tables/MARBLES_Exp_BatchAdj_CIBERSORT.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# EARLI File
exp_earli$probes <- rownames(exp_earli)
exp_earli <- merge(x = matched, y = exp_earli, by = "probes", all = FALSE)
exp_earli <- exp_earli[, !colnames(exp_earli) == "probes"]
table(exp_marbles$hgnc_symbol == exp_earli$hgnc_symbol) # All TRUE
write.table(exp_earli, file = "Tables/EARLI_Exp_noBatchAdj_CIBERSORT.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(ensembl, exp_earli, exp_marbles, genes, matched, probes)

# Pheno Data ####
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
rownames(cov_earli) <- cov_earli$sampleID
cov_earli <- cov_earli[,c("Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", "gest_age_deliv_wk", 
                          "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth", "MomAgeYr", "DadAgeYr", 
                          "Mat_BMI_PrePreg", "DeliveryMethod", "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", 
                          "MomEduBachelors", "DadEduBachelors", "HomeOwn")]
cov_earli <- as.matrix(cov_earli)
pheno <- list(marbles = list(data = cov_marbles),
              earli = list(data = cov_earli))
rm(cov_earli, cov_marbles)

# Run with CIBERSORT ####
# Ran with relative and absolute modes together
# LM22 Default Reference File
# Permutations = 100
# Disabled quantile normalization
ciber_marbles <- read.delim("Tables/MARBLES CIBERSORT Results.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ciber_marbles <- ciber_marbles[,!colnames(ciber_marbles) %in% c("P.value", "Pearson.Correlation", "RMSE"), ]
ciber_earli <- read.delim("Tables/EARLI CIBERSORT Results.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ciber_earli <- ciber_earli[,!colnames(ciber_earli) %in% c("P.value", "Pearson.Correlation", "RMSE"), ]
ciber_combined <- rbind(ciber_marbles, ciber_earli)
ciber_combined$Study <- c(rep("MARBLES", nrow(ciber_marbles)), rep("EARLI", nrow(ciber_earli)))
table(ciber_combined$Input.Sample == c(rownames(pheno$marbles$data), rownames(pheno$earli$data))) # All TRUE
ciber_combined$Diagnosis <- c(pheno$marbles$data[, "Dx_alg"], pheno$earli$data[, "Dx_alg"])
ciber_combined$Diagnosis[ciber_combined$Diagnosis == 1] <- "TD"
ciber_combined$Diagnosis[ciber_combined$Diagnosis == 2] <- "NonTD"
ciber_combined$Diagnosis[ciber_combined$Diagnosis == 3] <- "ASD"
ciber_combined <- melt(ciber_combined, id.vars = c("Input.Sample", "Study", "Diagnosis"))
colnames(ciber_combined) <- c("SampleID", "Study", "Diagnosis", "CellType", "Fraction")
ciber_combined$Study <- factor(ciber_combined$Study, levels = c("MARBLES", "EARLI"))
ciber_combined$Diagnosis <- factor(ciber_combined$Diagnosis, levels = c("TD", "NonTD", "ASD"))
overall_mean <- aggregate(Fraction ~ CellType, data = ciber_combined, FUN = mean)

# Get average proportion of each cell type by study with meta-analysis
ciber_mean <- aggregate(Fraction ~ CellType + Diagnosis + Study, data = ciber_combined, FUN = mean)
cellType <- as.character(ciber_mean$CellType) %>% unique
cellType[cellType == "T.cells.regulatory..Tregs."] <- "T.cells.regulatory"
cellType <- gsub(pattern = ".", replacement = " ", x = cellType, fixed = TRUE)
ciber_mean$CellType <- cellType # Will repeat to fill in
ciber_mean <- subset(ciber_mean, !CellType %in% c("T cells follicular helper", "Dendritic cells resting"))
ciber_mean$CellType <- factor(ciber_mean$CellType, levels = unique(ciber_mean$CellType))
ciber_meta <- aggregate(Fraction ~ CellType + Diagnosis, data = ciber_combined, FUN = mean)
cellType <- as.character(ciber_meta$CellType) %>% unique
cellType[cellType == "T.cells.regulatory..Tregs."] <- "T.cells.regulatory"
cellType <- gsub(pattern = ".", replacement = " ", x = cellType, fixed = TRUE)
ciber_meta$CellType <- cellType # Will repeat to fill in
ciber_meta <- subset(ciber_meta, !CellType %in% c("T cells follicular helper", "Dendritic cells resting"))
ciber_meta$CellType <- factor(ciber_meta$CellType, levels = unique(ciber_meta$CellType))
ciber_meta$Study <- "Meta"
ciber_meta <- rbind(ciber_mean, ciber_meta[, c("CellType", "Diagnosis", "Study", "Fraction")])
ciber_meta$Diagnosis <- as.character(ciber_meta$Diagnosis)
ciber_meta$Diagnosis[ciber_meta$Diagnosis == "NonTD"] <- "Non-TD"
ciber_meta$Diagnosis <- factor(ciber_meta$Diagnosis, levels = c("TD", "Non-TD", "ASD"))
gg <- ggplot(ciber_meta, aes(x = Diagnosis, y = Fraction*100, fill = CellType, color = CellType))
gg + 
        geom_bar(stat="identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(1.27, 0.5), legend.background = element_blank(), 
              legend.key.size = unit(0.5, "cm"), strip.text.x = element_text(size = 16), 
              axis.ticks = element_line(size = 1.1, color = "black"), 
              legend.text = element_text(size = 12, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "vertical", 
              panel.spacing.y = unit(0, "lines"), axis.title.y = element_text(size = 17, color = "black"), 
              plot.margin = unit(c(0.25, 14, 0.5, 0.75), "lines"), axis.title.x = element_blank(), 
              axis.text = element_text(size = 14, color = "black"), legend.title = element_blank(),
              axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) +
        ylab("Cell Populations (%)") +
        scale_fill_manual(values = wes_palette(n = 22, name = "FantasticFox1", type = "continuous"), 
                          guide = guide_legend(ncol = 1)) +
        scale_color_manual(values = wes_palette(n = 22, name = "FantasticFox1", type = "continuous"),
                           guide = guide_legend(ncol = 1)) +
        scale_y_continuous(expand = c(0.005, 0)) +
        facet_grid(cols = vars(Study), scales = "free")
ggsave("Figures/CIBERSORT Mean Relative Fractions by Diagnosis and Study with Meta Stacked Barplot.png", dpi = 600, width = 9, 
       height = 5, units = "in")

# Correlate fraction of each cell type with covariates ####
# Prep cell type data
table(ciber_marbles$Input.Sample == rownames(pheno$marbles$data)) # All TRUE
ciber_marbles <- ciber_marbles[, !colnames(ciber_marbles) == "Input.Sample"]
rownames(ciber_marbles) <- rownames(pheno$marbles$data)
table(ciber_earli$Input.Sample == rownames(pheno$earli$data)) # All TRUE
ciber_earli <- ciber_earli[, !colnames(ciber_earli) == "Input.Sample"]
rownames(ciber_earli) <- rownames(pheno$earli$data)
ciber <- list(marbles = list(data = as.matrix(ciber_marbles)), earli = list(data = as.matrix(ciber_earli)))

# Get Meta-Analysis Correlations
metaCor <- list()
for (t in 1:ncol(pheno$marbles$data)){
        metaCor[[t]] = metaAnalysis(ciber, mtd.subset(pheno, colIndex = t), useRankPvalue = FALSE, corFnc = bicor,
                                  corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs"), 
                                  getQvalues = TRUE)
}
names(metaCor) <- colnames(pheno$marbles$data)
metaStats <- cbind(colnames(ciber$marbles$data),
                   sapply(metaCor, function(x) x[, "Z.marbles"]),
                   sapply(metaCor, function(x) x[, "Z.earli"]),
                   sapply(metaCor, function(x) x[, "Z.RootDoFWeights"]),
                   sapply(metaCor, function(x) x[, "p.RootDoFWeights"]),
                   sapply(metaCor, function(x) x[, "q.RootDoFWeights"])) %>% as.data.frame
colnames(metaStats) <- c("CellType", paste(rep(c("Zscore_MARBLES", "Zscore_EARLI", "Zscore_Meta", "pvalue_Meta", 
                                                 "qvalue_Meta"), each = length(colnames(pheno$marbles$data))),
                                           rep(colnames(pheno$marbles$data), 5), sep = "_"))
write.table(metaStats, "tables/Revision CIBERSORT Meta Covariate Correlation Stats.txt", sep = "\t", quote = FALSE)

zscores <- sapply(metaCor, function(x) x[["Z.RootDoFWeights"]])
rownames(zscores) <- colnames(ciber$marbles$data)
pvalues <- sapply(metaCor, function(x) x[["p.RootDoFWeights"]])
rownames(pvalues) <- colnames(ciber$marbles$data)
qvalues <- sapply(metaCor, function(x) x[["q.RootDoFWeights"]])
rownames(qvalues) <- colnames(ciber$marbles$data)
zscores <- zscores[!rownames(zscores) %in% c("T.cells.follicular.helper", "Dendritic.cells.resting"),] # Zero fraction in both
cellType <- rownames(zscores)
cellType[cellType == "T.cells.regulatory..Tregs."] <- "T.cells.regulatory"
cellType <- gsub(pattern = ".", replacement = " ", x = cellType, fixed = TRUE)
rownames(zscores) <- cellType
qvalues <- qvalues[!rownames(qvalues) %in% c("T.cells.follicular.helper", "Dendritic.cells.resting"),] # Zero fraction in both

# Plot All Correlations for Meta Analysis (Z-scores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
xLabels <- c("Diagnosis", "ASD vs TD", "Non-TD vs TD", "ASD vs Non-TD", "Sex", "Batch", "Gestational Age", "Birthweight",
             "Hispanic Ethnicity", "Race", "Race and Ethnicity", "Mother Age", "Father Age", "Mother BMI", 
             "Cesarean Delivery", "Urine Cotinine", "Mother Smoking", "Mother Bachelor's", "Father Bachelor's",
             "Own Home")
pdf("Figures/Revisions CIBERSORT Meta Covariate Correlation Plot.pdf", width = 11, height = 7)
sizeGrWindow(width = 11, height = 7)
par(mar = c(8, 14, 1, 2))
labeledHeatmap(Matrix = zscores, xLabels = xLabels, yLabels = rownames(zscores), 
               ySymbols = rownames(zscores), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2.5, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab = 1.2)
dev.off()
rm(ciber_earli, ciber_marbles, metaCor, metaStats, pheno, pvalues, qvalues, star, zscores, t)

# Correlate Cell Type Fractions with Module Eigengenes ####
# Get MEs
exp_marbles <- readRDS("R Objects/MARBLES_Exp_BatchAdj.rds")
exp_earli <- readRDS("R Objects/EARLI_Exp_noBatchAdj.rds")
table(colnames(exp_marbles) == colnames(exp_earli)) # All TRUE
table(rownames(exp_marbles) == rownames(ciber$marbles$data)) # All TRUE
table(rownames(exp_earli) == rownames(ciber$earli$data)) # All TRUE
moduleMembership <- read.delim("Tables/Revision Consensus Modules MARBLES Probe Module Membership.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE)
consensusMods <- moduleMembership$Module
MEs_marbles <- moduleEigengenes(exp_marbles, colors = consensusMods)$eigengenes
rownames(MEs_marbles) <- rownames(exp_marbles)
MEs_earli <- moduleEigengenes(exp_earli, colors = consensusMods)$eigengenes
rownames(MEs_earli) <- rownames(exp_earli)
consensusMEs <- list(marbles = list(data = MEs_marbles), earli = list(data = MEs_earli))
consensusMEs <- orderMEs(consensusMEs)

# Get Meta-Analysis Correlations
metaCor <- list()
for (t in 1:ncol(consensusMEs$marbles$data)){
        metaCor[[t]] = metaAnalysis(ciber, mtd.subset(consensusMEs, colIndex = t), useRankPvalue = FALSE, corFnc = bicor,
                                    corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs"), 
                                    getQvalues = TRUE)
}
names(metaCor) <- colnames(consensusMEs$marbles$data)
metaStats <- cbind(colnames(ciber$marbles$data),
                   sapply(metaCor, function(x) x[, "Z.RootDoFWeights"]),
                   sapply(metaCor, function(x) x[, "p.RootDoFWeights"]),
                   sapply(metaCor, function(x) x[, "q.RootDoFWeights"])) %>% as.data.frame
colnames(metaStats) <- c("CellType", paste(rep(c("Zscore_Meta", "pvalue_Meta", "qvalue_Meta"), 
                                               each = length(colnames(consensusMEs$marbles$data))),
                                           rep(colnames(consensusMEs$marbles$data), 3), sep = "_"))
write.table(metaStats, "tables/Revision CIBERSORT Meta Module Eigengene Correlation Stats.txt", sep = "\t", quote = FALSE)

zscores <- sapply(metaCor, function(x) x[["Z.RootDoFWeights"]])
rownames(zscores) <- colnames(ciber$marbles$data)
pvalues <- sapply(metaCor, function(x) x[["p.RootDoFWeights"]])
rownames(pvalues) <- colnames(ciber$marbles$data)
qvalues <- sapply(metaCor, function(x) x[["q.RootDoFWeights"]])
rownames(qvalues) <- colnames(ciber$marbles$data)
zscores <- zscores[!rownames(zscores) %in% c("T.cells.follicular.helper", "Dendritic.cells.resting"),] # Zero fraction in both
cellType <- rownames(zscores)
cellType[cellType == "T.cells.regulatory..Tregs."] <- "T.cells.regulatory"
cellType <- gsub(pattern = ".", replacement = " ", x = cellType, fixed = TRUE)
rownames(zscores) <- cellType
zscores <- t(zscores)
qvalues <- qvalues[!rownames(qvalues) %in% c("T.cells.follicular.helper", "Dendritic.cells.resting"),] # Zero fraction in both
qvalues <- t(qvalues)

# Plot All Correlations for Meta Analysis (Z-scores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Figures/Revisions CIBERSORT Meta Module Eigengene Correlation Plot.pdf", width = 11, height = 15)
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
star <- apply(qvalues_sub, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Figures/Revisions CIBERSORT Meta Module Eigengene Correlation Plot Significant Only.pdf", width = 9, height = 13)
sizeGrWindow(width = 9, height = 13)
par(mar = c(9, 8, 1, 1))
labeledHeatmap(Matrix = zscores_sub, xLabels = colnames(zscores_sub), yLabels = rownames(zscores_sub), 
               ySymbols = gsub("ME", "", rownames(zscores_sub)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Grey60 Module and B-cells Plots ####
table(rownames(ciber$marbles$data) == rownames(consensusMEs$marbles$data) &
      rownames(ciber$marbles$data) == rownames(pheno$marbles$data) &
      rownames(consensusMEs$marbles$data) == rownames(pheno$marbles$data)) # All TRUE
table(rownames(ciber$earli$data) == rownames(consensusMEs$earli$data) &
              rownames(ciber$earli$data) == rownames(pheno$earli$data) &
              rownames(consensusMEs$earli$data) == rownames(pheno$earli$data)) # All TRUE
grey60_Bcell <- data.frame(sampleID = c(rownames(pheno$marbles$data), rownames(pheno$earli$data)),
                           Diagnosis = c(pheno$marbles$data[, "Dx_alg"], pheno$earli$data[, "Dx_alg"]),
                           DeliveryMethod = c(pheno$marbles$data[, "DeliveryMethod"], pheno$earli$data[, "DeliveryMethod"]),
                           Study = c(rep("MARBLES", nrow(pheno$marbles$data)), rep("EARLI", nrow(pheno$earli$data))),
                           MEgrey60 = c(consensusMEs$marbles$data$MEgrey60, consensusMEs$earli$data$MEgrey60),
                           Bcells = c(ciber$marbles$data[, "B.cells.naive"], ciber$earli$data[, "B.cells.naive"]))
grey60_Bcell$Diagnosis[grey60_Bcell$Diagnosis == 1] <- "TD"
grey60_Bcell$Diagnosis[grey60_Bcell$Diagnosis == 2] <- "NonTD"
grey60_Bcell$Diagnosis[grey60_Bcell$Diagnosis == 3] <- "ASD"
grey60_Bcell$DeliveryMethod[grey60_Bcell$DeliveryMethod == 1] <- "Vaginal"
grey60_Bcell$DeliveryMethod[grey60_Bcell$DeliveryMethod == 2] <- "Cesarean"
grey60_Bcell$Diagnosis <- factor(grey60_Bcell$Diagnosis, levels = c("TD", "NonTD", "ASD"))
grey60_Bcell$DeliveryMethod <- factor(grey60_Bcell$DeliveryMethod, levels = c("Vaginal", "Cesarean"))
grey60_Bcell$Study <- factor(grey60_Bcell$Study, levels = c("MARBLES", "EARLI"))
gg <- ggplot()
gg + 
        geom_point(data = grey60_Bcell, aes(MEgrey60, Bcells*100), size = 2.5, color = "#3366CC") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.87, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_text(size = 22),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0, 1, 1, 0.5), "lines"),
              axis.text = element_text(size = 22, color = "black")) +
        xlab("grey60 Module Eigengene") +
        ylab("Naive B Cells (%)") +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        facet_wrap(facets = vars(Study), nrow = 2)
ggsave("Figures/Revisions grey60 ME B cells Scatterplot.png", dpi = 600, height = 9, width = 8)

grey60_Bcell <- subset(grey60_Bcell, !is.na(grey60_Bcell$DeliveryMethod))
gg <- ggplot(data = grey60_Bcell)
gg <- gg +
        geom_boxplot(aes(x = DeliveryMethod, y = Bcells*100, fill = Diagnosis), size = 0.8, outlier.size = 0.8) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.87, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_text(size = 22),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0, 1, 1, 0.5), "lines"), axis.title.x = element_blank(), 
              axis.text = element_text(size = 22, color = "black")) +
        ylab("Naive B Cells (%)") +
        scale_fill_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        scale_color_manual(name = "Diagnosis", values = c("#3366CC", "#009933", "#FF3366")) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        coord_cartesian(ylim = c(0, 40)) +
        facet_wrap(facets = vars(Study), nrow = 2)
ggsave(filename = "Figures/Naive B Cells by Diagnosis and Delivery Method Boxplot.png", plot = gg, dpi = 600, width = 8, height = 9, units = "in")
