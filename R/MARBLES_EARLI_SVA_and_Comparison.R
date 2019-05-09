# MARBLES EARLI SVA and Comparison Analysis --------------------------------------
# Cord Expression Revisions
# Charles Mordaunt
# 4/21/19

# Packages ####
sapply(c("tidyverse", "scales", "variancePartition", "sm", "biomaRt", "reshape2", "wesanderson", "WGCNA", "matrixStats"), 
       require, character.only = TRUE)

# SVA-Covariate Association and Variance Partition Analysis -----------------------------------------------
# SVA-Covariate Association for MARBLES ####
# Data
cov <- read.delim("Tables/Revisions MARBLES Microarray Cord Sample Covariates Matching EARLI for Table 1.txt", 
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pheno <- read.delim("/Users/charles/Documents/Programming/MARBLES Cord Blood/Microarray/Tables/Pheno Data for MARBLES Dx_alg Samples B123.txt",
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sv <- read.delim("/Users/charles/Documents/Programming/MARBLES Cord Blood/Microarray/Tables/MARBLES Dx_alg with SV Design Matrix B123.txt",
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# SV Table
sv$Dx_alg <- "TD"
sv$Dx_alg[sv$Dx_alg_ASD == 1] <- "ASD"
sv$Dx_alg[sv$Dx_alg_NonTD == 1] <- "NonTD"
table(sv$Dx_alg == pheno$Dx_alg) # All TRUE
sv$IBC <- pheno$IBC
sv <- sv[match(cov$IBC, sv$IBC), ]
table(sv$IBC == cov$IBC) # All TRUE
rownames(sv) <- sv$IBC
sv <- sv[, grepl("Sv", colnames(sv), fixed = TRUE)] %>% as.matrix
colnames(sv) <- paste("SV", 1:ncol(sv), sep = "")
rm(pheno)

# Cov Table
colnames(cov)[colnames(cov) == "COI_GENDER"] <- "Sex"
cov$Dx_alg[cov$Dx_alg == "Non-TD"] <- "NonTD"
cov$ASDvsTD <- sapply(cov$Dx_alg, function(x){ifelse(x == "NonTD", NA, x)})
cov$NonTDvsTD <- sapply(cov$Dx_alg, function(x){ifelse(x == "ASD", NA, x)})
cov$ASDvsNonTD <- sapply(cov$Dx_alg, function(x){ifelse(x == "TD", NA, x)})
cov <- cov[,c("IBC", "Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", "gest_age_deliv_wk", 
                  "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth", "MomAgeYr", "DadAgeYr", 
                  "Mat_BMI_PrePreg", "DeliveryMethod", "DM1_EQ", "DM2_EQ", "GDM_EQ", "HTN_EQ", "PE_EQ", 
                  "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", "MomEduBachelors", "DadEduBachelors", "HomeOwn")]
cov <- sapply(cov, as.character) %>% as.data.frame(stringsAsFactors = FALSE)  # All columns are character vectors
contCols <- c("IBC", "Batch", "gest_age_deliv_wk", "birthweightkg", "ChildHisp", "MomAgeYr", "DadAgeYr", "Mat_BMI_PrePreg", 
              "Cotinine_Urine_Conc_ng_ml")
cov[,contCols] <- lapply(cov[,contCols], as.numeric)
cov$Dx_alg <- factor(cov$Dx_alg, levels = c("TD", "NonTD", "ASD")) %>% as.numeric
cov$ASDvsTD <- factor(cov$ASDvsTD, levels = c("TD", "ASD")) %>% as.numeric
cov$NonTDvsTD <- factor(cov$NonTDvsTD, levels = c("TD", "NonTD")) %>% as.numeric
cov$ASDvsNonTD <- factor(cov$ASDvsNonTD, levels = c("NonTD", "ASD")) %>% as.numeric
cov$Sex <- factor(cov$Sex, levels = c("M", "F")) %>% as.numeric
cov$ChildRace <- factor(cov$ChildRace, levels = c("1", "2", "4", "5", "9")) %>% as.numeric
cov$ChildRaceEth <- factor(cov$ChildRaceEth, levels = c("1", "2", "4", "5", "61", "62", "71", "72", "9")) %>% as.numeric
cov$DeliveryMethod <- factor(cov$DeliveryMethod, levels = c("Vaginal", "Cesarean")) %>% as.numeric
cov$DM1_EQ <- factor(cov$DM1_EQ, levels = c("N", "Y")) %>% as.numeric
cov$DM2_EQ <- factor(cov$DM2_EQ, levels = c("N", "Y")) %>% as.numeric
cov$GDM_EQ <- factor(cov$GDM_EQ, levels = c("N", "Y")) %>% as.numeric
cov$HTN_EQ <- factor(cov$HTN_EQ, levels = c("N", "Y")) %>% as.numeric
cov$PE_EQ <- factor(cov$PE_EQ, levels = c("N", "Y")) %>% as.numeric
cov$Cotinine_Urine_Smoker <- factor(cov$Cotinine_Urine_Smoker, levels = c("FALSE", "TRUE")) %>% as.numeric
cov$MomEduBachelors <- factor(cov$MomEduBachelors, levels = c("FALSE", "TRUE")) %>% as.numeric
cov$DadEduBachelors <- factor(cov$DadEduBachelors, levels = c("FALSE", "TRUE")) %>% as.numeric
cov$HomeOwn <- factor(cov$HomeOwn, levels = c("Renter", "Owner")) %>% as.numeric
rownames(cov) <- cov$IBC
cov <- cov[,c("Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", "gest_age_deliv_wk", 
                  "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth", "MomAgeYr", "DadAgeYr", 
                  "Mat_BMI_PrePreg", "DeliveryMethod", "DM1_EQ", "DM2_EQ", "GDM_EQ", "HTN_EQ", "PE_EQ", 
                  "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", "MomEduBachelors", "DadEduBachelors", "HomeOwn")]
cov <- as.matrix(cov)
rm(contCols)

# Associate Surrogate variables with Covariates
assoc <- NULL
for(i in 1:ncol(sv)){
        for(j in 1:ncol(cov)){
                temp <- c(colnames(sv)[i], colnames(cov)[j], summary(lm(cov[, j] ~ sv[, i]))$coefficients[2, ])
                assoc <- rbind(assoc, temp)
        }
}
rm(i, j, temp)
assoc <- as.data.frame(assoc, row.names = 1:nrow(assoc), stringsAsFactors = FALSE) 
colnames(assoc) <- c("SV", "Covariate", "Estimate", "StdError", "tvalue", "pvalue")
assoc$SV <- factor(assoc$SV, levels = unique(assoc$SV) %>% rev)
assoc$Covariate <- factor(assoc$Covariate, levels = unique(assoc$Covariate))
assoc[, c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(assoc[, c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)
assoc$qvalue <- p.adjust(assoc$pvalue, method = "fdr")
assoc$log_pvalue <- -log10(assoc$pvalue)
assoc$log_qvalue <- -log10(assoc$qvalue)
write.table(assoc, "Tables/Revisions MARBLES SVA-Covariate Association Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Make Plots of SV-Covariate Associations
gg <- ggplot(data = assoc)
gg +
        geom_tile(aes(x = Covariate, y = SV, fill = log_pvalue)) +
        scale_fill_gradientn("-log(p-value)", colors = c("Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                             limits = c(0, 4), breaks = pretty_breaks(n = 4)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = c(1.12, 0.83), legend.background = element_blank(), 
              plot.margin = unit(c(1, 8, 1, 1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 16, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title.y = element_text(size = 18), legend.title = element_text(size = 18),
              legend.text = element_text(size = 16), axis.title.x = element_blank()) +
        ylab("Surrogate Variable")
ggsave("Figures/Revisions MARBLES SVA-Covariate Association log pvalue Heatmap.png", dpi = 600, width = 10, height = 7, 
       units = "in")
rm(gg)

# Variance Partition Analysis for MARBLES ####
# Expression Data
sum <- readRDS("R Objects/child_sumData_noOutliers_filtered_B123.rds")
pheno <- pData(sum)
sampleNames(sum) <- pheno$IBC
exp <- exprs(sum)
exp <- exp[, match(rownames(cov), colnames(exp))]
table(colnames(exp) == rownames(cov)) # All TRUE
table(colnames(exp) == rownames(sv)) # All TRUE
rm(sum, pheno)

cov <- as.data.frame(cov, stringsAsFactors = FALSE)
sv <- as.data.frame(sv, stringsAsFactors = FALSE)

# Fit Model for SVs and Plot
varPart <- fitExtractVarPartModel(exp, formula = ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 +
                                          SV12 + SV13 + SV14 + SV15 + SV16 + SV17 + SV18 + SV19 + SV20 + SV21, data = sv)
varPart <- sortCols(varPart)
gg <- plotVarPart(varPart, label.angle = 90)
gg <- gg +
        scale_y_continuous(expand = c(0.01, 0), breaks = pretty_breaks(n = 5)) +
        coord_cartesian(ylim = c(0, 100)) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = "none",
              plot.margin = unit(c(1, 1, 1, 1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 14, color = "Black", hjust = 1, vjust = 0.5),
              axis.title.y = element_text(size = 16), axis.title.x = element_blank())
ggsave("Figures/Revisions MARBLES SVA Variance Partition Violin Plots.png", plot = gg, dpi = 600, width = 7, height = 5, 
       units = "in")
rm(gg, varPart)

# SVA-Covariate Association for EARLI ####
# Data
pheno <- read.csv("Tables/Revisions EARLI Array Sample Info.csv", header = TRUE, stringsAsFactors = FALSE)
cov <- read.delim("Tables/Revisions EARLI Microarray Cord Sample Covariates.txt", sep = "\t", header = TRUE, 
                  stringsAsFactors = FALSE)
race <- read.delim("Tables/EARLI Child Race and Ethnicity MARBLES Coding Array Samples.txt", sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)
load("R Objects/EARLI_SVs_ASD_031919.rda")
svASD <- as.data.frame(svobj$sv)
rownames(svASD) <- rownames(mod)
colnames(svASD) <- paste("SV", 1:ncol(svASD), sep = "")
rm(mod, svobj)
load("R Objects/EARLI_SVs_NonTD_031919.rda")
svNonTD <- as.data.frame(svobj$sv)
rownames(svNonTD) <- rownames(mod)
colnames(svNonTD) <- paste("SV", 1:ncol(svNonTD), sep = "")
rm(mod, svobj)

# Covariate Table
cov <- merge(x = cov, y = race[,c("Sample.ID", "ChildHisp", "ChildRace", "ChildRaceEth")], by.x = "sampleID", 
             by.y = "Sample.ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
cov <- merge(x = pheno, y = cov, by = "sampleID", all.x = TRUE, all.y = FALSE, sort = FALSE)
rm(pheno, race)
table(cov$sex == cov$COI_GENDER) # All TRUE
cov$Dx_alg[cov$Dx_alg == "Non-TD"] <- "NonTD"
table(cov$diagnosis == cov$Dx_alg) # All TRUE
table(cov$batch == cov$Batch) # All TRUE
cov$Cotinine_Urine_Smoker <- cov$Cotinine_Urine_Conc_ngml > 50
cov$svID <- strsplit(cov$CELfile, split = "-") %>% sapply(., function(x) x[2]) %>% gsub("_", replacement = "", x = .)
cov <- cov[,c("sampleID", "svID", "diagnosis", "sex", "batch", "gest_age_deliv_wk", "birthweightkg", "ChildHisp", "ChildRace", 
                  "ChildRaceEth", "MomAgeYr", "DadAgeYr", "Mat_BMI_PrePreg", "DeliveryMethod", "Cotinine_Urine_Conc_ngml", 
                  "Cotinine_Urine_Smoker", "MomEduBachelors", "DadEduBachelors", "home_ownership")]
colnames(cov)[colnames(cov) == "diagnosis"] <- "Dx_alg"
colnames(cov)[colnames(cov) == "sex"] <- "Sex"
colnames(cov)[colnames(cov) == "batch"] <- "Batch"
colnames(cov)[colnames(cov) == "Cotinine_Urine_Conc_ngml"] <- "Cotinine_Urine_Conc_ng_ml"
colnames(cov)[colnames(cov) == "home_ownership"] <- "HomeOwn"
cov$ASDvsTD <- sapply(cov$Dx_alg, function(x){ifelse(x == "NonTD", NA, x)})
cov$NonTDvsTD <- sapply(cov$Dx_alg, function(x){ifelse(x == "ASD", NA, x)})
cov$ASDvsNonTD <- sapply(cov$Dx_alg, function(x){ifelse(x == "TD", NA, x)})
cov <- cov[,c("sampleID", "svID", "Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", "gest_age_deliv_wk", 
                  "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth", "MomAgeYr", "DadAgeYr", 
                  "Mat_BMI_PrePreg", "DeliveryMethod", "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", 
                  "MomEduBachelors", "DadEduBachelors", "HomeOwn")]
cov <- sapply(cov, as.character) %>% as.data.frame(stringsAsFactors = FALSE) # All columns are character vectors
contCols <- c("sampleID", "Batch", "gest_age_deliv_wk", "birthweightkg", "ChildHisp", "MomAgeYr", "DadAgeYr", "Mat_BMI_PrePreg", 
              "Cotinine_Urine_Conc_ng_ml")
cov[,contCols] <- lapply(cov[,contCols], as.numeric)
cov$Dx_alg <- factor(cov$Dx_alg, levels = c("TD", "NonTD", "ASD")) %>% as.numeric
cov$ASDvsTD <- factor(cov$ASDvsTD, levels = c("TD", "ASD")) %>% as.numeric
cov$NonTDvsTD <- factor(cov$NonTDvsTD, levels = c("TD", "NonTD")) %>% as.numeric
cov$ASDvsNonTD <- factor(cov$ASDvsNonTD, levels = c("NonTD", "ASD")) %>% as.numeric
cov$Sex <- factor(cov$Sex, levels = c("M", "F")) %>% as.numeric
cov$ChildRace <- factor(cov$ChildRace, levels = c("1", "2", "4", "9")) %>% as.numeric
cov$ChildRaceEth <- factor(cov$ChildRaceEth, levels = c("1", "2", "4", "61", "62", "71", "72", "9")) %>% as.numeric
cov$DeliveryMethod <- factor(cov$DeliveryMethod, levels = c("Vaginal", "Cesarean")) %>% as.numeric
cov$Cotinine_Urine_Smoker <- factor(cov$Cotinine_Urine_Smoker, levels = c("FALSE", "TRUE")) %>% as.numeric
cov$MomEduBachelors <- factor(cov$MomEduBachelors, levels = c("FALSE", "TRUE")) %>% as.numeric
cov$DadEduBachelors <- factor(cov$DadEduBachelors, levels = c("FALSE", "TRUE")) %>% as.numeric
cov$HomeOwn <- factor(cov$HomeOwn, levels = c("Renter", "Owner")) %>% as.numeric
rownames(cov) <- cov$svID
covASD <- cov[match(rownames(svASD), rownames(cov)), c("Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", "gest_age_deliv_wk", 
                                                       "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth", "MomAgeYr", "DadAgeYr", 
                                                       "Mat_BMI_PrePreg", "DeliveryMethod", "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", 
                                                       "MomEduBachelors", "DadEduBachelors", "HomeOwn")] %>% as.matrix
table(rownames(covASD) == rownames(svASD)) # All TRUE
covNonTD <- cov[match(rownames(svNonTD), rownames(cov)), c("Dx_alg", "ASDvsTD", "NonTDvsTD", "ASDvsNonTD", "Sex", "Batch", "gest_age_deliv_wk", 
                                                           "birthweightkg", "ChildHisp", "ChildRace", "ChildRaceEth", "MomAgeYr", "DadAgeYr", 
                                                           "Mat_BMI_PrePreg", "DeliveryMethod", "Cotinine_Urine_Conc_ng_ml", "Cotinine_Urine_Smoker", 
                                                           "MomEduBachelors", "DadEduBachelors", "HomeOwn")] %>% as.matrix
table(rownames(covNonTD) == rownames(svNonTD)) # All TRUE
rm(contCols)

# Associate Surrogate variables with Covariates for ASD
assocASD <- NULL
for(i in 1:ncol(svASD)){
        for(j in 1:ncol(covASD)){
                temp <- c(colnames(svASD)[i], colnames(covASD)[j], summary(lm(covASD[, j] ~ svASD[, i]))$coefficients[2, ])
                assocASD <- rbind(assocASD, temp)
        }
}
rm(i, j, temp)
assocASD <- as.data.frame(assocASD, row.names = 1:nrow(assocASD), stringsAsFactors = FALSE) 
colnames(assocASD) <- c("SV", "Covariate", "Estimate", "StdError", "tvalue", "pvalue")
assocASD$SV <- factor(assocASD$SV, levels = unique(assocASD$SV) %>% rev)
assocASD$Covariate <- factor(assocASD$Covariate, levels = unique(assocASD$Covariate))
assocASD[, c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(assocASD[, c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)
assocASD$qvalue <- p.adjust(assocASD$pvalue, method = "fdr")
assocASD$log_pvalue <- -log10(assocASD$pvalue)
assocASD$log_qvalue <- -log10(assocASD$qvalue)
write.table(assocASD, "Tables/Revisions EARLI SVA-Covariate Association Stats ASD.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Associate Surrogate variables with Covariates for NonTD
assocNonTD <- NULL
for(i in 1:ncol(svNonTD)){
        for(j in 1:ncol(covNonTD)){
                temp <- c(colnames(svNonTD)[i], colnames(covNonTD)[j], summary(lm(covNonTD[, j] ~ svNonTD[, i]))$coefficients[2, ])
                assocNonTD <- rbind(assocNonTD, temp)
        }
}
rm(i, j, temp)
assocNonTD <- as.data.frame(assocNonTD, row.names = 1:nrow(assocNonTD), stringsAsFactors = FALSE) 
colnames(assocNonTD) <- c("SV", "Covariate", "Estimate", "StdError", "tvalue", "pvalue")
assocNonTD$SV <- factor(assocNonTD$SV, levels = unique(assocNonTD$SV) %>% rev)
assocNonTD$Covariate <- factor(assocNonTD$Covariate, levels = unique(assocNonTD$Covariate))
assocNonTD[, c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(assocNonTD[, c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)
assocNonTD$qvalue <- p.adjust(assocNonTD$pvalue, method = "fdr")
assocNonTD$log_pvalue <- -log10(assocNonTD$pvalue)
assocNonTD$log_qvalue <- -log10(assocNonTD$qvalue)
write.table(assocNonTD, "Tables/Revisions EARLI SVA-Covariate Association Stats NonTD.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Make Plots of SV-Covariate Associations for ASD
gg <- ggplot(data = assocASD)
gg +
        geom_tile(aes(x = Covariate, y = SV, fill = log_pvalue)) +
        scale_fill_gradientn("-log(p-value)", colors = c("Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                             limits = c(0, 4), breaks = pretty_breaks(n = 4)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = c(1.12, 0.83), legend.background = element_blank(), 
              plot.margin = unit(c(1, 8, 1, 1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 16, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title.y = element_text(size = 18), legend.title = element_text(size = 18),
              legend.text = element_text(size = 16), axis.title.x = element_blank()) +
        ylab("Surrogate Variable")
ggsave("Figures/Revisions EARLI SVA-Covariate Association log pvalue Heatmap ASD.png", dpi = 600, width = 10, height = 7, 
       units = "in")

# Make Plots of SV-Covariate Associations for NonTD
gg <- ggplot(data = assocNonTD)
gg +
        geom_tile(aes(x = Covariate, y = SV, fill = log_pvalue)) +
        scale_fill_gradientn("-log(p-value)", colors = c("Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                             limits = c(0, 4), breaks = pretty_breaks(n = 4)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = c(1.12, 0.83), legend.background = element_blank(), 
              plot.margin = unit(c(1, 8, 1, 1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 16, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title.y = element_text(size = 18), legend.title = element_text(size = 18),
              legend.text = element_text(size = 16), axis.title.x = element_blank()) +
        ylab("Surrogate Variable")
ggsave("Figures/Revisions EARLI SVA-Covariate Association log pvalue Heatmap NonTD.png", dpi = 600, width = 10, height = 7, 
       units = "in")
rm(gg, assocASD, assocNonTD)

# Variance Partition Analysis for EARLI ####
# Expression Data
exp <- readRDS("R Objects/EARLI_Exp_noBatchAdj_DiffExp.rds")
exp <- t(exp) %>% as.data.frame
table(cov$sampleID == colnames(exp)) # All TRUE
colnames(exp) <- cov$svID
expASD <- exp[, match(rownames(covASD), colnames(exp))]
table(colnames(expASD) == rownames(covASD)) # All TRUE
table(colnames(expASD) == rownames(svASD)) # All TRUE
expNonTD <- exp[, match(rownames(covNonTD), colnames(exp))]
table(colnames(expNonTD) == rownames(covNonTD)) # All TRUE
table(colnames(expNonTD) == rownames(svNonTD)) # All TRUE
covASD <- as.data.frame(covASD, stringsAsFactors = FALSE)
covNonTD <- as.data.frame(covNonTD, stringsAsFactors = FALSE)

# Fit Model for SVs and Plot for ASD
varPartASD <- fitExtractVarPartModel(expASD, formula = ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11, 
                                     data = svASD) %>% sortCols
gg <- plotVarPart(varPartASD, label.angle = 90)
gg <- gg +
        scale_y_continuous(expand = c(0.01, 0), breaks = pretty_breaks(n = 5)) +
        coord_cartesian(ylim = c(0, 100)) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = "none",
              plot.margin = unit(c(1, 1, 1, 1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 14, color = "Black", hjust = 1, vjust = 0.5),
              axis.title.y = element_text(size = 16), axis.title.x = element_blank())
ggsave("Figures/Revisions EARLI SVA Variance Partition Violin Plots ASD.png", plot = gg, dpi = 600, width = 7, height = 5, 
       units = "in")

# Fit Model for SVs and Plot for NonTD
varPartNonTD <- fitExtractVarPartModel(expNonTD, formula = ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12, 
                                     data = svNonTD) %>% sortCols
gg <- plotVarPart(varPartNonTD, label.angle = 90)
gg <- gg +
        scale_y_continuous(expand = c(0.01, 0), breaks = pretty_breaks(n = 5)) +
        coord_cartesian(ylim = c(0, 100)) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = "none",
              plot.margin = unit(c(1, 1, 1, 1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 14, color = "Black", hjust = 1, vjust = 0.5),
              axis.title.y = element_text(size = 16), axis.title.x = element_blank())
ggsave("Figures/Revisions EARLI SVA Variance Partition Violin Plots NonTD.png", plot = gg, dpi = 600, width = 7, height = 5, 
       units = "in")

# MARBLES EARLI Fold Change Scatterplots ------------------------------------------------
# Data ####
ASD <- read.delim("Tables/ASD Meta-Analysis with METAL All Genes B123.txt", sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
NonTD <- read.delim("Tables/NonTD Meta-Analysis with METAL All Genes B123.txt", sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
table(ASD$Probe == NonTD$Probe) # All TRUE

# Meta ASD vs NonTD ####
ASDvsNonTD <- data.frame(Probe = ASD$Probe, ASD_logFC = ASD$Meta_logFC, NonTD_logFC = NonTD$Meta_logFC)
ASDvsNonTD$Density <- sm.density(ASDvsNonTD[, c("ASD_logFC", "NonTD_logFC")], 
                                 eval.points = ASDvsNonTD[, c("ASD_logFC", "NonTD_logFC")], display = "none", 
                                 eval.grid = FALSE, panel = FALSE)$estimate
test <- cor.test(ASDvsNonTD$ASD_logFC, ASDvsNonTD$NonTD_logFC)
r <- round(test$estimate, 3)
p <- test$p.value #0
p <- 2.2e-16
gg <- ggplot()
gg + 
        geom_point(data = ASDvsNonTD, aes(ASD_logFC, NonTD_logFC, color = Density), size = 1.4) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "red") +
        annotate("text", x = -0.6, y = 1.2, label = paste("r =", r, sep = " "), size = 7, hjust = 0) +
        annotate("text", x = -0.6, y = 1.05, label = paste("p <", p, sep = " "), size = 7, hjust = 0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.1, color = "black"), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text = element_text(color = "black"), axis.title = element_text(size = 22),
              plot.margin = unit(c(1.5, 2, 1, 1), "lines")) +
        xlab(expr("Meta ASD log"[2]*"(Fold Change)")) +
        ylab(expr("Meta Non-TD log"[2]*"(Fold Change)")) +
        coord_cartesian(xlim = c(-0.6, 1.2), ylim = c(-0.6, 1.2)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5))
ggsave("Figures/Revisions ASD vs NonTD log Fold Change Scatterplot.png", dpi = 600, height = 6, width = 7)

# MARBLES ASD vs NonTD ####
ASDvsNonTD_marbles <- data.frame(Probe = ASD$Probe, ASD_logFC = ASD$marbles_logFC, NonTD_logFC = NonTD$marbles_logFC)
ASDvsNonTD_marbles$Density <- sm.density(ASDvsNonTD_marbles[, c("ASD_logFC", "NonTD_logFC")], 
                                 eval.points = ASDvsNonTD_marbles[, c("ASD_logFC", "NonTD_logFC")], display = "none", 
                                 eval.grid = FALSE, panel = FALSE)$estimate
test <- cor.test(ASDvsNonTD_marbles$ASD_logFC, ASDvsNonTD_marbles$NonTD_logFC)
r <- round(test$estimate, 3)
p <- test$p.value #0
p <- 2.2e-16
gg <- ggplot()
gg + 
        geom_point(data = ASDvsNonTD_marbles, aes(ASD_logFC, NonTD_logFC, color = Density), size = 1.4) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "red") +
        annotate("text", x = -0.6, y = 1.2, label = paste("r =", r, sep = " "), size = 7, hjust = 0) +
        annotate("text", x = -0.6, y = 1.05, label = paste("p <", p, sep = " "), size = 7, hjust = 0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.1, color = "black"), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text = element_text(color = "black"), axis.title = element_text(size = 20),
              plot.margin = unit(c(1.5, 2, 1, 1), "lines")) +
        xlab(expr("MARBLES ASD log"[2]*"(Fold Change)")) +
        ylab(expr("MARBLES Non-TD log"[2]*"(Fold Change)")) +
        coord_cartesian(xlim = c(-0.6, 1.2), ylim = c(-0.6, 1.2)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5))
ggsave("Figures/Revisions MARBLES ASD vs NonTD log Fold Change Scatterplot.png", dpi = 600, height = 6, width = 7)

# EARLI ASD vs NonTD ####
ASDvsNonTD_earli <- data.frame(Probe = ASD$Probe, ASD_logFC = ASD$earli_logFC, NonTD_logFC = NonTD$earli_logFC)
ASDvsNonTD_earli$Density <- sm.density(ASDvsNonTD_earli[, c("ASD_logFC", "NonTD_logFC")], 
                                         eval.points = ASDvsNonTD_earli[, c("ASD_logFC", "NonTD_logFC")], display = "none", 
                                         eval.grid = FALSE, panel = FALSE)$estimate
test <- cor.test(ASDvsNonTD_earli$ASD_logFC, ASDvsNonTD_earli$NonTD_logFC)
r <- round(test$estimate, 3)
p <- test$p.value #0
p <- 2.2e-16
gg <- ggplot()
gg + 
        geom_point(data = ASDvsNonTD_earli, aes(ASD_logFC, NonTD_logFC, color = Density), size = 1.4) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "red") +
        annotate("text", x = -0.6, y = 1.2, label = paste("r =", r, sep = " "), size = 7, hjust = 0) +
        annotate("text", x = -0.6, y = 1.05, label = paste("p <", p, sep = " "), size = 7, hjust = 0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.1, color = "black"), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text = element_text(color = "black"), axis.title = element_text(size = 22),
              plot.margin = unit(c(1.5, 2, 1, 1), "lines")) +
        xlab(expr("EARLI ASD log"[2]*"(Fold Change)")) +
        ylab(expr("EARLI Non-TD log"[2]*"(Fold Change)")) +
        coord_cartesian(xlim = c(-0.6, 1.2), ylim = c(-0.6, 1.2)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5))
ggsave("Figures/Revisions EARLI ASD vs NonTD log Fold Change Scatterplot.png", dpi = 600, height = 6, width = 7)

# ASD MARBLES vs EARLI ####
ASD_MvsE <- data.frame(Probe = ASD$Probe, MARBLES = ASD$marbles_logFC, EARLI = ASD$earli_logFC)
ASD_MvsE$Density <- sm.density(ASD_MvsE[, c("MARBLES", "EARLI")], 
                                 eval.points = ASD_MvsE[, c("MARBLES", "EARLI")], display = "none", 
                                 eval.grid = FALSE, panel = FALSE)$estimate
test <- cor.test(ASD_MvsE$MARBLES, ASD_MvsE$EARLI)
r <- round(test$estimate, 3)
p <- scientific(test$p.value, 2)
gg <- ggplot()
gg + 
        geom_point(data = ASD_MvsE, aes(MARBLES, EARLI, color = Density), size = 1.4) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "red") +
        annotate("text", x = -0.6, y = 1.2, label = paste("r =", r, sep = " "), size = 7, hjust = 0) +
        annotate("text", x = -0.6, y = 1.05, label = paste("p =", p, sep = " "), size = 7, hjust = 0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.1, color = "black"), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text = element_text(color = "black"), axis.title = element_text(size = 22),
              plot.margin = unit(c(1.5, 2, 1, 1), "lines")) +
        xlab(expr("MARBLES ASD log"[2]*"(Fold Change)")) +
        ylab(expr("EARLI ASD log"[2]*"(Fold Change)")) +
        coord_cartesian(xlim = c(-0.6, 1.2), ylim = c(-0.6, 1.2)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5))
ggsave("Figures/Revisions ASD MARBLES vs EARLI log Fold Change Scatterplot.png", dpi = 600, height = 6, width = 7)

# NonTD MARBLES vs EARLI ####
NonTD_MvsE <- data.frame(Probe = NonTD$Probe, MARBLES = NonTD$marbles_logFC, EARLI = NonTD$earli_logFC)
NonTD_MvsE$Density <- sm.density(NonTD_MvsE[, c("MARBLES", "EARLI")], 
                               eval.points = NonTD_MvsE[, c("MARBLES", "EARLI")], display = "none", 
                               eval.grid = FALSE, panel = FALSE)$estimate
test <- cor.test(NonTD_MvsE$MARBLES, NonTD_MvsE$EARLI)
r <- round(test$estimate, 3)
p <- round(test$p.value, 3)
gg <- ggplot()
gg + 
        geom_point(data = NonTD_MvsE, aes(MARBLES, EARLI, color = Density), size = 1.4) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "red") +
        annotate("text", x = -0.6, y = 1.23, label = paste("r =", r, sep = " "), size = 7, hjust = 0) +
        annotate("text", x = -0.6, y = 1.08, label = paste("p =", p, sep = " "), size = 7, hjust = 0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.1, color = "black"), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text = element_text(color = "black"), axis.title = element_text(size = 22),
              plot.margin = unit(c(1.5, 2, 1, 1), "lines")) +
        xlab(expr("MARBLES NonTD log"[2]*"(Fold Change)")) +
        ylab(expr("EARLI NonTD log"[2]*"(Fold Change)")) +
        coord_cartesian(xlim = c(-0.6, 1.2), ylim = c(-0.6, 1.2)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5))
ggsave("Figures/Revisions NonTD MARBLES vs EARLI log Fold Change Scatterplot.png", dpi = 600, height = 6, width = 7)

gg <- ggplot()
gg + 
        geom_point(data = NonTD_MvsE, aes(MARBLES, EARLI, color = Density), size = 1.4) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "red") +
        annotate("text", x = -0.6, y = 1.23, label = paste("r =", r, sep = " "), size = 7, hjust = 0) +
        annotate("text", x = -0.6, y = 1.08, label = paste("p =", p, sep = " "), size = 7, hjust = 0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.1, color = "black"), panel.grid.minor = element_blank(),
              legend.position = c(1.2, 0.8), axis.text = element_text(color = "black"), axis.title = element_text(size = 22),
              plot.margin = unit(c(1.5, 8, 1, 1), "lines"), legend.background = element_blank()) +
        xlab(expr("MARBLES NonTD log"[2]*"(Fold Change)")) +
        ylab(expr("EARLI NonTD log"[2]*"(Fold Change)")) +
        coord_cartesian(xlim = c(-0.6, 1.2), ylim = c(-0.6, 1.2)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5))
ggsave("Figures/Revisions NonTD MARBLES vs EARLI log Fold Change Scatterplot with legend.png", dpi = 600, height = 6, width = 7)
