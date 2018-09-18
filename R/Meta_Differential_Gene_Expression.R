# Meta-Analysis Differential Gene Expression ####
# MARBLES and EARLI
# Charles Mordaunt
# 6/15/18

setwd()

# Packages ####
library(ggplot2)
library(scales)
library(proto)
library(VennDiagram)
library(circlize)
library(biomaRt)
library(RColorBrewer)
library(WebGestaltR)

# Functions ####
# Detect and prevent collisions.
# Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
        # Determine height
        if (!is.null(height)) {
                # height set manually
                if (!(all(c("ymin", "ymax") %in% names(data)))) {
                        data$ymin <- data$y - height / 2
                        data$ymax <- data$y + height / 2
                }
        } else {
                if (!(all(c("ymin", "ymax") %in% names(data)))) {
                        data$ymin <- data$y
                        data$ymax <- data$y
                }
                
                # height determined from data, must be floating point constant
                heights <- unique(data$ymax - data$ymin)
                heights <- heights[!is.na(heights)]
                
                #   # Suppress warning message since it's not reliable
                #     if (!zero_range(range(heights))) {
                #       warning(name, " requires constant height: output may be incorrect",
                #         call. = FALSE)
                #     }
                height <- heights[1]
        }
        
        # Reorder by x position, relying on stable sort to preserve existing
        # ordering, which may be by group or order.
        data <- data[order(data$ymin), ]
        
        # Check for overlap
        intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
        intervals <- intervals[!is.na(intervals)]
        
        if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
                warning(name, " requires non-overlapping y intervals", call. = FALSE)
                # This is where the algorithm from [L. Wilkinson. Dot plots.
                # The American Statistician, 1999.] should be used
        }
        
        if (!is.null(data$xmax)) {
                plyr::ddply(data, "ymin", strategy, height = height)
        } else if (!is.null(data$x)) {
                data$xmax <- data$x
                data <- plyr::ddply(data, "ymin", strategy, height = height)
                data$x <- data$xmax
                data
        } else {
                stop("Neither x nor xmax defined")
        }
}

# Stack overlapping intervals.
# Assumes that each set has the same horizontal position
pos_stackv <- function(df, height) {
        if (nrow(df) == 1) return(df)
        
        n <- nrow(df) + 1
        x <- ifelse(is.na(df$x), 0, df$x)
        if (all(is.na(df$y))) {
                heights <- rep(NA, n)
        } else {
                heights <- c(0, cumsum(x))
        }
        
        df$xmin <- heights[-n]
        df$xmax <- heights[-1]
        df$x <- df$xmax
        df
}

# Stack overlapping intervals and set height to 1.
# Assumes that each set has the same horizontal position.
pos_fillv <- function(df, height) {
        stacked <- pos_stackv(df, height)
        stacked$xmin <- stacked$xmin / max(stacked$xmax)
        stacked$xmax <- stacked$xmax / max(stacked$xmax)
        stacked$x <- stacked$xmax
        stacked
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
        n <- length(unique(df$group))
        if (n == 1) return(df)
        
        if (!all(c("ymin", "ymax") %in% names(df))) {
                df$ymin <- df$y
                df$ymax <- df$y
        }
        
        d_height <- max(df$ymax - df$ymin)
        
        # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
        # ggplot(df, aes(n, div)) + geom_point()
        
        # Have a new group index from 1 to number of groups.
        # This might be needed if the group numbers in this set don't include all of 1:n
        groupidy <- match(df$group, sort(unique(df$group)))
        
        # Find the center for each group, then use that to calculate xmin and xmax
        df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
        df$ymin <- df$y - d_height / n / 2
        df$ymax <- df$y + d_height / n / 2
        
        df
}

position_dodgev <- function(height = NULL) {
        ggproto(NULL, PositionDodgeV, height = height)
}

PositionDodgeV <- ggproto(`_class` = "PositionDodgeV", `_inherit` = Position,
                          required_aes = "y",
                          height = NULL,
                          setup_params = function(self, data) {
                                  if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
                                          warning("height not defined. Set with `position_dodgev(height = ?)`",
                                                  call. = FALSE)
                                  }
                                  list(height = self$height)
                          },
                          
                          compute_panel = function(data, params, scales) {
                                  collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
                          }
)

# Data ####
marbles <- read.delim("Tables/MARBLES All Genes Stats Dx_alg SVAonly B123.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
earliASD <- read.csv("EARLI Data/ASD_TD limma output table 092217.csv")
earliNonTD <- read.csv("EARLI Data/Non_TD limma output table 092217.csv")

# Prep Tables for METAL Effect Size Analysis #### 
# Use this
# 95% CI.R = mean + 1.96*se
# se = (95% CI.R - mean)/1.96

# MARBLES ASD
METAL_marblesASD_eff <- marbles[,c("ASD_P.Value", "Probe")]
METAL_marblesASD_eff$Effect <- marbles$Dx_alg_ASD
METAL_marblesASD_eff$RefAllele <- rep("T", nrow(METAL_marblesASD_eff))
METAL_marblesASD_eff$OtherAllele <- rep("A", nrow(METAL_marblesASD_eff))
METAL_marblesASD_eff <- METAL_marblesASD_eff[order(METAL_marblesASD_eff$Probe),]
table(METAL_marblesASD_eff$Probe == marbles$Probe) #All TRUE
METAL_marblesASD_eff$CI.R <- marbles$ASD_CI.R
METAL_marblesASD_eff$se <- (METAL_marblesASD_eff$CI.R - METAL_marblesASD_eff$Effect)/1.96
METAL_marblesASD_eff <- METAL_marblesASD_eff[,c("Probe", "RefAllele", "OtherAllele", "Effect", "se")]
write.table(METAL_marblesASD_eff, "Tables/METAL_MARBLES_ASD_eff_B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# MARBLES NonTD
METAL_marblesNonTD_eff <- marbles[,c("NonTD_P.Value", "Probe")]
METAL_marblesNonTD_eff$Effect <- marbles$Dx_alg_NonTD
METAL_marblesNonTD_eff$RefAllele <- rep("T", nrow(METAL_marblesNonTD_eff))
METAL_marblesNonTD_eff$OtherAllele <- rep("A", nrow(METAL_marblesNonTD_eff))
METAL_marblesNonTD_eff <- METAL_marblesNonTD_eff[order(METAL_marblesNonTD_eff$Probe),]
table(METAL_marblesNonTD_eff$Probe == marbles$Probe) #All TRUE
METAL_marblesNonTD_eff$CI.R <- marbles$NonTD_CI.R
METAL_marblesNonTD_eff$se <- (METAL_marblesNonTD_eff$CI.R - METAL_marblesNonTD_eff$Effect)/1.96
METAL_marblesNonTD_eff <- METAL_marblesNonTD_eff[,c("Probe", "RefAllele", "OtherAllele", "Effect", "se")]
write.table(METAL_marblesNonTD_eff, "Tables/METAL_MARBLES_NonTD_eff_B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# earli ASD
earliASD <- earliASD[order(earliASD$X),]
METAL_earliASD_eff <- earliASD[,c("P.Value", "X")]
colnames(METAL_earliASD_eff) <- c("P.Value", "Probe")
METAL_earliASD_eff$Effect <- earliASD$logFC
METAL_earliASD_eff$RefAllele <- rep("T", nrow(METAL_earliASD_eff))
METAL_earliASD_eff$OtherAllele <- rep("A", nrow(METAL_earliASD_eff))
METAL_earliASD_eff <- METAL_earliASD_eff[order(METAL_earliASD_eff$Probe),]
table(METAL_earliASD_eff$Probe == earliASD$X) #All TRUE
METAL_earliASD_eff$CI.R <- earliASD$CI.R
METAL_earliASD_eff$se <- (METAL_earliASD_eff$CI.R - METAL_earliASD_eff$Effect)/1.96
METAL_earliASD_eff <- METAL_earliASD_eff[,c("Probe", "RefAllele", "OtherAllele", "Effect", "se")]
write.table(METAL_earliASD_eff, "Tables/METAL_earli_ASD_eff_B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# earli NonTD
earliNonTD <- earliNonTD[order(earliNonTD$X),]
METAL_earliNonTD_eff <- earliNonTD[,c("P.Value", "X")]
colnames(METAL_earliNonTD_eff) <- c("P.Value", "Probe")
METAL_earliNonTD_eff$Effect <- earliNonTD$logFC
METAL_earliNonTD_eff$RefAllele <- rep("T", nrow(METAL_earliNonTD_eff))
METAL_earliNonTD_eff$OtherAllele <- rep("A", nrow(METAL_earliNonTD_eff))
METAL_earliNonTD_eff <- METAL_earliNonTD_eff[order(METAL_earliNonTD_eff$Probe),]
table(METAL_earliNonTD_eff$Probe == earliNonTD$X) #All TRUE
METAL_earliNonTD_eff$CI.R <- earliNonTD$CI.R
METAL_earliNonTD_eff$se <- (METAL_earliNonTD_eff$CI.R - METAL_earliNonTD_eff$Effect)/1.96
METAL_earliNonTD_eff <- METAL_earliNonTD_eff[,c("Probe", "RefAllele", "OtherAllele", "Effect", "se")]
write.table(METAL_earliNonTD_eff, "Tables/METAL_earli_NonTD_eff_B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# Process with METAL

# Results from METAL Effect Size Analysis ####
# Effect size estimates weighted by inverse of standard error, so effect sizes with smaller standard errors are weighted more heavily

# ASD
marbles <- marbles[order(marbles$Probe),]
earliASD <- earliASD[order(earliASD$X),]
METAL_ASD <- read.delim("Tables/MetaAnalysis_ASD_eff_B123.tbl")
METAL_ASD <- METAL_ASD[order(METAL_ASD$MarkerName),]
METAL_ASD$Effect <- -METAL_ASD$Effect #METAL switched effect direction

ASDstats <- cbind(marbles[,c("Probe", "GeneID", "GeneName")], METAL_marblesASD_eff[,c("Effect", "se")], marbles[,c("ASD_P.Value", "ASD_adj.P.Val")], 
                  METAL_earliASD_eff[,c("Effect", "se")], earliASD[,c("P.Value", "adj.P.Val")], METAL_ASD[,c("Effect", "StdErr", "P.value", "Direction")])
colnames(ASDstats) <- c("Probe", "GeneID", "GeneName", "marbles_logFC", "marbles_StdErr", "marbles_pValue", "marbles_adj_pValue", 
                        "earli_logFC", "earli_StdErr", "earli_pValue", "earli_adj_pValue", "Meta_logFC", "Meta_StdErr", "Meta_pValue", "Meta_Direction")
ASDstats$Meta_adj_pValue <- p.adjust(ASDstats$Meta_pValue, "fdr")
ASDdiff <- subset(ASDstats, abs(Meta_logFC) > 0.1 & Meta_pValue < 0.01) #185 probes 
ASDdiff <- subset(ASDdiff, Meta_Direction == "++" | Meta_Direction == "--") #178 probes

ASDstats <- ASDstats[,c("Probe", "GeneID", "GeneName", "marbles_logFC", "marbles_StdErr", "marbles_pValue", "marbles_adj_pValue", 
                        "earli_logFC", "earli_StdErr", "earli_pValue", "earli_adj_pValue", "Meta_logFC", "Meta_StdErr", "Meta_pValue", 
                        "Meta_adj_pValue")]
write.table(ASDstats, "Tables/ASD Meta-Analysis with METAL All Genes B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

ASDdiff <- ASDdiff[,c("Probe", "GeneID", "GeneName", "marbles_logFC", "marbles_StdErr", "marbles_pValue", "marbles_adj_pValue", 
                      "earli_logFC", "earli_StdErr", "earli_pValue", "earli_adj_pValue", "Meta_logFC", "Meta_StdErr", "Meta_pValue", 
                      "Meta_adj_pValue")]
write.table(ASDdiff, "Tables/ASD Meta-Analysis with METAL Differential Genes B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# NonTD
marbles <- marbles[order(marbles$Probe),]
earliNonTD <- earliNonTD[order(earliNonTD$X),]
METAL_NonTD <- read.delim("Tables/MetaAnalysis_NonTD_eff_B123.tbl")
METAL_NonTD <- METAL_NonTD[order(METAL_NonTD$MarkerName),]
METAL_NonTD$Effect <- -METAL_NonTD$Effect #METAL switched effect direction

NonTDstats <- cbind(marbles[,c("Probe", "GeneID", "GeneName")], METAL_marblesNonTD_eff[,c("Effect", "se")], marbles[,c("NonTD_P.Value", "NonTD_adj.P.Val")], 
                  METAL_earliNonTD_eff[,c("Effect", "se")], earliNonTD[,c("P.Value", "adj.P.Val")], METAL_NonTD[,c("Effect", "StdErr", "P.value", "Direction")])
colnames(NonTDstats) <- c("Probe", "GeneID", "GeneName", "marbles_logFC", "marbles_StdErr", "marbles_pValue", "marbles_adj_pValue", 
                        "earli_logFC", "earli_StdErr", "earli_pValue", "earli_adj_pValue", "Meta_logFC", "Meta_StdErr", "Meta_pValue", "Meta_Direction")
NonTDstats$Meta_adj_pValue <- p.adjust(NonTDstats$Meta_pValue, "fdr")
NonTDdiff <- subset(NonTDstats, abs(Meta_logFC) > 0.1 & Meta_pValue < 0.01) #69 probes 
NonTDdiff <- subset(NonTDdiff, Meta_Direction == "++" | Meta_Direction == "--") #66 probes

NonTDstats <- NonTDstats[,c("Probe", "GeneID", "GeneName", "marbles_logFC", "marbles_StdErr", "marbles_pValue", "marbles_adj_pValue", 
                        "earli_logFC", "earli_StdErr", "earli_pValue", "earli_adj_pValue", "Meta_logFC", "Meta_StdErr", "Meta_pValue", 
                        "Meta_adj_pValue")]
write.table(NonTDstats, "Tables/NonTD Meta-Analysis with METAL All Genes B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

NonTDdiff <- NonTDdiff[,c("Probe", "GeneID", "GeneName", "marbles_logFC", "marbles_StdErr", "marbles_pValue", "marbles_adj_pValue", 
                      "earli_logFC", "earli_StdErr", "earli_pValue", "earli_adj_pValue", "Meta_logFC", "Meta_StdErr", "Meta_pValue", 
                      "Meta_adj_pValue")]
write.table(NonTDdiff, "Tables/NonTD Meta-Analysis with METAL Differential Genes B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# WebGestalt GSEA ####
GSEA_probe_ASD <- data.frame(probe = ASDstats$Probe, log2FC = ASDstats$Meta_logFC)
write.table(GSEA_probe_ASD, "Tables/GSEA Meta ASD Dx_alg log2FC probe SVAonly B123.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

GSEA_probe_NonTD <- data.frame(probe = NonTDstats$Probe, log2FC = NonTDstats$Meta_logFC)
write.table(GSEA_probe_NonTD, "Tables/GSEA Meta NonTD Dx_alg log2FC probe SVAonly B123.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Differential Genes ####
ASDdiffGenes <- sort(unique(as.character(ASDdiff$GeneName))) #172 genes
NonTDdiffGenes <- sort(unique(as.character(NonTDdiff$GeneName))) #66 genes
ASDandNonTDdiffGenes <- intersect(ASDdiffGenes, NonTDdiffGenes) #8 genes
ASDupGenes <- sort(unique(as.character(ASDdiff$GeneName[ASDdiff$Meta_logFC > 0]))) #87 genes
ASDdownGenes <- sort(unique(as.character(ASDdiff$GeneName[ASDdiff$Meta_logFC < 0]))) #85
NonTDupGenes <- sort(unique(as.character(NonTDdiff$GeneName[NonTDdiff$Meta_logFC > 0]))) #38
NonTDdownGenes <- sort(unique(as.character(NonTDdiff$GeneName[NonTDdiff$Meta_logFC < 0]))) #28

# ASD-NonTD Overlap 
intersect(ASDupGenes, NonTDupGenes) #[1] " IGLV1-40 "  " LRRC37A4P " " MIR1299 "   " PMCHL2 "    " TRBV11-2 " 
intersect(ASDdownGenes, NonTDdownGenes) #[1] " RNU4ATAC11P " " TPTE2P5 "     " TRIM74 "     

venn.diagram(list("ASD"=ASDdiffGenes, "Non-TD"=NonTDdiffGenes), 
             rotation.degree = 0, margin = 0.02, filename = "Figures/Meta-Analysis with METAL ASD vs Non-TD Venn B123.png", 
             cat.cex = 1, cex = 1.4, fill = c("lightblue", "lightpink"), ext.text=FALSE,
             resolution = 600, imagetype = "png", cat.pos = c(0,0), fontfamily = "sans", cat.fontfamily = "sans", 
             height = 2100, width = 2200, cat.dist = c(0.05, 0.05))
length(unique(ASDstats$GeneName)) #33897 unique genes
ASD_nonTD_diff <- matrix(c(33897-164-58-8, 58, 164, 8), nrow=2, byrow=TRUE, dimnames=list(ASDdiff=c("FALSE", "TRUE"), nonTDdiff=c("FALSE", "TRUE")))
fisher.test(ASD_nonTD_diff)
# p-value = 1.668e-09
# 95% CI = 11.48585 60.77793
# odds ratio = 28.30499 

# TD vs ASD NonTD Gene Forest Plot (Meta Analysis values)
ASDdiff$GeneName <- gsub(" ", "", as.character(ASDdiff$GeneName), fixed=TRUE)
NonTDdiff$GeneName <- gsub(" ", "", as.character(NonTDdiff$GeneName), fixed=TRUE)
TDvsASD_diff_NonTD <- merge(ASDdiff, NonTDdiff, by="Probe")
TDvsASD_diff_NonTD <- TDvsASD_diff_NonTD[,c("Probe", "GeneID.x", "GeneName.x", "Meta_logFC.x", "Meta_StdErr.x", "Meta_pValue.x", "Meta_adj_pValue.x",
                                            "GeneID.y", "GeneName.y", "Meta_logFC.y", "Meta_StdErr.y", "Meta_pValue.y", "Meta_adj_pValue.y")]
colnames(TDvsASD_diff_NonTD) <- c("Probe", "GeneID_ASD", "GeneName_ASD", "Meta_logFC_ASD", "Meta_StdErr_ASD", "Meta_pValue_ASD", "Meta_adj_pValue_ASD",
                                  "GeneID_NonTD", "GeneName_NonTD", "Meta_logFC_NonTD", "Meta_StdErr_NonTD", "Meta_pValue_NonTD", "Meta_adj_pValue_NonTD")
TDvsASD_diff_NonTD$Meta_ErrBar.L_ASD <- TDvsASD_diff_NonTD$Meta_logFC_ASD - TDvsASD_diff_NonTD$Meta_StdErr_ASD
TDvsASD_diff_NonTD$Meta_ErrBar.R_ASD <- TDvsASD_diff_NonTD$Meta_logFC_ASD + TDvsASD_diff_NonTD$Meta_StdErr_ASD
TDvsASD_diff_NonTD$Meta_ErrBar.L_NonTD <- TDvsASD_diff_NonTD$Meta_logFC_NonTD - TDvsASD_diff_NonTD$Meta_StdErr_NonTD
TDvsASD_diff_NonTD$Meta_ErrBar.R_NonTD <- TDvsASD_diff_NonTD$Meta_logFC_NonTD + TDvsASD_diff_NonTD$Meta_StdErr_NonTD
TDvsASD_forest_NonTD <- as.data.frame(rbind(as.matrix(TDvsASD_diff_NonTD[,c("GeneName_ASD", "Meta_logFC_ASD", "Meta_ErrBar.L_ASD", "Meta_ErrBar.R_ASD", "Meta_pValue_ASD")]),
                                            as.matrix(TDvsASD_diff_NonTD[,c("GeneName_NonTD", "Meta_logFC_NonTD", "Meta_ErrBar.L_NonTD", "Meta_ErrBar.R_NonTD", "Meta_pValue_NonTD")])), 
                                      stringsAsFactors=FALSE)
TDvsASD_forest_NonTD$Diagnosis <- factor(rep(c("ASD", "NonTD"), each=dim(TDvsASD_diff_NonTD)[1]), levels=c("ASD", "NonTD"), ordered=TRUE)
TDvsASD_forest_NonTD$GeneName_ASD <- factor(TDvsASD_forest_NonTD$GeneName_ASD, levels=TDvsASD_diff_NonTD$GeneName_ASD[order(TDvsASD_diff_NonTD$Meta_logFC_NonTD)], ordered=TRUE)
TDvsASD_forest_NonTD$Meta_logFC_ASD <- as.numeric(TDvsASD_forest_NonTD$Meta_logFC_ASD)
TDvsASD_forest_NonTD$Meta_ErrBar.L_ASD <- as.numeric(TDvsASD_forest_NonTD$Meta_ErrBar.L_ASD)
TDvsASD_forest_NonTD$Meta_ErrBar.R_ASD <- as.numeric(TDvsASD_forest_NonTD$Meta_ErrBar.R_ASD)
TDvsASD_forest_NonTD$Meta_pValue_ASD <- as.numeric(TDvsASD_forest_NonTD$Meta_pValue_ASD)
colnames(TDvsASD_forest_NonTD) <- c("Gene", "log2FC", "ErrBar.L", "ErrBar.R", "p-value", "Diagnosis")
TDvsASD_forest_NonTD$logp <- -log10(TDvsASD_forest_NonTD$`p-value`)

gg <- ggplot(data = TDvsASD_forest_NonTD)
gg +
        geom_vline(xintercept=0, color="black", lty=2, size=1.05) +
        geom_errorbarh(aes(x=log2FC, xmin=ErrBar.L, xmax=ErrBar.R, y=Gene, color=Diagnosis), position=position_dodgev(height=0.75), size = 1, height=0) +
        geom_point(aes(x = log2FC, y = Gene, color = Diagnosis), size=3, position=position_dodgev(height=0.75)) +
        theme_bw(base_size = 26) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.15, 1.05), 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "black"),
              legend.key.height = unit(2, "lines"), legend.title=element_blank(), legend.direction="horizontal",
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22, face = "italic"),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title.x = element_text(color = "black", size = 24), 
              axis.title.y = element_blank()) +
        scale_color_manual(breaks=c("NonTD", "ASD"), values = c("NonTD" = "#FF3366", "ASD" = "#3366CC")) +
        xlab(expression(log[2]*"(Fold Change)")) +
        coord_cartesian(xlim=c(-0.25, 0.3)) +
        scale_x_continuous(breaks=pretty_breaks(n=6))
ggsave("Figures/EARLI MARBLES Meta-analysis with METAL TD vs ASD NonTD B123.png", dpi = 600, width = 10, height = 7, units = "in")

# Volcano Plots ####
# Meta ASD
volcanoData <- ASDstats
volcanoData$Meta_log10p <- -log10(volcanoData$Meta_pValue)
volcanoData$Diff <- abs(volcanoData$Meta_logFC) > 0.1 & volcanoData$Meta_pValue < 0.01
volcanoData$Diff <- factor(volcanoData$Diff, levels = c("TRUE", "FALSE"), ordered = TRUE)

gg <- ggplot(data = volcanoData)
gg +
        geom_point(aes(x = Meta_logFC, y = Meta_log10p, color = Diff), size = 2) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = "none", 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black"),
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 22),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 24), 
              legend.title = element_text(color = "black", size = 24)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        xlab(expression(log[2]*"(Fold Change)")) +
        ylab(expression(-log[10]*"(p-value)")) +
        coord_cartesian(xlim=c(-0.4,0.8), ylim=c(0,5)) +
        scale_x_continuous(breaks=pretty_breaks(n=5)) +
        scale_y_continuous(expand=c(0.01,0))
ggsave("Figures/ASD Meta-Analysis with METAL Dx_alg Diff Volcano Plot B123.png", dpi = 600, width = 10, height = 10, units = "in")

# Meta NonTD
volcanoData <- NonTDstats
volcanoData$Meta_log10p <- -log10(volcanoData$Meta_pValue)
volcanoData$Diff <- abs(volcanoData$Meta_logFC) > 0.1 & volcanoData$Meta_pValue < 0.01
volcanoData$Diff <- factor(volcanoData$Diff, levels = c("TRUE", "FALSE"), ordered = TRUE)

gg <- ggplot(data = volcanoData)
gg +
        geom_point(aes(x = Meta_logFC, y = Meta_log10p, color = Diff), size = 2) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = "none", 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black"),
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 22),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 24), 
              legend.title = element_text(color = "black", size = 24)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        xlab(expression(log[2]*"(Fold Change)")) +
        ylab(expression(-log[10]*"(p-value)")) +
        coord_cartesian(xlim=c(-0.3,0.5), ylim=c(0,5)) +
        scale_y_continuous(expand=c(0.01,0))
ggsave("Figures/NonTD Meta-Analysis with METAL Dx_alg Diff Volcano Plot B123.png", dpi = 600, width = 10, height = 10, units = "in")

# marbles ASD
volcanoData <- ASDstats
volcanoData$marbles_log10p <- -log10(volcanoData$marbles_pValue)
volcanoData$Diff <- abs(volcanoData$marbles_logFC) > 0.1 & volcanoData$marbles_pValue < 0.01
volcanoData$Diff <- factor(volcanoData$Diff, levels = c("TRUE", "FALSE"), ordered = TRUE)
table(volcanoData$Diff)
# TRUE FALSE 
# 295 36164 
length(unique(volcanoData$GeneName[volcanoData$Diff == "TRUE"])) # 291 genes

summary(volcanoData$marbles_log10p) # Min 0, Max 4.45
summary(volcanoData$marbles_logFC) # Min -0.48, Max 0.78

gg <- ggplot(data = volcanoData)
gg +
        geom_point(aes(x = marbles_logFC, y = marbles_log10p, color = Diff), size = 2) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = "none", 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black"),
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 22),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 24), 
              legend.title = element_text(color = "black", size = 24)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        xlab(expression(log[2]*"(Fold Change)")) +
        ylab(expression(-log[10]*"(p-value)")) +
        coord_cartesian(xlim=c(-0.6,1), ylim=c(0,4.8)) +
        scale_y_continuous(breaks=pretty_breaks(n=5), expand=c(0.01,0))
ggsave("Figures/ASD MARBLES Dx_alg Diff Volcano Plot B123.png", dpi = 600, width = 10, height = 10, units = "in")

# marbles NonTD
volcanoData <- NonTDstats
volcanoData$marbles_log10p <- -log10(volcanoData$marbles_pValue)
volcanoData$Diff <- abs(volcanoData$marbles_logFC) > 0.1 & volcanoData$marbles_pValue < 0.01
volcanoData$Diff <- factor(volcanoData$Diff, levels = c("TRUE", "FALSE"), ordered = TRUE)
table(volcanoData$Diff)
# TRUE FALSE 
# 208 36251 

length(unique(volcanoData$GeneName[volcanoData$Diff == "TRUE"])) # 201 genes

summary(volcanoData$marbles_log10p) # Min 0, Max 4.67 
summary(volcanoData$marbles_logFC) # Min -0.49, Max 0.55

gg <- ggplot(data = volcanoData)
gg +
        geom_point(aes(x = marbles_logFC, y = marbles_log10p, color = Diff), size = 2) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = "none", 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black"),
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 22),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 24), 
              legend.title = element_text(color = "black", size = 24)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        xlab(expression(log[2]*"(Fold Change)")) +
        ylab(expression(-log[10]*"(p-value)")) +
        coord_cartesian(xlim=c(-0.6,1.2), ylim=c(0,4.7)) +
        scale_y_continuous(expand=c(0.01,0))
ggsave("Figures/NonTD MARBLES Dx_alg Diff Volcano Plot B123.png", dpi = 600, width = 10, height = 10, units = "in")

# earli ASD
volcanoData <- ASDstats
volcanoData$earli_log10p <- -log10(volcanoData$earli_pValue)
volcanoData$Diff <- abs(volcanoData$earli_logFC) > 0.1 & volcanoData$earli_pValue < 0.01
volcanoData$Diff <- factor(volcanoData$Diff, levels = c("TRUE", "FALSE"), ordered = TRUE)
table(volcanoData$Diff)
# TRUE FALSE 
#  392 36067 

length(unique(volcanoData$GeneName[volcanoData$Diff == "TRUE"])) # 386 genes

summary(volcanoData$earli_log10p) # Min 0, Max 4.75
summary(volcanoData$earli_logFC) # Min -0.57, Max 0.94

gg <- ggplot(data = volcanoData)
gg +
        geom_point(aes(x = earli_logFC, y = earli_log10p, color = Diff), size = 2) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = "none", 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black"),
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 22),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 24), 
              legend.title = element_text(color = "black", size = 24)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        xlab(expression(log[2]*"(Fold Change)")) +
        ylab(expression(-log[10]*"(p-value)")) +
        coord_cartesian(xlim=c(-0.6,1), ylim=c(0,4.8)) +
        scale_y_continuous(breaks=pretty_breaks(n=5),expand=c(0.01,0))
ggsave("Figures/ASD earli Dx_alg Diff Volcano Plot B123.png", dpi = 600, width = 10, height = 10, units = "in")

# earli NonTD
volcanoData <- NonTDstats
volcanoData$earli_log10p <- -log10(volcanoData$earli_pValue)
volcanoData$Diff <- abs(volcanoData$earli_logFC) > 0.1 & volcanoData$earli_pValue < 0.01
volcanoData$Diff <- factor(volcanoData$Diff, levels = c("TRUE", "FALSE"), ordered = TRUE)
table(volcanoData$Diff)
# TRUE FALSE 
#  271 36188 

length(unique(volcanoData$GeneName[volcanoData$Diff == "TRUE"])) # 260 genes

summary(volcanoData$earli_log10p) # Min 0, Max 4.44
summary(volcanoData$earli_logFC) # Min -0.53, Max 1.17

gg <- ggplot(data = volcanoData)
gg +
        geom_point(aes(x = earli_logFC, y = earli_log10p, color = Diff), size = 2) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = "none", 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black"),
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 22),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 24), 
              legend.title = element_text(color = "black", size = 24)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        xlab(expression(log[2]*"(Fold Change)")) +
        ylab(expression(-log[10]*"(p-value)")) +
        coord_cartesian(xlim=c(-0.6,1.2), ylim=c(0,4.7)) +
        scale_y_continuous(expand=c(0.01,0))
ggsave("Figures/NonTD earli Dx_alg Diff Volcano Plot B123.png", dpi = 600, width = 10, height = 10, units = "in")

rm(volcanoData, gg)

# Forest Plots Top 20 Log2FC ####
# TD vs ASD
TDvsASD_diffTopFC <- subset(ASDdiff, rank(-abs(ASDdiff$Meta_logFC)) <= 20)
TDvsASD_diffTopFC$marbles_ErrBar.L <- TDvsASD_diffTopFC$marbles_logFC - TDvsASD_diffTopFC$marbles_StdErr
TDvsASD_diffTopFC$marbles_ErrBar.R <- TDvsASD_diffTopFC$marbles_logFC + TDvsASD_diffTopFC$marbles_StdErr
TDvsASD_diffTopFC$earli_ErrBar.L <- TDvsASD_diffTopFC$earli_logFC - TDvsASD_diffTopFC$earli_StdErr
TDvsASD_diffTopFC$earli_ErrBar.R <- TDvsASD_diffTopFC$earli_logFC + TDvsASD_diffTopFC$earli_StdErr
TDvsASD_diffTopFC$Meta_ErrBar.L <- TDvsASD_diffTopFC$Meta_logFC - TDvsASD_diffTopFC$Meta_StdErr
TDvsASD_diffTopFC$Meta_ErrBar.R <- TDvsASD_diffTopFC$Meta_logFC + TDvsASD_diffTopFC$Meta_StdErr
TDvsASD_diffTopFC <- subset(TDvsASD_diffTopFC, !duplicated(TDvsASD_diffTopFC$GeneName))

TDvsASD_forestFC <- as.data.frame(rbind(as.matrix(TDvsASD_diffTopFC[,c("GeneName", "marbles_logFC", "marbles_ErrBar.L", "marbles_ErrBar.R", "marbles_pValue")]),
                                        as.matrix(TDvsASD_diffTopFC[,c("GeneName", "earli_logFC", "earli_ErrBar.L", "earli_ErrBar.R", "earli_pValue")]), 
                                        as.matrix(TDvsASD_diffTopFC[,c("GeneName", "Meta_logFC", "Meta_ErrBar.L", "Meta_ErrBar.R", "Meta_pValue")])),
                                  stringsAsFactors=FALSE)
TDvsASD_forestFC$Study <- factor(rep(c("MARBLES", "EARLI", "Meta"), each=dim(TDvsASD_diffTopFC)[1]), levels=c("Meta", "EARLI", "MARBLES"), ordered=TRUE)
TDvsASD_forestFC$GeneName <- factor(TDvsASD_forestFC$GeneName, levels=TDvsASD_diffTopFC$GeneName[order(TDvsASD_diffTopFC$Meta_logFC)], ordered=TRUE)
TDvsASD_forestFC$marbles_logFC <- as.numeric(TDvsASD_forestFC$marbles_logFC)
TDvsASD_forestFC$marbles_ErrBar.L <- as.numeric(TDvsASD_forestFC$marbles_ErrBar.L)
TDvsASD_forestFC$marbles_ErrBar.R <- as.numeric(TDvsASD_forestFC$marbles_ErrBar.R)
TDvsASD_forestFC$marbles_pValue <- as.numeric(TDvsASD_forestFC$marbles_pValue)
colnames(TDvsASD_forestFC) <- c("Gene", "log2FC", "ErrBar.L", "ErrBar.R", "p-value", "Study")
TDvsASD_forestFC$logp <- -log10(TDvsASD_forestFC$`p-value`)

gg <- ggplot(data = TDvsASD_forestFC)
gg +
        geom_vline(xintercept=0, color="black", lty=2, size=1.05) +
        geom_errorbarh(aes(x=log2FC, xmin=ErrBar.L, xmax=ErrBar.R, y=Gene, color=Study), position=position_dodgev(height=0.85), size = 1, height=0) +
        geom_point(aes(x = log2FC, y = Gene, color = Study, size=Study), position=position_dodgev(height=0.85)) +
        theme_bw(base_size = 26) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.28, 1.03), 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "black"),
              legend.key.height = unit(2, "lines"), legend.title=element_blank(), legend.direction="horizontal",
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22, face = "italic"),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title.x = element_text(color = "black", size = 24), 
              axis.title.y = element_blank()) +
        scale_color_manual(breaks=c("MARBLES", "EARLI", "Meta"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366", "Meta" = "#009966")) +
        scale_size_manual(breaks=c("MARBLES", "EARLI", "Meta"), values = c("MARBLES" = 2.5, "EARLI" = 2.5, "Meta" = 3.25)) +
        xlab(expression(log[2]*"(Fold Change)")) +
        coord_cartesian(xlim=c(-0.7, 0.6)) +
        scale_x_continuous(breaks=pretty_breaks(n=6))
ggsave("Figures/EARLI MARBLES Meta-analysis with METAL TD vs ASD B123 by FC.png", dpi = 600, width = 10, height = 10, units = "in")

# TD vs NonTD
TDvsNonTD_diffTopFC <- subset(NonTDdiff, rank(-abs(NonTDdiff$Meta_logFC)) <= 20)
TDvsNonTD_diffTopFC$marbles_ErrBar.L <- TDvsNonTD_diffTopFC$marbles_logFC - TDvsNonTD_diffTopFC$marbles_StdErr
TDvsNonTD_diffTopFC$marbles_ErrBar.R <- TDvsNonTD_diffTopFC$marbles_logFC + TDvsNonTD_diffTopFC$marbles_StdErr
TDvsNonTD_diffTopFC$earli_ErrBar.L <- TDvsNonTD_diffTopFC$earli_logFC - TDvsNonTD_diffTopFC$earli_StdErr
TDvsNonTD_diffTopFC$earli_ErrBar.R <- TDvsNonTD_diffTopFC$earli_logFC + TDvsNonTD_diffTopFC$earli_StdErr
TDvsNonTD_diffTopFC$Meta_ErrBar.L <- TDvsNonTD_diffTopFC$Meta_logFC - TDvsNonTD_diffTopFC$Meta_StdErr
TDvsNonTD_diffTopFC$Meta_ErrBar.R <- TDvsNonTD_diffTopFC$Meta_logFC + TDvsNonTD_diffTopFC$Meta_StdErr
TDvsNonTD_diffTopFC <- subset(TDvsNonTD_diffTopFC, !duplicated(TDvsNonTD_diffTopFC$GeneName))

TDvsNonTD_forestFC <- as.data.frame(rbind(as.matrix(TDvsNonTD_diffTopFC[,c("GeneName", "marbles_logFC", "marbles_ErrBar.L", "marbles_ErrBar.R", "marbles_pValue")]),
                                          as.matrix(TDvsNonTD_diffTopFC[,c("GeneName", "earli_logFC", "earli_ErrBar.L", "earli_ErrBar.R", "earli_pValue")]), 
                                          as.matrix(TDvsNonTD_diffTopFC[,c("GeneName", "Meta_logFC", "Meta_ErrBar.L", "Meta_ErrBar.R", "Meta_pValue")])),
                                    stringsAsFactors=FALSE)
TDvsNonTD_forestFC$Study <- factor(rep(c("MARBLES", "EARLI", "Meta"), each=dim(TDvsNonTD_diffTopFC)[1]), levels=c("Meta", "EARLI", "MARBLES"), ordered=TRUE)
TDvsNonTD_forestFC$GeneName <- factor(TDvsNonTD_forestFC$GeneName, levels=TDvsNonTD_diffTopFC$GeneName[order(TDvsNonTD_diffTopFC$Meta_logFC)], ordered=TRUE)
TDvsNonTD_forestFC$marbles_logFC <- as.numeric(TDvsNonTD_forestFC$marbles_logFC)
TDvsNonTD_forestFC$marbles_ErrBar.L <- as.numeric(TDvsNonTD_forestFC$marbles_ErrBar.L)
TDvsNonTD_forestFC$marbles_ErrBar.R <- as.numeric(TDvsNonTD_forestFC$marbles_ErrBar.R)
TDvsNonTD_forestFC$marbles_pValue <- as.numeric(TDvsNonTD_forestFC$marbles_pValue)
colnames(TDvsNonTD_forestFC) <- c("Gene", "log2FC", "ErrBar.L", "ErrBar.R", "p-value", "Study")
TDvsNonTD_forestFC$logp <- -log10(TDvsNonTD_forestFC$`p-value`)

gg <- ggplot(data = TDvsNonTD_forestFC)
gg +
        geom_vline(xintercept=0, color="black", lty=2, size=1.05) +
        geom_errorbarh(aes(x=log2FC, xmin=ErrBar.L, xmax=ErrBar.R, y=Gene, color=Study), position=position_dodgev(height=0.8), size = 1, height=0) +
        geom_point(aes(x = log2FC, y = Gene, color = Study, size=Study), position=position_dodgev(height=0.8)) +
        theme_bw(base_size = 26) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.28, 1.03), 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "black"),
              legend.key.height = unit(2, "lines"), legend.title=element_blank(), legend.direction="horizontal",
              plot.margin = unit(c(3,1,1,1), "lines"), 
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22, face = "italic"),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title.x = element_text(color = "black", size = 24), 
              axis.title.y = element_blank()) +
        scale_color_manual(breaks=c("MARBLES", "EARLI", "Meta"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366", "Meta" = "#009966")) +
        scale_size_manual(breaks=c("MARBLES", "EARLI", "Meta"), values = c("MARBLES" = 2.5, "EARLI" = 2.5, "Meta" = 3.25)) +
        xlab(expression(log[2]*"(Fold Change)")) +
        coord_cartesian(xlim=c(-0.45, 0.6)) +
        scale_x_continuous(breaks=pretty_breaks(n=6))
ggsave("Figures/EARLI MARBLES Meta-analysis with METAL TD vs NonTD B123 by FC.png", dpi = 600, width = 10, height = 10, units = "in")
