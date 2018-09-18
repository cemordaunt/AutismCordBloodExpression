# Meta-Analysis Functional Enrichment Analyses ####
# MARBLES and EARLI
# Charles Mordaunt
# 6/29/18

setwd()

# Packages ####
library(biomaRt)
library(RDAVIDWebService)
library(WebGestaltR)
library(reshape)
library(ggplot2)
library(scales)
library(GeneOverlap)
library(dplyr)

# Functions ####
webGestalt <- function(interestGene, referenceGene=NULL, testGeneSets){
        enrichResult <- NULL
        for(i in 1:length(testGeneSets)){
                cat("\n\nRunning WebGestalt ORA for", testGeneSets[i], "\n")
                temp <- NULL
                temp <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", enrichDatabase=testGeneSets[i], 
                                    interestGene=interestGene, interestGeneType="entrezgene", referenceGene=referenceGene, 
                                    referenceGeneType="entrezgene", minNum=5, maxNum=2000, fdrMethod="BH", 
                                    sigMethod="fdr", fdrThr=0.1, dNum=100, is.output=FALSE, methodType="R",
                                    dagColor="binary",hostName="http://www.webgestalt.org/")
                if(!is.null(temp) & !is.null(dim(temp))){
                        cat("\nSignificant gene set is identified based on FDR 0.1!\n")
                        if(testGeneSets[i] %in% c("chromosomalLocation_CytogenicBand", "network_miRNA_target", 
                                                  "network_Transcription_Factor_target", "network_CPTAC_Proteomics_OV",
                                                  "network_CPTAC_Proteomics_COADREAD", "community-contributed_Hallmark50")){
                                temp$description <- temp$geneset
                                temp <- temp[,c(1,dim(temp)[2], 2:(dim(temp)[2]-1))]
                        }
                        temp$Database <- rep(testGeneSets[i], dim(temp)[1])
                        enrichResult <- rbind(enrichResult, temp)
                }
        }
        enrichResult
}

# Load Data ####
ASD <- read.delim("Tables/ASD Meta-Analysis with METAL All Genes B123.txt", sep="\t", header=TRUE)
nonTD <- read.delim("Tables/NonTD Meta-Analysis with METAL All Genes B123.txt", sep="\t", header=TRUE)
all_probes <- as.character(ASD$Probe)
all_genes <- sort(unique(as.character(ASD$GeneName)))

# Subset Differential Genes ####
# Meta ASD
ASD_meta_diff <- subset(ASD, abs(ASD$Meta_logFC) > 0.1 & Meta_pValue < 0.01)
ASD_meta_diff <- subset(ASD_meta_diff, (marbles_logFC < 0 & earli_logFC < 0) | (marbles_logFC > 0 & earli_logFC > 0))
ASD_meta_up <- subset(ASD_meta_diff, Meta_logFC > 0)
ASD_meta_down <- subset(ASD_meta_diff, Meta_logFC < 0)
ASD_meta_diff_probes <- as.character(ASD_meta_diff$Probe)
ASD_meta_up_probes <- as.character(ASD_meta_up$Probe)
ASD_meta_down_probes <- as.character(ASD_meta_down$Probe)
ASD_meta_diff_genes <- sort(unique(as.character(ASD_meta_diff$GeneName)))
ASD_meta_up_genes <- sort(unique(as.character(ASD_meta_up$GeneName)))
ASD_meta_down_genes <- sort(unique(as.character(ASD_meta_down$GeneName)))

# Meta NonTD
nonTD_meta_diff <- subset(nonTD, abs(nonTD$Meta_logFC) > 0.1 & Meta_pValue < 0.01)
nonTD_meta_diff <- subset(nonTD_meta_diff, (marbles_logFC < 0 & earli_logFC < 0) | (marbles_logFC > 0 & earli_logFC > 0))
nonTD_meta_up <- subset(nonTD_meta_diff, Meta_logFC > 0)
nonTD_meta_down <- subset(nonTD_meta_diff, Meta_logFC < 0)
nonTD_meta_diff_probes <- as.character(nonTD_meta_diff$Probe)
nonTD_meta_up_probes <- as.character(nonTD_meta_up$Probe)
nonTD_meta_down_probes <- as.character(nonTD_meta_down$Probe)
nonTD_meta_diff_genes <- sort(unique(as.character(nonTD_meta_diff$GeneName)))
nonTD_meta_up_genes <- sort(unique(as.character(nonTD_meta_up$GeneName)))
nonTD_meta_down_genes <- sort(unique(as.character(nonTD_meta_down$GeneName)))
genes_array <-list("ASD_meta_diff_genes"=ASD_meta_diff_genes, "ASD_meta_up_genes"=ASD_meta_up_genes, "ASD_meta_down_genes"=ASD_meta_down_genes, 
                   "nonTD_meta_diff_genes"=nonTD_meta_diff_genes, "nonTD_meta_up_genes"=nonTD_meta_up_genes, 
                   "nonTD_meta_down_genes"=nonTD_meta_down_genes)
genes_array <- lapply(genes_array, function(x){gsub(" ", "", x, fixed = TRUE)})

# Get Gene IDs from Ensembl ####
probes <- list("ASD_meta_diff_probes"=ASD_meta_diff_probes, "ASD_meta_up_probes"=ASD_meta_up_probes, "ASD_meta_down_probes"=ASD_meta_down_probes, 
               "nonTD_meta_diff_probes"=nonTD_meta_diff_probes, "nonTD_meta_up_probes"=nonTD_meta_up_probes, 
               "nonTD_meta_down_probes"=nonTD_meta_down_probes, "all_probes"=all_probes)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="apr2018.archive.ensembl.org", dataset="hsapiens_gene_ensembl", ensemblRedirect = FALSE, verbose=TRUE, archive=FALSE)
genes <- sapply(probes, function(x){getBM(attributes="hgnc_symbol", filters="affy_hugene_2_0_st_v1", values=x, mart=ensembl, verbose=TRUE)})
#write.table(sort(unique(genes$ASD_meta_diff_probes.hgnc_symbol)),file="Tables/ASD_Meta_Differential_Genes_ensembl_B123.txt", sep="\t", quote=FALSE, row.names=FALSE, 
#            col.names=FALSE)
#write.table(sort(unique(genes$nonTD_meta_diff_probes.hgnc_symbol)),file="Tables/NonTD_Meta_Differential_Genes_ensembl_B123.txt", sep="\t", quote=FALSE, row.names=FALSE, 
#            col.names=FALSE)

entrezIDs <- sapply(probes, function(x){getBM(attributes="entrezgene", filters="affy_hugene_2_0_st_v1", values=x, mart=ensembl, verbose=TRUE)})
entrezIDs_background <- as.character(entrezIDs$all_probes.entrezgene)
entrezIDs <- entrezIDs[1:6]
entrezIDs <- sapply(entrezIDs, unique)
entrezIDs_background <- unique(entrezIDs_background)

rm(all_probes, ASD_meta_diff_probes, ASD_meta_down_probes, ASD_meta_up_probes, nonTD_meta_diff_probes,
   nonTD_meta_down_probes, nonTD_meta_up_probes, ensembl)

# Run WebGestalt ORA from entrezIDs ####
allGeneSets <- listGeneSet(organism="hsapiens")
testGeneSets <- allGeneSets[c(1:11,15:17,51:59)]

# Meta Analysis Genes
entrezIDs_meta <- entrezIDs
WebGestalt_Results_meta <- lapply(entrezIDs_meta, function(x){webGestalt(interestGene=x, referenceGene=entrezIDs_background, testGeneSets=testGeneSets)})
names(WebGestalt_Results_meta) <- sapply(strsplit(names(WebGestalt_Results_meta), ".entrezgene"), function(x) x[1])
WebGestalt_Results_meta2 <- WebGestalt_Results_meta[sapply(WebGestalt_Results_meta, length) > 0]

WebGestalt_Sum_meta <- NULL
for(i in 1:length(WebGestalt_Results_meta2)){
        temp <- as.data.frame(WebGestalt_Results_meta2[[i]])
        temp$GeneList <- as.character(rep(names(WebGestalt_Results_meta2)[i], dim(WebGestalt_Results_meta2[[i]])[1]))
        WebGestalt_Sum_meta <- rbind(WebGestalt_Sum_meta, temp)
}
write.table(WebGestalt_Sum_meta, "Enrichment Gene Lists/MARBLES-EARLI Meta-Analysis Enrichment WebGestalt ORA Results B123.txt", quote=FALSE, row.names=FALSE, 
            sep = "\t")
#Get ORA genes
WebGestalt_Sum_meta_ASD_sig <- subset(WebGestalt_Sum_meta, grepl("ASD", WebGestalt_Sum_meta$GeneList) & FDR < 0.05)
ora_entrezIDs <- strsplit(WebGestalt_Sum_meta_ASD_sig$overlapGene, ";")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
ora_symbols <- lapply(ora_entrezIDs, function(x){getBM(attributes="hgnc_symbol", filters="entrezgene", values=x, mart=ensembl)})
ora_symbols <- lapply(ora_symbols, function(x) {as.character(unlist(x))})
WebGestalt_Sum_meta_ASD_sig$hgnc_symbol <- as.character(ora_symbols)
write.table(WebGestalt_Sum_meta_ASD_sig, "Enrichment Gene Lists/MARBLES-EARLI Meta-Analysis Enrichment WebGestalt ORA Results ASD sig B123.txt", quote=FALSE, row.names=FALSE, 
            sep = "\t")

# Overlap SFARI with GeneOverlap ####
sfari <- read.csv("Enrichment Gene Lists/SFARI-Gene_genes_export18-10-2017.csv", stringsAsFactors=FALSE)
sfari_hyp <- sort(unique(as.character(sfari$gene.symbol[sfari$gene.score <= 5])))
sfari_min <- sort(unique(as.character(sfari$gene.symbol[sfari$gene.score <= 4])))
sfari_sugg <- sort(unique(as.character(sfari$gene.symbol[sfari$gene.score <= 3])))
sfari_strong <- sort(unique(as.character(sfari$gene.symbol[sfari$gene.score <= 2])))
sfari_high <- sort(unique(as.character(sfari$gene.symbol[sfari$gene.score <= 1])))
sfari_synd <- sort(unique(as.character(sfari$gene.symbol[sfari$syndromic == 1])))
sfari_genes <- list("Hypothesized"=sfari_hyp, "Minimal"=sfari_min, "Suggestive"=sfari_sugg, "Strong"=sfari_strong, "HighConf"=sfari_high, "Syndromic"=sfari_synd)

names(genes_array) <- c("ASD_Diff", "ASD_Up", "ASD_Down", "NonTD_Diff", "NonTD_Up", "NonTD_Down")
gom <- newGOM(genes_array, sfari_genes, genome.size = length(all_genes)) # genome.size = genes on the array by Affy annotation
oddsRatio <- getMatrix(gom, "odds.ratio") 
pValue <- getMatrix(gom, "pval")
fdr <- matrix(p.adjust(pValue, method="fdr"), nrow=nrow(pValue), ncol=ncol(pValue), byrow=FALSE)
dimnames(fdr) <- dimnames(pValue)
intersects <- getMatrix(gom, "intersection")
intersects_genes <- getNestedList(gom, "intersection")
overlapResults_sfari <- as.data.frame(cbind(intersects, oddsRatio, pValue, fdr))
colnames(overlapResults_sfari) <- paste(rep(c("intersect", "odds_ratio", "p_value", "q_value"), each=ncol(intersects)), colnames(overlapResults_sfari), sep="_")
write.table(overlapResults_sfari, "Enrichment Gene Lists/MARBLES-EARLI Meta-Analysis SFARI GeneOverlap Results.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# oddsRatio Heatmap
oddsRatio <- melt(oddsRatio)
colnames(oddsRatio) <- c("Meta_Gene", "SFARI_Gene", "oddsRatio")
oddsRatio$Meta_Gene <- factor(oddsRatio$Meta_Gene, levels=rev(names(genes_array)), ordered=TRUE)
oddsRatio$SFARI_Gene <- factor(oddsRatio$SFARI_Gene, levels=names(sfari_genes), ordered=TRUE)

intersects <- getMatrix(gom, "intersection")
intersects <- melt(intersects)
colnames(intersects) <- c("Meta_Gene", "SFARI_Gene", "intersects")
intersects$Meta_Gene <- factor(intersects$Meta_Gene, levels=rev(names(genes_array)), ordered=TRUE)
intersects$SFARI_Gene <- factor(intersects$SFARI_Gene, levels=names(sfari_genes), ordered=TRUE)
oddsRatio$intersects <- intersects$intersects

gg <- ggplot(data = oddsRatio)
gg +
        geom_tile(aes(x = SFARI_Gene, y = Meta_Gene, fill = oddsRatio)) +
        geom_text(aes(x = SFARI_Gene, y = Meta_Gene, label = intersects), color="white", size=7) +
        scale_fill_gradientn("Odds Ratio\n", colors = c("#0000FF", "Black", "#FF0000"), values = c(0,0.5,1), na.value = "Black", limits=c(0,2)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1), legend.key = element_blank(), legend.text=element_text(size=14),
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.8), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 14, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_text(size=16, color="Black"), 
              legend.title = element_text(size = 16)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        xlab("SFARI Gene") +
        ylab("Meta-Analysis Gene")
ggsave("Figures/EARLI MARBLES Meta Gene SFARI GeneOverlap OddsRatio Heatmap.png", dpi = 600, width = 8, height = 6, units = "in")

# Primate-specific WebGestalt Setup ####
primate <- read.csv("Tables/Zhang 2011 primate-specific genes.csv", header=TRUE, stringsAsFactors = FALSE)
primate_ENSG <- list("Frog"=sort(unique(as.character(primate$Gene[primate$Branch >= 1]))),
                     "Lizard"=sort(unique(as.character(primate$Gene[primate$Branch >= 2]))),
                     "Platypus"=sort(unique(as.character(primate$Gene[primate$Branch >= 3]))),
                     "Opossum"=sort(unique(as.character(primate$Gene[primate$Branch >= 4]))),
                     "Armadillo"=sort(unique(as.character(primate$Gene[primate$Branch >= 5]))),
                     "Dog"=sort(unique(as.character(primate$Gene[primate$Branch >= 6]))),
                     "Mouse"=sort(unique(as.character(primate$Gene[primate$Branch >= 7]))),
                     "Marmoset"=sort(unique(as.character(primate$Gene[primate$Branch >=8]))),
                     "Rhesus"=sort(unique(as.character(primate$Gene[primate$Branch >=9]))),
                     "Orangutan"=sort(unique(as.character(primate$Gene[primate$Branch >= 10]))),
                     "Chimp"=sort(unique(as.character(primate$Gene[primate$Branch >= 11]))),
                     "Human"=sort(unique(as.character(primate$Gene[primate$Branch >= 12]))))

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="apr2018.archive.ensembl.org", ensemblRedirect = FALSE, verbose=TRUE, archive=FALSE)
datasets <- listDatasets(ensembl)
attributes <- listAttributes(mart=ensembl)
ENSGs <- sapply(probes, function(x){getBM(attributes="ensembl_gene_id", filters="affy_hugene_2_0_st_v1", values=x, mart=ensembl, verbose=TRUE)})
ENSGs_background <- as.character(ENSGs$all_probes.ensembl_gene_id)
ENSGs <- ENSGs[1:6]
ENSGs <- sapply(ENSGs, unique)
ENSGs_background <- unique(ENSGs_background)
names(ENSGs) <- c("ASD_Diff", "ASD_Up", "ASD_Down", "NonTD_Diff", "NonTD_Up", "NonTD_Down")
primate_entrezID <- sapply(primate_ENSG, function(x){getBM(attributes="entrezgene", filters="ensembl_gene_id", values=x, mart=ensembl)})
primate_entrezID <- sapply(primate_entrezID, unique)
for(i in 1:length(primate_entrezID)){
        write.table(t(primate_entrezID[i]), file="Enrichment Gene Lists/primate_genes.gmt", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
}

# Primate-specific WebGestalt ORA ####
primate_WebGestalt_ORA <- lapply(entrezIDs, function(x) {WebGestaltR(enrichMethod="ORA", organism="hsapiens", enrichDatabase="others", enrichDatabaseFile = "Enrichment Gene Lists/primate_genes.gmt",
                                                                       enrichDatabaseType = "entrezgene", interestGene=x, interestGeneType="entrezgene", 
                                                                       referenceGene=entrezIDs_background, referenceGeneType="entrezgene", minNum=1, maxNum=10000, fdrMethod="BH", 
                                                                       sigMethod="top", topThr=15, dNum=100, is.output=FALSE, methodType="R",
                                                                       dagColor="binary",hostName="http://www.webgestalt.org/")})
primate_WebGestalt_ORA_sum <- NULL
for(i in 1:length(primate_WebGestalt_ORA)){
        temp <- as.data.frame(primate_WebGestalt_ORA[[i]])
        temp$GeneList <- as.character(rep(names(primate_WebGestalt_ORA)[i], dim(primate_WebGestalt_ORA[[i]])[1]))
        primate_WebGestalt_ORA_sum <- rbind(primate_WebGestalt_ORA_sum, temp)
}
write.table(primate_WebGestalt_ORA_sum, "Enrichment Gene Lists/MARBLES-EARLI Meta-Analysis Enrichment primate WebGestalt ORA Results B123.txt", quote=FALSE, row.names=FALSE, 
            sep = "\t")

# PGC GeneOverlap ####
pgc <- read.csv("Enrichment Gene Lists/PGC Autism.csv", stringsAsFactors = FALSE)
OR <- strsplit(pgc$ODDS.RATIO..95.CI., " ")
OR <- as.numeric(lapply(OR, function(x){x[1]}))
pgc$OR <- OR

# All Genes
pgc_genes_all <- as.character(pgc$GENES.IN.REGION)
pgc_genes_all <- strsplit(pgc_genes_all, " ")
pgc_genes_all <- melt(pgc_genes_all)
pgc_genes_all <- sort(unique(as.character(pgc_genes_all$value)))

# Up Genes OR > 1
pgc_genes_up <- as.character(pgc$GENES.IN.REGION[pgc$OR > 1])
pgc_genes_up <- strsplit(pgc_genes_up, " ")
pgc_genes_up <- melt(pgc_genes_up)
pgc_genes_up <- sort(unique(as.character(pgc_genes_up$value)))

# Down Genes OR < 1
pgc_genes_down <- as.character(pgc$GENES.IN.REGION[pgc$OR < 1])
pgc_genes_down <- strsplit(pgc_genes_down, " ")
pgc_genes_down <- melt(pgc_genes_down)
pgc_genes_down <- sort(unique(as.character(pgc_genes_down$value)))

# Entrez IDs
pgc_genes <- list("All"=pgc_genes_all, "Up"=pgc_genes_up, "Down"=pgc_genes_down)
entrezIDs_pgc <- sapply(pgc_genes, function(x){getBM(attributes="entrezgene", filters="hgnc_symbol", values=x, mart=ensembl)})
for(i in 1:length(entrezIDs_pgc)){
        write.table(t(entrezIDs_pgc[i]), file="Enrichment Gene Lists/pgc_genes.gmt", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
}

# ORA
WebGestalt_Results_pgc <- lapply(entrezIDs, function(x) {WebGestaltR(enrichMethod="ORA", organism="hsapiens", enrichDatabase="others", enrichDatabaseFile = "Enrichment Gene Lists/pgc_genes.gmt",
                                                                       enrichDatabaseType = "entrezgene", interestGene=x, interestGeneType="entrezgene", 
                                                                       referenceGene=entrezIDs_background, referenceGeneType="entrezgene", minNum=1, maxNum=2000, fdrMethod="BH", 
                                                                       sigMethod="fdr", fdrThr=1, dNum=100, is.output=FALSE, methodType="R",
                                                                       dagColor="binary",hostName="http://www.webgestalt.org/")})
# None enriched

pgc_and_diff_genes <- NULL
for(i in 1:length(genes_array)){
        for(j in 1:length(pgc_genes)){
                temp <- c(names(genes_array)[i],names(pgc_genes)[j],intersect(genes_array[[i]], pgc_genes[[j]]))
                pgc_and_diff_genes <- rbind(pgc_and_diff_genes, temp)
        }
}

# PGC Genes GeneOverlap ####
gom <- newGOM(genes_array, pgc_genes, genome.size = length(all_genes)) # genome.size = genes on the array by Affy annotation
oddsRatio <- getMatrix(gom, "odds.ratio") 
pValue <- getMatrix(gom, "pval")
fdr <- matrix(p.adjust(pValue, method="fdr"), nrow=nrow(pValue), ncol=ncol(pValue), byrow=FALSE)
dimnames(fdr) <- dimnames(pValue)
intersects <- getMatrix(gom, "intersection")
intersects_genes <- getNestedList(gom, "intersection")
overlapResults_pgc <- as.data.frame(cbind(intersects, oddsRatio, pValue, fdr))
colnames(overlapResults_pgc) <- paste(rep(c("intersect", "odds_ratio", "p_value", "q_value"), each=ncol(intersects)), colnames(overlapResults_pgc), sep="_")
write.table(overlapResults_pgc, "Enrichment Gene Lists/MARBLES-EARLI Meta-Analysis pgc GeneOverlap Results.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# oddsRatio Heatmap
oddsRatio <- melt(oddsRatio)
colnames(oddsRatio) <- c("Meta_Gene", "pgc_Gene", "oddsRatio")
oddsRatio$Meta_Gene <- factor(oddsRatio$Meta_Gene, levels=rev(names(genes_array)), ordered=TRUE)
oddsRatio$pgc_Gene <- factor(oddsRatio$pgc_Gene, levels=names(pgc_genes), ordered=TRUE)

intersects <- getMatrix(gom, "intersection")
intersects <- melt(intersects)
colnames(intersects) <- c("Meta_Gene", "pgc_Gene", "intersects")
intersects$Meta_Gene <- factor(intersects$Meta_Gene, levels=rev(names(genes_array)), ordered=TRUE)
intersects$pgc_Gene <- factor(intersects$pgc_Gene, levels=names(pgc_genes), ordered=TRUE)
oddsRatio$intersects <- intersects$intersects

gg <- ggplot(data = oddsRatio)
gg +
        geom_tile(aes(x = pgc_Gene, y = Meta_Gene, fill = oddsRatio)) +
        geom_text(aes(x = pgc_Gene, y = Meta_Gene, label = intersects), color="white", size=7) +
        scale_fill_gradientn("Odds Ratio\n", colors = c("#0000FF", "Black", "#FF0000"), values = c(0,0.13,1), na.value = "Black", limits=c(0,6)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1), legend.key = element_blank(), legend.text=element_text(size=14),
              panel.grid.minor = element_blank(), legend.position = c(1.23, 0.8), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 14, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_text(size=16, color="Black"), 
              legend.title = element_text(size = 16)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        xlab("PGC Gene") +
        ylab("Meta-Analysis Gene")
ggsave("Figures/EARLI MARBLES Meta Gene pgc GeneOverlap OddsRatio Heatmap.png", dpi = 600, width = 7, height = 6, units = "in")

# Diff Gene Expression from other studies in blood and brain ####
# Load Gene Lists
# Gupta 2014 Cortex RNAseq p < 0.01 (gene symbol)
cortex1_up <- as.character(unlist(read.delim("Enrichment Gene Lists/Gupta_2014_Cortex_Up_RNAseq_p01.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
cortex1_down <- as.character(unlist(read.delim("Enrichment Gene Lists/Gupta_2014_Cortex_Down_RNAseq_p01.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
cortex1_all <- sort(unique(as.character(c(cortex1_up, cortex1_down))))

# Parikshak 2016 Cortex RNAseq p < 0.01 (ENSG)
cortex2_up <- as.character(unlist(read.delim("Enrichment Gene Lists/Parikshak_2016_Cortex_RNAseq_Up_p01.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
cortex2_down <- as.character(unlist(read.delim("Enrichment Gene Lists/Parikshak_2016_Cortex_RNAseq_Down_p01.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
cortex2_all <- sort(unique(as.character(c(cortex2_up, cortex2_down))))

# Parikshak 2016 Cerebellum RNAseq p < 0.01 (ENSG)
cereb_up <- as.character(unlist(read.delim("Enrichment Gene Lists/Parikshak_2016_Cerebellum_RNAseq_Up_p01.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
cereb_down <- as.character(unlist(read.delim("Enrichment Gene Lists/Parikshak_2016_Cerebellum_RNAseq_Down_p01.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
cereb_all <- sort(unique(as.character(c(cereb_up, cereb_down))))

# Tylee 2017 Blood Array Mega-Analysis p < 0.01 (gene symbol)
blood_up <- as.character(unlist(read.delim("Enrichment Gene Lists/Tylee_2017_Blood_Array_Up_p01.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
blood_down <- as.character(unlist(read.delim("Enrichment Gene Lists/Tylee_2017_Blood_Array_Down_p01.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
blood_all <- sort(unique(as.character(c(blood_up, blood_down))))

# Tylee 2017 Lymphoblast cell line RNAseq p < 0.05 (gene symbol)
lcl_up <- as.character(unlist(read.delim("Enrichment Gene Lists/Tylee_2017_LCL_RNAseq_Up_p05.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
lcl_down <- as.character(unlist(read.delim("Enrichment Gene Lists/Tylee_2017_LCL_RNAseq_Down_p05.txt", sep="\n", header=FALSE, stringsAsFactors = FALSE)))
lcl_all <- sort(unique(as.character(c(lcl_up, lcl_down))))

brain_blood_genes <- list("Cortex1_All"=cortex1_all, "Cortex1_Up"=cortex1_up, "Cortex1_Down"=cortex1_down,
                          "Blood_All"=blood_all, "Blood_Up"=blood_up, "Blood_Down"=blood_down,
                          "LCL_All"=lcl_all, "LCL_Up"=lcl_up, "LCL_Down"=lcl_down)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="apr2018.archive.ensembl.org", dataset="hsapiens_gene_ensembl", ensemblRedirect = FALSE, verbose=TRUE, archive=FALSE)
brain_blood_ENSG <- sapply(brain_blood_genes, function(x){getBM(attributes="ensembl_gene_id", filters="external_gene_name", values=x, mart=ensembl)})
brain_blood_ENSG <- list("Cortex1_All"=brain_blood_ENSG[[1]], "Cortex1_Up"=brain_blood_ENSG[[2]], "Cortex1_Down"=brain_blood_ENSG[[3]], 
                         "Cortex2_All"=cortex2_all, "Cortex2_Up"=cortex2_up, "Cortex2_Down"=cortex2_down,
                         "Cerebellum_All"=cereb_all, "Cerebellum_Up"=cereb_up, "Cerebellum_Down"=cereb_down, 
                         "Blood_All"=brain_blood_ENSG[[4]], "Blood_Up"=brain_blood_ENSG[[5]], "Blood_Down"=brain_blood_ENSG[[6]], 
                         "LCL_All"=brain_blood_ENSG[[7]], "LCL_Up"=brain_blood_ENSG[[8]], "LCL_Down"=brain_blood_ENSG[[9]])
brain_blood_ENSG <- sapply(brain_blood_ENSG, unique)

# brain_blood Genes GeneOverlap ####
gom <- newGOM(ENSGs, brain_blood_ENSG, genome.size = length(ENSGs_background)) # genome.size = ENSGs from array probes
oddsRatio <- getMatrix(gom, "odds.ratio") 
pValue <- getMatrix(gom, "pval")
fdr <- matrix(p.adjust(pValue, method="fdr"), nrow=nrow(pValue), ncol=ncol(pValue), byrow=FALSE)
dimnames(fdr) <- dimnames(pValue)
intersects <- getMatrix(gom, "intersection")
intersects_genes <- getNestedList(gom, "intersection")
overlapResults_brain_blood <- as.data.frame(cbind(intersects, oddsRatio, pValue, fdr))
colnames(overlapResults_brain_blood) <- paste(rep(c("intersect", "odds_ratio", "p_value", "q_value"), each=ncol(intersects)), colnames(overlapResults_brain_blood), sep="_")
write.table(overlapResults_brain_blood, "Enrichment Gene Lists/MARBLES-EARLI Meta-Analysis brain_blood GeneOverlap Results.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# oddsRatio Heatmap
oddsRatio <- melt(oddsRatio)
colnames(oddsRatio) <- c("Meta_Gene", "brain_blood_Gene", "oddsRatio")
oddsRatio$Meta_Gene <- factor(oddsRatio$Meta_Gene, levels=rev(names(ENSGs)), ordered=TRUE)
oddsRatio$brain_blood_Gene <- factor(oddsRatio$brain_blood_Gene, levels=names(brain_blood_ENSG), ordered=TRUE)

intersects <- getMatrix(gom, "intersection")
intersects <- melt(intersects)
colnames(intersects) <- c("Meta_Gene", "brain_blood_Gene", "intersects")
intersects$Meta_Gene <- factor(intersects$Meta_Gene, levels=rev(names(ENSGs)), ordered=TRUE)
intersects$brain_blood_Gene <- factor(intersects$brain_blood_Gene, levels=names(brain_blood_ENSG), ordered=TRUE)
oddsRatio$intersects <- intersects$intersects

gg <- ggplot(data = oddsRatio)
gg +
        geom_tile(aes(x = brain_blood_Gene, y = Meta_Gene, fill = oddsRatio)) +
        geom_text(aes(x = brain_blood_Gene, y = Meta_Gene, label = intersects), color="white", size=7) +
        scale_fill_gradientn("Odds Ratio\n", colors = c("#0000FF", "Black", "#FF0000"), values = c(0,0.13,1), na.value = "Black", limits=c(0,6)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1), legend.key = element_blank(), legend.text=element_text(size=14),
              panel.grid.minor = element_blank(), legend.position = c(1.15, 0.77), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 14, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_text(size=16, color="Black"), 
              legend.title = element_text(size = 16)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        xlab("Differentially-Expressed Genes") +
        ylab("Meta-Analysis Gene")
ggsave("Figures/EARLI MARBLES Meta Gene brain_blood GeneOverlap OddsRatio Heatmap.png", dpi = 600, width = 9, height = 6, units = "in")

# Get Intersecting genes in ASD
brain_blood_genessymbols <- NULL
for(i in 1:length(intersects_genes)){
        cat(i, "\n")
        temp <- NULL
        gene_symbols <- NULL
        temp <- intersects_genes[[i]]
        gene_symbols <- sapply(temp, function(x){
                if(length(x) > 0){
                        as.character(unlist(getBM(attributes="external_gene_name", filters="ensembl_gene_id", values=x, mart=ensembl, verbose=TRUE)))
                } else {
                        NULL
                }
        })
        brain_blood_genessymbols[[i]] <- gene_symbols
}
names(brain_blood_genessymbols) <- names(intersects_genes)

# Blood Cell Type GeneOverlap ####
cellGenes <- read.csv("Tables/LM22 DEGs.csv", header=TRUE, stringsAsFactors = FALSE)

# Get Gene Lists
cellList <- lapply(cellGenes[,2:ncol(cellGenes)], function(x){
        temp <- NULL
        temp <- cellGenes$Gene[x == 1]
        temp %>% as.character %>% unique %>% sort
        return(temp)
})
rm(cellGenes)

# Get HGNC Gene Symbols
genes <- sapply(probes, function(x){getBM(attributes="hgnc_symbol", filters="affy_hugene_2_0_st_v1", values=x, mart=ensembl, 
                                          verbose=TRUE)})
genes <- sapply(genes, function(x) x %>% unique %>% sort)
names(genes) <- c("ASD_Diff", "ASD_Up", "ASD_Down", "nonTD_Diff", "nonTD_Up", "nonTD_Down")
all_genes <- getBM(attributes="hgnc_symbol", filters="affy_hugene_2_0_st_v1", values=all_probes, mart=ensembl, verbose=TRUE)
all_genes <- all_genes %>% unlist %>% as.character %>% unique %>% sort
rm(ensembl, probes, all_probes)

# Overlap Meta and Cell Type Genes
gom <- newGOM(genes, cellList, genome.size = length(all_genes)) # genome.size = genes on the array by Affy annotation
oddsRatio <- getMatrix(gom, "odds.ratio") 
pValue <- getMatrix(gom, "pval")
fdr <- matrix(p.adjust(pValue, method="fdr"), nrow=nrow(pValue), ncol=ncol(pValue), byrow=FALSE)
dimnames(fdr) <- dimnames(pValue)
intersects <- getMatrix(gom, "intersection")
intersects_genes <- getNestedList(gom, "intersection")
overlapResults <- as.data.frame(cbind(intersects, oddsRatio, pValue, fdr))
colnames(overlapResults) <- paste(rep(c("intersect", "odds_ratio", "p_value", "q_value"), each=ncol(intersects)), 
                                  colnames(overlapResults), sep="_")
write.table(overlapResults, "Enrichment Gene Lists/MARBLES-EARLI Meta-Analysis Cibersort CellType GeneOverlap Results.txt", sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)

# oddsRatio Heatmap
oddsRatio <- melt(oddsRatio)
colnames(oddsRatio) <- c("Meta_Gene", "cell_Gene", "oddsRatio")
oddsRatio$Meta_Gene <- factor(oddsRatio$Meta_Gene, levels=rev(names(genes)), ordered=TRUE)
oddsRatio$cell_Gene <- factor(oddsRatio$cell_Gene, levels=names(cellList), ordered=TRUE)
oddsRatio$intersects <- intersects$intersects

gg <- ggplot(data = oddsRatio)
gg +
        geom_tile(aes(x = cell_Gene, y = Meta_Gene, fill = oddsRatio)) +
        geom_text(aes(x = cell_Gene, y = Meta_Gene, label = intersects), color="white", size=6) +
        scale_fill_gradientn("Odds Ratio\n", colors = c("#0000FF", "Black", "#FF0000"), values = c(0,0.5,1), na.value = "Black", limits=c(0,2)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1), legend.key = element_blank(), legend.text=element_text(size=14),
              panel.grid.minor = element_blank(), legend.position = c(1.13, 0.74), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 12, color = "Black", angle = 60, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 14, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title.y = element_text(size=16, color="Black"), axis.title.x=element_blank(),
              legend.title = element_text(size = 16)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        ylab("Meta-Analysis Gene")
ggsave("Figures/EARLI MARBLES Meta Gene Cibersort Cell type GeneOverlap OddsRatio Heatmap.png", dpi = 600, width = 10, height = 6, units = "in")

# Get Cell-Type Genes
intersects_genes_df <- NULL
for(i in 1:length(intersects_genes)){
        for(j in 1:length(intersects_genes[[i]])){
                temp <- NULL
                temp <- data.frame("CellType"=rep(names(intersects_genes)[i], length(intersects_genes[[i]][[j]])),
                                   "Meta_Diff"=rep(names(intersects_genes[[i]])[j], length(intersects_genes[[i]][[j]])),
                                   "Gene"=intersects_genes[[i]][[j]])
                intersects_genes_df <- rbind(intersects_genes_df, temp)
        }
}
write.table(intersects_genes_df, "Enrichment Gene Lists/Cibersort Cell Type and Meta-Analysis Overlapping Genes.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Blood CellType GSEA ####
cellExp <- read.csv("Tables/LM22 DEG Exp.csv", header=TRUE, stringsAsFactors = FALSE)
cellGenes <- read.csv("Tables/LM22 DEGs.csv", header=TRUE, stringsAsFactors = FALSE)
table(colnames(cellExp) == colnames(cellGenes)) # All TRUE
table(cellExp$Gene == cellGenes$Gene) # All TRUE
cellTypes <- colnames(cellExp)[2:length(colnames(cellExp))] %>% as.character
diffTest <- list(NULL)
for(i in 1:length(cellTypes)){
        temp <- NULL
        cell <- NULL
        cell <- cellTypes[i]
        temp <- cellExp[cellGenes[,cell] == 1, cell] / rowMeans(cellExp[cellGenes[,cell] == 1, cellTypes[!cellTypes == cell]])
        diffTest[[i]] <- temp
}
lapply(diffTest, function(x) table(x >= 1))
cell_entrezIDs <- sapply(cellList, function(x){getBM(attributes="entrezgene", filters="hgnc_symbol", values=x, mart=ensembl)})
for(i in 1:length(cell_entrezIDs)){
        write.table(t(cell_entrezIDs[i]), file="Enrichment Gene Lists/cibersort_cell_genes.gmt", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
}
