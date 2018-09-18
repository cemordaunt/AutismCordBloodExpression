# Meta-Analysis Differentially Expressed Gene Expression Levels ####
# Charles Mordaunt
# 7/3/18

setwd()

# Packages ####
library(ggplot2)
library(scales)

# Data ####
# Load
earli_asd <- read.csv("EARLI Data/ASD_TD limma output table 092217.csv")
earli_nonTD <- read.csv("EARLI Data/Non_TD limma output table 092217.csv")
marbles <- read.delim("Tables/MARBLES All Genes Stats Dx_alg SVAonly B123.txt", sep="\t", header=TRUE)
meta_asd_diff_probes <- as.integer(as.character(unlist(read.delim("Tables/ASD Meta-Analysis with METAL Differential Probe List B123.txt", sep="\t", header=FALSE))))
meta_nonTD_diff_probes <- as.integer(as.character(unlist(read.delim("Tables/Non-TD Meta-Analysis with METAL Differential Probe List B123.txt", sep="\t", header=FALSE))))

# Meta Diff Tests
earli_asd$MetaASDdiff <- earli_asd$X %in% meta_asd_diff_probes
earli_asd$MetaASDdiff <- factor(earli_asd$MetaASDdiff, levels=c("TRUE", "FALSE"), ordered=TRUE)
earli_nonTD$MetaNonTDdiff <- earli_nonTD$X %in% meta_nonTD_diff_probes
earli_nonTD$MetaNonTDdiff <- factor(earli_nonTD$MetaNonTDdiff, levels=c("TRUE", "FALSE"), ordered=TRUE)

marbles$MetaASDdiff <- marbles$Probe %in% meta_asd_diff_probes
marbles$MetaASDdiff <- factor(marbles$MetaASDdiff, levels=c("TRUE", "FALSE"), ordered=TRUE)
marbles$MetaNonTDdiff <- marbles$Probe %in% meta_nonTD_diff_probes
marbles$MetaNonTDdiff <- factor(marbles$MetaNonTDdiff, levels=c("TRUE", "FALSE"), ordered=TRUE)

# Average Expression Density Plots ####
# EARLI TD vs ASD
gg <- ggplot()
gg +
        geom_density(data=earli_asd, aes(x = AveExpr, color = MetaASDdiff), size=1.25) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.height = unit(1.5, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.85, 0.85), 
              legend.background = element_blank(), legend.text = element_text(size = 18, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 18),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 20), 
              legend.title = element_text(color = "black", size = 20)) +
        xlab(expression(log[2]*"(Mean Expression)")) +
        ylab("Density") +
        scale_x_continuous(breaks=pretty_breaks(n=6)) +
        scale_y_continuous(expand=c(0,0), breaks=pretty_breaks(n=3)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        coord_cartesian(xlim=c(1,14), ylim=c(0,0.28))
ggsave("Figures/Meta ASD Diff in EARLI Dx_alg Mean Expression Density Plot B123.png", dpi = 600, width = 8, height = 8, units = "in")

# EARLI TD vs NonTD
gg <- ggplot()
gg +
        geom_density(data=earli_nonTD, aes(x = AveExpr, color = MetaNonTDdiff), size=1.25) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.height = unit(1.5, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.85, 0.85), 
              legend.background = element_blank(), legend.text = element_text(size = 18, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 18),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 20), 
              legend.title = element_text(color = "black", size = 20)) +
        xlab(expression(log[2]*"(Mean Expression)")) +
        ylab("Density") +
        scale_x_continuous(breaks=pretty_breaks(n=6)) +
        scale_y_continuous(expand=c(0,0), breaks=pretty_breaks(n=3)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        coord_cartesian(xlim=c(1,14), ylim=c(0,0.28))
ggsave("Figures/Meta NonTD Diff in EARLI Dx_alg Mean Expression Density Plot B123.png", dpi = 600, width = 8, height = 8, units = "in")

# MARBLES TD vs ASD
gg <- ggplot()
gg +
        geom_density(data=marbles, aes(x = AveExpr, color = MetaASDdiff), size=1.25) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.height = unit(1.5, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.85, 0.85), 
              legend.background = element_blank(), legend.text = element_text(size = 18, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 18),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 20), 
              legend.title = element_text(color = "black", size = 20)) +
        xlab(expression(log[2]*"(Mean Expression)")) +
        ylab("Density") +
        scale_x_continuous(breaks=pretty_breaks(n=6)) +
        scale_y_continuous(expand=c(0,0), breaks=pretty_breaks(n=3)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        coord_cartesian(xlim=c(1,14), ylim=c(0,0.28))
ggsave("Figures/Meta ASD Diff in MARBLES Dx_alg Mean Expression Density Plot B123.png", dpi = 600, width = 8, height = 8, units = "in")

# MARBLES TD vs NonTD
gg <- ggplot()
gg +
        geom_density(data=marbles, aes(x = AveExpr, color = MetaNonTDdiff), size=1.25) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.height = unit(1.5, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.85, 0.85), 
              legend.background = element_blank(), legend.text = element_text(size = 18, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 18),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 20), 
              legend.title = element_text(color = "black", size = 20)) +
        xlab(expression(log[2]*"(Mean Expression)")) +
        ylab("Density") +
        scale_x_continuous(breaks=pretty_breaks(n=6)) +
        scale_y_continuous(expand=c(0,0), breaks=pretty_breaks(n=3)) +
        scale_color_manual(name = "Differentially\nExpressed\n", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        coord_cartesian(xlim=c(1,14), ylim=c(0,0.28))
ggsave("Figures/Meta NonTD Diff in MARBLES Dx_alg Mean Expression Density Plot B123.png", dpi = 600, width = 8, height = 8, units = "in")

# Median Expression Stats ####
# EARLI TD vs ASD
ASDmedians <- NULL
for(i in 1:10000){
        temp <- sample(x=earli_asd$AveExpr, size=sum(earli_asd$MetaASDdiff == "TRUE"))
        ASDmedians <- c(ASDmedians, median(temp))
}
ASDmedians <- as.data.frame(ASDmedians)
ASDtrueMedian <- median(earli_asd$AveExpr[earli_asd$MetaASDdiff == "TRUE"]) #4.34

ASDp <- 2*(1-ecdf(as.numeric(unlist(ASDmedians$ASDmedians)))(ASDtrueMedian)) #0.516

gg <- ggplot()
gg +
        geom_histogram(data=ASDmedians, aes(x = ASDmedians), fill = "#3366CC", bins=50) +
        geom_vline(xintercept = ASDtrueMedian, color="#FF3366", size=1.25) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.9, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 16, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 16),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 18), 
              legend.title = element_text(color = "black", size = 18)) +
        xlab(expression("Median "*log[2]*"(Mean Expression)")) +
        ylab("Null Gene Sets") +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4), expand=c(0,0)) +
        coord_cartesian(xlim=c(3.5,5.5), ylim=c(0,805)) +
        annotate("text", x=4.4, y=765, label=paste(expression(paste(italic(p), " = 0.52"))), hjust=0, size=7, parse=TRUE)
ggsave("Figures/Meta ASD Diff in EARLI Dx_alg Sampled Median Mean Expression Histogram B123.png", dpi = 600, width = 8, height = 8, units = "in")

# EARLI TD vs NonTD
NonTDmedians <- NULL
for(i in 1:10000){
        temp <- sample(x=earli_nonTD$AveExpr, size=sum(earli_nonTD$MetaNonTDdiff == "TRUE"))
        NonTDmedians <- c(NonTDmedians, median(temp))
}
NonTDmedians <- as.data.frame(NonTDmedians)
NonTDtrueMedian <- median(earli_nonTD$AveExpr[earli_nonTD$MetaNonTDdiff == "TRUE"])

NonTDp <- 2*ecdf(as.numeric(unlist(NonTDmedians$NonTDmedians)))(NonTDtrueMedian) #0.8978

gg <- ggplot()
gg +
        geom_histogram(data=NonTDmedians, aes(x = NonTDmedians), fill = "#3366CC", bins=50) +
        geom_vline(xintercept = NonTDtrueMedian, color="#FF3366", size=1.25) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.9, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 16, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 16),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 18), 
              legend.title = element_text(color = "black", size = 18)) +
        xlab(expression("Median "*log[2]*"(Mean Expression)")) +
        ylab("Null Gene Sets") +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4), expand=c(0,0)) +
        coord_cartesian(xlim=c(3.5,5.5), ylim=c(0,805)) +
        annotate("text", x=4.23, y=765, label=paste(expression(paste(italic(p), " = 0.90"))), hjust=0, size=7, parse=TRUE)
ggsave("Figures/Meta NonTD Diff in EARLI Dx_alg Sampled Median Mean Expression Histogram B123.png", dpi = 600, width = 8, height = 8, units = "in")

# MARBLES TD vs ASD
ASDmedians <- NULL
for(i in 1:10000){
        temp <- sample(x=marbles$AveExpr, size=sum(marbles$MetaASDdiff == "TRUE"))
        ASDmedians <- c(ASDmedians, median(temp))
}
ASDmedians <- as.data.frame(ASDmedians)
ASDtrueMedian <- median(marbles$AveExpr[marbles$MetaASDdiff == "TRUE"])

ASDp <- 2*(1-ecdf(as.numeric(unlist(ASDmedians$ASDmedians)))(ASDtrueMedian)) #0.7372

gg <- ggplot()
gg +
        geom_histogram(data=ASDmedians, aes(x = ASDmedians), fill = "#3366CC", bins=50) +
        geom_vline(xintercept = ASDtrueMedian, color="#FF3366", size=1.25) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.9, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 16, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 16),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 18), 
              legend.title = element_text(color = "black", size = 18)) +
        xlab(expression("Median "*log[2]*"(Mean Expression)")) +
        ylab("Null Gene Sets") +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4), expand=c(0,0)) +
        coord_cartesian(xlim=c(3.5,5.5), ylim=c(0,805)) +
        annotate("text", x=4.75, y=765, label=paste(expression(paste(italic(p), " = 0.74"))), hjust=0, size=7, parse=TRUE)
ggsave("Figures/Meta ASD Diff in MARBLES Dx_alg Sampled Median Mean Expression Histogram B123.png", dpi = 600, width = 8, height = 8, units = "in")

# MARBLES TD vs NonTD
NonTDmedians <- NULL
for(i in 1:10000){
        temp <- sample(x=marbles$AveExpr, size=sum(marbles$MetaNonTDdiff == "TRUE"))
        NonTDmedians <- c(NonTDmedians, median(temp))
}
NonTDmedians <- as.data.frame(NonTDmedians)
NonTDtrueMedian <- median(marbles$AveExpr[marbles$MetaNonTDdiff == "TRUE"])

NonTDp <- 2*(ecdf(as.numeric(unlist(NonTDmedians$NonTDmedians)))(NonTDtrueMedian)) #0.6496

gg <- ggplot()
gg +
        geom_histogram(data=NonTDmedians, aes(x = NonTDmedians), fill = "#3366CC", bins=50) +
        geom_vline(xintercept = NonTDtrueMedian, color="#FF3366", size=1.25) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.9, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 16, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 16),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 18), 
              legend.title = element_text(color = "black", size = 18)) +
        xlab(expression("Median "*log[2]*"(Mean Expression)")) +
        ylab("Null Gene Sets") +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4), expand=c(0,0)) +
        coord_cartesian(xlim=c(3.5,5.5), ylim=c(0,805)) +
        annotate("text", x=4.55, y=765, label=paste(expression(paste(italic(p), " = 0.65"))), hjust=0, size=7, parse=TRUE)
ggsave("Figures/Meta NonTD Diff in MARBLES Dx_alg Sampled Median Mean Expression Histogram B123.png", dpi = 600, width = 8, height = 8, units = "in")
