# MARBLES Demographic Characteristics ####
# Charles Mordaunt
# 9/18/18

setwd()

# Data ####
cov <- read.delim("Sample Info/MARBLES Microarray Cord Sample Covariates.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
cc <- read.csv("Sample Info/Creatine Cotinine Data MARBLES.csv", header=TRUE, stringsAsFactors = FALSE)
cov <- merge(x=cov, y=cc, by="COI_ID", all.x=TRUE, all.y=FALSE, sort=FALSE)
cov <- subset(cov, Outlier == "N")
cov$ChildRaceWhite <- cov$ChildRace == 1
cov$ChildHispYN <- !cov$ChildHisp == "N"
cov$MomBirthPlaceOutsideUS <- cov$MomBirthPlace > 3
cov$DadBirthPlaceOutsideUS <- cov$DadBirthPlace > 3
cov$MomEduBachelors <- cov$MomEdu %in% 4:5
cov$DadEduBachelors <- cov$DadEdu %in% 4:5
contVars <- c("MomAgeYr", "DadAgeYr", "Mat_BMI_PrePreg","Mat_BMI_Deliv", "birthweight", "birthlength", "gest_age_deliv", "Hour_Birth", 
              "Creatine_Conc_mg_dl", "Cotinine_Extract_Conc_ng_ml", "Cotinine_Urine_Conc_ng_ml")
catVars <- c("COI_GENDER", "ChildRaceWhite", "ChildHispYN", "MomBirthPlaceOutsideUS", "DadBirthPlaceOutsideUS", "SmokeYN_PrePreg", 
             "SmokeYN_Pregnancy", "OthSmokeYN_PrePreg", "OthSmokeYN_Pregnancy", "DeliveryPayer", "HomeOwn", "MomEduBachelors", "DadEduBachelors", 
             "WantedPreg_5cat", "BirthSeason", "DeliveryMethod")
numVars <- length(contVars) + length(catVars) #27
cov_table1 <- cov[,c("COI_ID", "IBC", "Dx_alg", contVars, catVars)]
write.table(cov_table1, "Sample Info/MARBLES Microarray Cord Sample Covariates for Table 1.txt", sep="\t", quote=FALSE, col.names=TRUE, 
            row.names=FALSE)

# Continuous Variables (One-way ANOVA) ####
continuous <- cov[,c("COI_ID", "IBC", "Dx_alg", contVars)]
continuous$Dx_alg <- factor(continuous$Dx_alg)
diagnoses <- c("ASD", "non TD/ASD", "TD")
contStats <- NULL
for(i in 1:length(contVars)){
        means <- aggregate(continuous[,contVars[i]] ~ continuous$Dx_alg, FUN=mean)[,2]
        sds <- aggregate(continuous[,contVars[i]] ~ continuous$Dx_alg, FUN=sd)[,2]
        p <- summary(aov(continuous[,contVars[i]] ~ continuous$Dx_alg))[[1]][1, "Pr(>F)"]
        temp <- cbind(rep(contVars[i], 3), diagnoses, means, sds, rep(p, 3))
        contStats <- rbind(contStats, temp)
}
colnames(contStats) <- c("Variable", "Dx_alg", "Mean", "SD", "pvalue")
contStats <- as.data.frame(contStats, stringsAsFactors=FALSE)
contStats[,c("Variable", "Dx_alg")] <- lapply(contStats[,c("Variable", "Dx_alg")], as.factor)
contStats[,c("Mean", "SD", "pvalue")] <- lapply(contStats[,c("Mean", "SD", "pvalue")], as.numeric)
cont_pvalues <- aggregate(contStats$pvalue, by=list(contStats$Variable), FUN=mean)
cont_pvalues$qvalue <- p.adjust(cont_pvalues$x, method="fdr", n=numVars)
contStats$qvalue <- cont_pvalues$qvalue[match(contStats$Variable, cont_pvalues$Group.1)]
write.table(contStats, "Tables/Continuous Variables Table 1 MARBLES B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Categorical Variables (Fisher's Exact Test) ####
categorical <- cov[, c("COI_ID", "IBC", "Dx_alg", catVars)]
diagnoses <- c("ASD", "non TD/ASD", "TD")
categorical[,3:ncol(categorical)] <- lapply(categorical[,3:ncol(categorical)], function(x){as.factor(as.character(x))})
catStats <- NULL
for(i in 1:length(catVars)){
        counts <- table(categorical[,catVars[i]], categorical$Dx_alg)
        p <- fisher.test(categorical[,catVars[i]], categorical$Dx_alg, workspace=2e7)$p.value
        temp <- cbind(catVars[i], rownames(counts), counts, p)
        catStats <- rbind(catStats, temp)
}
colnames(catStats)[1:2] <- c("Variable", "Value")
catStats <- as.data.frame(catStats, stringsAsFactors = FALSE)
catStats[,c("Variable", "Value")] <- lapply(catStats[,c("Variable", "Value")], as.factor)
catStats[,c("ASD", "Non-TD", "TD")] <- lapply(catStats[,c("ASD", "Non-TD", "TD")], as.integer)
catStats$p <- as.numeric(catStats$p)
varSums <- aggregate(as.matrix(catStats[,c("ASD", "Non-TD", "TD")]) ~ Variable, data=catStats, FUN=sum)
catStats$per_ASD <- catStats$ASD / varSums$ASD[match(catStats$Variable, varSums$Variable)]
catStats$'per_Non-TD' <- catStats$'Non-TD' / varSums$'Non-TD'[match(catStats$Variable, varSums$Variable)]
catStats$per_TD <- catStats$TD / varSums$TD[match(catStats$Variable, varSums$Variable)]
cat_pvalues <- aggregate(catStats$p, by=list(catStats$Variable), FUN=mean)
cat_pvalues$qvalue <- p.adjust(cat_pvalues$x, method="fdr", n=numVars)
catStats$qvalue <- cat_pvalues$qvalue[match(catStats$Variable, cat_pvalues$Group.1)]
write.table(catStats, "Tables/Categorical Variables Table 1 MARBLES B123.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Missing Values
catMissing <- aggregate(categorical, by = list(categorical$Dx_alg), FUN = function(x) tryCatch(table(is.na(x))[["TRUE"]], error=function(x) NA))
catMissing$Group.1 <- as.character(catMissing$Group.1)
notEmpty <- sapply(catMissing, function(x) tryCatch(sum(x, na.rm=TRUE) > 0, error=function(x) TRUE))
catMissing <- catMissing[,notEmpty]
# Group.1 MomBirthPlaceOutsideUS DadBirthPlaceOutsideUS WantedPreg_5cat
#     ASD                      6                      7               1
#  Non-TD                     NA                     NA               1
#      TD                     NA                      1               1

contMissing <- aggregate(continuous, by = list(continuous$Dx_alg), FUN = function(x) tryCatch(table(is.na(x))[["TRUE"]], error=function(x) NA))
contMissing$Group.1 <- as.character(contMissing$Group.1)
notEmpty <- sapply(contMissing, function(x) tryCatch(sum(x, na.rm=TRUE) > 0, error=function(x) TRUE))
contMissing <- contMissing[,notEmpty]
# Group.1 MomAgeYr DadAgeYr Mat_BMI_PrePreg Mat_BMI_Deliv birthlength Creatine_Conc_mg_dl Cotinine_Extract_Conc_ng_ml Cotinine_Urine_Conc_ng_ml
#     ASD        1        2               1             1          26                  15                          18                        18
#  Non-TD       NA       NA              NA            NA          18                  22                          21                        21
#      TD       NA        1               1             1          28                  34                          42                        42


