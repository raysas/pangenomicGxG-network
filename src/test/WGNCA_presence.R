library(WGCNA)

Rtab_path="data/pangenomes/Campylobacter_coli/roary_pangenome/gene_presence_absence.Rtab"
df <- read.table(Rtab_path, header=T, row.names=1)


df <- t(df)
dim(df)

csv_file="data/presence_matrices/Campylobacter_coli_presence_absence.csv"
df <-read.csv(csv_file, header=T, row.names=1)

#removing all cols that are all 1s (no variation <=> st. dev.=0 <=> cor=NaN <=> error :/)
stdevs <- apply(df, 2, sd)
df <- df[,stdevs!=0]

#correlation matrix
cor_matrix <- cor(df)

#soft threshold
powers = c(1:10)
sft = pickSoftThreshold(df, powerVector = powers, verbose = 5)
# TOM = TOMsimilarityFromExpr(df, power = sft$power)

net = blockwiseModules(df, power = sft$power, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "TOM", verbose = 3)
