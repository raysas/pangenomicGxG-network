install.packages(c("BiocManager", "devtools", "tidyverse"))

#corelated genes: have similar gene expression patterns, they are functionally associated

BiocManager::install("WGCNA")
BiocManager::install("DESeq2")
BiocManager::install("GEOquery")
BiocManager::install("gridExtra")

install.packages("remotes")
remotes::install_github("kevinblighe/CorLevelPlot")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(CorLevelPlot)
library(gridExtra)
library(tidyverse)

#downlaod data from a link
accesion_link<-"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152418&format=file&file=GSE152418%5Fp20047%5FStudy1%5FRawCounts%2Etxt%2Egz"
download.file(accesion_link, "data/test/GSE152418_p20047_Study1_RawCounts.txt.gz")
system("gunzip data/test/GSE152418_p20047_Study1_RawCounts.txt.gz")


gse <- getGEO("GSE152418")
pheno_data <- pData(phenoData(gse[[1]]))
pheno_data <- pheno_data[,c(1,2,46:50)]

df <- read.delim("data/test/GSE152418_p20047_Study1_RawCounts.txt", header=T)

df %>%
    gather(key='samples', value='counts', -ENEMBLID) %>%
    mutate(samples=gsub('\\.', '-', samples)) %>%
    inner_join(., pheno_data, by=c('samples'='titles')) %>%
    select(1,3,4) %>%
    spread(key="geo_accession", value='counts') %>%
    column_to_rownames(var="ENEMBLID") -> df


View(df)


