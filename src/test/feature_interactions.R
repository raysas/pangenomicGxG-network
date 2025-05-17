install.packages(c('iml','rpart','MASS',"e1071"))

library(iml)
library(rpart)
library(MASS)
library(e1071)

PHENO_FILE_PATH <- "data/phenotypes/Campylobacter_coli_ciprofloxacin.csv"
PRESENC_FILE_PATH <- "data/pangenomes/Campylobacter_coli/roary_pangenome/gene_presence_absence.Rtab"

get_presence_df <-function(Rtab_path, pheno_path){
    df<- read.table(Rtab_path, header=TRUE, sep="\t")
    rownames(df) <- df$Gene
    df <- df[,-1]
    df <- t(df)

    labels <- read.csv(pheno_path)
    rownames(labels) <- labels$genome_id

    new_df <- cbind(labels, df)
    new_df <- new_df[,-1]
    new_df <- as.data.frame(new_df)

    new_df$SIR <- as.factor(new_df$SIR)

    return(new_df)
}

df <- get_presence_df(PRESENC_FILE_PATH, PHENO_FILE_PATH)

#################################
# IE


rf<- rpart(SIR ~ ., data=df, method="class")

model_obj <- Predictor$new(rf, data=df[-1,])

interaction_strength <- Interaction$new(model_obj)
plot(interaction_strength)

head(interaction_strength$results)


######################

#train a SVM with 5 fold cv
svm <- svm(SIR ~ ., data=df, kernel="radial", cost=10, gamma=0.1)

