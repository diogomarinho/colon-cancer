# model logistic regression where:
# phenotype = pc1 + pc2 + pc3 : logistic regression

# . 
sample_annotation <- read.csv('~/colon-cancer/MethylationEPIC_Sample_Sheet_2018.csv', stringsAsFactors = F)
base_dir <- '~/colon-cancer/genes/'
gene_files <- list.files(base_dir)

cases <- which(sample_annotation$Status != 'Normal')

# .
model_significance <- c()
for (i in 1:length(gene_files)){
    df <- read.csv(paste0(base_dir, gene_files[i]), head=F)
    df <- df[cases, ]
    colnames(df) <- c('TS1', 'TS2', 'TS3')
    df$Y <- as.factor(sample_annotation[cases, ]$Status)
    fit_glm <- summary(glm(Y~., data=df, family=binomial(link='logit')))
    # analgous to the ftest for the lm function: more details -> http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R7_LogisticRegression-Survival/R7_LogisticRegression-Survival3.html
    model_significance[i] <- 1 - pchisq(fit_glm$null.deviance - fit_glm$deviance, fit_glm$df.null - fit_glm$df.residual)
}
names(model_significance) <- gene_files
#example <- summary(glm(Y~., data=df, family=binomial(link='logit')))
#v <- c(example$null.deviance - example$deviance, example$df.null - example$df.residual)
