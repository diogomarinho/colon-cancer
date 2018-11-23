library(minfi)
library(ggplot2)
library(sva)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# .  
print('loading basic info')
sample_sheet <- read.csv('~/colon-cancer/patients_phenotype.csv', stringsAsFactors=F)
colnames(sample_sheet) <- c('Slide', 'Array', 'Sample_Name', 'Group')
sample_sheet$Status <- as.numeric(as.factor(sample_sheet$Group))
sample_sheet$Basename <- paste0('~/colon-cancer/idat/', sample_sheet$Slide, '_',sample_sheet$Array)
sample_sheet <- sample_sheet[order(sample_sheet$Status), ]

# pre-processing
print('reading methylation data')
RGset <- read.metharray.exp(targets = sample_sheet)
MSet <- preprocessRaw(RGset)
#ratioSet <- ratioConvert(MSet, what='both', keepCN = T)

# plotting quality control
print('plotting basic quality control')
qc <- getQC(MSet)
png('~/colon-cancer/plots/QC.png')
print(plotQC(qc))
dev.off()

#plotting density of methylation values coloured by group
png('~/colon-cancer/plots/methylation_density_distribution.png')
densityPlot(MSet, as.factor(sample_sheet$Group))
dev.off()

print('applying functional normalization')
# applying functional normalization 
norm.data <- preprocessFunnorm(RGset)
gset <- dropLociWithSnps(norm.data, snps=c("SBE","CpG"), maf=0)
# annotation <- getAnnotationObject(gset)

# Extracting beta values
beta_values <- getBeta(gset)
m_values <- getM(gset)
colnames(m_values) <- sample_sheet$Sample_Name
colnames(beta_values) <- sample_sheet$Sample_Name
annotation <- getAnnotation(gset)

# including registry so we correct for gender and age among the patients
registry <- read.csv('~/colon-cancer/registry.csv', stringsAsFactors=F)
row.names(registry) <- registry$Codigo_amostra
row.names(sample_sheet) <- sample_sheet$Sample_Name

#  . . . 
write.csv(beta_values,  '~/colon-cancer/processed_data/beta_values.csv', row.names = T, quote = F)
write.csv(m_values,  '~/colon-cancer/processed_data/m_values.csv', row.names = T, quote = F)
write.csv(annotation, '~/colon-cancer/processed_data/annotation.csv', row.names = T, quote = F)

# starting analysis treatment vs non-treatment:
common_samples <- intersect(registry$Codigo_amostra, sample_sheet$Sample_Name)
m_values_patients <- m_values[,common_samples] 
phenotype_patients <- sample_sheet[common_samples, ]
registry <- registry[common_samples, ]

# no NA CpGs 
m_values_patients <- m_values_patients[complete.cases(m_values_patients),]

# using dmp finder 
dmp <- dmpFinder(m_values_patients, pheno = as.factor(registry$RESPOSTA), type='categorical')
cpgs_df <- data.frame(chr=annotation$chr, pos=annotation$pos, strand=annotation$strand, name=annotation$Name, P=dmp$pval)
write.csv(cpgs_df, '~/colon-cancer/processed_data/dmp_pvalues.csv', quote=F, row.names=T)

# linear model corrected by age and sex 
print('compute logistic regression corrected by age and sex')
pvalues <- apply(m_values_patients, 1, function(x) {
    pval <- tryCatch({
        data_set <- data.frame(cpg=x,treatment=as.factor(registry$RESPOSTA), age=registry$Idade_ao_diagnostico, sex=as.factor(registry$Sexo))
         mod1 <- lm(cpg~., data_set)
         mod1_summary <- summary(mod1)
         return(mod1_summary$coefficients['treatmentpIR',4])
    }, error = function(e){
        return (NA)
    })
   return(pval)
})
# . 
# . 
cpgs_df <- data.frame(chr=annotation$chr, pos=annotation$pos, strand=annotation$strand, name=annotation$Name, P=pvalues)
# . 
write.csv(cpgs_df, '~/colon-cancer/processed_data/lm_pvalues.csv', quote=F, row.names=T)

#
#cpgs_df <- data.frame(chr=annotation$chr, pos=annotation$pos, strand=annotation$strand, name=annotation$name, p=pvalues)
#cpgs_df$bonferroni <- p.adjust(pvalues, method='bonferroni')
#cpgs_df$fdr <- p.adjust(pvalues, method='fdr')
#write.csv(cpgs_df, '~/colon-cancer/processed_data/pvalues.csv', quote=f, row.names=t)

# . 
# phenotype file 
#status <- as.factor(sample_sheet$Group)
#gender <- as.factor(getSex(norm.data)$predictedSex)
#Array <- as.factor(sample_sheet$Array)
#pheno <- data.frame(status, Array)
#print('CpG wise differential methylation finder')
#cpgs <- sample(1:nrow(beta_values), 1000)

# comment when the protyped scrip is over 
#beta_values <- beta_values[cpgs, ]
#annotation <- annotation[cpgs, ]
# . Correcting for batch effect 

#mod <- model.matrix(~as.factor(status), data=pheno)
#mod0 <- model.matrix(~1, data=pheno)
#sva.results <- sva(beta_values, mod, mod0)
#sv <- data.frame(sva.results$sv)
# . 
##for (i in 1:ncol(sv)) {
#    colnames(sv)[i] <- paste0('sv', i)
#}
# . 
#pheno <- data.frame(status, gender, sv)
## create a ANOVA analysis to understand which treatments a significant different
#pvalues <- apply(beta_values, 1, function(cpg) {
#    data_set <- data.frame(cpg, pheno)
#    reg <- lm(cpg~status+gender+sv1+sv2+sv3+sv4+sv5+sv6+sv7, data=data_set)
#    reg.anova <- anova(reg)
#    pval <- reg.anova[[5]][1]
#    
#    return(pval)
#})

# cpg finder
#dmp <- dmpFinder(beta_values, pheno = pheno$status, type='categorical')
#dmp$bonferroni <- p.adjust(dmp$pval, method='bonferroni')

# CpG test 
#colnames(beta_values) <- sample_sheet$Sample_Name
#cpg1 <- beta_values[1, ]
#data_test <- data.frame(cpg1, status, gender)

#mod1 <- lm(cpg1 ~ status + gender, data=data_test)
#mod1.anova <- anova(mod1)
#s.mod1 <- summary(mod1)
#pval <- mod1.anova[[5]][1]

# beta_values <- getBeta(RGset)
# applying functional normalization in in the dataset
#beta_values <- getBeta(norm.data)
#beta_values <- beta_values[complete.cases(beta_values),]
#write.csv('~/colon-cancer/processed_data/beta_values.csv', quote=F, row.names=T)
# to make the frame work faster, comment this line afterwards ...
#beta_values <- beta_values[sample(1:nrow(beta_values), 10000),]

# checking for batch effects
#t_beta_values <- t(beta_values)
#row.names(t_beta_values) <-  colnames(beta_values)
#colnames(t_beta_values) <- row.names(beta_values)
#pca <- prcomp(t_beta_values, center=T, scale=T)
# 2 principal components
#pcs <- data.frame(predict(pca, t_beta_values))
#pcs$Array <- sample_sheet$Array
# visualing PCS
#png('~/colon-cancer/plots/batch_1.png')
#print(ggplot(pcs, aes(x=PC1, y=PC2, col=Array)) + geom_point())
#dev.off()
#
#png('~/colon-cancer/plots/batch_2.png')
#print(ggplot(pcs, aes(x=PC1, y=PC3, col=Array)) + geom_point())
#dev.off()
#
#png('~/colon-cancer/plots/batch_3.png')
#print(ggplot(pcs, aes(x=PC1, y=PC4, col=Array)) + geom_point())
#dev.off()
# . 
# . Removing batch effect
# . applying batch effect correction 
#pheno <- sample_sheet[, c('Sample_Name', 'Status', 'Array')]
#batch <- pheno$Array
#mod <- model.matrix(~Status, data=pheno)
#combat = ComBat(dat=beta_values, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

# once the data is batch-free we can perform a linear regression on the model:
# . Methylation ~ Status


