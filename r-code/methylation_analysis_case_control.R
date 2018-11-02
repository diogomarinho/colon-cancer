library(minfi)
library(ggplot2)
library(sva)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# . 
# OBS exlude cases RET2, RET14 e RET23
# Age estimation needs to be performed for controls 
# As well gender of the patients

# formatting sample sheet 

print('loading basic info')
sample_sheet <- read.csv('~/colon-cancer/MethylationEPIC_Sample_Sheet_2018.csv', stringsAsFactors=F)
colnames(sample_sheet) <- c('Slide', 'Array', 'Sample_Name', 'Group')
sample_sheet$Status <- as.numeric(as.factor(sample_sheet$Group))
sample_sheet$Basename <- paste0('~/colon-cancer/idat/', sample_sheet$Slide, '_',sample_sheet$Array)
sample_sheet <- sample_sheet[order(sample_sheet$Status), ]
# . 
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
row.names(sample_sheet) <- sample_sheet$Sample_Name

# .  . . . 
write.csv(beta_values,  '~/colon-cancer/processed_data/beta_values.csv', row.names = T, quote = F)
write.csv(m_values,  '~/colon-cancer/processed_data/m_values.csv', row.names = T, quote = F)
write.csv(annotation, '~/colon-cancer/processed_data/annotation.csv', row.names = T, quote = F)

# no NA CpGs
m_values <- m_values[complete.cases(m_values), ]

# .
# creating t-test for control Vs pIR
pIR <- which(sample_sheet$Group == 'pIR')
control <- which(sample_sheet$Group == 'Normal')
# . 
pvals <- apply(m_values, 1, function(x) {
    p <- tryCatch({
        return(t.test(x[control], x[pIR])$p.val)
    }, error = function(e){
        return(NA)
    })
})
# . 
pvals.fdr <- p.adjust(pvals, method='fdr')
pIR_cpgs <- names(pvals.fdr[which(pvals.fdr <= 0.05)])

# creating t-test for control pCR
pCR <- which(sample_sheet$Group == 'pCR')
control <- which(sample_sheet$Group == 'Normal')
# . 
pvals <- apply(m_values, 1, function(x) {
    p <- tryCatch({
        return(t.test(x[control], x[pCR])$p.val)
    }, error = function(e){
        return(NA)
    })
})
# . 
pvals.fdr <- p.adjust(pvals, method='fdr')
pCR_cpgs <- names(pvals.fdr[which(pvals.fdr <= 0.05)])

unique_cpgs <- unique(c(pCR_cpgs, pIR_cpgs))

registry <- read.csv('~/colon-cancer/registry.csv', stringsAsFactors=F)
row.names(registry) <- registry$Codigo_amostra
row.names(sample_sheet) <- sample_sheet$Sample_Name


common_samples <- intersect(registry$Codigo_amostra, sample_sheet$Sample_Name)
m_values_patients <- m_values[unique_cpgs, common_samples]
phenotype_patients <- sample_sheet[common_samples, ]
registry <- registry[common_samples, ]


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
cpgs_df <- data.frame(chr=annotation$chr, pos=annotation$pos, strand=annotation$strand, name=annotation$Name, P=pvalues)
write.csv(cpgs_df, '~/colon-cancer/processed_data/lm_pvalues_control_filtered.csv', row.names=T, quote=F)
# using dmp finder

