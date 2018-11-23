library(foreach)
library(minfi)
library(ggplot2)
library(sva)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(doParallel)

# . pre-processing  
print('loading basic info')
sample_sheet <- read.csv('~/colon-cancer/patients_phenotype.csv', stringsAsFactors=F)
colnames(sample_sheet) <- c('Slide', 'Array', 'Sample_Name', 'Group')
sample_sheet$Status <- as.numeric(as.factor(sample_sheet$Group))
sample_sheet$Basename <- paste0('~/colon-cancer/idat/', sample_sheet$Slide, '_',sample_sheet$Array)
sample_sheet <- sample_sheet[order(sample_sheet$Status), ]

# . pre-processing
print('reading methylation data')
RGset <- read.metharray.exp(targets = sample_sheet)
MSet <- preprocessRaw(RGset)
#ratioSet <- ratioConvert(MSet, what='both', keepCN = T)

# . plotting quality control
print('plotting basic quality control')
qc <- getQC(MSet)
png('~/colon-cancer/plots/QC.png')
print(plotQC(qc))
dev.off()

# . plotting density of methylation values coloured by group
png('~/colon-cancer/plots/methylation_density_distribution.png')
densityPlot(MSet, as.factor(sample_sheet$Group))
dev.off()

print('applying functional normalization')
# applying functional normalization 
norm.data <- preprocessFunnorm(RGset)
gset <- dropLociWithSnps(norm.data, snps=c("SBE","CpG"), maf=0)
# annotation <- getAnnotationObject(gset)

# . Extracting beta values
beta_values <- getBeta(gset)
m_values <- getM(gset)
colnames(m_values) <- sample_sheet$Sample_Name
colnames(beta_values) <- sample_sheet$Sample_Name
annotation <- getAnnotation(gset)

# Including registry so we correct for gender and age among the patients
registry <- read.csv('~/colon-cancer/registry.csv', stringsAsFactors=F)
row.names(registry) <- registry$Codigo_amostra
row.names(sample_sheet) <- sample_sheet$Sample_Name

# .  
cl <- makeCluster(10)
registerDoParallel(cl)
pheno <- sample_sheet$Status
designMatrix <- model.matrix(~pheno)
dmrs <- bumphunter(gset, design = designMatrix, cutoff = 0.2, B=1000, type="Beta")
# 
