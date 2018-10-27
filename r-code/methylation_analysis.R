library(minfi)
library(ggplot2)
library(sva)

#OBS exlude cases RET2, RET14 e RET23
# Age estimation needs to be performed for controls 
# As well gender of the patients

# formatting sample sheet 
sample_sheet <- read.csv('~/cancer/MethylationEPIC_Sample_Sheet_21.02.17.csv', stringsAsFactors=F)
colnames(sample_sheet) <- c('Slide', 'Array', 'Sample_Name', 'Group')
sample_sheet$Status <- as.numeric(as.factor(sample_sheet$Group))
sample_sheet$Basename <- paste0(getwd(),'/idat/', sample_sheet$Slide, '_',sample_sheet$Array)

# . 
registry <- read.csv('~/cancer/registry.csv')

# 
RGset <- read.metharray.exp(targets = sample_sheet)
beta_values <- getBeta(RGset)
beta_values <- beta_values[complete.cases(beta_values),]


# checking for batch effects
t_beta_values <- t(beta_values)
row.names(t_beta_values) <-  colnames(beta_values)
colnames(t_beta_values) <- row.names(beta_values)
pca <- prcomp(t_beta_values, center=T, scale=T)
# 2 principal components
pcs <- data.frame(predict(pca, t_beta_values))
pcs$Array <- sample_sheet$Array
# visualing PCS
png('~/cancer/plots/batch_1.png')
print(ggplot(pcs, aes(x=PC1, y=PC2, col=Array)) + geom_point())
dev.off()

png('~/cancer/plots/batch_2.png')
print(ggplot(pcs, aes(x=PC1, y=PC3, col=Array)) + geom_point())
dev.off()

png('~/cancer/plots/batch_3.png')
print(ggplot(pcs, aes(x=PC1, y=PC4, col=Array)) + geom_point())
dev.off()

# .
combat_edata = ComBat(dat=beta_values, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)









