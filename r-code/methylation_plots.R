library(qqman)
library(ggplot2)
#  
# . pvalues from anova analysis using the mode methylation = treatment1 + treatment2 + gender + surrogatevariables1..7
# since age is not available for the control samples it is not used as covariate
# . loading data 
pvalues <- read.csv('~/colon-cancer/processed_data/dmp_pvalues.csv', row.names=1)
m_values <- read.csv('~/colon-cancer/processed_data/m_values.csv', row.names=1)
phenotype <- read.csv('~/colon-cancer/MethylationEPIC_Sample_Sheet_21.02.17.csv', stringsAsFactors = F)
registry <- read.csv('~/colon-cancer/registry.csv', stringsAsFactors = F)
row.names(phenotype) <- phenotype$Sample_Name
row.names(registry) <- registry$Codigo_amostra

#  gettin only patients for now
common_samples <- intersect(registry$Codigo_amostra, phenotype$Sample_Name)
m_values_patients <- m_values[, common_samples]
phenotype_patients <- phenotype[common_samples, ]
registry_patients <- registry[common_samples, ]

# . 
colnames(m_values_patients) <- phenotype_patients$Status
pCR <- which(colnames(m_values_patients) == 'pCR')
pIR <- which(colnames(m_values_patients) == 'pIR')
# . 
mean_difference <- apply(m_values_patients, 1, function(x){ 
    mean_difference <- mean(x[pCR]) - mean(x[pIR])
    return(mean_difference)
})

# bonferroni bounderies for plotting with the original p-values
significant <- 0.05/nrow(pvalues)
suggestive <- 0.1/nrow(pvalues)

# volcano plot of the methylation difference from treatmen vs non-treatment vs anova p-value
volcano_plot_data <- data.frame(fold_change=mean_difference, P=pvalues$P)
row.names(volcano_plot_data) <- row.names(m_values_patients)

#
plot <- ggplot(volcano_plot_data, aes(fold_change, -log(P))) + geom_point() + 
    geom_hline(yintercept = -log(significant), linetype='dashed', col='red', size = 0.2) + 
    geom_hline(yintercept = -log(suggestive), linetype='dashed', col='blue', size = 0.2)  

png('~/colon-cancer/plots/v_plot.png')
print(plot)
dev.off()

png('~/colon-cancer/plots/qq.png')
print(qq(volcano_plot_data$P))
dev.off()
#
