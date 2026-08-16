library('DESeq2')
library('RColorBrewer')
library('ggplot2')
rawCountTable <- read.table('READemption_analysis/output/salmonella_gene_quanti_combined/gene_wise_quantifications_combined.csv', skip=1, sep='\t', quote='', comment.char='', colClasses=c(rep('character',10), rep('numeric',4)))
countTable <- round(rawCountTable[,11:length(names(rawCountTable))])
colnames(countTable) <- c('InSPI2_R1','InSPI2_R2','LSP_R1','LSP_R2')
# Select only the libraries of this species
countTable <- countTable[, c('InSPI2_R1','InSPI2_R2','LSP_R1','LSP_R2')]
libs <- c('InSPI2_R1','InSPI2_R2','LSP_R1','LSP_R2')
conds <- c('InSPI2', 'InSPI2', 'LSP', 'LSP')
reps <- c('1', '2', '1', '2')
samples <- data.frame(row.names=libs, condition=conds, lib=libs, replicate=reps)
dds <- DESeqDataSetFromMatrix(countData=countTable, colData=samples, design=~condition)
dds <- DESeq(dds, betaPrior=TRUE)

# PCA plot
pdf('READemption_analysis/output/salmonella_deseq/deseq_raw/sample_comparison_pca_heatmap.pdf')
rld <- rlog(dds)
pcaData <- plotPCA(rld, 'condition', intgroup=c('condition', 'replicate'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, 'percentVar'))
print(ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
geom_point(size=3) +
xlab(paste0('PC1: ',percentVar[1],'% variance')) +
ylab(paste0('PC2: ',percentVar[2],'% variance')) +
coord_fixed())
# Heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- with(colData(dds), paste(lib, sep=' : '))
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap(mat, col = rev(hmcol), margins=c(13, 13))
dev.off()
comp0 <- results(dds, contrast=c('condition','InSPI2', 'LSP'))
write.table(comp0, file='READemption_analysis/output/salmonella_deseq/deseq_raw/deseq_comp_InSPI2_vs_LSP.csv', quote=FALSE, sep='\t')
comp1 <- results(dds, contrast=c('condition','LSP', 'InSPI2'))
write.table(comp1, file='READemption_analysis/output/salmonella_deseq/deseq_raw/deseq_comp_LSP_vs_InSPI2.csv', quote=FALSE, sep='\t')
