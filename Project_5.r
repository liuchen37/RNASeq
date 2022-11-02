
# 4
# 4.1 Preparations

library(DESeq2)
library(tidyverse)
library(ggfortify)
library(gplots)
library(RColorBrewer)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(fgsea)
library(labeling)
library(ggplot2)
library(affy)
library(scales)
library(Hmisc)

file.choose()
setwd("/Users/chenliu/Documents/Tmp workspace/P5_DAI")

counts_raw <- read.table("readouts.tsv",
                         header = T,
                         sep = "\t",
                         fill = T,
                         check.names = F)
#remove decimal
counts_raw$Name <- substr(counts_raw$Name, 1, 18)
rownames(counts_raw) <- counts_raw$Name
counts_raw <- counts_raw[,-c(1)]

#remove non-expressing genes
countdata <- counts_raw[rowSums(counts_raw) > 10,]

#make a sample table
samples <- data.frame(
  sampleID = c(1:9),
  sample = c("WT", "KO", "WT", "WT", "WT", "WT", "KO", "KO", "KO"))
samples$sample <- factor(samples$sample, levels = c("WT", "KO"))
rownames(samples) = samples$sampleID

# 4.2 QC plots
#define colour
colour <- colorRampPalette(c("blue4", "white", "red"))(256)

#log2
logcounts <- log2(countdata + 1)

# 4.2.1 Boxplot reflecting basic stats cross samples
boxplot(logcounts,
        ylab = "log2(expression)", xlab = "Sample ID",
        col = "green3")

# 4.2.2 Density plot for read counts
plotDensity(logcounts, 
            col =1:9,
            lty = c(1:ncol(countdata)), 
            xlab = "Log2(expression)",
            ylab = "Density",
            main = "Expression Distribution")
#add legend
par(cex = 1)
legend("right",
       names(countdata),
       lty = c(1:ncol(countdata)),
       col = 1:9)
abline(v = -1.5, lwd = 1, col = "red", lty = 2)

# 4.2.3 Heatmap (Pearson correlation matrix)
cormat <- round(cor(logcounts, method = c("pearson")),2)
#calculate the p value for each correlation coefficient
res <- rcorr(as.matrix(cormat))
#check the value
res[["P"]]

ha_column = HeatmapAnnotation(df = data.frame(status = c("WT", "KO", "WT", "WT", "WT", "WT", "KO", "KO", "KO")),
                              col = list(status = c("WT" =  "green3",
                                                    "KO" = "red3")))
Heatmap(cormat, name = "PCC",
        column_title = "Pearson Correlation Matrix",
        top_annotation = ha_column)

# 4.2.4 PCA plot
#convert to matrix, then dds
countdata <- as.matrix(countdata)
ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = samples,
                              design = ~sample)
#DESeq2
ddsObj <- DESeq(ddsObj.raw)
results <- results(ddsObj, alpha = 0.5)
#variance stabalisation, to allow further normalisation and plot
vsd <- vst(ddsObj, blind = F)
plotPCA(vsd, intgroup = c("sample"))
#return coordinate and put coordinate to a table
plotPCA(vsd, intgroup = c("sample"), returnData = T)
PCAreturn <- "
coordinates here
"
PCAreturn <- read.table(header = T,
                        text = PCAreturn)
#convert to factor
PCAreturn$sample = factor(PCAreturn$sample,
                          levels = c("WT", "KO"))
#plot with ggplot2
ggplot(PCAreturn) +
  geom_point(aes(PC1, PC2, color = sample), size = 3) +
  labs(x = "PC1: 71% variance", y = "PC2: 11% variance", face = "bold") +
  theme_classic(base_size = 16)

# 4.3 Analysis

# 4.3.1 DESeq and normalisation
dds <- DESeq(ddsObj.raw)
sizeFactors(dds)
#normalisation with DESeq2 and export
normalised_counts <- counts(dds, normalized=TRUE)
write.table(normalised_counts, file="normalised_counts.txt",
            sep="\t",
            quote=F,
            col.names=NA)

# 5 check up/down-regulation

#5.1 in R
normalised_counts <- as.data.frame(normalised_counts)
normalised_counts$logFC_KO_vs_WT <-
  log2((rowSums(normalised_counts[c(2, 7, 8, 9)])) / 
         (rowSums(normalised_counts[c(1, 3, 5, 6)])))

#extract WT vs KO from dds and export as data frame
resultsNames(dds)
KO_vs_WT <- results(dds, name = "sample_KO_vs_WT")
KO_vs_WT <- as.data.frame(KO_vs_WT)

#check the existence of N/A, if yes, replace by 1
KO_vs_WT$padj[is.na(KO_vs_WT$padj)] = 1

#combine
normalised_counts <- merge(normalised_counts, KO_vs_WT, by = "row.names")
rownames(normalised_counts) <- normalised_counts$Row.names
normalised_counts <- normalised_counts[,-1]

write.table(normalised_counts, file="dge.txt",
            sep="\t",
            quote=F,
            col.names=NA)

# 5.2 P value histogram
hist(normalised_counts$pvalue, 
     col = "green3",
     xlab = "p value",
     main = "P value distribution")

#what are these genes?
info <- read.csv("gene_info.csv",
                         header = T,
                         fill = T,
                         check.names = F)
info$EnsemblID <- substr(info$EnsemblID, 1, 18)
rownames(info) <- info$EnsemblID
info <- info[,-c(1)]
dge_results <- merge(info, KO_vs_WT, by = "row.names")
#add a significant filter
dge_results_sig <- subset(dge_results, abs(dge_results$log2FoldChange) >1 &
                            dge_results$padj <0.05)
write.table(dge_results_sig, file="dge_results_sig.txt",
            sep="\t",
            quote=F,
            col.names=NA)
