## using edgeR to estimate differential gene expression (DGE)
# using edgeR from a count table

rm(list =ls())

library(dplyr)
library(readr)
library(magrittr)
library(edgeR)

# load count table and phenotype file
file.choose()
count_dataframe <- readr::read_tsv("file")
glimpse(count_dataframe)
pheno_data <- readr::read_table("file")


# extract gene from the dataframe and replace with "NULL"
genes <- count_dataframe[['gene']]
count_dataframe[['gene']] <- NULL

# convert into matrix, genes as row names
count_matrix <- as.matrix(count_dataframe)
rownames(count_matrix) <- genes

# specify experiment of interest
experiments_of_interest <- c("L1Larvae", "L2Larvae")
columns_of_interest <- which( pheno_data[['stage']] %in% experiments_of_interest )

# form grouping factor
grouping <- pheno_data[['stage']][columns_of_interest] %>% 
  forcats::as_factor()

# form subset count data
counts_of_interest <-  count_matrix[,columns_of_interest]

# build DGE object
count_dge <- edgeR::DGEList(counts = counts_of_interest, group = grouping)

# perform DGE
design <- model.matrix(~ grouping)
eset_dge <- edgeR::estimateDisp(count_dge, design)

# fit into a generalized linear model and quasi-likelihood F-test
fit <- edgeR::glmQLFit(eset_dge, design)
result <- edgeR::glmQLFTest(fit, coef=2)

# present result
topTags(result)
write.table(result, file = "file",
            sep = '\t', quote = F)
