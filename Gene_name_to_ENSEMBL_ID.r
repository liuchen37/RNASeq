library(tidyverse)
library(biomaRt)

df <- read.csv("diffexpr.csv", check.names = TRUE)
gene_symbols <- df$Gene

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Map gene symbols to Ensembl IDs
mapping <- getBM(
  attributes = c('mgi_symbol', 'ensembl_gene_id'),
  filters = 'mgi_symbol',
  values = gene_symbols,
  mart = ensembl
)

# Merge Ensembl IDs with your dataframe
df_merged <- left_join(df, mapping, by = c("Gene" = "mgi_symbol"))

write.csv(df_merged, "file.csv", row.names = FALSE)
