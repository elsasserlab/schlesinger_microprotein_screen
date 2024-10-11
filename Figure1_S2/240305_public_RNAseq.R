library(ggplot2)
library(tidyverse)
library(EnsDb.Hsapiens.v86)

setwd("~/work/experiments/240530_Public_RNAseq/")
gene_df <- read.table("./gene_summaries_final_v3.csv", sep = ",", header =TRUE, comment.char = "")

# stupid thing, gene_summaries HAS the gene_id column I just didn't see it
# pep_to_gene <- read.table("./A375_TRAC_expression_per_peptide_no_ties_only_overlapping_peptides.tab", sep = "\t", header = TRUE)

# GSM5330783: RNAseq_df
# GSM2343135 RNAseq_K562
expr_df <- read.table("GSM5330783_ENCFF435PHM_gene_quantifications_GRCh38.tsv", sep = "\t", header = TRUE)
expr_df <- expr_df %>% mutate(gene_id_nov = gsub("\\.[0-9]+$", "", gene_id))
ens2symbol <- ensembldb::select(EnsDb.Hsapiens.v86, keys= expr_df$gene_id_nov , keytype = "GENEID", columns = c("SYMBOL","GENEID"))
expr_df <- expr_df %>% left_join(ens2symbol, by = c("gene_id_nov" = "GENEID"))

expr_df_tpm <- expr_df %>% 
  dplyr::select("gene_id", "gene_id_nov", "SYMBOL", "TPM")

# If more than one row corresponds to the same gene, we just exclude these
# (majority of values are zero, or very conflicting (see MALAT1 extremely high in one, zero in another))
# About 220 gene names that appear more than once, covering almost 2000 rows
# of the table (out of 59k)
gene_symbol_counts <- table(expr_df_tpm$SYMBOL)
tpms_conflict <- expr_df_tpm %>% dplyr::filter(SYMBOL %in% names(gene_symbol_counts[gene_symbol_counts > 1])) %>% arrange(SYMBOL)

gene_df_plus_rnaseq <- gene_df %>% 
  dplyr::left_join(
    expr_df_tpm %>% 
      dplyr::select("gene_id", "SYMBOL", "TPM") %>%   
      dplyr::filter(!SYMBOL %in% unique(tpms_conflict$SYMBOL) & ! is.na(SYMBOL)) %>% 
      dplyr::rename_with( ~ paste0("HCT116_GSM5330783_", .x) ),
    by = c("gene_id" = "HCT116_GSM5330783_SYMBOL")
  )

# GSM2343135 RNAseq_K562
expr_df <- read.table("GSM2343135_ENCFF477XRV_gene_quantifications_GRCh38.tsv", sep = "\t", header = TRUE)
expr_df <- expr_df %>% mutate(gene_id_nov = gsub("\\.[0-9]+$", "", gene_id))
ens2symbol <- ensembldb::select(EnsDb.Hsapiens.v86, keys= expr_df$gene_id_nov , keytype = "GENEID", columns = c("SYMBOL","GENEID"))
expr_df <- expr_df %>% left_join(ens2symbol, by = c("gene_id_nov" = "GENEID"))

expr_df_tpm <- expr_df %>% 
  dplyr::select("gene_id", "gene_id_nov", "SYMBOL", "TPM")

# If more than one row corresponds to the same gene, we just exclude these
gene_symbol_counts <- table(expr_df_tpm$SYMBOL)
tpms_conflict <- expr_df_tpm %>% dplyr::filter(SYMBOL %in% names(gene_symbol_counts[gene_symbol_counts > 1])) %>% arrange(SYMBOL)

gene_df_plus_rnaseq_k562 <- gene_df_plus_rnaseq %>% 
  dplyr::left_join(
    expr_df_tpm %>% 
      dplyr::select("gene_id", "SYMBOL", "TPM") %>% 
      dplyr::filter(!SYMBOL %in% unique(tpms_conflict$SYMBOL) & ! is.na(SYMBOL)) %>% 
      dplyr::rename_with( ~ paste0("K562_GSM5330783_", .x) ),
    by = c("gene_id" = "K562_GSM5330783_SYMBOL")
  )

write.table(gene_df_plus_rnaseq_k562, file = "./gene_summaries_final_v4.csv", col.names = TRUE, row.names = FALSE, sep = ",")
