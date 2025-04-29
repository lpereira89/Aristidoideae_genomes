library(readr)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("~/Desktop/Sheff/Aristidoideae-genomes/correct-losses-denovo/gene-expression/")

################ Maize dataset from Li et al., 2010) ########################
# Load gene expression data and convert between different gene ID versions
Converter_GeneID <- read_table("B73v3_to_B73v5.tsv", col_names = FALSE)
colnames(Converter_GeneID) <- c("Old_geneID", "New_geneID")

expression_BS <- read_excel("41588_2010_BFng703_MOESM5_ESM.xls", sheet = "BS ")
expression_M <- read_excel("41588_2010_BFng703_MOESM5_ESM.xls", sheet = "ME")

BS <- merge(Converter_GeneID, expression_BS, by.x = "Old_geneID", by.y = "geneID", all.y = TRUE)
BS <- BS %>%
  separate(New_geneID, into = c("New_geneID1", "New_geneID2", "New_geneID3", "New_geneID4",
                                "New_geneID5", "New_geneID6", "New_geneID7", "New_geneID8",
                                "New_geneID9", "New_geneID10"),
           sep = ",")
M <- merge(Converter_GeneID, expression_M, by.x = "Old_geneID", by.y = "geneID", all.y = TRUE)
M <- M %>%
  separate(New_geneID, into = c("New_geneID1", "New_geneID2", "New_geneID3", "New_geneID4",
                                "New_geneID5", "New_geneID6", "New_geneID7"),
           sep = ",")

# Load gene candidate lists and process them to get a table with OG, gene ID and type of candidate
maize_dups <- read.delim("maize_dups.txt", header=FALSE)
colnames(maize_dups) <- c("OG", "Gene")
maize_dups <- maize_dups %>% mutate(type='DUP') %>% filter(Gene != '')

maize_de_novo <- read.delim("maize_de_novo.txt", header=FALSE)
colnames(maize_de_novo) <- c("OG", "Gene")
maize_de_novo <- maize_de_novo %>% mutate(type='DN')

maize_losses <- read.delim("maize_gene_loss.txt", header=FALSE)
colnames(maize_losses) <- c("OG", "Gene")
maize_losses <- maize_losses %>% mutate(type='LOSS')

maize_candidates <- bind_rows(maize_dups, maize_de_novo, maize_losses)
maize_candidates <- maize_candidates %>% filter(Gene != '')

# There are a few genes for which there was one gramene ID (old ID) but the gene model was split
# in new annotation. In those cases, since the expression data comes from the old ID, we only use
# one new gene model.
expression_BS <- merge(BS, maize_candidates, by.x = "New_geneID1", by.y = "Gene")
expression_BS <- expression_BS %>% select(-"New_geneID4", -"New_geneID5", -"New_geneID6", -"New_geneID7",
                                          -"New_geneID8", -"New_geneID9", -"New_geneID10")

expression_M <- merge(M, maize_candidates, by.x = "New_geneID1", by.y = "Gene")
expression_M <- expression_M %>% select(-"New_geneID4", -"New_geneID5", -"New_geneID6", -"New_geneID7")

expression_maize <- bind_rows(expression_BS, expression_M)
expression_maize <- expression_maize %>% mutate(foldChange = `RPM-BS`/`RPM-ME`)

write_delim(expression_maize, "results_expression_maize.txt")

plot1 <- ggplot(expression_maize) + geom_point(aes(x = log(`RPM-ME`), y = log(`RPM-BS`),
                                             color = foldChange, size = foldChange, shape = type)) +
  scale_color_gradient(low = "#ffb09c", high = "#900000", name = "Fold Change") +
  theme_classic() + xlab('Log(Mesophyll_RPM)') + ylab('Log(BundleSheath_RPM)')
pdf('Expression_plot_maize.pdf', width = 6, height = 6)
print(plot1)
dev.off()

# Focus on the DUP:
expression_maize_onlydup <- expression_maize %>% filter(type == 'DUP')
OG_withmaizegene <- levels(as.factor(maize_dups$OG))
OG_maize_DEGs <- levels(as.factor(expression_maize_onlydup$OG))
ME_maize_genes <- expression_maize_onlydup %>% filter(foldChange < 1)
BS_maize_genes <- expression_maize_onlydup %>% filter(foldChange > 1)

##################Sorghum dataset from Doring et al., 2016, only Illumina data ########################

# Load gene expression dataset
sorghum_expression <- read_excel("Supplementary_table_S1.xlsx",
                                    col_types = c("text", "text", "text",
                                                  "text", "text", "text", "text", "text",
                                                  "text", "text", "text", "numeric",
                                                  "numeric", "numeric", "numeric",
                                                  "numeric", "numeric", "numeric",
                                                  "numeric", "text", "numeric", "numeric",
                                                  "text", "text"))
sorghum_expr_DEG <- sorghum_expression %>% filter(p.valueBonferronicorrected < 0.05)

# Load sorghum gene candidate lists
sorghum_dups <- read.delim("sorghum_dups.txt", header=FALSE)
colnames(sorghum_dups) <- c("OG", "Gene")
sorghum_dups <- sorghum_dups %>% mutate(type='DUP')
sorghum_dups <- sorghum_dups %>% filter(Gene != '')

sorghum_de_novo <- read.delim("sorghum_de_novo.txt", header=FALSE)
colnames(sorghum_de_novo) <- c("OG", "Gene")
sorghum_de_novo <- sorghum_de_novo %>% mutate(type='DN')

sorghum_losses <- read.delim("sorghum_gene_loss.txt", header=FALSE)
colnames(sorghum_losses) <- c("OG", "Gene")
sorghum_losses <- sorghum_losses %>% mutate(type='LOSS')

sorghum_candidates <- bind_rows(sorghum_dups, sorghum_de_novo, sorghum_losses)
sorghum_candidates <- sorghum_candidates %>% filter(Gene != '')

expression_sorghum <- merge(sorghum_expr_DEG, sorghum_candidates, by.x = "locusName", by.y = "Gene")
expression_sorghum <- expression_sorghum %>% mutate(foldChange = B_rpkm/M_rpkm)

write_delim(expression_sorghum, "results_expression_sorghum.txt")

expression_sorghum <- expression_sorghum %>% mutate(foldChange = B_rpkm/M_rpkm)

plot2 <- ggplot(expression_sorghum) + geom_point(aes(x = log(M_rpkm), y = log(B_rpkm),
                                             color = foldChange, size = foldChange, shape = type)) +
  scale_color_gradient(low = "#ffb09c", high = "#900000", name = "Fold Change") +
  theme_classic() + xlab('Log(Mesophyll_RKPM)') + ylab('Log(BundleSheath_RKPM)')
pdf('Expression_plot_sorghum.pdf', width = 6, height = 6)
print(plot2)
dev.off()

# Focus on the DUP:
expression_sorghum_onlydup <- expression_sorghum %>% filter(type == 'DUP')
OG_withsorghumgene <- levels(as.factor(sorghum_dups$OG))
OG_sorghum_DEGs <- levels(as.factor(expression_sorghum_onlydup$OG))
ME_sorghum_genes <- expression_sorghum_onlydup %>% filter(foldChange < 1)
BS_sorghum_genes <- expression_sorghum_onlydup %>% filter(foldChange > 1)

## Adding both together

pdf("Gene_expr_both.pdf", height = 6, width = 10)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

diffexpr_both <- merge(expression_maize, expression_sorghum, by="OG")
OG_DEGs_both <- levels(as.factor(diffexpr_both$OG))
write.table(diffexpr_both, 'Common_maize_sorghum_DEG.txt', row.names = FALSE, quote = FALSE, sep = '\t')

#################### Rice dataset from Hua et al., 2021 ##############################
# Load rice expression dataset after deleting vein data (not relevant here)
rice_expression <- read_excel("tpj15292-sup-0004-tables3.xlsx")
# Keep only genes that are expressed in at least one replicate
rice_expression <- rice_expression %>% filter(!(TPM.M1 == 0 & TPM.M2 == 0 & TPM.M3 == 0 & TPM.M4 == 0 & TPM.M5 == 0 &
                                                  TPM.BS1 == 0 & TPM.BS2 == 0& TPM.BS3 == 0 & TPM.BS4 == 0& TPM.BS5 == 0))
rice_expr_DEG <- rice_expression %>% filter(`DESeq2 adjP BSvsM` < 0.05)

# Load candidate genes
rice_dups <- read.delim("rice_dups.txt", header = FALSE)
colnames(rice_dups) <- c("OG", "Gene")

# Merge dataset with candidates
expression_rice <- merge(rice_expr_DEG, rice_dups, by.x = "GENEID", by.y = "Gene")
write_delim(expression_rice, "results_expression_rice.txt")

################# Significant differences in % of DEGS? ################################

# First, when looking at all duplicated, we don't see much difference when comparing expected and observed
total_genome_sorghum <- 34211
total_genome_maize <- 32540
DEG_sorghum <- 1705
DEG_maize <- 3441
total_dup_maize <- 917
total_dup_sorghum <- 1035
expected_dup_DEG_maize <- DEG_maize * total_dup_maize / total_genome_maize
expected_dup_DEG_sorghum <- DEG_sorghum * total_dup_sorghum / total_genome_sorghum
observed_dup_DEG_maize <- 71
observed_dup_DEG_sorghum <- 54

dup_ret_OGs <- c("OG0000101", "OG0000180", "OG0000379", "OG0000391", "OG0000461", "OG0000505",
                 "OG0000548", "OG0000809", "OG0000861", "OG0001190", "OG0002111", "OG0002585",
                 "OG0002653", "OG0002789", "OG0003194")

maize_dups_ret <- maize_dups %>% filter(OG %in% dup_ret_OGs)
sorghum_dup_ret <- sorghum_dups %>% filter(OG %in% dup_ret_OGs)
total_dupret <- 40
expected_dupret_DEG_maize <- DEG_maize * total_dupret / total_genome_maize
expected_dupret_DEG_sorghum <- DEG_sorghum * total_dupret / total_genome_sorghum
dup_ret_expr <- expression_maize %>% filter(OG %in% dup_ret_OGs)
dup_ret_expr <- expression_sorghum %>% filter(OG %in% dup_ret_OGs)

observed_dupret_DEG <- 5

# contingency tables for chi-squared test
maize <- read_excel("contingency-table.xlsx", sheet = 'Maize-dup')
sorghum <- read_excel("contingency-table.xlsx", sheet = 'Sorghum-dup')

matrix_dup_maize <- maize[c(1,2),2:3]
matrix_ret_maize <- maize[c(1,3),2:3]

maize_dup_results <- chisq.test(matrix_dup_maize)
maize_ret_results <- chisq.test(matrix_ret_maize)

matrix_dup_sorghum <- sorghum[c(1,2),2:3]
matrix_ret_sorghum <- sorghum[c(1,3),2:3]

sorghum_dup_results <- chisq.test(matrix_dup_sorghum)
sorghum_ret_results <- chisq.test(matrix_ret_sorghum)


# calculate proportions for Z-test

p.genome = DEG_maize / total_genome_maize
p.dupret = observed_dupret_DEG / total_dupret
# Pooled proportion
p_combined <- (DEG_maize + observed_dupret_DEG) / (total_genome_maize + total_dupret)
# Standard error of the difference in proportions
se <- sqrt(p_combined * (1 - p_combined) * (1 / total_genome_maize + 1 / total_dupret))
# Z-statistic
z_stat <- (p.genome - p.dupret) / se
# P-value for two-tailed test
p_value <- 2 * pnorm(-abs(z_stat))
# Output results
z_stat
p_value

# Repeat for sorghum
p.genome = DEG_sorghum / total_genome_sorghum
p_combined <- (DEG_sorghum + observed_dupret_DEG) / (total_genome_sorghum + total_dupret)
se <- sqrt(p_combined * (1 - p_combined) * (1 / total_genome_sorghum + 1 / total_dupret))
z_stat <- (p.genome - p.dupret) / se
p_value <- 2 * pnorm(-abs(z_stat))
z_stat
p_value

# Repeat z-test for the duplicated genes (not retained)
p.genome = DEG_sorghum / total_genome_sorghum
p.dup_sorghum = observed_dup_DEG_sorghum / total_dup_sorghum
p_combined <- (DEG_sorghum + observed_dup_DEG_sorghum) / (total_genome_sorghum + total_dup_sorghum)
se <- sqrt(p_combined * (1 - p_combined) * (1 / total_genome_sorghum + 1 / total_dup_sorghum))
z_stat <- (p.genome - p.dup_sorghum) / se
p_value <- 2 * pnorm(-abs(z_stat))
z_stat
p_value

p.genome = DEG_maize / total_genome_maize
p.dup_maize = observed_dup_DEG_maize / total_dup_maize
p_combined <- (DEG_maize + observed_dup_DEG_maize) / (total_genome_maize + total_dup_maize)
se <- sqrt(p_combined * (1 - p_combined) * (1 / total_genome_maize + 1 / total_dup_maize))
z_stat <- (p.genome - p.dup_maize) / se
p_value <- 2 * pnorm(-abs(z_stat))
z_stat
p_value
