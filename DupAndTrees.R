library(dplyr)
library(readr)
library("ggtree")
library("ape")
library(treeio)
library(ggplot2)

setwd("~/Desktop/Sheff/Aristidoideae-genomes/Aristidoideae-resultsNov/")

# Upload all the duplication events identified by Orthofinder
duplications <- read_delim("Duplications.tsv",
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE)

# Filter by support values >0.5 and non-terminal duplications
filtered_dup <- duplications %>% filter(Type == 'Non-Terminal') %>%
                filter(Support >= 0.5)

# Get the number of duplications per node
filtered_dup$`Species Tree Node` <- as.factor(filtered_dup$`Species Tree Node`)

N7_dup <- filtered_dup %>% filter(`Species Tree Node` == 'N7')

results <- data.frame()

species <- c("Tillandsia_fasciculata", "Ananas_comosus", "Juncus_effusus",
          "Carex_scoparia", "Pharus_latifolius", "Streptochaeta_angustifolia",
          "Phragmites_australis", "Chasmanthium_laxum", "Dichantelium_oligosantes",
          "ASEM_RSA5", "Oryza_sativa", "Leersia_perrieri",
          "Zizania_latifolia", "Dendrocalamus_sinicus", "Ampelocalamus_luodianensis",
          "Olyra_latifolia", "Brachypodium_distachyon", "Hordeum_vulgareMorex",
          "Hordeum_spontaneum", "Stipagrostis_hirtigluma", "Aristida_adscensionis", "Eragrostis_curvula",
          "Eleusine_coracana", "Sorghum_bicolor", "Zea_mays", "Echinocloa_haploclada",
          "ASEM_AUS1", "Urochloa_fusta", "Cenchrus_americanus", "Setaria_viridis",
          "Panicum_virgatum", "Panicum_hallii", "Digitaria_exilis", "Digitaria_radicosa")

for (node in levels(filtered_dup$`Species Tree Node`)) {
  node_dup <- filtered_dup %>% filter(`Species Tree Node` == node)
  num_node_dup <- nrow(node_dup)
  results <- rbind(results, c(node, num_node_dup))
}

colnames(results) <- c('Node', '#Dup')

# Import species tree into R and transform it in an object
tree <- treeio::read.tree("SpeciesTree_Gene_Duplications_0.5_Support.txt")
#object_tree <- as.phylo(tree)
#tree$root.edge <- 0.02
photosynthetic_type <- list(`non-C4`=c("Tillandsia_fasciculata", "Ananas_comosus", "Juncus_effusus",
                                     "Carex_scoparia", "Pharus_latifolius", "Streptochaeta_angustifolia",
                                     "Phragmites_australis", "Chasmanthium_laxum", "Dichantelium_oligosantes",
                                     "Alloteropsis_semialata_RSA5", "Oryza_sativa", "Leersia_perrieri",
                                     "Zizania_latifolia", "Dendrocalamus_sinicus", "Ampelocalamus_luodianensis",
                                     "Olyra_latifolia", "Brachypodium_distachyon", "Hordeum_vulgareMorex",
                                     "Hordeum_spontaneum"),
                            C4=c("Stipagrostis_hirtigluma", "Aristida_adscensionis", "Eragrostis_curvula",
                                 "Eleusine_coracana", "Sorghum_bicolor", "Zea_mays", "Echinocloa_haploclada",
                                 "Alloteropsis_semialata_AUS1", "Urochloa_fusta", "Cenchrus_americanus", "Setaria_viridis",
                                 "Panicum_virgatum", "Panicum_hallii", "Digitaria_exilis", "Digitaria_radicosa"))
tree <- groupOTU(tree, .node=photosynthetic_type)
plot <- ggtree(tree, size=0.8) + geom_tippoint(aes(color=group), size = 3, alpha=0.6) +
  scale_color_manual(values=c("firebrick", "black")) + theme(legend.position='none')
tiff(file='species_tree-v2.tiff', width = 14, height = 14, res=300, units='cm')
print(plot)
dev.off()

###################### Check differences in gene copy number in PACMAD and BOP #########################

library(readxl)
library(tidyr)

SelectedOG_GeneCount <- read_excel("functional-analyses/SelectedOG.GeneCount.xlsx", sheet = "ToR")
GeneCount_long <- pivot_longer(SelectedOG_GeneCount, cols = colnames(SelectedOG_GeneCount)[4:18],
                               names_to = 'OG')

GeneCount_long <- GeneCount_long %>% dplyr::mutate(rel_genecopy = value / Ploidy)
GeneCount_long$Group <- factor(GeneCount_long$Group, levels = c('PACMAD', 'BOP', 'OUT'))


plot <- ggplot(GeneCount_long) + geom_boxplot(aes(x = Group, y = rel_genecopy), outlier.shape = NA) +
  geom_jitter(aes(x = Group, y = rel_genecopy, color = Group), alpha = 0.6) + scale_color_manual(values=c("firebrick", "black", "#556b2f")) +
  facet_wrap(~ OG, scales = 'free_y') + labs(y='Gene copy number / ploidy level') +
  theme_classic() +
  theme(legend.position = 'inside', legend.position.inside = c(0.9, 0.1), axis.title.x = element_blank(), axis.text.x = element_text(size=8))
tiff('GeneCounts_plot.tiff', units='cm', width = 20, height = 14, res = 200)
print(plot)
dev.off()

Orthogroups <- levels(as.factor(GeneCount_long$OG))
for (Orth in Orthogroups) {
  print(Orth)
  data <- GeneCount_long %>% dplyr::filter(OG == Orth) %>% dplyr::filter(Group != 'OUT')
  PACMAD <- data %>% dplyr::filter(Group == 'PACMAD')
  BOP <- data %>% dplyr::filter(Group == 'BOP')
  print(shapiro.test(PACMAD$rel_genecopy))
  print(shapiro.test(BOP$rel_genecopy))
  print(wilcox.test(data$rel_genecopy ~ data$Group))
  print(t.test(data$rel_genecopy ~ data$Group))
}


########## put everything in table and apply multiple testing correction ###########

library(dplyr)

# Initialize a data frame to store results
results <- data.frame(Orthogroup = character(),
                      Shapiro_PACMAD = numeric(),
                      Shapiro_BOP = numeric(),
                      Wilcox_P = numeric(),
                      Ttest_P = numeric(),
                      stringsAsFactors = FALSE)

# Get the list of Orthogroups
Orthogroups <- levels(as.factor(GeneCount_long$OG))

# Loop over each orthogroup
for (Orth in Orthogroups) {

  print(Orth)  # Print current orthogroup

  # Subset the data for the current orthogroup and exclude 'OUT' group
  data <- GeneCount_long %>%
    dplyr::filter(OG == Orth) %>%
    dplyr::filter(Group != 'OUT')

  # Subset data for each group
  PACMAD <- data %>% dplyr::filter(Group == 'PACMAD')
  BOP <- data %>% dplyr::filter(Group == 'BOP')

  # Perform Shapiro-Wilk test (normality test) for each group
  shapiro_PACMAD <- shapiro.test(PACMAD$rel_genecopy)$p.value
  shapiro_BOP <- shapiro.test(BOP$rel_genecopy)$p.value

  # Perform Wilcoxon test (non-parametric test) for differences between groups
  wilcox_p <- wilcox.test(rel_genecopy ~ Group, data = data)$p.value

  # Perform t-test (parametric test) for differences between groups
  ttest_p <- t.test(rel_genecopy ~ Group, data = data)$p.value

  # Add the results to the data frame
  results <- rbind(results, data.frame(
    Orthogroup = Orth,
    Shapiro_PACMAD = shapiro_PACMAD,
    Shapiro_BOP = shapiro_BOP,
    Wilcox_P = wilcox_p,
    Ttest_P = ttest_p
  ))
}

# Apply Benjamini-Hochberg correction for multiple testing (adjust p-values)
results$Shapiro_PACMAD_adj <- p.adjust(results$Shapiro_PACMAD, method = "BH")
results$Shapiro_BOP_adj <- p.adjust(results$Shapiro_BOP, method = "BH")
results$Wilcox_P_adj <- p.adjust(results$Wilcox_P, method = "BH")
results$Ttest_P_adj <- p.adjust(results$Ttest_P, method = "BH")

# Print the results table with adjusted p-values
print(results)
# Write the results to a text file (replace "results.txt" with your desired file name)
write.table(results, file = "gene_copy_number_testresults.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
