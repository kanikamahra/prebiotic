setwd("D:/prebiotic.exp/miseq/R.figures")

bacn <- qza_to_phyloseq(
  features = "D:/prebiotic.exp/miseq/R.figures/table_final.qza", 
  tree = "D:/prebiotic.exp/miseq/R.figures/rooted-tree.qza", 
  taxonomy = "D:/prebiotic.exp/miseq/R.figures/taxonomy_final.qza", 
  metadata = "D:/prebiotic.exp/miseq/R.figures/metadata.txt"
)

sample_names(bacn)
rank_names(bacn)
sample_variables(bacn)

#remove no abundance taxa, noise reduction,focus dominant taxa only 
bacteria <- prune_taxa(taxa_sums(bacn)>1, bacn)

# Convert treatment to character
treatment_char <- as.character(sample_data(bacteria)$treatment)

# taxonomy_barplots.R

# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define custom color palette (extend or adjust as needed)
custom_palette <- c("yellow","#1f77b4", "#ff7f0e","darkblue", "#2ca02c", "#d62728", "#9467bd","#999933","#D55E00",
                    "#882255","#009E73","#8c564b", "#bcbd22","#98df8a","#e377c2", "#17becf","darkgreen","#0072B2",
                    "#DDCC77","#332288","#aec7e8", "#ffbb78", "#ff9896", "cyan","purple", "#c49c94","aquamarine","black")

# List of samples to remove
samples_to_remove <- c("A1_1", "A1_2", "A1_3", "A2_1", "A2_2", "A2_3",
                       "A3_1", "A3_2", "A3_3", "S1_1", "S1_2", "S1_3",
                       "S2_1", "S2_2", "S2_3", "S3_1", "S3_2", "S3_3")

# Define the desired order for treatments
desired_order <- c("Blank_12hr_1", "Blank_12hr_2", "Blank_12hr_3",
                   "Pectin_1", "Pectin_2", "Pectin_3",
                   "Basil_1", "Basil_2", "Basil_3",
                   "Chia_1", "Chia_2", "Chia_3",
                   "Fenugreek_1", "Fenugreek_2", "Fenugreek_3",
                   "Flaxseed_1", "Flaxseed_2", "Flaxseed_3",
                   "Mustard_1", "Mustard_2", "Mustard_3",
                   "Mix1%_1", "Mix1%_2", "Mix1%_3")

#-------------------------------- PHYLUM LEVEL ------------------------------------------------#

bacnew_rel <- transform_sample_counts(bacteria, function(x) x / sum(x) * 100)
glom_phylum <- tax_glom(bacnew_rel, taxrank = "Phylum", NArm = FALSE)
bacnew_melt <- psmelt(glom_phylum)
bacnew_melt$Phylum <- as.character(bacnew_melt$Phylum)

bacnew_melt_filtered <- bacnew_melt %>%
  filter(!is.na(Phylum)) %>%
  filter(!(Sample %in% samples_to_remove)) %>%
  mutate(treatment = factor(treatment, levels = desired_order))

taxonomy_barplot_phylum <- ggplot(bacnew_melt_filtered, aes(x = treatment, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_palette) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum") +
  facet_wrap(~ group, scales = "free_x", ncol = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))
print(taxonomy_barplot_phylum)

#-------------------------------- FAMILY LEVEL ------------------------------------------------#

glom_family <- tax_glom(bacnew_rel, taxrank = "Family", NArm = FALSE)
bacnew_melt_family <- psmelt(glom_family)
bacnew_melt_family$Family <- as.character(bacnew_melt_family$Family)

top_25_families <- bacnew_melt_family %>%
  group_by(Family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  top_n(25, total_abundance) %>%
  pull(Family)

top_25_bacnew_melt_family <- bacnew_melt_family %>%
  filter(Family %in% top_25_families) %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

top_25_bacnew_melt_family_filtered <- top_25_bacnew_melt_family %>%
  filter(!(Sample %in% samples_to_remove)) %>%
  mutate(treatment = factor(treatment, levels = desired_order))

family_barplot_top_25_relative <- ggplot(top_25_bacnew_melt_family_filtered, aes(x = treatment, y = Relative_Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_palette) +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Family") +
  facet_wrap(~ group, scales = "free_x", ncol = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))
print(family_barplot_top_25_relative)

#-------------------------------- GENUS LEVEL ------------------------------------------------#

glom_genus <- tax_glom(bacnew_rel, taxrank = "Genus", NArm = FALSE)
bacnew_melt_genus <- psmelt(glom_genus)
bacnew_melt_genus$Genus <- as.character(bacnew_melt_genus$Genus)

top_25_genera <- bacnew_melt_genus %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  top_n(25, total_abundance) %>%
  pull(Genus)

top_25_bacnew_melt_genus <- bacnew_melt_genus %>%
  filter(Genus %in% top_25_genera) %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

top_25_bacnew_melt_genus_filtered <- top_25_bacnew_melt_genus %>%
  filter(!(Sample %in% samples_to_remove)) %>%
  mutate(treatment = factor(treatment, levels = desired_order))

genus_barplot_top_25_relative <- ggplot(top_25_bacnew_melt_genus_filtered, aes(x = treatment, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_palette) +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Genus") +
  facet_wrap(~ group, scales = "free_x", ncol = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))
print(genus_barplot_top_25_relative)
