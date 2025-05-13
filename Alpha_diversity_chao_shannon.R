setwd("D:/prebiotic.exp/miseq/R.figures")

library(ggplot2)
library(ggpubr)
library(ggsignif)

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

#calculates various alpha diversity measures for dataset pcos
tab <- microbiome::alpha(bacteria, index = "all")
#extract metadata from the microbiome dataset stored in pcos.meta
bac.meta <-meta(bacteria)

bac.meta$Shannon <- tab$diversity_shannon
bac.meta$Chao1 <- tab$chao1

#anova test
stat.test1 <- compare_means(Shannon ~ group,bac.meta, method = "anova")
stat.test2 <- compare_means(Chao1 ~ group, bac.meta, method = "anova")

sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.04999){"*"}
  else{NA}}



# Custom function to clean treatment names
clean_treatment_names <- function(df) {
  df$treatment <- gsub("Basil_.*", "Basil", df$treatment)
  df$treatment <- gsub("Flaxseed_.*", "Flaxseed", df$treatment)
  df$treatment <- gsub("Mustard_.*", "Mustard", df$treatment)
  df$treatment <- gsub("Fenugreek_.*", "Fenugreek", df$treatment)
  df$treatment <- gsub("Chia_.*", "Chia", df$treatment)
  df$treatment <- gsub("Mix1%_.*", "Mix1%", df$treatment)
  df$treatment <- gsub("Pectin_.*", "Pectin", df$treatment)
  df$treatment <- gsub("Blank_0hr.*", "Blank_0hr", df$treatment)
  df$treatment <- gsub("Blank_12hr.*", "Blank_12hr", df$treatment)
  
  df$treatment <- factor(df$treatment, levels = c(
    "Blank_0hr", "Blank_12hr", "Pectin", "Basil", "Flaxseed", 
    "Mustard", "Fenugreek", "Chia", "Mix1%"
  ))
  return(df)
}

# Color palette for treatment
treatment_colors <- c(
  "Blank_0hr" = "#808080",
  "Blank_12hr" = "#41424c",
  "Basil" = "#ED7014",
  "Flaxseed" = "#990F02",
  "Mustard" = "#FEC20C",
  "Fenugreek" = "#6F2DA8",
  "Chia" = "#3CB043",
  "Pectin" = "cornflowerblue",
  "Mix1%" = "#622A0F"
)

# Loop over each group
for (grp in unique(bac.meta$group)) {
  
  cat("Processing group:", grp, "\n")
  
  # Subset and clean
  sub_df <- bac.meta[bac.meta$group == grp, ]
  sub_df <- clean_treatment_names(sub_df)
  
  # ANOVA + Tukey test (optional)
  if (length(unique(sub_df$treatment)) > 1) {
    aov_result <- aov(Chao1 ~ treatment, data = sub_df)
    tukey_result <- TukeyHSD(aov_result)
    print(tukey_result)
  }
  
  # Create boxplot
  boxplot <- ggplot(sub_df, aes(x = treatment, y = Chao1, fill = treatment)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
    scale_fill_manual(values = treatment_colors) +
    theme_bw() +
    labs(x = NULL, y = "Chao1", title = paste("Chao1 -", grp)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  # Save boxplot
  ggsave(filename = paste0("Chao1_boxplot_", grp, ".png"),
         plot = boxplot, width = 8, height = 6)
  
  # Create violin plot
  violin <- ggplot(sub_df, aes(x = treatment, y = Chao1, fill = treatment)) +
    geom_violin(trim = FALSE, alpha = 0.9) +
    geom_boxplot(width = 0.2, fill = "white") +
    scale_fill_manual(values = treatment_colors) +
    theme_bw() +
    labs(x = NULL, y = "Chao1", title = paste("Chao1 -", grp)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  # Save violin plot
  ggsave(filename = paste0("Chao1_violin_", grp, ".png"),
         plot = violin, width = 8, height = 6)
}




######################################################################################

# shannon 
# Loop over each group for Shannon
for (grp in unique(bac.meta$group)) {
  
  cat("Processing Shannon for group:", grp, "\n")
  
  # Subset and clean
  sub_df <- bac.meta[bac.meta$group == grp, ]
  sub_df <- clean_treatment_names(sub_df)
  
  # ANOVA + Tukey test (optional)
  if (length(unique(sub_df$treatment)) > 1) {
    aov_result <- aov(Shannon ~ treatment, data = sub_df)
    tukey_result <- TukeyHSD(aov_result)
    print(tukey_result)
  }
  
  # Create boxplot
  boxplot <- ggplot(sub_df, aes(x = treatment, y = Shannon, fill = treatment)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
    scale_fill_manual(values = treatment_colors) +
    theme_bw() +
    labs(x = NULL, y = "Shannon", title = paste("Shannon -", grp)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(filename = paste0("Shannon_boxplot_", grp, ".png"),
         plot = boxplot, width = 8, height = 6)
  
  # Create violin plot
  violin <- ggplot(sub_df, aes(x = treatment, y = Shannon, fill = treatment)) +
    geom_violin(trim = FALSE, alpha = 0.9) +
    geom_boxplot(width = 0.2, fill = "white") +
    scale_fill_manual(values = treatment_colors) +
    theme_bw() +
    labs(x = NULL, y = "Shannon", title = paste("Shannon -", grp)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(filename = paste0("Shannon_violin_", grp, ".png"),
         plot = violin, width = 8, height = 6)
}
