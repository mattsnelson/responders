library(tidyverse)
library(Maaslin2)
library(mixOmics)
library(conflicted)

conflict_prefer("select","dplyr")
conflict_prefer("filter", "dplyr")

load("data_wrangled/hamsab_baseline_data.RData")

# str(groups.hamsab_start)
groups.hamsab_start$improved_bp <- as.factor(groups.hamsab_start$improved_bp)
# str(groups.hamsab_start.metadata)
groups.hamsab_start.metadata$improved_bp <- as.factor(groups.hamsab_start.metadata$improved_bp)


# ===== Unadjusted Model =====
folder_path <- "maaslin2/l6_full_taxa_unadjusted"

if (dir.exists(folder_path)) {
  cat("The folder exists. Getting results file")
  fit_data.full.unadj.results <- read.delim(file.path(folder_path, "all_results.tsv"), sep = "\t", header = TRUE)
} else {
  cat("The folder does not exist. Running Maaslin2")
  fit_data.full.unadj = Maaslin2(
    input_data = taxa.full.hamsab_start, 
    input_metadata = groups.hamsab_start, 
    output = folder_path, 
    fixed_effects = c("improved_bp"),
    reference = c("improved_bp,No"),
    transform = "NONE",      
    normalization = "TMM",  
    analysis_method = "NEGBIN", 
    min_abundance =0.1, min_prevalence = 0.25)
  
  fit_data.full.unadj.results <- fit_data.full.unadj$results
}

# PLOT
plot1.unadj <- ggplot(fit_data.full.unadj.results %>% filter(qval < 0.05),     # Filter for significant taxa
               aes(x = coef, y = reorder(feature, coef), colour = coef > 0)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = coef - stderr, xmax = coef + stderr), height = 0.2) +
  scale_colour_manual(values = c("FALSE" = "lightpink", "TRUE" = "lightblue"), 
                      labels = c("FALSE" = "Enriched in Responders", "TRUE" = "Enriched in Non-Responders")) +
  scale_x_continuous(limits = c(-10, 10)) +  # Set x-axis limits
  theme_bw() + 
  labs(x = "Coefficient", y = "Genus", colour = "Legend", title = "Differential Taxa, Unadjusted Model") +
  theme(axis.text.y = element_text(face = "italic"))  # Make y labels italic

plot1.unadj 


# ===== Age-adjusted Model =====
folder_path <- "maaslin2/l6_full_taxa_age"

if (dir.exists(folder_path)) {
  cat("The folder exists. Getting results file")
  fit_data.results <- read.delim(file.path(folder_path, "all_results.tsv"), sep = "\t", header = TRUE)
} else {

fit_data = Maaslin2(
  input_data = taxa.full.hamsab_start, 
  input_metadata = groups.hamsab_start.metadata, 
  output = folder_path, 
  fixed_effects = c("improved_bp", "Age"),
  reference = c("improved_bp,No"),
  transform = "NONE",
  normalization = "TMM",  
  analysis_method = "NEGBIN",
  min_abundance =0.1, min_prevalence = 0.25)

fit_data.results <- fit_data$results
}
## 2. Re-calc Q values
fit_data.results<- fit_data.results %>% filter(metadata == 'improved_bp') 
fit_data.results$qval_updated <- p.adjust(fit_data.results$pval, method = 'BH') # FDR correction using 'BH'

# rename clearer for final plot
fit_data.full.adj_age_results <- fit_data.results

## 3. Plot the model:
plot2_adj_model <- ggplot(fit_data.full.adj_age_results %>% filter(qval_updated < 0.05),   # Filter for significant taxa
               aes(x = coef, y = reorder(feature, coef), colour = coef > 0)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = coef - stderr, xmax = coef + stderr), height = 0.2) +
  scale_colour_manual(values = c("FALSE" = "lightpink", "TRUE" = "lightblue"), 
                      labels = c("FALSE" = "Decreased in Responders", "TRUE" = "Enriched in Responders")) +
  #scale_x_continuous(limits = c(-5, 5)) +  # Set x-axis limits
  theme_bw() + 
  labs(x = "Coefficient", y = "Genus", colour = "Legend", title = "Differential Taxa, Adjusted for Age") +
  theme(axis.text.y = element_text(face = "italic"))  # Make y labels italic

plot2_adj_model

# ===== Age and BMI-adjusted Model =====
folder_path <- "maaslin2/l6_full_taxa_age_bmi"

if (dir.exists(folder_path)) {
  cat("The folder exists. Getting results file")
  fit_data.results <- read.delim(file.path(folder_path, "all_results.tsv"), sep = "\t", header = TRUE)
} else {

fit_data = Maaslin2(
  input_data = taxa.full.hamsab_start, 
  input_metadata = groups.hamsab_start.metadata, 
  output = folder_path, 
  fixed_effects = c("improved_bp", "Age", "BMI"),
  reference = c("improved_bp,No"),
  transform = "NONE",
  normalization = "TMM",
  analysis_method = "NEGBIN",
  min_abundance =0.1, min_prevalence = 0.25)

fit_data.results <- fit_data$results
}
## 2. Re-calc Q values
fit_data.results<- fit_data.results %>% filter(metadata == 'improved_bp')
fit_data.results$qval_updated <- p.adjust(fit_data.results$pval, method = 'BH') # FDR correction using 'BH'

# rename clearer for final plot
fit_data.full.adj_age_bmi_results <- fit_data.results

## 3. Plot the model:
plot3_adj_model <- ggplot(fit_data.full.adj_age_bmi_results %>% filter(qval_updated < 0.05), # Filter for significant taxa
               aes(x = coef, y = reorder(feature, coef), colour = coef > 0)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = coef - stderr, xmax = coef + stderr), height = 0.2) +
  scale_colour_manual(values = c("FALSE" = "lightpink", "TRUE" = "lightblue"), 
                      labels = c("FALSE" = "Decreased in Responders", "TRUE" = "Enriched in Responders")) +
  #scale_x_continuous(limits = c(-5, 5)) +  # Set x-axis limits
  theme_bw() + 
  labs(x = "Coefficient", y = "Genus", colour = "Legend", title = "Differential Taxa, Adjusted for Age and BMI") +
  theme(axis.text.y = element_text(face = "italic"))  # Make y labels italic

plot3_adj_model


# ===== Combined Plot =====
df.unadj <- fit_data.full.unadj.results %>% 
        mutate(qval_updated = qval) %>% # make qval column name same for all
        mutate(adj = "Unadjusted") 
        
df.unadj <- df.unadj %>% 
        select(feature, metadata, value, coef, stderr, qval_updated, adj) # for rbind

df.adj_age <- fit_data.full.adj_age_results %>% 
        mutate(adj = "Adjusted for Age") %>%
        select(feature, metadata, value, coef, stderr, qval_updated, adj) # for rbind
    
df.adj_age_bmi <- fit_data.full.adj_age_bmi_results %>% 
        mutate(adj = "Adjusted for Age & BMI") %>%
        select(feature, metadata, value, coef, stderr, qval_updated, adj) # for rbind

# rbind them to one dataframe:
fit_combined_df <- rbind(df.unadj,
                         df.adj_age,
                         df.adj_age_bmi)

# tidy taxa names for neater graph:
fit_combined_df_filt <- fit_combined_df_filt %>%
  mutate(feature = feature %>%
           str_replace_all("\\.", " ") %>%         # Replace "." with a space
           str_replace_all("  ", " ") %>%         # Replace all doulbe spaces with a signle space
           str_replace_all("__", "") %>%          # Remove all "___"
           str_replace_all("Unclassifed", "Unclassified") %>% # Correct spelling
           str_trim()  # Remove spaces at the end of the name
            )

fit_combined_df_filt$feature

# ===== FILTER =====
fit_combined_df.sig <- fit_combined_df_filt  %>% 
  filter(qval_updated < 0.10) # Filter for significant taxa at Q=0.10

str(fit_combined_df.sig)

# ===== Combined Plot Graph =====
maaslin2_plot <- ggplot(fit_combined_df.sig, aes(x = coef, y = reorder(feature, coef), 
                                color = adj, shape = adj)) +  # colour and shape for degree of adjustment
    # Add background colors based on coef
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf), 
            fill = "#ebf2fb", alpha = 0.2) +  # Background for coef < 0
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf), 
            fill = "#fef3ea", alpha = 0.2) +  # Background for coef > 0
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = coef - stderr, xmax = coef + stderr), height = 0.2) +
  theme_bw() + 
  labs(x = "Coefficient of Enrichment", y = "Genus") +
  scale_shape_manual(values = c(0,2, 15)) +  # diff shape for each adjustment 
  scale_color_manual(values = c("#0b74b3", "#e1892a", "black")) + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        legend.position = "bottom",
        panel.grid.major.y = element_line(size = 0.5))

maaslin2_plot 

ggsave("fig/baseline_taxa2.png", 
       plot = maaslin2_plot, 
       width = 12,
       height = 8,
       dpi = 600)