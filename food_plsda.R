library(mixOmics)
library(tidyverse)
library(conflicted)
library(janitor)

conflict_prefer("select","dplyr")
conflict_prefer("filter", "dplyr")

load("data_wrangled/hamsab_baseline_data.RData")

# Food Level Analysis
diet.foods <- diet_data.food_serves.hamsab_start # rename for ease of use

# ===== Food serves PCA =====
diet.foods.pca = pca(diet.foods, ncomp = 5)

plot(diet.foods.pca)

plotIndiv(diet.foods.pca, 
           group = groups.hamsab_start$improved_bp, legend = T, legend.title = "Responder",
           ellipse = T, 
           title = "Baseline Diet Food Groups PCA Plot", ind.names = F)
 
plotVar(diet.foods.pca, cex =3, pch = 2, var.names = T, 
        title = 'Foods, PCA, comp 1 - 2') 

# ===== Food Serves PLS-DA =====
diet_plsda <- plsda(X = diet.foods,   # food data
       Y = groups.hamsab_start$improved_bp,   # grouping data
        ncomp = 10)

plotIndiv(diet_plsda , comp = 1:2, 
          group = groups.hamsab_start$improved_bp, ind.names = FALSE,  # colour points by responder status
          ellipse = TRUE,
          legend = TRUE)


# Plot the variables (features) on the first two components
plotVar(diet_plsda, comp = c(1, 2),
        title = "Variable Plot")

# quick check at cutoff =0.5. Wholegrains seem discriminatory
plotVar(diet_plsda, comp = c(1, 2),
        title = "Variable Plot, cutoff=0.5", cutoff = 0.5)

# ===== CHECK PERFORMANCE =====

# Use LOO cross validation
perf.plsda <- perf(diet_plsda, validation = "loo", progressBar = TRUE)

# Plot the classification error rate and optimal ncomp to keep for final
plot(perf.plsda, sd = TRUE, legend.position = 'horizontal')
perf.plsda$choice.ncomp

# rerun with optimal ncomp (which is 2 based on max.dist BER)
diet_plsda <- plsda(X = diet.foods,   # food data
                    Y = groups.hamsab_start$improved_bp,   # grouping data
                    ncomp = 2)

perf.plsda <- perf(diet_plsda, validation = "loo", progressBar = TRUE)
plot(perf.plsda, sd = TRUE, legend.position = 'horizontal')
perf.plsda$error.rate$BER 

# ===== Fig PLSDA =====

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(diet_plsda, comp.predicted=2, dist = "max.dist")

png("fig/Food_Groups_PLS_DA.png", width = 3000, height = 2000, res = 600)

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(diet_plsda, comp = 1:2, 
          group = groups.hamsab_start$improved_bp, ind.names = FALSE,  # colour points by responder status
          background = background, # include prediction background for each class
          legend = TRUE, title = "Food Groups PLS-DA")

dev.off()

# Create Figure for plotLoadings
plotLoadings(diet_plsda)

# Make a nicer coloured plotLoadings graph
comp1_loadings <- as.data.frame(diet_plsda$loadings$X)
comp1_loadings <- data.frame(rownames(comp1_loadings),
                             comp1_loadings$comp1)
colnames(comp1_loadings) <- c("Item", "Loading")

# make nicer names for plotLoadings graph
comp1_loadings <- comp1_loadings %>%
  mutate(
    Nice_Item = case_when(
      Item == "Poultry_serve" ~ "Poultry",
      Item == "Wholegrains_serve" ~ "Wholegrains",
      Item == "Processed_meats_serve" ~ "Processed Meats",
      Item == "Red_meats_serve" ~ "Red Meats",
      Item == "Other_fruit_serve" ~ "Other Fruits",
      Item == "Fruit_juice_serve" ~ "Fruit Juice",
      Item == "Other_starchy_veg_serve" ~ "Other Starchy Vegetables",
      Item == "Other_vegetables_serve" ~ "Other Vegetables",
      Item == "Milk_alternatives_serve" ~ "Milk Alternatives",
      Item == "food_low_LC_N3_serve" ~ "Low-LC N3 Foods",
      Item == "Organ_meats_serve" ~ "Organ Meats",
      Item == "food_high_LC_N3_serve" ~ "High-LC N3 Foods",
      Item == "Yoghurt_serve" ~ "Yoghurt",
      Item == "Soy_products_serve" ~ "Soy Products",
      Item == "red_orange_veg_serve" ~ "Red and Orange Vegetables",
      Item == "Cheese_serve" ~ "Cheese",
      Item == "Eggs_serve" ~ "Eggs",
      Item == "Legumes_veg_serve" ~ "Legumes (Vegetable)",
      Item == "Legumes_protein_serve" ~ "Legumes (Protein)",
      Item == "green_vegetables_serve" ~ "Green Vegetables",
      Item == "Potatoes_serve" ~ "Potatoes",
      Item == "Nuts_seeds_serve" ~ "Nuts and Seeds",
      Item == "Tomatoes_serve" ~ "Tomatoes",
      Item == "Milk_serve" ~ "Milk",
      Item == "Refined_serve" ~ "Refined Grains",
      Item == "Dark_green_vegetables_serve" ~ "Dark Green Vegetables",      
      TRUE ~ Item # Keep original if no match found
    )
  )

#select only the top 10 components from food groups to report on
comp1_loadings_top10 <- comp1_loadings %>%
  arrange(desc(abs(Loading))) %>%
  slice(1:10)

comp1_loadings_top10

# nice plot, top 10 food groups
comp_loadings <- ggplot(comp1_loadings_top10, aes(x = Loading, y = reorder(Nice_Item, -abs(Loading)), fill = Loading > 0)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#ebf2fb", "FALSE" = "#fef3ea")) +
  labs(x = "", y = "Food Group (Serves Per Day)", ) +
  xlim(-0.4, 0.4) +  
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.y = element_text(size=30),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24) ) 

ggsave(comp_loadings, file = "fig/diet_plsda_comploadings.png",
       width = 10, height = 8, dpi = 600)
