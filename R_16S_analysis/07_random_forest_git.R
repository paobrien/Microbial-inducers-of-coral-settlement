# Random Forests analysis

# Using random forests machine learning method to identify settlement inducing microbes

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(randomForest)

## import data ----

# filtered by 0.1% relative abundance and 10% prevalence
asv_one_sp.l <- readRDS(
  "Data/16S/formatted_data/networks/one_species/asv_one_species_list_f.rds"
  )

# metadata
meta_one_sp.l <- readRDS(
  "Data/16S/formatted_data/networks/one_species/meta_one_species_list.rds"
  )

# taxonomy
tax <- readRDS(
  "Data/16S/formatted_data/taxonomy_filtered.rds"
  )

## format data ----

# attach metadata
for (i in seq_along(asv_one_sp.l)) {
  asv_one_sp.l[[i]] <- left_join(asv_one_sp.l[[i]],
                                 meta_one_sp.l[[i]],
                                 join_by(Sample_ID == Sample_id))
}

# get row names
asv_one_sp.l <- lapply(asv_one_sp.l, function(x) {
  rownames(x) <- x$Sample_ID
  return(x)
})

# make model df - remove unwanted metadata

## classification model
df_class.l <- lapply(asv_one_sp.l, function(x) {
  x <- select(x, -c(Sample_ID, Coral, Treatment, Settlement, Percent_settled))
})

## regression model
df_reg.l <- lapply(asv_one_sp.l, function(x) {
  x <- select(x, -c(Sample_ID, Coral, Treatment, Settlement, Settlement_2))
})


## build the Random Forest model ----

set.seed(1107)

# build classification model (discrete data)
# get outcome (i.e., response) variable
out_var_class <- lapply(df_class.l, function(x) {
  x <- as.factor(x$Settlement_2)
})

# remove outcome from main df
df_class.l <- lapply(df_class.l, function(x) {
  x <- select(x, -Settlement_2)
})

# run model
rf_class <- list()
for (i in seq_along(df_class.l)) {
  rf_class[[i]] <- randomForest(x = df_class.l[[i]],
                                y = out_var_class[[i]])
  names(rf_class)[i] <- names(df_class.l)[i]
}

## build regression model (continuous data)
# get outcome variable
out_var_reg <- lapply(df_reg.l, function(x) {
  x <- x$Percent_settled
})

# remove outcome from main df
df_reg.l <- lapply(df_reg.l, function(x) {
  x <- select(x, -Percent_settled)
})

# run model
rf_reg <- list()
for (i in seq_along(df_reg.l)) {
  rf_reg[[i]] <- randomForest(x = df_reg.l[[i]],
                              y = out_var_reg[[i]])
  names(rf_reg)[i] <- names(df_reg.l)[i]
}

# feature importance (higher values = better predictive ability) 
head(importance(rf_class$psin))
head(importance(rf_reg$psin))

# prediction accuracy on training data
for (i in seq_along(rf_class)) {
  print(names(rf_class)[i])
  print(rf_class[[i]])
}

for (i in seq_along(rf_reg)) {
  print(names(rf_reg)[i])
  print(rf_reg[[i]])
}

# visualize the first tree in the Random Forest
plot(rf_class$psin$forest[[1]])
plot(rf_reg$dfav)

## save model ----

saveRDS(rf_class, "Results/16S/random_forests/rf_classification_list.rds")
saveRDS(rf_reg, "Results/16S/random_forests/rf_regression_list.rds")

## format output ----

# get asv predictors
rf_importance_class <- lapply(rf_class, function(x) {
  importance(x) %>% as.data.frame
}) 

rf_importance_reg <- lapply(rf_reg, function(x) {
  importance(x) %>% as.data.frame
}) 

# make feature id column
rf_importance_class <- lapply(rf_importance_class, function(x) {
  x$Feature_ID <- rownames(x)
  return(x)
})

rf_importance_reg <- lapply(rf_importance_reg, function(x) {
  x$Feature_ID <- rownames(x)
  return(x)
})

# combine with taxonomy
rf_importance_class <- lapply(rf_importance_class, function(x) {
  x <- left_join(x, tax, by = "Feature_ID")
})

rf_importance_reg <- lapply(rf_importance_reg, function(x) {
  x <- left_join(x, tax, by = "Feature_ID")
})

# sort highest to lowest
rf_importance_class <- lapply(rf_importance_class, function(x) {
  x <- x[order(x$MeanDecreaseGini, decreasing = T),]
})

rf_importance_reg <- lapply(rf_importance_reg, function(x) {
  x <- x[order(x$IncNodePurity, decreasing = T),]
})

## save importance ----

saveRDS(rf_importance_class, 
        "Results/16S/random_forests/model_importance_classification.rds")

saveRDS(rf_importance_reg, 
        "Results/16S/random_forests/model_importance_regression.rds")

# as table
for (i in seq_along(rf_importance_reg)) {
  write_tsv(as.data.frame(rf_importance_reg[[i]]),
            paste0("Results/16S/random_forests/", 
                   names(rf_importance_reg)[i],
                   "_importance_regression.tsv"))
}

## save confusion matrix ----

for (i in seq_along(rf_class)) {
  write_tsv(as.data.frame(rf_class[[i]]$confusion),
            paste0("Results/16S/random_forests/", 
                   names(rf_class)[i],
                   "_confusion_matrix.tsv"))
}




