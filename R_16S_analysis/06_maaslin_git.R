# MaAsLin2 - Microbiome Multivariable Association with Linear Models

# Using linear models to identify settlement inducing microbes

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(Maaslin2)

## import data ----

asv <- readRDS("Data/16S/formatted_data/asv_filtered.rds")
tax <- readRDS("Data/16S/formatted_data/taxonomy_filtered.rds")
meta <- read.table("Data/16S/formatted_data/metadata_files/metadata.tsv", 
                   header = T,
                   sep = "\t",
                   strip.white = T)

meta <- select(meta, Sample_id, Coral, Treatment, Tank, 
               Settlement, Percent_settled, Settlement_2)
meta_2ml <- filter(meta, Treatment == "Light_2M") 

## create data frames ----

# subset to one species dfs (note: easp and plob contain other data, subset to L2M)
#metadata
meta_psin <- filter(meta, Coral == "Psin")
meta_dfav <- filter(meta, Coral == "Dfav")
meta_easp <- filter(meta_2ml, Coral == "Easp")
meta_plob <- filter(meta_2ml, Coral == "Plob")

# asvs
asv_psin <- select(asv, Feature_ID, meta_psin$Sample_id) 
asv_dfav <- select(asv, Feature_ID, meta_dfav$Sample_id) 
asv_easp <- select(asv, Feature_ID, meta_easp$Sample_id) 
asv_plob <- select(asv, Feature_ID, meta_plob$Sample_id)

# group into list
# get names
one_sp.n <- c("psin", "dfav", "easp", "plob")

#list meta
meta_one_sp.l <- list(
  meta_psin,
  meta_dfav,
  meta_easp,
  meta_plob
)
names(meta_one_sp.l) <- one_sp.n

#list asvs
asv_one_sp.l <- list(
  asv_psin,
  asv_dfav,
  asv_easp,
  asv_plob,
)
names(asv_one_sp.l) <- one_sp.n

# remove asvs no longer in df
# function - filter by value (n = lowest number of ASVs) and then transpose
asv_filter <- function(df, n) {
  rownames(df) <- df$Feature_ID
  df <- df[,-1]
  df <- df[rowSums(df !=0) > n,] # edit n different filtering
  df <- t(df) %>% as.data.frame
  df$Sample_ID <- rownames(df)
  df <- select(df, Sample_ID, everything())
  return(df)
}

# run
asv_one_sp.l <- lapply(asv_one_sp.l, function(x) {
  asv_filter(x, 0)
})

# check number of asvs
lapply(asv_one_sp.l, dim)

## filter out low abundance and prevalence asvs
## using function filter_taxa

asv_filter_ra_pres <- function(df, n, p) {
  # change to rel anbundnace
  tmp <- df[,-1] / rowSums(df[,-1])*100
  # transpose so columns are samples
  tmp <- t(tmp)
  # run filter function
  tmp <- filter_taxa(tmp, n, p)
  # select filtered taxa from data from
  df <- select(df, Sample_ID, rownames(tmp))
  return(df)
}

# 0.1% rel abundance + 10% prevalence
asv_one_sp.l <- lapply(asv_one_sp.l, function(x) {
  asv_filter_ra_pres(x, 0.1, 10)
})

# check asv number
lapply(asv_one_sp.l, dim)

# check lowest occurence for asv
lapply(asv_one_sp.l, function(x) {
  print(min(colSums(x[,-1])))
})

# add sample columns
asv_one_sp.l <- lapply(asv_one_sp.l, function(x) {
  rownames(x) <- x[,1] 
  select(x, Sample_ID, everything())
})

## save ----

# save asv lists
# filtered rel abund and pres
saveRDS(asv_one_sp.l, "Data/16S/formatted_data/networks/one_species/asv_one_species_list_f.rds")

# save metadata lists
saveRDS(meta_one_sp.l, "Data/16S/formatted_data/networks/one_species/meta_one_species_list.rds")


## import filtered data 
asv_one_sp.l <- readRDS("Data/16S/formatted_data/networks/one_species/asv_one_species_list_f.rds")
meta_one_sp.l <- readRDS("Data/16S/formatted_data/networks/one_species/meta_one_species_list.rds")


## format data ----

# Format taxonomy table
tax <- separate(
  tax, Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
  sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn"
  )

# replace NAs with Unknown
tax[is.na(tax)] <- "Unknown"

# remove underscores from tax string
tax <- data.frame(lapply(tax, function(x) gsub(".+__", "", x)))

# make new tax id column
tax <- unite(tax, "ASV_tax", 
                     c("Phylum", "Family", "Genus", "ASV_id_new"), 
                     sep = "; ", remove = F)

# make the new tax ids asv row names
# first remove sample id column, transpose and get Feature ID column
asv_one_sp.l <- lapply(asv_one_sp.l, function(x) {
  x <- select(x, -Sample_ID)
  x <- t(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x$Feature_ID <- rownames(x)
  return(x)
})

# join with taxonomy by Feature ID column
asv_one_sp.l <- lapply(asv_one_sp.l, function(x) {
  x <- left_join(x, select(tax, Feature_ID, ASV_tax), by = "Feature_ID")
})

# make new tax ids row names and remove unnecessary cols
asv_one_sp.l <- lapply(asv_one_sp.l, function(x) {
  rownames(x) <- x$ASV_tax
  x <- select(x, -c(Feature_ID, ASV_tax))
  return(x)
})

# transpose back so taxa are columns
asv_one_sp.l <- lapply(asv_one_sp.l, t)

# format metadata
# change row names to sample id and remove column
meta_one_sp.l <- lapply(meta_one_sp.l, function(x) {
  rownames(x) <- x$Sample_id
  x <- select(x, -Sample_id)
})


## run maaslin ----

# continuous variable
fit_data <- Maaslin2(input_data = asv_one_sp.l$psin,
                     input_metadata = meta_one_sp.l$psin,
                     output = "Results/16S/maaslin/psin_tss_lm_log",
                     normalization = "TSS",
                     transform = "LOG",
                     analysis_method = "LM",
                     random_effects = "Tank",
                     fixed_effects = "Percent_settled",
                     min_abundance = 0.01,
                     min_prevalence = 0.1
                     )

# categorical variable
fit_data <- Maaslin2(input_data = asv_one_sp.l$psin,
                     input_metadata = meta_one_sp.l$psin,
                     output = "Results/16S/maaslin/psin_tmm_negbin_cat",
                     normalization = "TSS",
                     transform = "NONE",
                     analysis_method = "NEGBIN",
                     random_effects = c("Tank", "Coral"),
                     fixed_effects = "Settlement_2",
                     reference = c("Settlement_2,Low"),
                     min_abundance = 0.01,
                     min_prevalence = 0.1
                     )


## investigate results ----

# get significant results
sig_results <- filter(fit_data$results, qval < 0.25)

# format feature id
sig_results$feature <- gsub("^X", "", sig_results$feature)

# add taxonomy
sig_results <- left_join(sig_results, tax, join_by(feature == Feature_ID))
sig_results[,c("feature", "Taxon")]




