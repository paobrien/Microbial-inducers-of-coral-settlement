## combine results from multiple analysis to identify settlement inducers

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)

## import data ----

## indicator analysis ##

# get file path
fp <- "Data/16S/formatted_data/indicator_species"

# get file names
f_one.s <- list.files(paste0(fp, "/one_species"))

# subset settlement group2
f_one.s <- f_one.s[grepl(".*_g2\\.tsv$", f_one.s)]
f_one.s

# load files
ind_df_one.l <- list()
for (i in seq_along(f_one.s)) {
  ind_df_one.l[[i]] <- read.table(paste0(fp, "/one_species/", f_one.s[i]),
                              header = T,
                              sep = "\t")
  names(ind_df_one.l)[i] <- gsub("indval_16S_|_g2.tsv", "", f_one.s[i])
}


## maaslin ##

# get file path
fp <- "Results/16S/maaslin/"

# get file list
f_maas <- list.files(fp)

# import significant results from each folder
maas_df.l <- list()
for (i in seq_along(f_maas)) {
  maas_df.l[[i]] <- read.table(paste0(fp, f_maas[i], "/significant_results.tsv"),
                               header = T,
                               sep = "\t")
  names(maas_df.l)[i] <- gsub("_tss_lm_log", "", f_maas[i])
}

## random forrests ##

rf_df.l <- readRDS(
  "Results/16S/random_forests/model_importance_regression.rds"
  )

## taxonomy ##

tax <- readRDS(
  "Data/16S/formatted_data/taxonomy_filtered.rds"
)

## format data ----

## indicator analysis ##

# keep significant indicators (AB > 0.5) and remove medium settlement
ind_df_one.l <- lapply(ind_df_one.l, function(x) {
  filter(x, A > 0.50 & B > 0.50 & Group != "Med")
})

# extract feature id from taxonomy
ind_df_one.l <- lapply(ind_df_one.l, function(x) {
  separate(
    x, Taxonomy, c("Feature_ID", "Taxonomy"), sep = "_", extra = "drop"
  ) %>%
    select(-Taxonomy)
})


## maaslin ##

# separate asv_id from taxonomy
# regex to capture character string from the last '..' delimeter
maas_df.l <- lapply(maas_df.l, function(x) {
  extract(x, feature, into = "asv_id", regex = ".*\\.\\.(.*)$")
})

# add feature id to df
maas_df.l <- lapply(maas_df.l, function(x) {
  left_join(x, select(tax, Feature_ID, ASV_id_new), 
            by = join_by(asv_id == ASV_id_new)) %>%
    select(-c(metadata, value, N, N.not.0, asv_id))
})


## random forest ##

# take the top n features of importance and remove taxonomy column (will add in to combined df)
rf_df.l <- lapply(rf_df.l, function(x) {
  x[1:20,] %>%
    select(-c(Taxon, ASV_id_new))
})

## combine results ----

# set df order
df_order <- c("psin", "dfav", "easp", "plob")

# subset to dataframes that are present in each analysis
ind_df.l <- ind_df.l[df_order]
maas_df.l <- maas_df.l[df_order]
rf_df.l <- rf_df.l[df_order]

# merge results from different analyses

# inner join to only keep features that were identified across all analysis
combined_df.l <- list()
for (i in seq_along(ind_df.l)) {
  combined_df.l[[i]] <- inner_join(
    inner_join(
      ind_df.l[[i]], maas_df.l[[i]], by = "Feature_ID"
    ), rf_df.l[[i]], by = "Feature_ID")
  names(combined_df.l)[i] <- df_order[i]
}

# add taxonomy (to dfs that have combined analysis results)
combined_df.l <- lapply(combined_df.l, function(x) {
  if (nrow(x) > 0) {
    x <- left_join(x, tax, by = "Feature_ID")
  }
})

# full join to keep features that were identified in each analysis
combined_df_full.l <- list()
for (i in seq_along(ind_df.l)) {
  combined_df_full.l[[i]] <- full_join(
    full_join(
      ind_df.l[[i]], maas_df.l[[i]], by = "Feature_ID"
    ), rf_df.l[[i]], by = "Feature_ID") %>% 
    left_join(tax, by = "Feature_ID")
  names(combined_df_full.l)[i] <- df_order[i]
}

# reorder rows based on the presence of NA values
combined_df_full.l <- lapply(combined_df_full.l, function(x) {
  x <- x[order(rowSums(is.na(x))), ]
})

# now by stat results
combined_df_full.l <- lapply(combined_df_full.l, function(x) {
  x <- x[order(x$Group, desc(x$coef), desc(x$IncNodePurity)), ]
})

# remove columns to simplify df
combined_df_full.l <- lapply(combined_df_full.l, function(x) {
  select(x, -c(stat, p.value, stderr, pval, qval, Feature_ID))
})

# rearrange columns
combined_df_full.l <- lapply(combined_df_full.l, function(x) {
  select(x, ASV_id_new, Taxon, everything())
})

# rename columns
combined_df_full.l <- lapply(combined_df_full.l, function(x) {
  x <- x %>%
    rename(asv_id = ASV_id_new,
           lm_coef = coef,
           rf_importance = IncNodePurity)
})


## save data ----

for (i in seq_along(combined_df_full.l)) {
  write_tsv(combined_df_full.l[[i]],
            paste0("Results/16S/combined_analyses/", 
                   names(combined_df_full.l)[i],
                   "_combined.tsv"))
}


## check which features are present across all dfs ----

# initialize an empty vector to store the shared values
shared_values <- character()

# iterate over each dataframe in the list
for (df in combined_df_full.l) {
  # get the unique values in the current dataframe's column
  unique_values <- unique(df[["asv_id"]])
  
  # check if the unique values are shared with the previous dataframes
  shared_values <- intersect(shared_values, unique_values)
  
  # if it's the first dataframe, store the unique values as shared values
  if (length(shared_values) == 0) {
    shared_values <- unique_values
  }
}
shared_values

# filter the original dataframes based on the shared values
filtered_dfs <- lapply(combined_df_full.l, function(df) {
  df[df[["asv_id"]] %in% shared_values, ]
})



