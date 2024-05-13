## Format 16S data from QIIME2

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(vegan)

## import data ----

# to read in orginal data
# ASV table
asv <- read.table("Data/16S/feature-table.tsv", 
                  header = T, 
                  sep = "\t",
                  strip.white = T)

# taxonomy
tax <- read.table("Data/16S/taxonomy.tsv", 
                  header = T,
                  sep = "\t",
                  strip.white = T)

# metadata
meta <- read.table("Data/16S/metadata.tsv", 
                   header = T,
                   sep = "\t",
                   strip.white = T)

# tree
tree <- read.tree("Data/16S/rooted-tree.nwk")


## format data ----

# combine asv with taxa
df <- full_join(asv, tax[,1:2], by = "Feature_ID")

## filter Chloroplast, Mitochondria, Eukaryota
df <- filter(df, !grepl('Eukaryota|Chloroplast|Mitochondria', 
                        Taxon, ignore.case = T))

## filter out repeat samples (samples that were repeated due to poor quality)
# get repeat sample ids
repeats <- data.frame(table(meta$Sample)) # makes a table of the sample_ids and frequencies
repeats <- repeats[repeats$Freq > 1,]
repeats <- repeats$Var1

# filter meta to repeats
meta_rep <- filter(meta, meta$Sample %in% repeats)

# subset repeats to sample with lowest sequences (to remove)
meta_rm <- meta_rep %>% 
  distinct %>% 
  group_by(Sample) %>% 
  top_n(-1, filtered) %>%
  as.data.frame

# get sample ids to remove
to_rm <- meta_rm$Sample_id

# filter from dataframe
df <- select(df, -all_of(to_rm))

# filter from metadata
meta <- filter(meta, !meta$Sample_id %in% to_rm)

## Re-run to remove additional repeats (some sequenced 3 times)

# get repeat sample ids
repeats <- data.frame(table(meta$Sample)) 
repeats <- repeats[repeats$Freq > 1,]
repeats <- repeats$Var1

# filter meta to repeats
meta_rep <- filter(meta, meta$Sample %in% repeats)

# subset repeats to sample with lowest sequences (to remove)
meta_rm <- meta_rep %>% 
  distinct %>% 
  group_by(Sample) %>% 
  top_n(-1, filtered) %>%
  as.data.frame

# get sample ids to remove
to_rm <- meta_rm$Sample_id

# filter from dataframe
df <- select(df, -all_of(to_rm))

# filter from metadata
meta <- filter(meta, !meta$Sample_id %in% to_rm)

# filter from taxonomy file
tax <- select(df, Feature_ID, Taxon)

# filter from asv file
asv <- select(df, !c(Taxon))

## Filter any asvs that are no longer present in the data
# get feature_IDS
to_rm <- asv$Feature_ID[rowSums(asv[,-1]) == 0]

# filter from asv table
asv <- filter(asv, !asv$Feature_ID %in% to_rm)

# filter from tax 
tax <- filter(tax, !tax$Feature_ID %in% to_rm)

# add new unique asv_id for readability
tax$ASV_id_new <- rep(paste0("asv_", 1:nrow(tax)))

# add new settlement group to metadata
settlement <- function(df) {
  tmp <- NULL
  for (i in seq_along(df$Percent_settled)) {
    if (df$Percent_settled[i] < 35) {
      tmp[i] <- "Low"
    } else if (df$Percent_settled[i] >= 35 && df$Percent_settled[i] <70) {
      tmp[i] <- "Med"
    } else {
      tmp[i] <- "High"
    }
  }
  return(tmp)
}
meta$Settlement_2 <- settlement(meta)

# save files
saveRDS(asv, file = "Data/16S/formatted_data/asv_filtered.rds")
saveRDS(tax, file = "Data/16S/formatted_data/taxonomy_filtered.rds")
write.table(meta, file = "Data/16S/formatted_data/metadata.tsv", sep = "\t")


## create rarefied dataset ----

# edit table
asv_r <- asv[,-1]
rownames(asv_r) <- asv$Feature_ID
asv_r <- t(asv_r)

# drop sample <4000 reads
asv_r <- asv_r[rowSums(asv_r) > 4000,]

# rarefy to depth
asv_r <- rrarefy(asv_r, sample = 4000)

# save rarefied table
saveRDS(asv_r, "Data/16S/formatted_data/asv_rarefied.rds")






