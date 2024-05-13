# Indicator species analysis using 16S data

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(indicspecies)

## import and format data ----

# import data
asv <- readRDS("Data/16S/formatted_data/asv_filtered.rds")
tax <- readRDS("Data/16S/formatted_data/taxonomy_filtered.rds")
meta <- read.table("Data/16S/formatted_data/metadata_files/metadata.tsv", 
                   header = T,
                   sep = "\t",
                   strip.white = T)

## remove blanks and controls
# get sample ids
to_rm <- filter(meta, meta$Treatment %in% 
                  c("Blank", "Tile_ctrl", "Rubble_ctrl", 
                    "tank_water", "room_water", "Larvae")) %>% .$Sample_id
# remove from data
asv <- dplyr::select(asv, -all_of(to_rm))
meta <- filter(meta, !meta$Sample_id %in% to_rm)

# remove asvs no longer in df
rownames(asv) <- asv$Feature_ID
asv <- asv[,-1]
asv <- asv[rowSums(asv != 0) > 0,]

# transform counts to relative abundances
Feature_ID <- rownames(asv) # save feature IDS (rownames removed below - somtimes?)
asv <- t(asv)
asv <- asv / rowSums(asv) * 100

# check transformation
rowSums(asv)[1:10]

# add featuer ids
asv <- t(asv) %>% as.data.frame
asv$Feature_ID <- Feature_ID

# edit tax names
tax$Taxon <- gsub("d__Bacteria; ", "", tax$Taxon)
tax$Taxon <- gsub(" ", "", tax$Taxon)
tax$Taxon <- gsub("__", "_", tax$Taxon)

# combine asv with tax
asv <- right_join(tax, asv, by = "Feature_ID")

# join tax with feature id
asv$Feature_taxon <- paste(asv$Feature_ID, asv$Taxon, sep = "_")
asv <- dplyr::select(asv, Feature_taxon, everything())

# get row names
rownames(asv) <- asv$Feature_taxon

# transpose (and remove tax columns)
asv <- t(asv[,-c(1:3)]) %>% as.data.frame()

# check for NAs
asv[is.na(asv)]


## create dataframes for testing ----

# subset to L2M data
# metadata
meta.l2m <- filter(meta, Treatment == "Light_2M")

# subset to one species (note: easp and plob contains other treatments, subset from L2M data)
#metadata
meta.psin <- filter(meta, Coral == "Psin")
meta.dfav <- filter(meta, Coral == "Dfav")
meta.easp <- filter(meta.l2m, Coral == "Easp")
meta.plob <- filter(meta.l2m, Coral == "Plob")

# asvs
asv.psin <- filter(asv, rownames(asv) %in% meta.psin$Sample_id) 
asv.dfav <- filter(asv, rownames(asv) %in% meta.dfav$Sample_id) 
asv.easp <- filter(asv, rownames(asv) %in% meta.easp$Sample_id) 
asv.plob <- filter(asv, rownames(asv) %in% meta.plob$Sample_id)

# make into df list
# names
one.sp_names <- c("psin", "dfav", "easp", "plob")

# meta
one.sp_meta <- list(meta.psin, meta.dfav, meta.easp, meta.plob)
names(one.sp_meta) <- one.sp_names

# asvs
one.sp_asv <- list(asv.psin, asv.dfav, asv.easp, asv.plob)
names(one.sp_asv) <- one.sp_names


## run indicator sp analysis ----

# across species list
indval.one.sp.l <- list()
for (i in seq_along(one.sp_asv)) {
  print(paste0("running analysis ", one.sp_names[i]))
  indval.one.sp.l[[i]] <- multipatt(one.sp_asv[[i]], 
                                    one.sp_meta[[i]]$Settlement_2, 
                                    func = 'IndVal.g', 
                                    duleg = TRUE,
                                    control=how(nperm=9999))
}
names(indval.one.sp.l) <- one.sp_names

# save
saveRDS(indval.one.sp.l, "Results/16S/indicator_species/one_species/indval_one_sp_list_g2.rds")


## summarise indicator analysis ----

# note: A = probability that the surveyed site belongs to the target site (in this case settlement group)
#       B = probability of finding the species in sites belonging to the site group

summary(indval.one.sp.l$psin, indvalcomp = TRUE)
summary(indval.one.sp.l$dfav, indvalcomp = TRUE)
summary(indval.one.sp.l$easp, indvalcomp = TRUE)
summary(indval.one.sp.l$plob, indvalcomp = TRUE)


## capture output as file ----

# check console width
width <- getOption("width")
width

# change to prevent text wrapping
options(width = 10000)

# one species
for (i in seq_along(indval.one.sp.l)) {
  capture.output(summary(indval.one.sp.l[[i]],
                         indvalcomp = TRUE),
                 file = paste0("Results/16S/indicator_species/indval_16S_", 
                               names(indval.one.sp.l[i]),
                               "_g2.txt"))
}

# to convert output from txt to tsv file - run the following in terminal
# awk -v OFS="\t" '{$1=$1; print}' input.txt > output.txt

# formatted output file moved to "/Data/16S/formatted_data/indicator_species"

# to change width back to default
getOption(width = width)

## add asv column ----

# after formatting, add asv column for plotting

# import tax file (new asv ids added after analysis)
tax <- readRDS("Data/16S/formatted_data/taxonomy_filtered.rds")

# import indval data
# get file path
fp <- "Data/16S/formatted_data/indicator_species"

# get file names with paths
# note: use recursive to include subdirectories, full names to include file paths
fp <- list.files(path = fp, pattern = "\\.tsv$", full.names = T, recursive = T)

# get file names
f.names <- c()
for (i in seq_along(fp)) {
  f.names[i] <- basename(fp[i])
}

# make dataframe list
df.l <- list()
for (i in seq_along(fp)) {
  df.l[[i]] <- read.table(fp[[i]], header = T, sep = "\t")
}
names(df.l) <- f.names

# split feature id and tax into two columns
for (i in seq_along(df.l)) {
  df.l[[i]] <- separate(df.l[[i]],
                        Taxonomy, 
                        c("Feature_ID", "Taxonomy"),
                        "_",
                        remove = T,
                        extra = "merge")
}

# add new asv id colums
for (i in seq_along(df.l)) {
  df.l[[i]] <- left_join(df.l[[i]], tax, by = "Feature_ID")
}

# clean columns
for (i in seq_along(df.l)) {
  df.l[[i]] <- select(df.l[[i]], 
                      Feature_ID,
                      ASV_id_new,
                      Taxon,
                      A,
                      B,
                      stat,
                      p.value,
                      Group)
}

# save tables
for (i in seq_along(df.l)) {
  write.table(df.l[[i]], 
              file = paste0("Results/16S/indicator_species/formatted_output_with_new_ids/", 
                            names(df.l[i])),
              sep = "\t",
              row.names = F)
}




