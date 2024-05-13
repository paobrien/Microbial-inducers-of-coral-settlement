# Create venn diagram to see shared ASV between groups

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(ggvenn)

## import data ----

## asv, tax, meta
asv <- readRDS("Data/16S/formatted_data/asv_filtered.rds")
tax <- readRDS("Data/16S/formatted_data/taxonomy_filtered.rds")
meta <- read.table("Data/16S/formatted_data/metadata_files/metadata.tsv", 
                   header = T,
                   sep = "\t",
                   strip.white = T)


## combined analysis files

# get file path and files
fp <- "Results/16S/combined_analyses/"
fl <- list.files(fp)[c(1,2,5,7)]

# import files
df.l <- list()
for (i in seq_along(fl)) {
  df.l[[i]] <- read.table(paste0(fp, fl[i]), 
                          header = T, 
                          sep = "\t")
  names(df.l)[i] <- gsub("_combined.tsv", "", fl[i])
}
head(df.l$dfav)

## format data ----

# transform counts to relative abundances
# function
rel_abun <- function(df) {
  Feature_ID <- df$Feature_ID # save feature IDS
  df <- t(df[,-1]) # transpose - does not work with colsums??
  df <- df / rowSums(df) * 100 # convert to rel abund
  df <- t(df) %>% as.data.frame # transpose back
  df$Feature_ID <- Feature_ID # add feature ids
  return(df)
}

# apply
asv <- rel_abun(asv)

# check transformation
colSums(asv[,1:10])

# subset to l2m samples
# meta
meta <- filter(
  meta, meta$Treatment == "Light_2M"
)

# asvs
asv <- dplyr::select(asv, meta$Sample_id, Feature_ID)

# remove asvs no longer in df
rownames(asv) <- asv$Feature_ID
asv <- dplyr::select(asv, -c("Feature_ID"))
asv <- asv[rowSums(asv != 0) > 0,]

# remove uneccessary columns from combined df
df_s.l <- lapply(df.l, function (x) {
  dplyr::select(x, -c(Taxon, A, B, isolate_id, percent_identity))
})

# identify duplicated rows
dups <- lapply(df_s.l, duplicated)

# remove duplicated rows
for (i in seq_along(df_s.l)) {
  df_s.l[[i]] <- df_s.l[[i]][!dups[[i]],]
}

# remove unneeded meta columns
meta <- dplyr::select(
  meta, Sample_id, Coral, Treatment, Tank, Settlement_2, Percent_settled
)

# taxonomy
tax <- separate(
  tax, Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
  sep = ";", remove = FALSE, convert = FALSE, extra = "warn", fill = "warn"
)

# remove __ and characters preceeding
tax <- tax %>%
  mutate_at(vars(-Taxon), ~ gsub(".{1}__", "", .))

# remove trailing white spaces 
tax <- tax %>%
  mutate_all(str_trim)

# change NAs to unclassified
tax <- tax %>%
  mutate_all(~ifelse(is.na(.), "Unclassified", .))

# add label columns
tax <- unite(tax, "label", c(Phylum, Family), sep = "; ", remove = F)


## get asvs to plot ----

# subset indicators function
subset_indicators <- function(df, group, lm, coral) {
  if (lm == ">") {
    filtered_df <- df %>%
      filter(Group == group | lm_coef > 0)
    
    n_asv_ids <- nrow(filtered_df) 
    
    filtered_df %>%
      dplyr::select(asv_id) %>%
      mutate(coral = rep(coral, n_asv_ids),
             indicator = rep(group, n_asv_ids)) %>%
      return()
    
  } else if (lm == "<") {
    filtered_df <- df %>%
      filter(Group == group | lm_coef < 0)
    
    n_asv_ids <- nrow(filtered_df) 
    
    filtered_df %>%
      dplyr::select(asv_id) %>%
      mutate(coral = rep(coral, n_asv_ids),
             indicator = rep(group, n_asv_ids)) %>%
      return()
  } else {
    stop("lm argument must be either '>' or '<'")
  }
}

## get high inidicators
asvs_h.l <- list()
for(i in seq_along(df_s.l)) {
  asvs_h.l[[i]] <- subset_indicators(
    df = df_s.l[[i]], 
    group = "High", 
    lm = ">", 
    coral = names(df_s.l)[i]
  )
  names(asvs_h.l)[i] <- names(df_s.l)[i]
} 

## get low inidicators
asvs_l.l <- list()
for(i in seq_along(df_s.l)) {
  asvs_l.l[[i]] <- subset_indicators(
    df = df_s.l[[i]], 
    group = "Low", 
    lm = "<", 
    coral = names(df_s.l)[i]
  )
  names(asvs_l.l)[i] <- names(df_s.l)[i]
} 

# total
lapply(asvs_h.l, nrow)
lapply(asvs_l.l, nrow)


## plot venn diagram ----

# high asvs
ven_h.l <- list(
  asvs_h.l$dfav$asv_id,
  asvs_h.l$easp$asv_id,
  asvs_h.l$plob$asv_id,
  asvs_h.l$psin$asv_id
)
names(ven_h.l) <- c("dfav", "easp", "plob", "psin")

# compute the intersections
p <- ggvenn(ven_h.l, 
            stroke_size = 0.5, 
            set_name_size = 6)
p

ggsave("venn_high_ASVs.svg", device = "svg", width = 15, height = 15, units = "cm", dpi = 300)


# low asvs
ven_l.l <- list(
  asvs_l.l$dfav$asv_id,
  asvs_l.l$easp$asv_id,
  asvs_l.l$plob$asv_id,
  asvs_l.l$psin$asv_id
)
names(ven_l.l) <- c("dfav", "easp", "plob", "psin")

# compute the intersections
p <- ggvenn(ven_l.l, 
       stroke_size = 0.5, 
       set_name_size = 6)
p

ggsave("venn_low_ASVs.svg", device = "svg", width = 15, height = 15, units = "cm", dpi = 300)

# high v low
high_all <- c(ven_h.l$dfav,
              ven_h.l$easp,
              ven_h.l$plob,
              ven_h.l$psin)

low_all <- c(ven_l.l$dfav,
             ven_l.l$easp,
             ven_l.l$plob,
             ven_l.l$psin)

ven_hl.l <- list("High" = high_all, 
                 "Low" = low_all)

# compute the intersections
p <- ggvenn(ven_hl.l, 
       stroke_size = 0.5, 
       set_name_size = 6)
p

ggsave("venn_high_low_ASVs.svg", device = "svg", width = 15, height = 15, units = "cm", dpi = 300)


## identify which ASVs are shared between corals
intersect(ven_h.l$easp, ven_h.l$plob)
intersect(ven_l.l$plob, ven_l.l$dfav)

# between high v low
intersect(ven_hl.l$High, ven_hl.l$Low)

lapply(ven_h.l, function(x) {
  intersect("asv_167", x)
})







