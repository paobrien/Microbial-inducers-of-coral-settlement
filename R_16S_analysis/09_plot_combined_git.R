## Plot combined results from indval, maaslin and random forests

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(viridis)
library(scales)
library(ggpubr)

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

# add feature id
# high
asvs_h.l <- lapply(asvs_h.l, function(x) {
  df <- left_join(
    x, dplyr::select(tax, Feature_ID, ASV_id_new),
    join_by(asv_id == ASV_id_new)
  )
  return(df)
})

# low
asvs_l.l <- lapply(asvs_l.l, function(x) {
  df <- left_join(
    x, dplyr::select(tax, Feature_ID, ASV_id_new),
    join_by(asv_id == ASV_id_new)
  )
  return(df)
})

# subset asv df to high/low asvs
asv$Feature_ID <- row.names(asv) # first fet feature id

# high
high_df.l <- lapply(asvs_h.l, function(x) {
  filter(asv, Feature_ID %in% x$Feature_ID)
})

# low
low_df.l <- lapply(asvs_l.l, function(x) {
  filter(asv, Feature_ID %in% x$Feature_ID)
})

# group by family/phylum
# first add family/phylum to df
high_df.l <- lapply(high_df.l, function(x) {
  x <- right_join(dplyr::select(tax, c(Feature_ID, label, Phylum)), 
                  x, by = "Feature_ID")
})

low_df.l <- lapply(low_df.l, function(x) {
  x <- right_join(dplyr::select(tax, c(Feature_ID, label, Phylum)), 
                  x, by = "Feature_ID")
})

# group by family
fam_h.l <- lapply(high_df.l, function(x) {
  x <- x %>%
    dplyr::select(!c(Feature_ID, Phylum)) %>%
    group_by(label) %>%
    summarise_all(list(sum)) %>% 
    as.data.frame()
    
})

fam_l.l <- lapply(low_df.l, function(x) {
  x <- x %>%
    dplyr::select(!c(Feature_ID, Phylum)) %>%
    group_by(label) %>%
    summarise_all(list(sum)) %>% 
    as.data.frame()
  
})

# group by phylum
phy_h.l <- lapply(high_df.l, function(x) {
  x <- x %>%
    dplyr::select(!c(Feature_ID, label)) %>%
    group_by(Phylum) %>%
    summarise_all(list(sum)) %>% 
    as.data.frame()
  
})

phy_l.l <- lapply(low_df.l, function(x) {
  x <- x %>%
    dplyr::select(!c(Feature_ID, label)) %>%
    group_by(Phylum) %>%
    summarise_all(list(sum)) %>% 
    as.data.frame()
  
})


# transpose so samples are rows
# family
fam_h.l <- lapply(fam_h.l, function(x) {
  rownames(x) <- x$label
  dplyr::select(x, -label) %>% 
    t %>% as.data.frame
})

fam_l.l <- lapply(fam_l.l, function(x) {
  rownames(x) <- x$label
  dplyr::select(x, -label) %>% 
    t %>% as.data.frame
})

# phylum
phy_h.l <- lapply(phy_h.l, function(x) {
  rownames(x) <- x$Phylum
  dplyr::select(x, -Phylum) %>% 
    t %>% as.data.frame
})

phy_l.l <- lapply(phy_l.l, function(x) {
  rownames(x) <- x$Phylum
  dplyr::select(x, -Phylum) %>% 
    t %>% as.data.frame
})

# add sample column
#family
fam_h.l <- lapply(fam_h.l, function(x) {
  x$Sample_id <- rownames(x)
  return(x)
})

fam_l.l <- lapply(fam_l.l, function(x) {
  x$Sample_id <- rownames(x)
  return(x)
})

#phylum
phy_h.l <- lapply(phy_h.l, function(x) {
  x$Sample_id <- rownames(x)
  return(x)
})

phy_l.l <- lapply(phy_l.l, function(x) {
  x$Sample_id <- rownames(x)
  return(x)
})

# add meta
#family
fam_h.l <- lapply(fam_h.l, function(x) {
  df <- right_join(meta, x, by = "Sample_id")
  return(df)
})

fam_l.l <- lapply(fam_l.l, function(x) {
  df <- right_join(meta, x, by = "Sample_id")
  return(df)
})

#phylum
phy_h.l <- lapply(phy_h.l, function(x) {
  df <- right_join(meta, x, by = "Sample_id")
  return(df)
})

phy_l.l <- lapply(phy_l.l, function(x) {
  df <- right_join(meta, x, by = "Sample_id")
  return(df)
})

# filter samples by coral

# dfav
fam_h.l$dfav <- filter(
  fam_h.l$dfav, Coral == "Dfav"
)

fam_l.l$dfav <- filter(
  fam_l.l$dfav, Coral == "Dfav"
)

phy_h.l$dfav <- filter(
  phy_h.l$dfav, Coral == "Dfav"
)

phy_l.l$dfav <- filter(
  phy_l.l$dfav, Coral == "Dfav"
)

# easp
fam_h.l$easp <- filter(
  fam_h.l$easp, Coral == "Easp"
)

fam_l.l$easp <- filter(
  fam_l.l$easp, Coral == "Easp"
)

phy_h.l$easp <- filter(
  phy_h.l$easp, Coral == "Easp"
)

phy_l.l$easp <- filter(
  phy_l.l$easp, Coral == "Easp"
)

# plob
fam_h.l$plob <- filter(
  fam_h.l$plob, Coral == "Plob"
)

fam_l.l$plob <- filter(
  fam_l.l$plob, Coral == "Plob"
)

phy_h.l$plob <- filter(
  phy_h.l$plob, Coral == "Plob"
)

phy_l.l$plob <- filter(
  phy_l.l$plob, Coral == "Plob"
)

# psin
fam_h.l$psin <- filter(
  fam_h.l$psin, Coral == "Psin"
)

fam_l.l$psin <- filter(
  fam_l.l$psin, Coral == "Psin"
)

phy_h.l$psin <- filter(
  phy_h.l$psin, Coral == "Psin"
)

phy_l.l$psin <- filter(
  phy_l.l$psin, Coral == "Psin"
)

# relevel samples by settlement then abundance
# first get abundnace
for (i in seq_along(fam_h.l)) {
  fam_h.l[[i]]$Total <- rowSums(fam_h.l[[i]] %>%
                                    dplyr::select(!c(Sample_id,
                                                     Coral,
                                                     Treatment,
                                                     Tank,
                                                     Settlement_2,
                                                     Percent_settled
                                    )))
}

# relevel using function (orders by percent settled then total abundance)
reorder_sample_id <- function(x, percent_settled, total) {
  # create a custom order based on Percent_settled and Total
  sample_order <- order(-percent_settled, -total)
  # reorder factor levels
  factor(x, levels = x[sample_order])
}

# apply function to reorder Sample_id factor in each dataframe
fam_h.l <- lapply(fam_h.l, function(df) {
  df$Sample_id <- reorder_sample_id(df$Sample_id, df$Percent_settled, df$Total)
  return(df)
})

# relevel remaining dfs to same order
# fam low
for (i in seq_along(fam_l.l)) {
  fam_l.l[[i]]$Sample_id <- factor(
    fam_l.l[[i]]$Sample_id, levels = fam_l.l[[i]]$Sample_id[
      order(fam_h.l[[i]]$Sample_id, decreasing = FALSE)
    ]
  )
}

# phy high
for (i in seq_along(phy_h.l)) {
  phy_h.l[[i]]$Sample_id <- factor(
    phy_h.l[[i]]$Sample_id, levels = phy_h.l[[i]]$Sample_id[
      order(fam_h.l[[i]]$Sample_id, decreasing = FALSE)
    ]
  )
}

# phy low
for (i in seq_along(phy_l.l)) {
  phy_l.l[[i]]$Sample_id <- factor(
    phy_l.l[[i]]$Sample_id, levels = phy_l.l[[i]]$Sample_id[
      order(fam_h.l[[i]]$Sample_id, decreasing = FALSE)
    ]
  )
}

# change to long format
# family
fam_h.l <- lapply(fam_h.l, function(x) {
  x <- dplyr::select(x, !Total) # check this works if running again
  df <- pivot_longer(
    x,
    cols = names(x)[-c(1:6)], # otherwise need to remove here
    names_to = "Family",
    values_to = "Relative Abundance"
  )
  return(df)
})

fam_l.l <- lapply(fam_l.l, function(x) {
  df <- pivot_longer(
    x,
    cols = names(x)[-c(1:6)],
    names_to = "Family",
    values_to = "Relative Abundance"
  )
  return(df)
})

# phylum
phy_h.l <- lapply(phy_h.l, function(x) {
  df <- pivot_longer(
    x,
    cols = names(x)[-c(1:6)], 
    names_to = "Phylum",
    values_to = "Relative Abundance"
  )
  return(df)
})

phy_l.l <- lapply(phy_l.l, function(x) {
  df <- pivot_longer(
    x,
    cols = names(x)[-c(1:6)],
    names_to = "Phylum",
    values_to = "Relative Abundance"
  )
  return(df)
})


# reorder family/phylum from highest to lowest rel abund
# first get total
# family
totals_fh <- lapply(fam_h.l, function (x) {
  aggregate(
    `Relative Abundance` ~ Family, x, sum
  )
})

totals_fl <- lapply(fam_l.l, function (x) {
  aggregate(
    `Relative Abundance` ~ Family, x, sum
  )
})

#phylum
totals_ph <- lapply(phy_h.l, function (x) {
  aggregate(
    `Relative Abundance` ~ Phylum, x, sum
  )
})

totals_pl <- lapply(phy_l.l, function (x) {
  aggregate(
    `Relative Abundance` ~ Phylum, x, sum
  )
})


# reorder levels of family/phylum based on the total values
# family
for (i in seq_along(fam_h.l)) {
  fam_h.l[[i]]$Family <- factor(
    fam_h.l[[i]]$Family, levels = totals_fh[[i]]$Family[
      order(totals_fh[[i]]$`Relative Abundance`, decreasing = TRUE)
    ]
  )
}

for (i in seq_along(fam_l.l)) {
  fam_l.l[[i]]$Family <- factor(
    fam_l.l[[i]]$Family, levels = totals_fl[[i]]$Family[
      order(totals_fl[[i]]$`Relative Abundance`, decreasing = TRUE)
    ]
  )
}

# phylum
for (i in seq_along(phy_h.l)) {
  phy_h.l[[i]]$Phylum <- factor(
    phy_h.l[[i]]$Phylum, levels = totals_ph[[i]]$Phylum[
      order(totals_ph[[i]]$`Relative Abundance`, decreasing = TRUE)
    ]
  )
}

for (i in seq_along(phy_l.l)) {
  phy_l.l[[i]]$Phylum <- factor(
    phy_l.l[[i]]$Phylum, levels = totals_pl[[i]]$Phylum[
      order(totals_pl[[i]]$`Relative Abundance`, decreasing = TRUE)
    ]
  )
}


# relevel settlement groups
# family
for (i in seq_along(fam_h.l)) {
  fam_h.l[[i]]$Settlement_2 <- factor(
    fam_h.l[[i]]$Settlement_2, levels = c("High", "Med", "Low")
  )
}

for (i in seq_along(fam_l.l)) {
  fam_l.l[[i]]$Settlement_2 <- factor(
    fam_l.l[[i]]$Settlement_2, levels = c("High", "Med", "Low")
  )
}

# family
for (i in seq_along(phy_h.l)) {
  phy_h.l[[i]]$Settlement_2 <- factor(
    phy_h.l[[i]]$Settlement_2, levels = c("High", "Med", "Low")
  )
}

for (i in seq_along(phy_l.l)) {
  phy_l.l[[i]]$Settlement_2 <- factor(
    phy_l.l[[i]]$Settlement_2, levels = c("High", "Med", "Low")
  )
}

## plot bar graph ----

# family high
plots_fh <- list()
for (i in seq_along(fam_h.l)) {
  plots_fh[[i]] <- ggplot(
    fam_h.l[[i]], aes(fill=Family, 
                        colour=Family, 
                        y=`Relative Abundance`, 
                        x=Sample_id)) + 
    geom_bar(stat="identity", width=1) +
    scale_color_viridis(option = "C", discrete = T, direction = 1) +
    scale_fill_viridis(option = "C", discrete = T, direction = 1) +
    labs(fill = "High Settlement (Inducers)", 
         colour = "High Settlement (Inducers)") +
    theme_bw() +
    theme(#axis.text.x = element_text(angle = 90),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 13),
      legend.position='right',
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 13)) +
    facet_grid(~Settlement_2, 
               scales = 'free',
               space = 'free')
  names(plots_fh)[i] <- names(fam_h.l)[i]
}

#view
plots_fh$psin

# save
for (i in seq_along(plots_fh)) {
  ggsave(paste0(names(plots_fh)[i], "_family_high_indicators_taxa.svg"), 
         plot = plots_fh[[i]],
         device = "svg",
         width = 25,
         height = 15,
         units = "cm",
         dpi = 300)
}

# family low
plots_fl <- list()
for (i in seq_along(fam_l.l)) {
  plots_fl[[i]] <- ggplot(
    fam_l.l[[i]], aes(fill=Family, 
                      colour=Family, 
                      y=`Relative Abundance`, 
                      x=Sample_id)) + 
    geom_bar(stat="identity", width=1) +
    scale_color_viridis(option = "E", discrete = T, direction = 1) +
    scale_fill_viridis(option = "E", discrete = T, direction = 1) +
    scale_y_continuous(trans = "reverse") +
    labs(fill = "Low Settlement (Inhibitors)", 
         colour = "Low Settlement (Inhibitors)",
         x = "Biofilm Sample") +
    guides(fill = guide_legend(reverse = TRUE),
           colour = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(#axis.text.x = element_text(angle = 90),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 13),
      legend.position='right',
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      #strip.text = element_text(size = 13),
      strip.text = element_blank()) +
    facet_grid(~Settlement_2, 
               scales = 'free',
               space = 'free')
  names(plots_fl)[i] <- names(fam_l.l)[i]
}

#view
plots_fl$psin

# save
for (i in seq_along(plots_fl)) {
  ggsave(paste0(names(plots_fl)[i], "_family_low_indicators_taxa.svg"), 
         plot = plots_fl[[i]],
         device = "svg",
         width = 25,
         height = 15,
         units = "cm",
         dpi = 300)
}

# phylum high
plots_ph <- list()
for (i in seq_along(phy_h.l)) {
  plots_ph[[i]] <- ggplot(
    phy_h.l[[i]], aes(fill=Phylum, 
                      colour=Phylum, 
                      y=`Relative Abundance`, 
                      x=Sample_id)) + 
    geom_bar(stat="identity", width=1) +
    scale_color_viridis(option = "C", discrete = T, direction = 1) +
    scale_fill_viridis(option = "C", discrete = T, direction = 1) +
    labs(fill = "High Settlement (Inducers)", 
         colour = "High Settlement (Inducers)") +
    theme_bw() +
    theme(#axis.text.x = element_text(angle = 90),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 13),
      legend.position='right',
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 13)) +
    facet_grid(~Settlement_2, 
               scales = 'free',
               space = 'free')
  names(plots_ph)[i] <- names(phy_h.l)[i]
}

#view
plots_ph$plob

# save
for (i in seq_along(plots_ph)) {
  ggsave(paste0(names(plots_ph)[i], "_phylum_high_indicators_taxa.svg"), 
         plot = plots_ph[[i]],
         device = "svg",
         width = 25,
         height = 15,
         units = "cm",
         dpi = 300)
}

# phylum low
plots_pl <- list()
for (i in seq_along(phy_l.l)) {
  plots_pl[[i]] <- ggplot(
    phy_l.l[[i]], aes(fill=Phylum, 
                      colour=Phylum, 
                      y=`Relative Abundance`, 
                      x=Sample_id)) + 
    geom_bar(stat="identity", width=1) +
    scale_color_viridis(option = "E", discrete = T, direction = 1) +
    scale_fill_viridis(option = "E", discrete = T, direction = 1) +
    scale_y_continuous(trans = "reverse") +
    labs(fill = "Low Settlement (Inhibitors)", 
         colour = "Low Settlement (Inhibitors)") +
    guides(fill = guide_legend(reverse = TRUE),
           colour = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(#axis.text.x = element_text(angle = 90),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 13),
      legend.position='right',
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      #strip.text = element_text(size = 13),
      strip.text = element_blank()) +
    facet_grid(~Settlement_2, 
               scales = 'free',
               space = 'free')
  names(plots_pl)[i] <- names(phy_l.l)[i]
}

#view
plots_pl$plob

# save
for (i in seq_along(plots_pl)) {
  ggsave(paste0(names(plots_pl)[i], "_phylum_low_indicators_taxa.svg"), 
         plot = plots_pl[[i]],
         device = "svg",
         width = 25,
         height = 15,
         units = "cm",
         dpi = 300)
}


## plot multipanel ----

# psin
p1 <- plots_fh$psin
p2 <- plots_fl$psin

ggarrange(p1 + 
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(5, "mm"),
                  axis.title.y = element_text(size = 14)),
          p2 + 
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(5, "mm"),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(size = 14)),
          ncol = 1, align = "v")

ggsave("psin_all_ordered_abun.svg", device = "svg", 
       width = 22, height = 20, 
       units = "cm", dpi = 300)

# dfav
p1 <- plots_fh$dfav
p2 <- plots_fl$dfav

ggarrange(p1 + 
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(5, "mm"),
                  axis.title.y = element_text(size = 14)),
          p2 + 
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(5, "mm"),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(size = 14)),
          ncol = 1, align = "v")

ggsave("dfav_all_ordered_abun.svg", device = "svg", 
       width = 22, height = 20, 
       units = "cm", dpi = 300)

# easp
p1 <- plots_fh$easp
p2 <- plots_fl$easp

ggarrange(p1 + 
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(5, "mm"),
                  axis.title.y = element_text(size = 14)),
          p2 + 
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(5, "mm"),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(size = 14)),
          ncol = 1, align = "v")

ggsave("easp_all_ordered_abun.svg", device = "svg", 
       width = 22, height = 20, 
       units = "cm", dpi = 300)

# plot
p1 <- plots_ph$plob
p2 <- plots_fl$plob

ggarrange(p1 + 
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(5, "mm"),
                  axis.title.y = element_text(size = 14)),
          p2 + 
            theme(legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(5, "mm"),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(size = 14)),
          ncol = 1, align = "v")

ggsave("plob_all_ordered_abun.svg", device = "svg", 
       width = 22, height = 20, 
       units = "cm", dpi = 300)


## get indicator ASV metrics ----

# check numbers of features
# per species
lapply(high_df.l, function(x) {
  length(unique(x$Feature_ID))
})

lapply(low_df.l, function(x) {
  length(unique(x$Feature_ID))
})

# check number of unique phyla
# per species
lapply(high_df.l, function(x) {
  print(unique(x$Phylum))
})

lapply(low_df.l, function(x) {
  print(unique(x$Phylum))
})

# total
phyla_unique <- unique(unlist(lapply(high_df.l, function(df) df$Phylum)))
length(phyla_unique)

phyla_unique <- unique(unlist(lapply(low_df.l, function(df) df$Phylum)))
length(phyla_unique)

# check number of unique families
# per species
lapply(high_df.l, function(x) {
  print(unique(x$label))
})

lapply(low_df.l, function(x) {
  print(unique(x$label))
})

# total
family_unique <- unique(unlist(lapply(high_df.l, function(df) df$label)))
length(family_unique)

family_unique <- unique(unlist(lapply(low_df.l, function(df) df$label)))
length(family_unique)

# get family names
fam_df.l <- lapply(high_df.l, function(x) {
  tmp <- unique(x$label)
  return(tmp)
})

# combine vectors into dataframe
fam_df <- do.call(rbind, lapply(names(fam_df.l), function(coral_name) {
  data.frame(family = fam_df.l[[coral_name]], coral = coral_name)
}))
fam_df






