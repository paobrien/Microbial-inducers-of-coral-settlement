# Beta diversity analysis

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(vegan)
library(viridis)
library(scales)
library(ggpubr)

## import files ----

# data
asv <- readRDS("Data/16S/formatted_data/asv_filtered.rds")
tax <- readRDS("Data/16S/formatted_data/taxonomy_filtered.rds")
meta <- read.table("Data/16S/formatted_data/metadata_files/metadata.tsv", 
                   header = T,
                   sep = "\t",
                   strip.white = T)

## create filtered datasets for testing ----

# make df
rownames(asv) <- asv$Feature_ID
asv <- t(asv[,-1]) %>% as.data.frame
asv$Sample_id <- rownames(asv)

# all data
df_full <- full_join(meta, asv, by = "Sample_id")
row.names(df_full) <- df_full$Sample

# remove poorly sequenced samples (due to low number of reads)
df_full <- filter(df_full, !df_full$Sample %in% c("SE3704", "SE3707"))

# remove rubble treatment
df_full <- filter(df_full, !df_full$Treatment %in% c("Rubble", "Rubble_ctrl"))

# remove asv no longer in dataframe after subsetting
df_full <- df_full[, colSums(df_full != 0) > 0] 

# remove NA's from merge
df_full[is.na(df_full)] <- 0

# remove samples (blanks/ctrls) with no seqs
to_rm <- df_full$Sample[df_full$non.chimeric == 0]
df_full <- filter(df_full, !df_full$Sample %in% to_rm)

# no blanks, larave or water
df_nblw <- filter(df_full, !df_full$Tank %in% c("Blank", "Larvae",
                                                 "room_water", "tank_water"))

# only 2M light
df_2ML <- filter(df_full, df_full$Treatment %in% "Light_2M")
df_2ML <- df_2ML[, colSums(df_2ML != 0) > 0]

# make df list
df.l <- list(
  "full" = df_full, 
  "no_blanks_larvae_water" = df_nblw,
  "L2M" = df_2ML
)

## Calculate and plot NMDS ----

# function for nmds plot
plot_nmds <- function(df, distance, trymax, autotransform, k, plot, metadata) {
  # run nmds
  dist.nmds <- metaMDS(df, 
                       distance = distance, 
                       trymax = trymax, 
                       autotransform = autotransform, 
                       k = k, 
                       plot = plot)
  # Extract scores from metaMDS and return as data matrix
  nmds.scores <- scores(dist.nmds)
  nmds.scores.site <- as.data.frame(nmds.scores$sites)
  # add metadata
  nmds.scores.site$Sample <- row.names(nmds.scores.site)
  nmds.df <- left_join(nmds.scores.site, metadata, by = "Sample")
  return(nmds.df)
}

# creates matrices
mat.l <- list()
for (i in seq_along(df.l)) {
  mat.l[[i]] <- select(df.l[[i]], !colnames(meta))
}
names(mat.l) <- names(df.l)

# run function over matrix list
nmds.df.l <- list() 
for (i in seq_along(mat.l)) {
  nmds.df.l[[i]] <- plot_nmds(df = mat.l[[i]], 
                              distance = "bray",
                              trymax = 500,
                              autotransform = T,
                              k = 2,
                              plot = T,
                              metadata = meta)
}
names(nmds.df.l) <- names(df.l)

# save nmds
saveRDS(nmds.df.l, "Results/16S/beta_diversity/nmds_list.rds")

## plot nmds
# reorder treatment
# full
nmds.df.l$full$Treatment <- factor(nmds.df.l$full$Treatment, levels = c(
  "Tile_ctrl", "Dark_2M", "Light_1M", "Light_2M", "tank_water", "room_water", "Larvae", "Blank"
))

# nblw
nmds.df.l$no_blanks_larvae_water$Treatment <- factor(nmds.df.l$no_blanks_larvae_water$Treatment, levels = c(
  "Tile_ctrl", "Dark_2M", "Light_1M", "Light_2M"
))

# get cols (for light dark - to match previous plots)
# nblw
clrs <- viridis(4, direction = 1, option = "turbo")
clrs <- clrs[c(1:3)]
clrs <- c("#C0C0C0", clrs)
show_col(clrs)

# full
clrs2 <- c(clrs, "#3E9BFEFF", "#7A0403FF", "#1BD0D5FF", "#D2E935FF")
show_col(clrs2)

# get labels
labs_full <- c("Control", "Dark 2M", "Light 1M", "Light 2M", "Tank Water", "Plate Water", "Larvae", "Blank")
labs_nblw <- c("Control", "Dark 2M", "Light 1M", "Light 2M")
labs_l2m <- c("Tank 1", "Tank 2", "Tank 3")

# nblw - shape = tank
p <- ggplot(nmds.df.l$no_blanks_larvae_water, aes(x=NMDS1, 
                         y=NMDS2, 
                         fill = Treatment,
                         shape = Tank)) +
  geom_point(size = 4, alpha = 0.7) + 
  scale_fill_manual(values = clrs,
                    labels = labs_tmt.tc) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = labs_tank.tc) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)) 
p

ggsave("16S_nmds_light_dark_ctrl_shape_corder.svg", device = "svg", width = 19, height = 15, 
       dpi = 300, units = "cm")

# nblw - colour = settlement
p <- ggplot(nmds.df.l$no_blanks_larvae_water, aes(x=NMDS1, 
                                                  y=NMDS2, 
                                                  colour = Percent_settled,
                                                  fill = Percent_settled,
                                                  shape = Treatment)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = labs_nblw) +
  scale_colour_viridis(discrete = F, option = "cividis", direction = 1) +
  scale_fill_viridis(discrete = F, option = "cividis", direction = 1) +
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 15),
        legend.position = "left") +
  labs(colour = "Percent Settled",
       fill = "Percent Settled",
       shape = "Treatment") +
  guides(colour = FALSE, fill = FALSE) # to remove percent settled legend for combined plot below 
p

ggsave("16S_nmds_light_dark_ctrl_colour_border.svg", device = "svg", width = 20, height = 15, 
       dpi = 300, units = "cm")

# full - colour = treatment
p <- ggplot(nmds.df.l$full, aes(x=NMDS1, 
                         y=NMDS2, 
                         fill = Treatment)) +
  geom_point(size = 4, alpha = 0.7, shape = 21) +
  scale_fill_manual(values = clrs2,
                      labels = labs_full) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)) 
p

ggsave("16S_nmds_full_colour_treatment.svg", device = "svg", width = 20, height = 15, 
       dpi = 300, units = "cm")


# l2m only - settlement + tank
p2 <- ggplot(nmds.df.l$L2M, aes(x=NMDS1, 
                         y=NMDS2, 
                         colour = Percent_settled,
                         fill = Percent_settled,
                         shape = Tank)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22, 24),
                     labels = labs_l2m) +
  scale_colour_viridis(discrete = F, option = "cividis", direction = 1) +
  scale_fill_viridis(discrete = F, option = "cividis", direction = 1) +
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 15),
        legend.position = "right") +
  labs(colour = "Percent Settled",
       fill = "Percent Settled",
       shape = "Tank")
p2

# add conditioning month
nmds.df.l$L2M$Month <- ifelse(nmds.df.l$L2M$Coral %in% c("Psin", "Dfav", "Easp"), "Oct", "Nov")

# 2ml only - settlement + conditioing month
p3 <- ggplot(nmds.df.l$L2M, aes(x=NMDS1, 
                                  y=NMDS2, 
                                  colour = Percent_settled,
                                  fill = Percent_settled,
                                  shape = Tank,
                                  linetype = Month)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(size = 0.4) +
  scale_shape_manual(values = c(21, 22, 24),
                     labels = labs_l2m) +
  scale_colour_viridis(discrete = F, option = "cividis", direction = 1) +
  scale_fill_viridis(discrete = F, option = "cividis", direction = 1) +
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 15),
        legend.position = "right") +
  labs(colour = "Percent Settled",
       fill = "Percent Settled",
       shape = "Tank")
p3

ggsave("16S_nmds_conditioning_month.svg", device = "svg", width = 20, height = 15, 
       dpi = 300, units = "cm")

# combine all treatments with l2m plot
ggarrange(p + 
            ggtitle("A) All Conditioning Treatments") + 
            theme(plot.title = element_text(size = 16)), 
          p2 + 
            ggtitle("B) Light 2-Month Conditioned Biofilms") +
            theme(plot.title = element_text(size = 16)),
          widths = c(0.92, 1))

ggsave("treatment_and_l2m_nmds_comb_2.svg", device = "svg", width = 32, height = 13, 
       dpi = 300, units = "cm")

## Calculate multivariate statistics using PERMANOVA ----

# note: only necessary to do it on treatments excluding blanks, larvae and sw

## create distance matrix 
# get matrix
mat.l.perm <- list(
  "no_blanks_larvae_water" = mat.l$no_blanks_larvae_water, 
  "L2M" = mat.l$L2M 
)

# standardise matrix
for (i in seq_along(mat.l.perm)) {
  mat.l.perm[[i]] <- wisconsin(mat.l.perm[[i]])
}

# get distance
bray_dist.l <- list()
for (i in seq_along(mat.l.perm)) {
  bray_dist.l[[i]] <- vegdist(mat.l.perm[[i]])
}
names(bray_dist.l) <- names(mat.l.perm)

# join df with metadata
df.perm.l <- list()
for (i in seq_along(mat.l.perm)) {
  mat.l.perm[[i]] <- as.data.frame(mat.l.perm[[i]])
  mat.l.perm[[i]]$Sample <- row.names(mat.l.perm[[i]])
  df.perm.l[[i]] <- right_join(meta, mat.l.perm[[i]], by = "Sample")
}
names(df.perm.l) <- names(mat.l.perm)

## Run PERMANOVAs

# all treatments and corals
perm.treatment <- adonis2(
  bray_dist.l$no_blanks_larvae_water~Treatment*Settlement_2*Tank*Coral,
  df.perm.l$no_blanks_larvae_water,
  permutations = 9999
)
perm.treatment

# l2m only
perm.l2m <- adonis2(
  bray_dist.l$L2M~Tank*Settlement_2*Coral,
  df.perm.l$L2M,
  permutations = 9999
)
perm.l2m

# combine to list
permanova_results.l <- list(
  perm.treatment, perm.l2m
)
names(permanova_results.l) <- c("all", "light2M")

# save results
saveRDS(permanova_results.l, "Results/16S/beta_diversity/permanova_results_list.rds")





