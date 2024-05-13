# 16S taxonomic profile

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(viridis)
library(scales)
library(RColorBrewer)


## import and format files ----

# data
asv <- readRDS("Data/16S/formatted_data/asv_filtered.rds")
tax <- readRDS("Data/16S/formatted_data/taxonomy_filtered.rds")
meta <- read.table("Data/16S/formatted_data/metadata_files/metadata.tsv", 
                   header = T,
                   sep = "\t",
                   strip.white = T)

# format taxonomy table
tax <- separate(tax, Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")

# replace NAs with Unknown
tax[is.na(tax)] <- "Unknown"

# join asv with tax data
tax_asv <-  inner_join(tax, asv, by = "Feature_ID")

## create taxonomic profile ----

# collapse by tax level of interest
phylum <- dplyr::select(
  tax_asv, !c(Feature_ID, Kingdom, Class, Order, Family, Genus, Species, ASV_id_new)) %>% 
  group_by(Phylum) %>%
  summarise_all(list(sum)) %>% 
  as.data.frame()

family <- dplyr::select(
  tax_asv, !c(Feature_ID, Kingdom, Phylum, Class, Order, Genus, Species, ASV_id_new)) %>% 
  group_by(Family) %>%
  summarise_all(list(sum)) %>% 
  as.data.frame()

## transpose and join with metadata
# clean tax strings
phylum$Phylum <- gsub(" p__", "", phylum$Phylum)
family$Family <- gsub("f__", "", family$Family)

# get row names
rownames(phylum) <- phylum$Phylum
rownames(family) <- family$Family

# remove tax column and transpose
phylum <- t(phylum[,-1]) %>% as.data.frame()
family <- t(family[,-1]) %>% as.data.frame()

# change to relative abundance
phylum <- phylum / rowSums(phylum) * 100
family <- family / rowSums(family) * 100

# attach metadata
phylum$Sample_id <- rownames(phylum) 
phylum <- inner_join(meta, phylum, by = "Sample_id")

family$Sample_id <- rownames(family)
family <- inner_join(meta, family, by = "Sample_id")

# remove samples (blanks/ctrls) with no seqs
# 1 easp ctrl removed, 2 blanks
to_rm <- phylum$Sample_id[phylum$non.chimeric == 0]
phylum <- filter(phylum, !phylum$Sample_id %in% to_rm)
family <- filter(family, !family$Sample_id %in% to_rm)

# filter out rubble (not using for this data for this MS)
phylum <- filter(phylum, !phylum$Treatment %in% c("Rubble", "Rubble_ctrl"))
family <- filter(family, !family$Treatment %in% c("Rubble", "Rubble_ctrl"))

# filter out blanks, larave, water (for no controls dataset)
phylum_nb <- filter(phylum, !phylum$Tank %in% c("Blank", "Larvae",
                                             "room_water", "tank_water"))

family_nb <- filter(family, !family$Tank %in% c("Blank", "Larvae",
                                             "room_water", "tank_water"))

# create df list
tax_profile.l <- list(phylum, phylum_nb, family, family_nb)
names(tax_profile.l) <- c("phylum", "phylum_nb", "family", "family_nb")

# remove taxa no longer occuring in df after filtering
lapply(tax_profile.l, dim)
tax_profile.l <- lapply(tax_profile.l, function(x) {
  x <- x[, colSums(x != 0) > 0]
  return(x)
})
lapply(tax_profile.l, dim)

# change to long format
# first remove unncessary columns
tax_profile.l <- lapply(tax_profile.l, function(x) {
  x <- dplyr::select(x, !c("Sample_id", "extraction_ID",
                           "input", "filtered",
                           "percentage.of.input.passed.filter",
                           "denoised", "merged",
                           "percentage.of.input.merged",
                           "non.chimeric",
                           "percentage.of.input.non.chimeric"))
  return(x)
})

# pivot longer
for (i in seq_along(tax_profile.l)) {
  tax_profile.l[[i]] <- pivot_longer(
    tax_profile.l[[i]], !c(Sample,  
                           Coral, 
                           Treatment, 
                           Tank, 
                           Settlement,
                           Percent_settled,
                           Settlement_2),
    names_to = names(tax_profile.l)[i], 
    values_to = "Relative Abundance")
}

# group by taxonony and summarise
for (i in seq_along(tax_profile.l)) {
  tax_profile.l[[i]] <- group_by(tax_profile.l[[i]],
                                 names(tax_profile.l[i]))
}

# reorder treatment
for (i in seq_along(tax_profile.l)) {
  tax_profile.l[[i]]$Treatment <- factor(tax_profile.l[[i]]$Treatment, 
                                         levels = c("Tile_ctrl", "Light_1M", 
                                                    "Light_2M", "Dark_2M",
                                                    "Rubble", "Larvae", 
                                                    "tank_water", "room_water", 
                                                    "Rubble_ctrl", "Blank"))
}

# reorder taxa
tax_profile.l$phylum <- mutate(tax_profile.l$phylum, 
                               Phylum = fct_relevel(phylum, 
                                                    "Unknown", after = Inf))
tax_profile.l$phylum_nb <- mutate(tax_profile.l$phylum_nb, 
                               phylum_nb = fct_relevel(phylum_nb, 
                                                       "Unknown", after = Inf))
tax_profile.l$family <- mutate(tax_profile.l$family, 
                                  family = fct_relevel(family, 
                                                       "Unknown", after = Inf))
tax_profile.l$family_nb <- mutate(tax_profile.l$family_nb, 
                                  family_nb = fct_relevel(family_nb, 
                                                          "Unknown", after = Inf))

## plot relative abundance ----

## plot as heatmap
# get range
range(tax_profile.l$phylum$`Relative Abundance`)
r_val.p <- c(0, 20, 100)

range(tax_profile.l$phylum_nb$`Relative Abundance`)
r_val.pnb <- c(0, 20, 98)

range(tax_profile.l$family$`Relative Abundance`)
r_val.f <- c(0, 10, 79)

range(tax_profile.l$family_nb$`Relative Abundance`)
r_val.fnb <- c(0, 10, 79)


## plot relative abundance
# get colours
heat_col <- viridis(n = 9, 
                    alpha = 0.9, 
                    begin = 0, 
                    end = 0.95, 
                    option = "turbo", 
                    direction = 1)
show_col(heat_col)

# get facet labels
#all
facet_names.all <- c(
  `Light_1M` = "Light 1Month",
  `Light_2M` = "Light 2Month",
  `Dark_2M` = "Dark 2Month",
  `Larvae` = "Lv",
  `tank_water` = "TW",
  `room_water` = "RW",
  `Tile_ctrl` = "TC",
  `Rubble_ctrl` = "RC",
  `Blank` = "Bl"
)

# for 'facet wrap'
facet_names.allfw <- c(
  `Light_1M` = "Light 1Month",
  `Light_2M` = "Light 2Month",
  `Dark_2M` = "Dark 2Month",
  `Larvae` = "Larvae",
  `tank_water` = "Tank Water",
  `room_water` = "Room Water",
  `Tile_ctrl` = "Tile Control",
  `Rubble_ctrl` = "Rubble Control",
  `Blank` = "Blank"
)

#nb
facet_names.nb <- c(
  `Light_1M` = "Light 1Month",
  `Light_2M` = "Light 2Month",
  `Dark_2M` = "Dark 2Month"
)

## plot heatmap

# phylum
p <- ggplot(tax_profile.l$phylum, aes(x = Sample, y = phylum, 
                                      fill = `Relative Abundance`)) +
  geom_raster() +
  scale_fill_gradientn(name = "Relative Abundance", 
                       colours = heat_col,
                       values = rescale(r_val.p)) + 
  scale_y_discrete(position = "right") +
  xlab("Biofilm Sample") +
  ylab("Microbial Phylum") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "left",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)) 
# facet wrap
p + facet_wrap("Treatment", 
               nrow = 1, 
               ncol = 10, 
               scales = "free_x",
               labeller = as_labeller(facet_names.allfw)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", linewidth = 1))

# save
ggsave("16S_rel_abun_phylum_heat_all.svg", device = "svg", width = 30, 
       height = 20, dpi = 300, units = "cm")


# phylum no-blanks
p <- ggplot(tax_profile.l$phylum_nb, aes(x = Sample, y = phylum_nb, 
                                         fill = `Relative Abundance`)) +
  geom_raster() +
  scale_fill_gradientn(name = "Relative Abundance", 
                       colours = heat_col,
                       values = rescale(r_val.pnb)) + 
  xlab("Biofilm Sample") +
  ylab("Microbial Phylum") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)) 
# facet grid
p + facet_grid(~Treatment, 
               scales = 'free_x', 
               space = 'free', 
               labeller = as_labeller(facet_names.all)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", linewidth = 1))

# save
ggsave("16S_rel_abun_phylum_heat_nb.svg", device = "svg", width = 25, 
       height = 20, dpi = 300, units = "cm")


# family
p <- ggplot(tax_profile.l$family, aes(x = Sample, y = family, 
                                      fill = `Relative Abundance`)) +
  geom_raster() +
  scale_fill_gradientn(name = "Relative Abundance", 
                       colours = heat_col,
                       values = rescale(r_val.f)) + 
  xlab("Biofilm Sample") +
  ylab("Microbial Family") +
  theme(axis.text = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)) 
# facet wrap
p + facet_wrap("Treatment", 
               nrow = 1, 
               ncol = 10, 
               scales = "free_x",
               labeller = as_labeller(facet_names.allfw)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", linewidth = 1))

# save
ggsave("16S_rel_abun_family_heat_all.svg", device = "svg", width = 30, 
       height = 18, dpi = 300, units = "cm")

# family no-blanks
p <- ggplot(tax_profile.l$family_nb, aes(x = Sample, y = family_nb, 
                                      fill = `Relative Abundance`)) +
  geom_raster() +
  scale_fill_gradientn(name = "Relative Abundance", 
                       colours = heat_col,
                       values = rescale(r_val.fnb)) + 
  xlab("Biofilm Sample") +
  ylab("Microbial Family") +
  theme(axis.text = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)) 
# facet grid
p + facet_grid(~Treatment, 
               scales = 'free_x', 
               space = 'free', 
               labeller = as_labeller(facet_names.all)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", linewidth = 1))

# save
ggsave("16S_rel_abun_family_heat_nb.svg", device = "svg", width = 25, 
       height = 18, dpi = 300, units = "cm")

## summary statistics ----

# Function get mean and SE of for each group for a particular variable
# take a dataframe, grouping variables and response variable as input
# can take 1 or 2 groups
sum_stats <- function(df, group1, group2 = NULL, variable) {
  # format  columns (convery string to symbol)
  group1_col <- sym(group1)
  variable_col <- sym(variable)
  
  # if no group2 (NULL)
  if (is.null(group2)) {
    new_df <- df %>%
      group_by(!!group1_col) %>% 
      summarise(mean = mean(!!variable_col, na.rm = TRUE), 
                stdev = sd(!!variable_col, na.rm = TRUE), 
                sterr = sd(!!variable_col, na.rm = TRUE)/sqrt(n()),
                min = min(!!variable_col, na.rm = TRUE),
                max = max(!!variable_col, na.rm = TRUE))
  } else {
    group2_col <- sym(group2)
  # if group2   
    new_df <- df %>%
      group_by(!!group1_col, !!group2_col) %>% 
      summarise(mean = mean(!!variable_col, na.rm = TRUE), 
                stdev = sd(!!variable_col, na.rm = TRUE), 
                sterr = sd(!!variable_col, na.rm = TRUE)/sqrt(n()),
                min = min(!!variable_col, na.rm = TRUE),
                max = max(!!variable_col, na.rm = TRUE))
  }
  
  return(new_df)
}

# by treatment and phylum
summary_treatment_p <- sum_stats(tax_profile.l$phylum_nb, 
                                 "Treatment", 
                                 "phylum_nb", 
                                 variable = "Relative Abundance")

# by phylum
summary_phylum <- sum_stats(tax_profile.l$phylum_nb, 
                               "phylum_nb", 
                               variable = "Relative Abundance")

# by treatment and family
summary_treatment_f <- sum_stats(tax_profile.l$family_nb, 
                                 "Treatment", 
                                 "family_nb", 
                                 variable = "Relative Abundance")

# by family
summary_family <- sum_stats(tax_profile.l$family_nb, 
                            "family_nb", 
                            variable = "Relative Abundance")


## create top 20 profile ----

# collapse by tax level of interest
phylum <- dplyr::select(
  tax_asv, !c(Feature_ID, Kingdom, Class, Order, Family, Genus, Species, ASV_id_new)) %>% 
  group_by(Phylum) %>%
  summarise_all(list(sum)) %>% 
  as.data.frame()

family <- dplyr::select(
  tax_asv, !c(Feature_ID, Kingdom, Phylum, Class, Order, Genus, Species, ASV_id_new)) %>% 
  group_by(Family) %>%
  summarise_all(list(sum)) %>% 
  as.data.frame()

# subset to top 20
# get tax totals
phylum$total <- rowSums(phylum[,-1])
family$total <- rowSums(family[,-1])

# re-order rows by total
phylum <- phylum[order(-phylum[,which(colnames(phylum) == 'total')]),]
family <- family[order(-family[,which(colnames(family) == 'total')]),]  

# group families outside top 20 into 'other' category (using 19 to include other as 20)
Other <- colSums(phylum[-c(1:19), -1]) 
phylum <- rbind(c("Other", Other), phylum)

Other <- colSums(family[-c(1:19), -1])
family <- rbind(c("Other", Other), family)

# remove total column 
phylum <- dplyr::select(phylum, !total)
family <- dplyr::select(family, !total)

# subset to top 20 (including other)
phylum <- phylum[1:20,]
family <- family[1:20,]

# change to numeric (rbind makes it character)
tmp <- as.data.frame(sapply(phylum[,-1], as.numeric))
Phylum <- phylum$Phylum
phylum <- cbind(Phylum, tmp)

tmp <- as.data.frame(sapply(family[,-1], as.numeric))
Family <- family$Family
family <- cbind(Family, tmp)

## transpose and join with metadata
# clean tax strings
phylum$Phylum <- gsub(" p__", "", phylum$Phylum)
family$Family <- gsub("f__", "", family$Family)

# get row names
rownames(phylum) <- phylum$Phylum
rownames(family) <- family$Family

# remove tax column and transpose
phylum <- t(phylum[,-1]) %>% as.data.frame()
family <- t(family[,-1]) %>% as.data.frame()

# change to relative abundance
phylum <- phylum / rowSums(phylum) * 100
family <- family / rowSums(family) * 100

# attach metadata
phylum$Sample_id <- rownames(phylum) 
phylum <- inner_join(meta, phylum, by = "Sample_id")

family$Sample_id <- rownames(family)
family <- inner_join(meta, family, by = "Sample_id")

# remove samples (blanks/ctrls) with no seqs
to_rm <- phylum$Sample_id[phylum$non.chimeric == 0]
phylum <- filter(phylum, !phylum$Sample_id %in% to_rm)
family <- filter(family, !family$Sample_id %in% to_rm)

# filter out rubble (not using for this data for this MS)
phylum <- filter(phylum, !phylum$Treatment %in% c("Rubble", "Rubble_ctrl"))
family <- filter(family, !family$Treatment %in% c("Rubble", "Rubble_ctrl"))

# filter out blanks, controls, larave, water (for no controls dataset)
phylum_nb <- filter(phylum, !phylum$Tank %in% c("Blank", "Larvae", 
                                                "room_water", "tank_water"))

family_nb <- filter(family, !family$Tank %in% c("Blank", "Larvae", 
                                                "room_water", "tank_water"))

# create df list
tax_profile.l <- list(phylum, phylum_nb, family, family_nb)
names(tax_profile.l) <- c("phylum", "phylum_nb", "family", "family_nb")

# change to long format
# first remove unncessary columns
tax_profile.l <- lapply(tax_profile.l, function(x) {
  x <- dplyr::select(x, !c("Sample_id", "extraction_ID",
                           "input", "filtered",
                           "percentage.of.input.passed.filter",
                           "denoised", "merged",
                           "percentage.of.input.merged",
                           "non.chimeric",
                           "percentage.of.input.non.chimeric"))
  return(x)
})

# pivot longer
for (i in seq_along(tax_profile.l)) {
  tax_profile.l[[i]] <- pivot_longer(
    tax_profile.l[[i]], !c(Sample,  
                           Coral, 
                           Treatment, 
                           Tank, 
                           Settlement,
                           Percent_settled,
                           Settlement_2),
    names_to = names(tax_profile.l)[i], 
    values_to = "Relative Abundance")
}

# group by taxonomy and summarise
for (i in seq_along(tax_profile.l)) {
  tax_profile.l[[i]] <- group_by(tax_profile.l[[i]],
                                 names(tax_profile.l[i]))
}

# reorder treatment
for (i in seq_along(tax_profile.l)) {
  tax_profile.l[[i]]$Treatment <- factor(tax_profile.l[[i]]$Treatment, 
                                         levels = c("Tile_ctrl", "Dark_2M", 
                                                    "Light_1M", "Light_2M", 
                                                    "Rubble", "Larvae", 
                                                    "tank_water", "room_water", 
                                                    "Rubble_ctrl", "Blank"))
}

# reorder taxa (from most abundant to least)
# function
reorder_taxa <- function(df, column_name) {
  
  # calculate total relative abundance for each phylum
  tmp <- df %>%
    group_by(!!sym(column_name)) %>%
    summarise(total_abundance = sum(`Relative Abundance`))
  
  # reorder phyla based on total abundance
  tmp <- tmp %>%
    arrange(desc(total_abundance)) %>%
    pull(!!sym(column_name))
  
  # move "Other" and "Unknown" phyla to the beginning
  tmp <- c("Other", "Unknown", tmp[tmp != "Other" & tmp != "Unknown"])
  
  # apply the new order to the phylum column
  df <- df %>%
    mutate(!!sym(column_name) := fct_relevel(!!sym(column_name), tmp))
  
  return(df)
}

# apply function to dfs
tax_profile.l$phylum <- reorder_taxa(tax_profile.l$phylum, "phylum")
tax_profile.l$phylum_nb <- reorder_taxa(tax_profile.l$phylum_nb, "phylum_nb")
tax_profile.l$family <- reorder_taxa(tax_profile.l$family, "family")
tax_profile.l$family_nb <- reorder_taxa(tax_profile.l$family_nb, "family_nb")

## plot top 20 ----

# colour brewer palette
cols <- colorRampPalette(brewer.pal(12, "Paired"))(20)
show_col(cols)

# get facet labels
#all
facet_names.all <- c(
  `Light_1M` = "Light 1Month",
  `Light_2M` = "Light 2Month",
  `Dark_2M` = "Dark 2Month",
  `Larvae` = "Lv",
  `tank_water` = "TW",
  `room_water` = "RW",
  `Tile_ctrl` = "TC",
  `Blank` = "Bl"
)

# for 'facet wrap'
facet_names.allfw <- c(
  `Light_1M` = "Light 1Month",
  `Light_2M` = "Light 2Month",
  `Dark_2M` = "Dark 2Month",
  `Rubble` = "Rubble",
  `Larvae` = "Larvae",
  `tank_water` = "Tank Water",
  `room_water` = "Plate Water",
  `Tile_ctrl` = "Tile Control",
  `Blank` = "Blank"
)

#nb
facet_names.nb <- c(
  `Light_1M` = "Light 1Month",
  `Light_2M` = "Light 2Month",
  `Dark_2M` = "Dark 2Month",
  `Tile_ctrl` = "Ctrl"
)

# phylum
p <- ggplot(tax_profile.l$phylum, aes(fill=phylum, colour=phylum, 
                                      y=`Relative Abundance`, 
                                      x=Sample)) + 
  geom_bar(position="fill", stat="identity", width=1) +
  scale_colour_manual(values = cols, name = "Phylum") +
  scale_fill_manual(values = cols, name = "Phylum") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position='bottom') 
# facet wrap
p + facet_wrap("Treatment", 
               nrow = 1, 
               ncol = 10, 
               scales = "free_x",
               labeller = as_labeller(facet_names.allfw)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", linewidth = 1))

# ggsave
ggsave("16S_rel_abun_phylum_bar_top20_all.svg", 
       device = "svg", width = 23, height = 16, 
       dpi = 300, units = "cm")

#phylum no blanks
p <- ggplot(tax_profile.l$phylum_nb, aes(fill=phylum_nb, colour=phylum_nb, 
                                      y=`Relative Abundance`, 
                                      x=Sample)) + 
  geom_bar(position="fill", stat="identity", width=1) +
  scale_colour_manual(values = cols, name = "Phylum") +
  scale_fill_manual(values = cols, name = "Phylum") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position='bottom') 
# facet
p + facet_grid(~Treatment, 
               scales = 'free_x', 
               space = 'free', 
               labeller = as_labeller(facet_names.nb)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        strip.text = element_text(size = 12))

# ggsave
ggsave("16S_rel_abun_phylum_bar_top20_nb.svg", 
       device = "svg", width = 25, height = 16, 
       dpi = 300, units = "cm")


#family
p <- ggplot(tax_profile.l$family, aes(fill=family, colour=family, 
                                         y=`Relative Abundance`, 
                                         x=Sample)) + 
  geom_bar(position="fill", stat="identity", width=1) +
  scale_colour_manual(values = cols, name = "Family") +
  scale_fill_manual(values = cols, name = "Family") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'bottom') 

# facet wrap
p + facet_wrap("Treatment", 
               nrow = 1, 
               ncol = 10, 
               scales = "free_x",
               labeller = as_labeller(facet_names.allfw)) + 
  theme(panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", linewidth = 1))

# ggsave
ggsave("16S_rel_abun_family_bar_top20_all.svg", 
       device = "svg", width = 23, height = 16, 
       dpi = 300, units = "cm")

#family nb
p <- ggplot(tax_profile.l$family_nb, aes(fill=family_nb, colour=family_nb, 
                        y=`Relative Abundance`, 
                        x=Sample)) + 
  geom_bar(position="fill", stat="identity", width=1) +
  scale_colour_manual(values = cols, name = "Family") +
  scale_fill_manual(values = cols, name = "Family") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'bottom') 

# facet grid
p + facet_grid(~Treatment, 
                   scales = 'free_x', 
                   space = 'free', 
                   labeller = as_labeller(facet_names.nb)) + 
  theme(panel.spacing = unit(0.1, "lines"))

# ggsave
ggsave("16S_rel_abun_family_bar_top20_nb.svg", 
       device = "svg", width = 23, height = 16, 
       dpi = 300, units = "cm")



