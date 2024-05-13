# Get alpha diversity of 16S data

setwd("~/Documents/R/Microbial_inducers")

library(tidyverse)
library(vegan)
library(viridis)
library(scales)
library(car)
library(multcomp)
library(ggpubr)

## import and format data ----

#asv <- readRDS("Data/16S/formatted_data/asv_filtered.rds")
asv_r <- readRDS("Data/16S/formatted_data/asv_rarefied.rds")
tax <- readRDS("Data/16S/formatted_data/taxonomy_filtered.rds")
meta <- read.table("Data/16S/formatted_data/metadata_files/metadata.tsv", 
                   header = T,
                   sep = "\t",
                   strip.white = T)

## calculate asv richness and diversity ----

# note: better to use rarefied data for richness

# get richness 
asv.rich <- apply(asv_r > 0, 1, sum) %>% 
  as.data.frame %>%
  dplyr::rename(richness = 1)

# calculate shannons
asv.shan <- diversity(asv_r, index = "shannon") %>% 
  as.data.frame %>%
  dplyr::rename(shannon = 1) 

# get sample_id column
asv.rich$Sample_id <- rownames(asv_r)
asv.shan$Sample_id <- rownames(asv_r)

# join with metadata
asv.div <- right_join(meta, asv.rich, by = "Sample_id")
asv.div <- right_join(asv.div, asv.shan, by = "Sample_id")

# remove rubble
asv.div <- filter(asv.div, !asv.div$Treatment %in% c("Rubble", "Rubble_ctrl"))

# remove controls (note: blanks and larvae removed during rarefaction due to low number of sequences)
# no blanks, water or rubble
asv.div.nb <- filter(asv.div, !asv.div$Tank %in% c(
  "Rubble", "Blank", "room_water", "tank_water"
))

# reorder treatment
#asv
asv.div$Treatment <- factor(asv.div$Treatment, levels = c(
  "Tile_ctrl", "Dark_2M", "Light_1M", "Light_2M", "tank_water", "room_water"
))

#nb
asv.div.nb$Treatment <- factor(asv.div.nb$Treatment, levels = c(
  "Tile_ctrl", "Dark_2M", "Light_1M", "Light_2M"
))

## plot richness and diversity ----

# get cols
clrs <- viridis(4, direction = 1, option = "turbo")
show_col(clrs)

clrs <- clrs[c(1:3)]
clrs <- c("#C0C0C0", clrs) # adding grey for controls 

# cols 2
clrs2 <- c(clrs, "#3E9BFEFF", "#7A0403FF")
show_col(clrs2)

# get labels
x.labs <- c( "Control", "Dark 2M", "Light 1M", "Light 2M", "Tank Water", "Plate Water")
x.labs.nb <- c("Control", "Dark 2M", "Light 1M", "Light 2M")

# plot richness
# all
p <- ggplot(asv.div, aes(x=Treatment, y=richness, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylab("ASV richness") +
  scale_fill_manual(values = clrs2, labels = x.labs) +
  #scale_fill_viridis(option = "turbo", direction = 1, 
  #                   discrete = T, labels = x.labs) +
  scale_x_discrete(labels = x.labs) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=11, angle = 45, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p

# to add significance values (code below)
p + geom_text(data=asv.rich.groups, 
              aes(x = Treatment, y = Ymax+100, label = letters),
              vjust=0, size = 4)
# save
ggsave("asv_richness_rarefied.svg", device = "svg", width = 18, height = 15, 
       units = "cm", dpi = 300)

# no blanks
p <- ggplot(asv.div.nb, aes(x=Treatment, y=richness, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylab("ASV richness") +
  scale_fill_manual(values = clrs, labels = x.labs.nb) +
  scale_x_discrete(labels = x.labs.nb) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p
# to add significance values (code below)
p1 <- p + geom_text(data=asv.rich.nb.groups, 
              aes(x = Treatment, y = Ymax+75, label = letters),
              vjust=0, size = 5)
p1

# save
ggsave("asv_nb_richness_rarefied_w_ctrl.svg", device = "svg", width = 15, height = 12, 
       units = "cm", dpi = 300)


# plot shannon
# all
p <- ggplot(asv.div, aes(x=Treatment, y=shannon, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylim(0, NA) +
  ylab("Shannon Diversity Index") +
  scale_fill_manual(values = clrs2, labels = x.labs) +
  #scale_fill_viridis(option = "turbo", direction = 1, 
  #                   discrete = T, labels = x.labs) +
  scale_x_discrete(labels = x.labs) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=11, angle = 45, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p

# to add significance values (code below)
p + geom_text(data=asv.shan.groups, 
              aes(x = Treatment, y = Ymax+0.25, label = letters),
              vjust=0, size = 4)
# save
ggsave("asv_shannon_rarefied_ylim.svg", device = "svg", width = 18, height = 15, 
       units = "cm", dpi = 300)

# asv.nb
p <- ggplot(asv.div.nb, aes(x=Treatment, y=shannon, fill=Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  ylim(0, NA) +
  ylab("Shannon Diversity Index") +
  scale_fill_manual(values = clrs, labels = x.labs.nb) +
  scale_x_discrete(labels = x.labs.nb) +
  guides(fill = guide_legend(title = "Treatment")) +
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 13),
        legend.position = "none") 
p
# to add significance values (code below)
p2 <- p + geom_text(data=asv.shan.nb.groups, 
              aes(x = Treatment, y = Ymax+0.25, label = letters),
              vjust=0, size = 5)
p2

# save
ggsave("asv_nb_shannon_rarefied_ylim_w_ctrl.svg", device = "svg", width = 15, height = 12, 
       units = "cm", dpi = 300)


# make into multi panel figure
ggarrange(p1 +
            theme(axis.text=element_text(size=13),
                  axis.title.y=element_text(size = 15)), 
          p2 +
            theme(axis.text=element_text(size=13),
                  axis.title.y=element_text(size = 15)),
          widths = c(1, 0.98))

ggsave("asv_nb_richness_and_shannon.svg", device = "svg", width = 25, height = 10, 
       units = "cm", dpi = 300)


## calculate significance values ----

# all treatments
# Transform data
asv.div$rich_fourth <- (asv.div$richness ^ (1/4))
asv.div$shan_fourth <- (asv.div$shannon ^ (1/4))

# Check assumptions 
#richness
residualPlots(lm(asv.div$richness~asv.div$Treatment))
residualPlots(lm(asv.div$rich_fourth~asv.div$Treatment)) # better to use transformed data as less variance between treatments

#shannon
residualPlots(lm(asv.div$shannon~asv.div$Treatment))
residualPlots(lm(asv.div$shan_fourth~asv.div$Treatment)) # better to use transformed data as less variance between treatments

# get lm
#richness
asv.rich.lm <- lm(rich_fourth ~ Treatment, asv.div)
plot(asv.rich.lm)

summary(asv.rich.lm) # summarise linear model. Note: treatment names are retained when factor levels are not re-ordered
anova(asv.rich.lm) # for anova

#shannon
asv.shan.lm <- lm(shan_fourth ~ Treatment, asv.div)
plot(asv.shan.lm)

summary(asv.shan.lm)
anova(asv.shan.lm) 

#post hoc test
# richness
glht.asv.rich <- summary(glht(asv.rich.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
asv.rich.groups <- cld(glht.asv.rich)
asv.rich.groups <- fortify(asv.rich.groups)
colnames(asv.rich.groups) <- c("Treatment", "letters")

ymax <- tapply(asv.div$richness, asv.div$Treatment, max)
asv.rich.groups$Ymax <- ymax # add to plot above

# shannon
glht.asv.shan <- summary(glht(asv.shan.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
asv.shan.groups <- cld(glht.asv.shan)
asv.shan.groups <- fortify(asv.shan.groups)
colnames(asv.shan.groups) <- c("Treatment", "letters")

ymax <- tapply(asv.div$shannon, asv.div$Treatment, max)
asv.shan.groups$Ymax <- ymax # add to plot above


# treatments no blanks/controls
# transform data
asv.div.nb$rich_fourth <- (asv.div.nb$richness ^ (1/4))
asv.div.nb$shan_fourth <- (asv.div.nb$shannon ^ (1/4))

# Check assumptions 
#richness
residualPlots(lm(asv.div.nb$richness~asv.div.nb$Treatment))
residualPlots(lm(asv.div.nb$rich_fourth~asv.div.nb$Treatment)) # better to use transformed data as less variance between treatments

#shannon
residualPlots(lm(asv.div.nb$shannon~asv.div.nb$Treatment))
residualPlots(lm(asv.div.nb$shan_fourth~asv.div.nb$Treatment))

# get lm
#richness
asv.rich.nb.lm <- lm(rich_fourth ~ Treatment, asv.div.nb)
plot(asv.rich.nb.lm)

summary(asv.rich.nb.lm) # summarise linear model. Note: treatment names are retained when factor levels are not re-ordered
anova(asv.rich.nb.lm) # for anova

#shannon
asv.shan.nb.lm <- lm(shan_fourth ~ Treatment, asv.div.nb)
plot(asv.shan.nb.lm)

summary(asv.shan.nb.lm)
anova(asv.shan.nb.lm) 

#post hoc test
# richness
glht.asv.rich.nb <- summary(glht(asv.rich.nb.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
asv.rich.nb.groups <- cld(glht.asv.rich.nb)
asv.rich.nb.groups <- fortify(asv.rich.nb.groups)
colnames(asv.rich.nb.groups) <- c("Treatment", "letters")

ymax <- tapply(asv.div.nb$richness, asv.div.nb$Treatment, max)
asv.rich.nb.groups$Ymax <- ymax # add to plot above

# shannon
glht.asv.shan.nb <- summary(glht(asv.shan.nb.lm, linfct = mcp(Treatment = "Tukey")), test = adjusted("bonferroni"))
asv.shan.nb.groups <- cld(glht.asv.shan.nb)
asv.shan.nb.groups <- fortify(asv.shan.nb.groups)
colnames(asv.shan.nb.groups) <- c("Treatment", "letters")

ymax <- tapply(asv.div.nb$shannon, asv.div.nb$Treatment, max)
asv.shan.nb.groups$Ymax <- ymax # add to plot above


## get group averages ----

# function for std error
std <- function(x) sd(x)/sqrt(length(x))

# richness
rich.mean <- asv.div %>%
  group_by(Treatment) %>%
  dplyr::summarise(Mean = mean(richness, na.rm=TRUE),
                   SD = sd(richness, na.rm=TRUE),
                   SE = std(richness),
                   min = min(richness, na.rm = TRUE),
                   max = max(richness, na.rm = TRUE))

write.table(rich.mean, file = "summary_stats_richness.tsv", sep = "\t", row.names = F)

# shannon
shan.mean <- asv.div %>%
  group_by(Treatment) %>%
  dplyr::summarise(Mean = mean(shannon, na.rm=TRUE),
                   SD = sd(shannon, na.rm=TRUE),
                   SE = std(shannon),
                   min = min(shannon, na.rm = TRUE),
                   max = max(shannon, na.rm = TRUE))

write.table(shan.mean, file = "summary_stats_shannon.tsv", sep = "\t", row.names = F)




