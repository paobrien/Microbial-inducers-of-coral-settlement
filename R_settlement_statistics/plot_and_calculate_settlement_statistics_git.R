## Summary statistics for settlement data

setwd("~/Documents/R/Microbial_inducers") 

library(tidyverse)
library(effects)
library(lme4)
library(emmeans)
library(multcomp)
library(ggpubr)


## load and format data ----

# Settlement dataframes
df <- readRDS("Data/Settlement_data/R_files/settlement_data_df_list") # biofilm assays
df_chem <- readRDS("Data/Settlement_data/R_files/settlement_data_chem_lists") # chem assays

# subset to relevant df

# main exp
df.l <- list(
  "lcor" = df$nov_l_corymbosa, 
  "plob" = df$nov_p_lobata, 
  "dfav" = df$oct_d_favus, 
  "easp" = df$oct_e_aspera, 
  "psin" = df$oct_p_sinensis
)

# check dataframes
lapply(df.l, head)

# rename cols to keep consistent
for(i in seq_along(df.l)) {
  colnames(df.l[[i]])[which(names(df.l[[i]]) == "Not_settled")] <- "Total_not_settled"
}

for(i in seq_along(df_chem)) {
  colnames(df_chem[[i]])[which(names(df_chem[[i]]) == "Not_settled")] <- "Total_not_settled"
}

#change percent settled to proportion (0-1)
for(i in seq_along(df.l)) {
  df.l[[i]]$Percent_settled <- df.l[[i]]$Percent_settled/100
}

for(i in seq_along(df_chem)) {
  df_chem[[i]]$Percent_settled <- df_chem[[i]]$Percent_settled/100
}

## Fit GLMM with random effects (main exp) ----

# fit a binomial model with tank and observation-level random effects
df.glmer <- list()
for(i in seq_along(df.l)) {
  df.glmer[[i]] <- glmer(cbind(Total_settled, Total_not_settled) ~ Treatment + (1 | Tank/Sample_ID), 
                         data = df.l[[i]], 
                         family = binomial)
  names(df.glmer)[i] <- names(df.l)[i]
}

## explore diagnostics
par(mfrow = c(1, 2))

# plot the random effects against the predicted values from the fixed effect component of the model and check for no trend:
random_effects <- function(model) {
  Xmat <- model.matrix(model)
  fit <- Xmat %*% fixef(model)
  ran <- ranef(model, drop = T)$Sample_ID
  return(plot(fit, ran, pch = 19, las = 1, cex = 1.4))
}

random_effects(df.glmer$easp)
abline(0, 0, lty = 1)

# also check for approximate normality of random effects:
qqnorm(ran, pch = 19, las = 1, cex = 1.4)
qqline(ran)

# check overdispersion
sum(resid(df.glmer$psin, type = "pearson")^2) / (nrow(df.l$psin) - length(coef(df.glmer$psin))) 
library(RVAideMemoire) #GLMM overdispersion test
overdisp.glmer(md1)


## explore model results
lapply(df.glmer, summary)

for(i in seq_along(df.glmer)) {
  print(names(df.glmer)[i])
  print(exp(fixef(df.glmer[[i]])))
}

## save fixed effects table

# extract table from model
fixef_table <- list()
for(i in seq_along(df.glmer)) {
  fixef_table[[i]] <- coef(summary(df.glmer[[i]]))
  names(fixef_table)[i] <- names(df.glmer)[i]
} 
fixef_table

# remove treatment from row names
for(i in seq_along(fixef_table)) {
  row.names(fixef_table[[i]]) <- gsub("Treatment", "", row.names(fixef_table[[i]]))
}

# write table
for(i in seq_along(fixef_table)) {
  write.csv(fixef_table[[i]], paste0(names(fixef_table)[i],"_fixed_effects.csv"))
}


## Fitting and predicting ----

# Plot - means +/- SE ----
plot_mean <- function(df, model) {
  # first extract mean and SE
  newdata <- data.frame(Treatment = levels(df$Treatment))  #creates a df with every treatment combination to predict
  Xmat <- model.matrix(~Treatment, data = newdata) # build model matrix
  coefs <- fixef(model)
  fit <- as.vector(coefs %*% t(Xmat))
  se <- sqrt(diag(Xmat %*% vcov(model) %*% t(Xmat)))
  Q <- qt(0.975, df = nrow(model@frame) - length(coefs) - 2)
  newdata <- cbind(newdata, fit = binomial()$linkinv(fit), 
                   lower = binomial()$linkinv(fit - Q * se), #95% CI lower
                   upper = binomial()$linkinv(fit + Q * se)) #95% CI upper
  # now to plot
  p <- ggplot(newdata, aes(y = fit, x = Treatment)) + 
    geom_blank() + 
    geom_pointrange(aes(ymin = lower, ymax = upper, x = as.numeric(Treatment))) + 
    scale_y_continuous("Proportion of larvae settled") +
    scale_x_discrete("Treatment") + 
    theme_classic() + 
    theme(axis.line.x = element_line(), 
          axis.line.y = element_line(), 
          axis.title.x = element_text(margin = margin(t = 2, unit = "lines")), 
          axis.title.y = element_text(margin = margin(r = 2, unit = "lines")))
  return(p)
}

plot_mean(df.l$lcor, df.glmer$lcor)


## Post-hoc test ----

# Once a significant factor effect has been established (i.e., treatment effect above), 
# explore which treatments differ using posthoc tets

# to compare only against the control 
# dunnet's test where control is the name of your control treatment
df.dun <- list()
for(i in seq_along(df.glmer)) {
  df.dun[[i]] <- contrast(emmeans(df.glmer[[i]],  ~Treatment), "trt.vs.ctrl1", ref = 'Control')
  names(df.dun)[i] <- names(df.glmer)[i]
}
df.dun
  
#Compare all combinations - Tukeys test
df.glht <- list()
for(i in seq_along(df.glmer)) {
  df.glht[[i]] <- summary(glht(df.glmer[[i]], mcp(Treatment="Tukey")), test = adjusted("bonferroni"))
  names(df.glht)[i] <- names(df.glmer)[i]
}


# get significance letters for plot
df.groups <- list()
ymax <- list()
for(i in seq_along(df.glht)) {
  #extract letters
  df.groups[[i]] <- cld(df.glht[[i]])
  # make columns
  df.groups[[i]] <- fortify(df.groups[[i]])
  colnames(df.groups[[i]]) <- c("Treatment", "letters")
  # plot coordinates
  ymax[[i]] <- tapply(df[[i]]$Percent_settled, df[[i]]$Treatment, max)
  df.groups[[i]]$Ymax <- ymax[[i]] 
  names(df.groups)[i] <- names(df.glht)[i]
}
df.groups

## Plot boxplots main experiment ----

# y = Ymax+0.1 for just above boxplot

# get colours
# get cols
clrs <- viridis(4, direction = 1, option = "turbo")
show_col(clrs)

clrs <- c(
  "#C0C0C0",
  "#7A0403FF",
  "#30123BFF",
  "#1AE4B6FF",
  "#FABA39FF"
)

# legend labels
l_labs <- c("Control", "Dark 1M", "Dark 2M", "Light 1M", "Light 2M")

#lcor
p_lc <- ggplot(df$lcor, aes(x = Treatment, y = Percent_settled, fill = Treatment)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.15, shape = 21) +
  scale_fill_manual(values = clrs, aesthetics = "fill", labels = l_labs) +
  scale_x_discrete(labels = c("Control \n n=72", 
                              "Dark 1M \n n=90", 
                              "Dark 2M \n n=90", 
                              "Light 1M \n n=90", 
                              "Light 2M \n n=90")) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0, 1.1)) +
  ylab("Proportion Settled") +
  theme_bw() + 
  geom_text(data=df.groups$lcor, aes(x = Treatment, y = 1.075, label = letters),vjust=0) # for significance letters
p_lc

ggsave("l_cor_settlement_significance.svg", device = "svg", width = 8, height = 6, units = "in", dpi = 300)

#plob
p_pl <- ggplot(df$plob, aes(x = Treatment, y = Percent_settled, fill = Treatment)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.15, shape = 21) +
  scale_fill_manual(values = clrs, aesthetics = "fill", labels = l_labs) +
  scale_x_discrete(labels = c("Control \n n=72", 
                              "Dark 1M \n n=90", 
                              "Dark 2M \n n=90", 
                              "Light 1M \n n=90", 
                              "Light 2M \n n=90")) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0, 1.1)) +
  ylab("Proportion Settled") +
  theme_bw() + 
  geom_text(data=df.groups$plob, aes(x = Treatment, y = 1.075, label = letters),vjust=0) # for significance letters
p_pl

ggsave("p_lob_settlement_significance.svg", device = "svg", width = 8, height = 6, units = "in", dpi = 300)

#dfav
p_df <- ggplot(df$dfav, aes(x = Treatment, y = Percent_settled, fill = Treatment)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.15, shape = 21) +
  scale_fill_manual(values = clrs, aesthetics = "fill", labels = l_labs) +
  scale_x_discrete(labels = c("Control \n n=72", 
                              "Dark 1M \n n=90", 
                              "Dark 2M \n n=90", 
                              "Light 1M \n n=90", 
                              "Light 2M \n n=90")) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0, 1.1)) +
  ylab("Proportion Settled") +
  theme_bw() + 
  geom_text(data=df.groups$dfav, aes(x = Treatment, y = 1.075, label = letters),vjust=0) # for significance letters
p_df

ggsave("d_fav_settlement_significance.svg", device = "svg", width = 8, height = 6, units = "in", dpi = 300)

#easp
p_ea <- ggplot(df$easp, aes(x = Treatment, y = Percent_settled, fill = Treatment)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.15, shape = 21) +
  scale_fill_manual(values = clrs, aesthetics = "fill", labels = l_labs) +
  scale_x_discrete(labels = c("Control \n n=72", 
                              "Dark 1M \n n=90", 
                              "Dark 2M \n n=90", 
                              "Light 1M \n n=90", 
                              "Light 2M \n n=90")) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0, 1.1)) +
  ylab("Proportion Settled") +
  theme_bw() + 
  geom_text(data=df.groups$easp, aes(x = Treatment, y = 1.075, label = letters),vjust=0) # for significance letters
p_ea

ggsave("e_asp_settlement_significance.svg", device = "svg", width = 8, height = 6, units = "in", dpi = 300)

#psin
p_ps <- ggplot(df$psin, aes(x = Treatment, y = Percent_settled, fill = Treatment)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.15, shape = 21) +
  scale_fill_manual(values = clrs, aesthetics = "fill", labels = l_labs) +
  scale_x_discrete(labels = c("Control \n n=72", 
                              "Dark 1M \n n=90", 
                              "Dark 2M \n n=90", 
                              "Light 1M \n n=90", 
                              "Light 2M \n n=90")) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0, 1.1)) +
  ylab("Proportion Settled") +
  theme_bw() + 
  geom_text(data=df.groups$psin, aes(x = Treatment, y = 1.075, label = letters),vjust=0) # for significance letters
p_ps

ggsave("p_sin_settlement_significance.svg", device = "svg", width = 8, height = 6, units = "in", dpi = 300)

#combine
ggarrange(p_ps +
            ylab("") + ggtitle("P. sinensis") +
            theme(axis.title.x = element_blank(), 
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  plot.title = element_text(face="italic"),
                  legend.text=element_text(size=11), legend.title=element_text(size=13)), 
          p_df + 
            ylab("")  + ggtitle("D. favus") +
            theme(axis.title.x = element_blank(), 
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  plot.title = element_text(face="italic")), 
          p_ea + 
            ggtitle("E. aspera") +
            theme(axis.title.x = element_blank(), 
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  plot.title = element_text(face="italic")),
          p_pl + 
            ylab("") + ggtitle("P. lobata") +
            theme(axis.title.x = element_blank(), 
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  plot.title = element_text(face="italic")),
          p_lc + 
            ylab("") + ggtitle("L. corymbosa") +
            theme(axis.title.x = element_blank(), plot.title = element_text(face="italic")),
          heights = c(0.9, 0.9, 1),
          #widths = c(1, 0.98),
          nrow = 3, ncol = 2,  common.legend = T, legend = "right")

ggsave("settlement_plot_all_coloured.svg", device = "svg", width = 20, height = 20, units = "cm", dpi = 300)




## Chemical trials ----

# first summarise data to plot using function below

# SE function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# summarise data
chem_sum <- list()
for(i in seq_along(df_chem)) {
  chem_sum[[i]] <- summarySE(df_chem[[i]], measurevar="Percent_settled", groupvars=c("Tank", "Treatment"))
  chem_sum[[i]]$Chemical <- c(rep("DCM", 6), rep("EtOH", 6), rep("DCM", 6), rep("EtOH", 6)) 
  names(chem_sum)[i] <- names(df_chem)[i]
}

# change factor levels
chem_levels <- c(
  "DCM 0", "DCM 5", "DCM 10", "DCM 25", "DCM 50", "DCM 100",
  "EtoH 0", "EtoH 5", "EtoH 10", "EtoH 25", "EtoH 50", "EtoH 100"
)

chem_sum <- lapply(chem_sum, function(x) {
  x$Treatment <- factor(x$Treatment, levels = chem_levels)
  return(x)
})

# dodge overlapping values
pd <- position_dodge(0.2) # move them .05 to the left and right

p1 <- ggplot(chem_sum$p_lob, aes(x=Treatment, y=Percent_settled, fill=Tank)) + 
  geom_errorbar(aes(ymin=Percent_settled-se, ymax=Percent_settled+se), width=.1, position = pd) +
  geom_point(size=4, shape=21, position = pd) +
  scale_fill_viridis(discrete = T, name = "Treatment",
                     labels = c("2-Month Light", "Unconditioned")) +
  scale_x_discrete(labels = c("0", "5", "10", "25", "50", "100")) +
  ylim(0, 1) +
  ylab("Proportion Settled") +
  xlab("Volume (µl)") +
  theme_bw() +
  facet_wrap("Chemical", scales = "free_x")
p1

ggsave("p_lob_chem_trials.svg", device = "svg", width = 8, height = 6, units = "in", dpi = 300)

p2 <- ggplot(chem_sum$l_cor, aes(x=Treatment, y=Percent_settled, fill=Tank)) + 
  geom_errorbar(aes(ymin=Percent_settled-se, ymax=Percent_settled+se), width=.1, position = pd) +
  geom_point(size=4, shape=21, position = pd) +
  scale_fill_viridis(discrete = T, name = "Treatment",
                     labels = c("2-Month Light", "Unconditioned")) +
  scale_x_discrete(labels = c("0", "5", "10", "25", "50", "100")) +
  ylim(0, 1) +
  ylab("Proportion Settled") +
  xlab("Volume (µl)") +
  theme_bw() +
  facet_wrap("Chemical", scales = "free_x")
p2

ggsave("l_cor_chem_trials.svg", device = "svg", width = 8, height = 6, units = "in", dpi = 300)

# combine
ggarrange(p1 +
            ylab("") + ggtitle("P. lobata") +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_text(size = 12),
                  plot.title = element_text(face="italic"),
                  legend.text=element_text(size=11), legend.title=element_text(size=13)), 
          p2 + 
            ggtitle("L. corymbosa") +
            theme(axis.title.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  plot.title = element_text(face="italic")), 
          heights = c(0.88, 1.0),
          nrow = 2, ncol = 1,  common.legend = T, legend = "right")

ggsave("settlement_trials_chem_noAten.svg", device = "svg", width = 20, height = 18, units = "cm", dpi = 300)


## summary stats (mean settlement & SD) ----

# Function get mean and SE of for each group 
sum_stats <- function(df) {
  new_df <- group_by(df, Treatment) %>% 
    summarise(mean = mean(Percent_settled, na.rm = T), 
              stdev = sd(Percent_settled, na.rm = T), 
              sterr = sd(Percent_settled, na.rm = T)/sqrt(n()),
              min = min(Percent_settled, na.rm = T),
              max = max(Percent_settled, na.rm = T)) %>% as.data.frame
  return(new_df)
}

# main exp
sum_df <- lapply(df, sum_stats)
names(sum_df) <- names(df)
sum_df

# write table
for(i in seq_along(sum_df)) {
  write.csv(sum_df[[i]], paste0(names(sum_df)[i],"_summary_stats.csv"))
}

# pos ctrl
sum_pos <- lapply(df_pos, sum_stats)
names(sum_pos) <- names(df_pos)
sum_pos

# pos ctrl
sum_dark_rub <- lapply(df_dark_rub, sum_stats)
names(sum_dark_rub) <- names(df_dark_rub)
sum_dark_rub

# pos ctrl
sum_chem <- lapply(df_chem, sum_stats)
names(sum_chem) <- names(df_chem)
sum_chem










