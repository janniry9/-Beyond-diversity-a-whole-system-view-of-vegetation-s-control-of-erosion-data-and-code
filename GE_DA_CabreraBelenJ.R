# "Beyond diversity: a whole system view of vegetation’s control of erosion"
#
#Belen, J.C; Livsey, D.N; Martin, C.B; Groh, T.A; Cardinale, B.J
#
#
####LIBRARIES####
install.packages(c("glmmTMB", "readr", "dplyr", "purrr", "tibble", "ggplot2", "stringr", "broom", "lavaan", "semPlot", "piecewiseSEM", "lme4", "lmerTest", "DHARMa", "mgcv", "ggnewscale","ggeffects", "patchwork"))

#call libraries
library(glmmTMB)
library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(stringr)
library(broom)
library(lavaan)
library(semPlot)
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(DHARMa)
library(mgcv)
library(ggnewscale)
library(patchwork)
library(ggeffects)

####DATA SETUP####
cleaned_new <- read_csv("~/Desktop/Chapter 2_GE_2024/summarized data/mastersheet.csv")
View(cleaned_new)

#removing controls#
cleaned_new_filtered <- cleaned_new %>% 
  filter(treatment != "n")
View(cleaned_new_filtered)

#check data structure #
str(cleaned_new_filtered)

#changing table and treatment to factor#
cleaned_new_filtered <- cleaned_new_filtered %>%
  mutate(
    treatment = factor(treatment),
    table = factor(table),
    richness = as.numeric(richness)  
  )
View(cleaned_new_filtered)

#summary stats#
summary(cleaned_new_filtered)

#NA check for all columns#
colSums(is.na(cleaned_new_filtered))

#problem checking
readr::problems(cleaned_new_filtered) 

#log transforming right skewed root traits, including controls (no plants) #
cleaned_new <- cleaned_new %>%
  mutate(
    log_root_biomass = log(total_root_biomass),
    log_root_length = log(total_root_length),
    log_root_volume = log(total_root_volume),
    log_root_tips = log(total_root_tips),
  )

view(cleaned_new)

cleaned_new[sapply(cleaned_new, is.infinite)] <- 0 #removes INF's and turns them into 0's #

#log transformation of dataset with no controls #
cleaned_new_filtered <- cleaned_new_filtered %>%
  mutate(
    log_root_biomass = log(total_root_biomass),
    log_root_length = log(total_root_length),
    log_root_volume = log(total_root_volume),
    log_root_tips = log(total_root_tips)
  )

View(cleaned_new_filtered)

#Checking for collinearity of root variables #
root_vars <- cleaned_new_filtered %>%
  dplyr::select(
    log_root_biomass,
    log_root_length,
    log_root_volume,
    log_root_tips,
    mean_root_diameter,
  )

cor_mat <- cor(root_vars, use = "complete.obs")
cor_mat 



######t-test####
#T-test on erosion rate in vegetated vs unvegetated (control) mesocosms #
#grouping by vegetation or unvegetated so we can then compare#
erosion_dat <- cleaned_new %>%
  filter(!is.na(erosion_rate), !is.na(richness)) %>%
  mutate(veg_status = factor(ifelse(richness == 0, "Unvegetated", "Vegetated"),
                             levels = c("Unvegetated", "Vegetated")))

#t test
t_result <- t.test(erosion_rate ~ veg_status, data = erosion_dat)
t_result



####GLMM's####
######Richness x Root traits######
#x = richness, y = logged root biomass
richnessxbiomass <- glmmTMB(
  log_root_biomass ~ richness + (1 | table),
  data = cleaned_new_filtered, 
  family = gaussian
)
summary(richnessxbiomass) 
#check residuals#
richness_residuals <- simulateResiduals(richnessxbiomass) 
plot(richness_residuals)#looks good


#x = richness, y = logged root length
richnessxlength <- glmmTMB(
  log_root_length ~ richness + (1 | table),
  data = cleaned_new_filtered, 
  family = gaussian
)
summary(richnessxlength) 
#check residuals#
richness_residuals <- simulateResiduals(richnessxlength) 
plot(richness_residuals)#looks good


# x = richness, y = logged root volume #
richnessxvolume <- glmmTMB(
  log_root_volume ~ richness + (1 | table),
  data = cleaned_new_filtered, 
  family = gaussian
)
summary(richnessxvolume) 
#check residuals#
richness_residuals <- simulateResiduals(richnessxvolume) 
plot(richness_residuals) #looks good


# x = richness, y = mean root diameter #
richnessxdiameter <- glmmTMB(
  mean_root_diameter ~ richness + (1 | table),
  data = cleaned_new_filtered, 
  family = gaussian
)
summary(richnessxdiameter) #significant
#check residuals#
richness_residuals <- simulateResiduals(richnessxdiameter) 
plot(richness_residuals) #looks good



#####Richness x Erosion rate######
# x = richness, y = erosion rate #
richnessxerosion <- glmmTMB(
  erosion_rate ~ richness + (1 | table),
  data = cleaned_new_filtered, 
  family = gaussian
)
summary(richnessxerosion) 
#check residuals#
richness_residuals <- simulateResiduals(richnessxerosion) 
plot(richness_residuals) #looks good



######Figure 2######
# Raw datapoints plotted against species richness, with the fitted GLMM line and its 95% CI.
# Solid line = significant richness effect (p < 0.05); dashed = not significant.

#specs for all the figures
theme_fig2 <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black", size = 12, family = "sans"),
    axis.title = element_text(color = "black", size = 14, face = "bold", family = "sans"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0, family = "sans"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5, unit = "mm")
  )

# Fits the GLMM for one response, then draws points + fitted line + 95% CI
richness_panel <- function(data, yvar, ylab, tag, xlab = "Plant species richness") {
  
  dat <- data %>% filter(!is.na(.data[[yvar]]), !is.na(richness))
  
  mod <- glmmTMB(as.formula(paste(yvar, "~ richness + (1 | table)")),
                 data = dat, family = gaussian())
  
  coefs <- summary(mod)$coefficients$cond
  beta  <- coefs["richness", "Estimate"]
  pval  <- coefs["richness", "Pr(>|z|)"]
  
  pred <- as.data.frame(ggpredict(mod, terms = "richness [1:6 by=0.1]"))
  
  beta_txt <- formatC(beta, digits = 3, format = "fg")
  p_txt    <- ifelse(pval < 0.001, "0.001", formatC(pval, digits = 3, format = "f"))
  p_op     <- ifelse(pval < 0.001, "<", "==")
  stat_lab <- sprintf("beta == '%s'*','~~italic(p) %s '%s'", beta_txt, p_op, p_txt)
  
  p <- ggplot(dat, aes(x = richness, y = .data[[yvar]])) +
    geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high),
                inherit.aes = FALSE, fill = "gray70", alpha = 0.4) +
    geom_point(position = position_jitter(width = 0.12, height = 0, seed = 123),
               alpha = 0.5, size = 2, shape = 16, color = "gray40") +
    geom_line(data = pred, aes(x = x, y = predicted), inherit.aes = FALSE,
              color = "black", linewidth = 1,
              linetype = ifelse(pval < 0.05, "solid", "dashed")) +
    annotate("text", x = Inf, y = Inf, label = stat_lab, parse = TRUE,
             hjust = 1.05, vjust = 1.3, size = 3.8, family = "sans") +
    scale_x_continuous(breaks = c(1, 2, 4, 6), limits = c(0.7, 6.3)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = xlab, y = ylab, title = tag) +
    theme_fig2
  
  message(tag, ": ", yvar, " ~ richness  beta = ", beta_txt, ", p = ", formatC(pval, digits = 4, format = "f"), "  (n = ", nrow(dat), ")")
  p
}

figA <- richness_panel(cleaned_new_filtered, "log_root_biomass",   "Root biomass(log)", "A", xlab = "")
figB <- richness_panel(cleaned_new_filtered, "log_root_volume",    "Root volume (log))",  "B", xlab = "")
figC <- richness_panel(cleaned_new_filtered, "log_root_length",    "Root length (log)",  "C", xlab = "")
figD <- richness_panel(cleaned_new_filtered, "mean_root_diameter", "Mean root diameter (mm)", "D")
figE <- richness_panel(cleaned_new_filtered, "erosion_rate",       "Erosion rate (mm/min)",   "E")

# Combined figure: 3 panels on top, 2 below, sized for a full-width 
fig2 <- (figA | figB | figC) / (figD | figE | plot_spacer())
fig2



######Roots x Erosion####
#x = root biomass, y = erosion rate
erosionxbiomass <- glmmTMB(
  erosion_rate ~ log_root_biomass + (1 | table),
  data = cleaned_new_filtered, 
  family = gaussian(),
  dispformula = ~ log_root_biomass
)
summary(erosionxbiomass)
#check residuals
erosion_residuals <- simulateResiduals(erosionxbiomass) 
plot(erosion_residuals) #quantile deviations#
# rechecking for normality due to quantile deviations
shapiro.test(residuals(erosionxbiomass)) #not significant so data is normal#


#x = root volume, y = erosion rate
erosionxvolume <- glmmTMB(
  erosion_rate ~ log_root_volume + (1 | table),
  data = cleaned_new_filtered, 
  family = gaussian
)
summary(erosionxvolume) 
#check residuals
erosion_residuals <- simulateResiduals(erosionxvolume) 
plot(erosion_residuals) #looks good



#####Interaction x Erosion#####
# RYT < 1 = competition, RYT > 1 = facilitation, RYT = 1 = monoculture (our reference)
#polycultures only b/c RYT measures interactions between species, so needs to be 2+
ryt_dat <- cleaned_new_filtered %>%
  filter(richness >= 2, !is.na(RYT_mean), !is.na(erosion_rate))

#x = RYT, y = Erosion rate #
erosionxRYT <- glmmTMB(
  erosion_rate ~ RYT_mean + (1 | table),
  data = ryt_dat,
  family = gaussian()
)
summary(erosionxRYT)
#checking residuals#
interaction_residuals <- simulateResiduals(erosionxRYT)
plot(interaction_residuals) #quantile deviations detected
#double check for normality 
shapiro.test(residuals(erosionxRYT)) #not significant so we are good



#####Interaction x Soil ####
#remove 1's (monocultures) from the dataset#
cleaned_new_RYT <- cleaned_new_filtered %>% 
  filter(RYT_mean != 1)

#x = RYT, y = soil moisture#
interactionxmoist <- glmmTMB(
  moisture_content ~ RYT_mean + (1 | table),
  data = cleaned_new_RYT, 
  family = beta_family(link = "logit")
)
summary(interactionxmoist) #significant
#check residuals
ryt_residuals <- simulateResiduals(interactionxmoist)
plot(ryt_residuals) #looks good



######Figure 3#####
# Shared theme for all three panels
theme_fig3 <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 16, face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

#(A) Root biomass and erosion
biomass_dat <- cleaned_new_filtered %>% filter(!is.na(log_root_biomass), !is.na(erosion_rate))
erosionxbiomass <- glmmTMB(erosion_rate ~ log_root_biomass + (1 | table),
                           data = biomass_dat, family = gaussian())
pred_biomass <- as.data.frame(ggpredict(erosionxbiomass, terms = "log_root_biomass [all]"))
bio_p <- summary(erosionxbiomass)$coefficients$cond["log_root_biomass", "Pr(>|z|)"]

fig3A <- ggplot(biomass_dat, aes(x = log_root_biomass, y = erosion_rate)) +
  geom_ribbon(data = pred_biomass, aes(x = x, ymin = conf.low, ymax = conf.high),
              inherit.aes = FALSE, fill = "gray80", alpha = 0.5) +
  geom_point(alpha = 0.7, size = 3.5, shape = 16, color = "black") +
  geom_line(data = pred_biomass, aes(x = x, y = predicted), inherit.aes = FALSE,
            color = "black", linewidth = 1.6, lineend = "butt",
            linetype = ifelse(bio_p < 0.05, "solid", "42")) +
  scale_x_continuous(breaks = seq(0, 6, 2)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_fig3 +
  labs(x = "Root biomass (log)", y = "Erosion rate (mm/min)", title = "A")
fig3A

#(B) Root volume and erosion
volume_dat <- cleaned_new_filtered %>% filter(!is.na(log_root_volume), !is.na(erosion_rate))
erosionxvolume <- glmmTMB(erosion_rate ~ log_root_volume + (1 | table),
                          data = volume_dat, family = gaussian())
pred_volume <- as.data.frame(ggpredict(erosionxvolume, terms = "log_root_volume [all]"))
vol_p <- summary(erosionxvolume)$coefficients$cond["log_root_volume", "Pr(>|z|)"]

fig3B <- ggplot(volume_dat, aes(x = log_root_volume, y = erosion_rate)) +
  geom_ribbon(data = pred_volume, aes(x = x, ymin = conf.low, ymax = conf.high),
              inherit.aes = FALSE, fill = "gray80", alpha = 0.5) +
  geom_point(alpha = 0.7, size = 3.5, shape = 16, color = "black") +
  geom_line(data = pred_volume, aes(x = x, y = predicted), inherit.aes = FALSE,
            color = "black", linewidth = 1.6, lineend = "butt",
            linetype = ifelse(vol_p < 0.05, "solid", "42")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_fig3 +
  labs(x = "Root volume (log)", y = "", title = "B")
fig3B

#(C) RYT and erosion (polycultures only)
pred_ryt <- as.data.frame(ggpredict(erosionxRYT, terms = "RYT_mean [all]"))
ryt_p <- summary(erosionxRYT)$coefficients$cond["RYT_mean", "Pr(>|z|)"]

figRYT <- ggplot(ryt_dat, aes(x = RYT_mean, y = erosion_rate)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray50", linewidth = 0.6) +
  geom_ribbon(data = pred_ryt, aes(x = x, ymin = conf.low, ymax = conf.high),
              inherit.aes = FALSE, fill = "gray80", alpha = 0.5) +
  geom_point(alpha = 0.7, size = 3.5, shape = 16, color = "black") +
  geom_line(data = pred_ryt, aes(x = x, y = predicted), inherit.aes = FALSE,
            color = "black", linewidth = 1.6, lineend = "butt",
            linetype = ifelse(ryt_p < 0.05, "solid", "42")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  coord_cartesian(xlim = c(0, 2)) +
  theme_fig3 +
  labs(x = "Relative yield total (RYT)", y = "Erosion rate (mm/min)", title = "C")
figRYT



######Figure 4 #####
#(A): moisture content based on color, with slope lines for dry and wet experimental units and no grid lines#
cleaned_new_filtered <- cleaned_new_filtered %>%
  mutate(moisture_group = ifelse(moisture_content < 0.12, "dry", "wet"))

ggplot(cleaned_new_filtered, aes(x = log_root_length, y = erosion_rate)) +
  
  geom_point(aes(color = moisture_content), alpha = 1, shape = 16, size = 3.5) +
  scale_color_gradient(name = "Moisture content", low = "lightblue", high = "darkblue") +
  ggnewscale::new_scale_color() +
  geom_smooth(
    data = cleaned_new_filtered, 
    mapping = aes(x = log_root_length, y = erosion_rate, color = moisture_group),
    method = "lm",
    se = FALSE,
    linewidth = 1.1,
    inherit.aes = FALSE 
  ) +
  scale_color_manual(name = "Group", values = c("dry" = "brown", "wet" = "navy")) +
  scale_x_continuous(breaks = seq(11, 14, 1), expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = "Root length (log)",
    y = "Erosion rate (mm/min)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.text = element_text(colour = "black", size = 12, family = "sans"),
    axis.title = element_text(colour = "black", size = 16, face = "bold", family = "sans"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "right"
  )

#stats for the wet and dry figure#
# Calculate regression models for each moisture group
dry_model <- lm(erosion_rate ~ log_root_length, data = cleaned_new_filtered %>% filter(moisture_group == "dry"))
wet_model <- lm(erosion_rate ~ log_root_length, data = cleaned_new_filtered %>% filter(moisture_group == "wet"))

# Extract statistics
dry_slope <- coef(dry_model)[2]
dry_intercept <- coef(dry_model)[1]
dry_r2 <- summary(dry_model)$r.squared

wet_slope <- coef(wet_model)[2]
wet_intercept <- coef(wet_model)[1]
wet_r2 <- summary(wet_model)$r.squared

# After extracting the stats, use this code to print to console
dry_label <- paste0("Dry: y = ", round(dry_slope, 3), "x + ", round(dry_intercept, 3), "\nR² = ", round(dry_r2, 3))
wet_label <- paste0("Wet: y = ", round(wet_slope, 3), "x + ", round(wet_intercept, 3), "\nR² = ", round(wet_r2, 3))

cat(dry_label, "\n\n", wet_label, "\n")


#B: RYT and moisture content with fitted beta GLMM line + 95% CI#
pred_moist <- as.data.frame(ggpredict(interactionxmoist, terms = "RYT_mean [all]"))
m_p <- summary(interactionxmoist)$coefficients$cond["RYT_mean", "Pr(>|z|)"]

ggplot(cleaned_new_RYT, aes(x = RYT_mean, y = moisture_content)) +
  # Reference line: RYT = 1 (no net interaction)
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray50", linewidth = 0.6) +
  # 95% confidence band from the model
  geom_ribbon(
    data = pred_moist,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    fill = "#08519c",
    alpha = 0.15
  ) +
  geom_point(
    shape = 21,                
    fill = "#000000",          
    color = "#000000",         
    size = 3.5,                
    stroke = 0,                
    alpha = 0.7                
  ) +
  # Fitted line from the model (solid if p < 0.05, dashed if not)
  geom_line(
    data = pred_moist,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    color = "#08519c",
    linewidth = 1.2,
    linetype = ifelse(m_p < 0.05, "solid", "42")
  ) +
  coord_cartesian(xlim = c(0, 2)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.text = element_text(color = "black", size = 12, family = "sans"),
    axis.title = element_text(color = "black", size = 16, face = "bold", family = "sans"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0, vjust = 2, family = "sans"),
    plot.subtitle = element_text(size = 11, hjust = 0, color = "gray40", family = "sans"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "mm"),
    legend.position = "none"
  ) +
  labs(
    x = "Relative Yield Total (RYT)",
    y = "Moisture Content",
    title = "",
    subtitle = ""
  ) +
  scale_x_continuous(
    breaks = seq(0, 2, 0.5),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1))
  )



####PSEM's####
# make sure table is a factor
cleaned_new_filtered$table <- factor(cleaned_new_filtered$table)

######m1:root length######
#First model: species presence/ absence -> root length
mod_root_length <- lmer(
  log_root_length ~ has_barnyard + has_ryegrass + has_switchgrass + (1 | table),
  data = cleaned_new_filtered
)
summary(mod_root_length)
#checking residuals#
model1L_residuals <- simulateResiduals(mod_root_length) 
plot(model1L_residuals) #within group deviations
#check because of within group deviation flag
mf <- model.frame(mod_root_length)
mf$combo <- interaction(mf$has_barnyard, mf$has_ryegrass, mf$has_switchgrass, sep = "_", drop = TRUE)
cat_test <- testCategorical(model1L_residuals, catPred = mf$combo, plot = FALSE)
cat_test$homogeneity# Levene not significant so we are good


# Second model: root length -> moisture content
mod_moist_length <- lm(
  moisture_content ~ log_root_length,
  data = cleaned_new_filtered
)
summary(mod_moist_length)
#checking residuals#
model2L_residuals <- simulateResiduals(mod_moist_length) 
plot(model2L_residuals)# looks good #


# Third model: moisture content + root length -> erosion rate
mod_erosion_moist_len <- lmer(
  erosion_rate ~ moisture_content + log_root_length + (1 | table),
  data = cleaned_new_filtered
)
summary(mod_erosion_moist_len)
#checking residuals#
model3L_residuals <- simulateResiduals(mod_erosion_moist_len) 
plot(model3L_residuals) #looks good


#PSEM for length#
sem_model_length <- psem(
  mod_root_length,
  mod_moist_length,
  mod_erosion_moist_len,
  data = cleaned_new_filtered
)
# summary
summary(sem_model_length)
coefs(sem_model_length, standardize = "scale") # standardized coefficients
fisherC(sem_model_length) # overall model fit
rsquared(sem_model_length) # R2 values



######m2:root biomass#######
# First model: species presence/absence -> root biomass
mod_root_biomass <- lmer(
  log_root_biomass ~ has_barnyard + has_ryegrass + has_switchgrass + has_clover + (1 | table),
  data = cleaned_new_filtered
)
summary(mod_root_biomass)
#checking residuals#
model1B_residuals <- simulateResiduals(mod_root_biomass) 
plot(model1B_residuals) #quantile deviations
##double check for normality 
shapiro.test(residuals(mod_root_biomass)) #not significant so we are good


# Second model: root biomass -> moisture content
mod_moist_biomass <- lm(
  moisture_content ~ log_root_biomass,
  data = cleaned_new_filtered
)
summary(mod_moist_biomass)
#check residuals
model2B_residuals <- simulateResiduals(mod_moist_biomass) 
plot(model2B_residuals) #looks good#


# Third model: moisture content + root biomass -> erosion rate
mod_erosion_moist_bio <- lmer(
  erosion_rate ~ moisture_content + log_root_biomass + (1 | table),
  data = cleaned_new_filtered
)
summary(mod_erosion_moist_bio)
#check residuals 
model3B_residuals <- simulateResiduals(mod_erosion_moist_bio) 
plot(model3B_residuals) #looks good


#PSEM model for biomass#
sem_model_biomass <- psem(
  mod_root_biomass,
  mod_moist_biomass,
  mod_erosion_moist_bio,
  data = cleaned_new_filtered
)

# summary
summary(sem_model_biomass)

coefs(sem_model_biomass, standardize = "scale") # standardized coefficients

fisherC(sem_model_biomass) # overall model fit

rsquared(sem_model_biomass) # R2 values



######m3:root volume#####
#first model: species presence/absence -> root volume
mod_root_volume <- lmer(
  log_root_volume ~ has_barnyard + has_ryegrass + has_virginia + (1 | table),
  data = cleaned_new_filtered
)
summary(mod_root_volume)
#residuals#
model1V_residuals <- simulateResiduals(mod_root_volume) 
plot(model1V_residuals) #looks good


# Second model: roots + virginia wild rye -> moisture content
mod_moist_volume <- lm(
  moisture_content ~ log_root_volume + has_virginia,
  data = cleaned_new_filtered
)
summary(mod_moist_volume)
#residuals#
model2V_residuals <- simulateResiduals(mod_moist_volume) 
plot(model2V_residuals) #looks good


# Third model: moisture content + root volume -> erosion rate
mod_erosion_moist_vol <- lmer(
  erosion_rate ~ moisture_content + log_root_volume + (1 | table),
  data = cleaned_new_filtered
)
summary(mod_erosion_moist_vol)
#residuals#
model3V_residuals <- simulateResiduals(mod_erosion_moist_vol) 
plot(model3V_residuals) #looks good


# Full SEM
sem_model_volume <- psem(
  mod_root_volume,
  mod_moist_volume,
  mod_erosion_moist_vol,
  data = cleaned_new_filtered
)

#summary
summary(sem_model_volume)

fisherC(sem_model_volume)# overall model fit

coefs(sem_model_volume, standardize = "scale") # standardized coefficients

rsquared(sem_model_volume)



######m4: root diameter#####
#First model: species presence/absence -> mean root diameter
mod_root_diameter <- lm(
  mean_root_diameter ~ has_barnyard + has_virginia,
  data = cleaned_new_filtered
)
summary(mod_root_diameter)
#residuals#
model1D_residuals <- simulateResiduals(mod_root_diameter) 
plot(model1D_residuals) #looks good#


# Second model: root diameter -> moisture content
mod_moist_diameter <- lm(
  moisture_content ~ mean_root_diameter + has_virginia + has_switchgrass,
  data = cleaned_new_filtered
)
summary(mod_moist_diameter)
#residuals#
model2D_residuals <- simulateResiduals(mod_moist_diameter) 
plot(model2D_residuals) #looks good


# Third model: moisture content * root volume -> erosion rate
mod_erosion_moist_diam <- lmer(
  erosion_rate ~ moisture_content + mean_root_diameter + (1 | table),
  data = cleaned_new_filtered
)
summary(mod_erosion_moist_diam)
#residuals#
model3D_residuals <- simulateResiduals(mod_erosion_moist_diam) 
plot(model3D_residuals)# looks good


# Full SEM
sem_model_diameter <- psem(
  mod_root_diameter,
  mod_moist_diameter,
  mod_erosion_moist_diam,
  data = cleaned_new_filtered
)

#summary
summary(sem_model_diameter)

fisherC(sem_model_diameter)# overall model fit

coefs(sem_model_diameter, standardize = "scale") # standardized coefficients

rsquared(sem_model_diameter)


####### AIC ####
AIC(sem_model_length, sem_model_biomass, sem_model_volume, sem_model_diameter,
    AIC.type = "dsep", aicc = TRUE)

