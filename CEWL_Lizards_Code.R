################################################################################
# Urban vs Non-Urban Lizard CEWL Project
# Cleaned Review Round 2 Analysis Script
#################################################################################

# ------------------------------------------------------------------------------
# 0.SETUP
# ------------------------------------------------------------------------------

setwd("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards")

# Core packages
library(tidyverse)
library(ggplot2)
library(patchwork)
library(emmeans)
library(car)
library(lme4)
library(lmerTest)
library(broom)
library(writexl)
library(ggrepel)
library(raster)
library(sp)
library(Cairo)

#Consistent colors
habitat_colors <- c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")

#Rename labels for publication plots without changing raw data
habitat_labs <- c("wall" = "Urban", "nonwall" = "Nonurban")

#Output folders
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("outputs/figures")) dir.create("outputs/figures", recursive = TRUE)
if (!dir.exists("outputs/tables")) dir.create("outputs/tables", recursive = TRUE)


# ------------------------------------------------------------------------------
# 1.LOAD DATA AND CREATE CEWL MASTER DATASET
# ------------------------------------------------------------------------------

dat <- read.csv("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/df3_final.csv")
location_dat1 <- read.csv("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/data-location.csv")

#Use only individuals with CEWL 
cewl.dat <- dat %>%
  filter(!is.na(CEWL)) %>%
  mutate(
    LogSVL = ifelse(!is.na(SVL), log(SVL), NA),
    es = 0.611 * exp((17.502 * ambient_temp) / (ambient_temp + 240.97)),
    ea = es * ambient_percent_rh / 100,
    VPD = es - ea
  )

if (nrow(cewl.dat) == 0) {
  stop("cewl.dat is empty. No CEWL values found. Check input data.")
}

#Limit site-level dataset to sites represented in CEWL data
location_dat <- location_dat1 %>%
  filter(Site %in% unique(cewl.dat$Site))

message("CEWL master dataset created: ", nrow(cewl.dat), " individuals across ",
        length(unique(cewl.dat$Site)), " sites.")

rm(dat)


# ------------------------------------------------------------------------------
# 2.EXTRACT / PREPARE MACROCLIMATE VARIABLES
# ------------------------------------------------------------------------------
# This section keeps only the raster extraction needed to build macroclimate PCA.

coords <- data.frame(lat = location_dat$Lat, long = location_dat$Long)
coordinates(coords) <- ~long + lat

# ---- precipitation ----
precip_files <- list.files(
  "~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/",
  pattern = "wc2.1_30s_prec_.*\\.tif$",
  full.names = TRUE
)
precip_stack <- stack(precip_files)
crs(coords) <- crs(precip_stack)
precip_values <- extract(precip_stack, coords, method = "bilinear")
location_dat$annual_precip <- rowSums(precip_values, na.rm = TRUE)
location_dat$monthly_avg_precip <- rowMeans(precip_values, na.rm = TRUE)

# ---- solar radiation ----
srad_files <- list.files(
  "~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/",
  pattern = "wc2.1_30s_srad_.*\\.tif$",
  full.names = TRUE
)
srad_stack <- stack(srad_files)
crs(coords) <- crs(srad_stack)
srad_values <- extract(srad_stack, coords, method = "bilinear")
location_dat[paste0("srad_", 1:12)] <- srad_values
location_dat$monthly_avg_srad <- rowMeans(srad_values, na.rm = TRUE)

# ---- water vapor pressure ----
vapr_files <- list.files(
  "~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/",
  pattern = "wc2.1_30s_vapr_.*\\.tif$",
  full.names = TRUE
)
vapr_stack <- stack(vapr_files)
crs(coords) <- crs(vapr_stack)
vapr_values <- extract(vapr_stack, coords, method = "bilinear")
location_dat[paste0("vapr_", 1:12)] <- vapr_values
location_dat$monthly_avg_vapr <- rowMeans(vapr_values, na.rm = TRUE)

# ---- maximum temperature ----
tmax_files <- list.files(
  "~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/",
  pattern = "wc2.1_30s_tmax_.*\\.tif$",
  full.names = TRUE
)
tmax_stack <- stack(tmax_files)
crs(coords) <- crs(tmax_stack)
tmax_values <- extract(tmax_stack, coords, method = "bilinear")
location_dat[paste0("tmax_", 1:12)] <- tmax_values
location_dat$monthly_avg_tmax <- rowMeans(tmax_values, na.rm = TRUE)

# ---- minimum temperature ----
tmin_files <- list.files(
  "~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/",
  pattern = "wc2.1_30s_tmin_.*\\.tif$",
  full.names = TRUE
)
tmin_stack <- stack(tmin_files)
crs(coords) <- crs(tmin_stack)
tmin_values <- extract(tmin_stack, coords, method = "bilinear")
location_dat[paste0("tmin_", 1:12)] <- tmin_values
location_dat$monthly_avg_tmin <- rowMeans(tmin_values, na.rm = TRUE)

# ---- average temperature ----
tavg_files <- list.files(
  "~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/",
  pattern = "wc2.1_30s_tavg_.*\\.tif$",
  full.names = TRUE
)
tavg_stack <- stack(tavg_files)
crs(coords) <- crs(tavg_stack)
tavg_values <- extract(tavg_stack, coords, method = "bilinear")
location_dat[paste0("tavg_", 1:12)] <- tavg_values
location_dat$monthly_avg_tavg <- rowMeans(tavg_values, na.rm = TRUE)

# ---- macroclimate VPD ----
# VPD is calculated from monthly tmax and vapor pressure.
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  vapr_col <- paste0("vapr_", i)
  vpd_col  <- paste0("VPD_", i)
  location_dat[[vpd_col]] <-
    0.611 * exp((17.502 * location_dat[[tmax_col]]) /
                  (location_dat[[tmax_col]] + 240.97)) -
    location_dat[[vapr_col]]
}
location_dat$VPD_avg <- rowMeans(location_dat[paste0("VPD_", 1:12)], na.rm = TRUE)


# ------------------------------------------------------------------------------
# 3.Does CEWL differ between urban and nonurban lizards?
# ------------------------------------------------------------------------------
# A priori model. Habitat is the predictor of interest; Sex and LogSVL are biologically motivated covariates. Site was tested as a random intercept, but in our results the random-effect variance was estimated as zero, so the simpler linear model is retained for reporting.

CEWL_mod_lmer <- lmer(
  CEWL ~ Habitat + Sex + LogSVL + (1 | Site),
  data = cewl.dat
)
summary(CEWL_mod_lmer)
isSingular(CEWL_mod_lmer)

CEWL_mod_lm <- lm(
  CEWL ~ Habitat + Sex + LogSVL,
  data = cewl.dat
)
summary(CEWL_mod_lm)
anova(CEWL_mod_lm)

# Estimated marginal means for plotting and reporting
CEWL_emm <- emmeans(CEWL_mod_lm, ~ Habitat)
CEWL_emm_df <- as.data.frame(summary(CEWL_emm))
CEWL_pairwise <- as.data.frame(pairs(CEWL_emm))

# Supplementary table for Question A
CEWL_table <- broom::tidy(CEWL_mod_lm) %>%
  mutate(model = "CEWL ~ Habitat + Sex + LogSVL")
write.csv(CEWL_table, "outputs/tables/Table_CEWL_habitat_model.csv", row.names = FALSE)
write.csv(CEWL_pairwise, "outputs/tables/Table_CEWL_emmeans_pairwise.csv", row.names = FALSE)

# Figure: CEWL by habitat
CEWL_plot <- ggplot(
  cewl.dat %>% mutate(Habitat = factor(Habitat, levels = c("wall", "nonwall"))),
  aes(x = Habitat, y = CEWL, color = Habitat)
) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.75) +
  geom_errorbar(
    data = CEWL_emm_df,
    aes(x = Habitat, ymin = lower.CL, ymax = upper.CL),
    width = 0.2,
    linewidth = 0.8,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = CEWL_emm_df,
    aes(x = Habitat, y = emmean),
    color = "black",
    size = 4,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = habitat_colors, labels = habitat_labs) +
  scale_x_discrete(labels = habitat_labs) +
  labs(
    x = "Habitat",
    y = expression(CEWL~(g/m^{2}~h))
  ) +
  theme_classic(base_size = 14) +
  guides(color = "none")

print(CEWL_plot)
ggsave("outputs/figures/Figure_CEWL_by_habitat.pdf", CEWL_plot,
       width = 5, height = 4, units = "in")


# ------------------------------------------------------------------------------
# 4. Does microclimate differ between habitats after accounting for macroclimate?
# ------------------------------------------------------------------------------
# Reviewer-requested framework:
#   MicroPC1 ~ Habitat * MacroPC1
# PC1 is used at both scales because it captures the dominant environmental axis.

#Differences in microclimatic variables
temp_mod <- lm(ambient_temp ~ Habitat, data = cewl.dat)
summary(temp_mod)
anova(temp_mod)

vpd_mod <- lm(VPD ~ Habitat, data = cewl.dat)
summary(vpd_mod)
anova(vpd_mod)

veg_mod <- lm(percent_veg_cover ~ Habitat, data = cewl.dat)
summary(veg_mod)
anova(veg_mod)

# ---- microclimate PCA ----
micro_pca_input <- cewl.dat %>%
  dplyr::select(ambient_temp, percent_veg_cover, VPD) %>%
  drop_na()

micro_pca <- prcomp(micro_pca_input, center = TRUE, scale. = TRUE)
summary(micro_pca)
micro_pca$rotation

# Add microclimate PC scores to matching rows
# This assumes the three microclimate variables have no missing values in cewl.dat.
# If missing values exist, use the complete-case object below instead.
cewl.dat <- cewl.dat %>%
  mutate(
    MicroPC1 = micro_pca$x[, 1],
    MicroPC2 = micro_pca$x[, 2],
    MicroPC3 = micro_pca$x[, 3]
  )

# ---- macroclimate PCA ----
macro_pca_input <- location_dat %>%
  dplyr::select(
    NDVI,
    NDMI,
    annual_precip,
    monthly_avg_srad,
    monthly_avg_precip,
    monthly_avg_vapr,
    monthly_avg_tmax,
    monthly_avg_tmin,
    monthly_avg_tavg,
    VPD_avg
  ) %>%
  drop_na()

macro_pca <- prcomp(macro_pca_input, center = TRUE, scale. = TRUE)
summary(macro_pca)
macro_pca$rotation

location_dat <- location_dat %>%
  mutate(
    MacroPC1 = macro_pca$x[, 1],
    MacroPC2 = macro_pca$x[, 2]
  )

# Merge MacroPC1 into the individual CEWL dataset
# Remove an old MacroPC1 if the script is rerun in the same R session.
cewl.dat <- cewl.dat %>%
  dplyr::select(-any_of(c("MacroPC1", "MacroPC2"))) %>%
  left_join(location_dat %>% dplyr::select(Site, MacroPC1, MacroPC2), by = "Site")

# Main model for Question B
MicroModel <- lm(
  MicroPC1 ~ Habitat * MacroPC1,
  data = cewl.dat
)
summary(MicroModel)
anova(MicroModel)

# Supplementary PCA and model tables
micro_pca_loadings <- as.data.frame(micro_pca$rotation) %>%
  rownames_to_column("Variable")
macro_pca_loadings <- as.data.frame(macro_pca$rotation) %>%
  rownames_to_column("Variable")

micro_pca_variance <- as.data.frame(t(summary(micro_pca)$importance)) %>%
  rownames_to_column("PC")
macro_pca_variance <- as.data.frame(t(summary(macro_pca)$importance)) %>%
  rownames_to_column("PC")

MicroModel_table <- broom::tidy(MicroModel) %>%
  mutate(model = "MicroPC1 ~ Habitat * MacroPC1")
MicroModel_anova <- as.data.frame(anova(MicroModel)) %>%
  rownames_to_column("Term")

write.csv(micro_pca_loadings, "outputs/tables/Table_microclimate_PCA_loadings.csv", row.names = FALSE)
write.csv(macro_pca_loadings, "outputs/tables/Table_macroclimate_PCA_loadings.csv", row.names = FALSE)
write.csv(micro_pca_variance, "outputs/tables/Table_microclimate_PCA_variance.csv", row.names = FALSE)
write.csv(macro_pca_variance, "outputs/tables/Table_macroclimate_PCA_variance.csv", row.names = FALSE)
write.csv(MicroModel_table, "outputs/tables/Table_microPC1_macroPC1_model.csv", row.names = FALSE)
write.csv(MicroModel_anova, "outputs/tables/Table_microPC1_macroPC1_ANOVA.csv", row.names = FALSE)

# Figure B1: Macroclimate PCA space by site
macro_scores <- as.data.frame(macro_pca$x) %>%
  mutate(
    Site = location_dat$Site,
    Habitat = location_dat$Habitat
  )

point_style <- list(
  geom_point(size = 2.5, alpha = 0.75)
)

ellipse_style <- list(
  stat_ellipse(level = 0.95, type = "norm", linewidth = 0.8)
)

MacroPCA_plot <- ggplot(macro_scores, aes(PC1, PC2, color = Habitat)) +
  point_style +
  ellipse_style +
  geom_text_repel(aes(label = Site), size = 3, show.legend = FALSE) +
  scale_color_manual(values = habitat_colors, labels = habitat_labs) +
  labs(
    x = paste0("Macroclimate PC1 (", round(summary(macro_pca)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("Macroclimate PC2 (", round(summary(macro_pca)$importance[2, 2] * 100, 1), "% variance)"),
    color = "Habitat"
  ) +
  theme_classic(base_size = 14)

print(MacroPCA_plot)
ggsave("outputs/figures/Figure_macroclimate_PCA.pdf", MacroPCA_plot,
       width = 6, height = 5, units = "in")

# Figure B2: Microclimate PCA space by individual lizard
micro_scores <- as.data.frame(micro_pca$x) %>%
  mutate(Habitat = cewl.dat$Habitat)

MicroPCA_plot <- ggplot(micro_scores, aes(PC1, PC2, color = Habitat)) +
  point_style +
  ellipse_style +
  scale_color_manual(values = habitat_colors, labels = habitat_labs) +
  labs(
    x = paste0("Microclimate PC1 (", round(summary(micro_pca)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("Microclimate PC2 (", round(summary(micro_pca)$importance[2, 2] * 100, 1), "% variance)"),
    color = "Habitat"
  ) +
  theme_classic(base_size = 14)

print(MicroPCA_plot)
ggsave("outputs/figures/Figure_microclimate_PCA.pdf", MicroPCA_plot,
       width = 6, height = 5, units = "in")

# Figure B3: reviewer-requested macro-PC1 by micro-PC1 plot
MacroMicro_plot <- ggplot(cewl.dat, aes(x = MacroPC1, y = MicroPC1, color = Habitat)) +
  point_style +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = habitat_colors, labels = habitat_labs) +
  labs(
    x = paste0("Macroclimate PC1 (", round(summary(macro_pca)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("Microclimate PC1 (", round(summary(micro_pca)$importance[2, 1] * 100, 1), "% variance)"),
    color = "Habitat"
  ) +
  theme_classic(base_size = 14)

print(MacroMicro_plot)
ggsave("outputs/figures/Figure_macroPC1_vs_microPC1.pdf", MacroMicro_plot,
       width = 6, height = 5, units = "in")

# Optional combined figure for manuscript
PCA_combined_plot <- MicroPCA_plot + MacroPCA_plot + MacroMicro_plot +
  plot_layout(ncol = 3)

ggsave("outputs/figures/Figure_combined_PCA_panels.pdf", PCA_combined_plot,
       width = 15, height = 5, units = "in")


#emtrends
micro_slopes <- emtrends(
  MicroModel,
  specs = "Habitat",
  var = "MacroPC1"
)

summary(micro_slopes)
pairs(micro_slopes)


# ------------------------------------------------------------------------------
# 5. Does microclimate explain CEWL across all individuals?
# ------------------------------------------------------------------------------

CEWL_env_lmer <- lmer(
  CEWL ~ VPD + ambient_temp + percent_veg_cover + (1 | Site),
  data = cewl.dat
)
summary(CEWL_env_lmer)
isSingular(CEWL_env_lmer)



CEWL_env_lm <- lm(
 CEWL ~ VPD + ambient_temp + percent_veg_cover,
  data = cewl.dat
)
summary(CEWL_env_lm)
anova(CEWL_env_lm)
car::vif(CEWL_env_lm)

# Diagnostic model only: shows why Habitat and VPD should not be in the same explanatory model. This is not intended as a main analysis because it produces high VIF for Habitat and VPD.
CEWL_env_habitat_diagnostic <- lm(
  CEWL ~ Habitat + VPD + ambient_temp + percent_veg_cover,
  data = cewl.dat
)
summary(CEWL_env_habitat_diagnostic)
car::vif(CEWL_env_habitat_diagnostic)

# Supplementary table for Question C
CEWL_env_table <- broom::tidy(CEWL_env_lm) %>%
  mutate(model = "CEWL ~ VPD + ambient_temp + percent_veg_cover")
CEWL_env_anova <- as.data.frame(anova(CEWL_env_lm)) %>%
  rownames_to_column("Term")
CEWL_env_vif <- data.frame(
  Variable = names(car::vif(CEWL_env_lm)),
  VIF = as.numeric(car::vif(CEWL_env_lm))
)

write.csv(CEWL_env_table, "outputs/tables/Table_CEWL_microclimate_model.csv", row.names = FALSE)
write.csv(CEWL_env_anova, "outputs/tables/Table_CEWL_microclimate_ANOVA.csv", row.names = FALSE)
write.csv(CEWL_env_vif, "outputs/tables/Table_CEWL_microclimate_VIF.csv", row.names = FALSE)

# Figure: CEWL vs VPD across all individuals
CEWL_VPD_plot <- ggplot(cewl.dat, aes(x = VPD, y = CEWL, color = Habitat)) +
  point_style +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  scale_color_manual(values = habitat_colors, labels = habitat_labs) +
  labs(
    x = "VPD at capture (kPa)",
    y = expression(CEWL~(g/m^{2}~h)),
    color = "Habitat"
  ) +
  theme_classic(base_size = 14)

print(CEWL_VPD_plot)
ggsave("outputs/figures/Figure_CEWL_vs_VPD.pdf", CEWL_VPD_plot,
       width = 5.5, height = 4.5, units = "in")


# ------------------------------------------------------------------------------
# 6. OPTIONAL DESCRIPTIVE MICROCLIMATE FIGURES
# ------------------------------------------------------------------------------


Temp_density <- ggplot(cewl.dat, aes(x = ambient_temp, fill = Habitat)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = habitat_colors, labels = habitat_labs) +
  labs(
    x = expression("Ambient temperature at capture (" * C * degree * ")"),
    y = "Density",
    fill = "Habitat"
  ) +
  theme_classic(base_size = 14)

VPD_density <- ggplot(cewl.dat, aes(x = VPD, fill = Habitat)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = habitat_colors, labels = habitat_labs) +
  labs(x = "VPD at capture (kPa)", y = "Density", fill = "Habitat") +
  theme_classic(base_size = 14)

Veg_density <- ggplot(cewl.dat, aes(x = percent_veg_cover, fill = Habitat)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = habitat_colors, labels = habitat_labs) +
  labs(x = "Vegetation cover (%)", y = "Density", fill = "Habitat") +
  theme_classic(base_size = 14)

Microclimate_density_plot <- Temp_density + VPD_density + Veg_density +
  plot_layout(ncol = 3, guides = "collect")

ggsave("outputs/figures/Figure_descriptive_microclimate_densities.pdf",
       Microclimate_density_plot, width = 12, height = 4, units = "in")


# ------------------------------------------------------------------------------
# 7. OPTIONAL ENVIRONMENTAL MAP FIGURE
# ------------------------------------------------------------------------------


make_environmental_map <- TRUE

if (make_environmental_map) {
  total_precip <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/total_precipitation.tif")

  climcol <- colorRampPalette(c("purple", "blue", "skyblue", "green",
                                "lightgreen", "yellow", "orange", "red", "darkred"))
  climcol_rev <- colorRampPalette(rev(c("purple", "blue", "skyblue", "green",
                                        "lightgreen", "yellow", "orange", "red", "darkred")))

  xlim_use <- c(25.33, 25.65)
  ylim_use <- c(36.9, 37.21)
  vpd_extent <- extent(25.31, 25.65, 36.9, 37.21)

  add_sites <- function() {
    points(location_dat$Long[location_dat$Habitat == "nonwall"],
           location_dat$Lat[location_dat$Habitat == "nonwall"],
           pch = 21, bg = habitat_colors["nonwall"], col = "black", cex = 1.4)
    points(location_dat$Long[location_dat$Habitat == "wall"],
           location_dat$Lat[location_dat$Habitat == "wall"],
           pch = 24, bg = habitat_colors["wall"], col = "black", cex = 1.4)
  }

  panel_label <- function(text) {
    mtext(text, side = 3, line = 0.15, adj = 0, font = 2, cex = 0.8)
  }

  tmax_5 <- raster(tmax_files[5])
  tmax_6 <- raster(tmax_files[6])
  vapr_5 <- raster(vapr_files[5])
  vapr_6 <- raster(vapr_files[6])

  tmax5_crop <- crop(tmax_5, vpd_extent)
  vapr5_crop <- crop(vapr_5, vpd_extent)
  VPD_May <- (0.611 * exp((17.502 * tmax5_crop) / (tmax5_crop + 240.97))) - vapr5_crop

  tmax6_crop <- crop(tmax_6, vpd_extent)
  vapr6_crop <- crop(vapr_6, vpd_extent)
  VPD_June <- (0.611 * exp((17.502 * tmax6_crop) / (tmax6_crop + 240.97))) - vapr6_crop

  upsample_factor <- 4
  VPD_May_hi <- disaggregate(VPD_May, fact = upsample_factor, method = "bilinear")
  VPD_June_hi <- disaggregate(VPD_June, fact = upsample_factor, method = "bilinear")
  tmax_5_hi <- disaggregate(crop(tmax_5, vpd_extent), fact = upsample_factor, method = "bilinear")
  tmax_6_hi <- disaggregate(crop(tmax_6, vpd_extent), fact = upsample_factor, method = "bilinear")
  precip_hi <- disaggregate(crop(total_precip, vpd_extent), fact = upsample_factor, method = "bilinear")

  Cairo::CairoPDF("outputs/figures/Figure_environmental_maps.pdf", width = 7.5, height = 10)

  layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE),
         widths = c(1, 1), heights = c(1, 1, 1))

  mar_left <- c(0.4, 0.05, 1.4, 0.12)
  mar_right <- c(0.4, 0.05, 1.4, 1.6)

  par(mar = mar_left)
  plot(VPD_May_hi, xlim = xlim_use, ylim = ylim_use, zlim = c(1, 2),
       col = climcol(100), axes = FALSE, ann = FALSE, box = FALSE,
       legend = FALSE, useRaster = TRUE, interpolate = FALSE)
  add_sites()
  panel_label("A. May VPD")

  par(mar = mar_right)
  plot(VPD_June_hi, xlim = xlim_use, ylim = ylim_use, zlim = c(1, 2),
       col = climcol(100), axes = FALSE, ann = FALSE, box = FALSE,
       legend = TRUE,
       legend.args = list(text = "VPD (kPa)", side = 4, line = 1.1, cex = 0.85),
       useRaster = TRUE, interpolate = FALSE)
  add_sites()
  panel_label("B. June VPD")

  par(mar = mar_left)
  plot(tmax_5_hi, xlim = xlim_use, ylim = ylim_use, zlim = c(19, 28),
       col = climcol(100), axes = FALSE, ann = FALSE, box = FALSE,
       legend = FALSE, useRaster = TRUE, interpolate = FALSE)
  add_sites()
  panel_label("C. May Maximum Temperature")

  par(mar = mar_right)
  plot(tmax_6_hi, xlim = xlim_use, ylim = ylim_use, zlim = c(19, 28),
       col = climcol(100), axes = FALSE, ann = FALSE, box = FALSE,
       legend = TRUE,
       legend.args = list(text = "Maximum Temperature (C)", side = 4, line = 1.1, cex = 0.85),
       useRaster = TRUE, interpolate = FALSE)
  add_sites()
  panel_label("D. June Maximum Temperature")

  par(mar = mar_left)
  plot(precip_hi, xlim = xlim_use, ylim = ylim_use,
       col = climcol_rev(100), axes = FALSE, ann = FALSE, box = FALSE,
       legend = FALSE, useRaster = TRUE, interpolate = FALSE)
  add_sites()
  panel_label("E. Annual Precipitation")

  par(mar = c(0, 0, 0, 0))
  plot.new()

  dev.off()
}


# ------------------------------------------------------------------------------
# 8. EXPORT A COMBINED SUPPLEMENTARY EXCEL FILE
# ------------------------------------------------------------------------------

supplement_tables <- list(
  CEWL_habitat_model = CEWL_table,
  CEWL_emmeans_pairwise = CEWL_pairwise,
  microclimate_model = MicroModel_table,
  microclimate_model_ANOVA = MicroModel_anova,
  CEWL_microclimate_model = CEWL_env_table,
  CEWL_microclimate_ANOVA = CEWL_env_anova,
  CEWL_microclimate_VIF = CEWL_env_vif,
  micro_PCA_loadings = micro_pca_loadings,
  macro_PCA_loadings = macro_pca_loadings,
  micro_PCA_variance = micro_pca_variance,
  macro_PCA_variance = macro_pca_variance
)

writexl::write_xlsx(supplement_tables, "outputs/tables/Supplementary_analysis_tables.xlsx")

message("Cleaned review-round analysis complete. Outputs written to outputs/figures and outputs/tables.")

