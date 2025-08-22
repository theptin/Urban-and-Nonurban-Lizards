#Urban VS Non-Urban Lizard Project

#Set WD
setwd("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards")

#Load Packages
library('scales')
library('lubridate')
library('factoextra')
library('RColorBrewer')
library('dplyr')
library('tidyverse')
library('ggpubr') 
library('rstatix')
library('broom')
library('vegan')
library('ggplot2')
library('nlme')
library('knitr')
library("writexl")
library("MuMIn")
library('skimr')
library("patchwork")
library('car')
library('emmeans')
library('ggvegan')
library('forcats')
library('raster')
library('sp')
library('rworldmap')
library('terra')
library('Cairo')
library('lme4')
library('lmerTest')
library('gridExtra')
library('sf')
library('ggrepel')
library('tibble')


##################################################################################
#Data setup
#Load Data
dat <- read.csv("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/df3_final.csv")
dat


#log transformation
dat$LogMass <- log(dat$Mass)

dat$LogSVL<- log(dat$SVL)

dat$LogHW<- log(dat$HeadWidth)

dat$LogHD<- log(dat$HeadDepth)

dat$LogHL<- log(dat$HeadLength)

dat$LogFL<- log(dat$FemurLength)

dat$LogTL<- log(dat$TibiaLength)

dat$LogBL<- log(dat$BicepLength)

dat$LogFAL<- log(dat$ForearmLength)

head(dat)

#(FemurLength + TibiaLength)/(BicepLength + ForearmLength)
dat$LimbRatio <- ((dat$FemurLength + dat$TibiaLength)/(dat$BicepLength + dat$ForearmLength))
dat$LogLR <- log(dat$LimbRatio)
head(dat)

#Donihue 2016- Standardized SVL to a mean of 0 and used linear mixed models
#Scaling SVL to a mean of 0 and st. deviation of 1
dat$SVLZERO <- scale(dat$SVL, center = TRUE, scale = TRUE)
head(dat)

#To Analyze Groups Seperately
#Make groups to test: Wall and Nonwall; Male Wall and Nonwall, Female Wall and Nonwall
datWall <- subset(dat, Habitat == "wall")
datNonWall <- subset(dat, Habitat == "nonwall")

datMWall <- subset(datWall, Sex == "M")
datFWall <- subset(datWall, Sex == "F")

datMNonWall <- subset(datNonWall, Sex == "M")
datFNonWall <- subset(datNonWall, Sex == "F")

datM <- subset(dat, Sex == "M")
datF <- subset(dat, Sex == "F")

#subset Data for CEWL
cewl.dat <- subset(dat, !is.na(CEWL))
head(cewl.dat)

cewl.datM <- subset(cewl.dat, Sex == "M")
cewl.datF <- subset(cewl.dat, Sex == "F")

#VPD
head(cewl.dat)
# Calculate VPD
cewl.dat$es <- 0.611 * exp((17.502 * cewl.dat$ambient_temp) / (cewl.dat$ambient_temp + 240.97))
cewl.dat$ea <- cewl.dat$es * cewl.dat$ambient_percent_rh / 100
cewl.dat$VPD <- cewl.dat$es - cewl.dat$ea

# View the updated data
head(cewl.dat)

#Caluclate VPD in dat
dat$es <- 0.611 * exp((17.502 * dat$ambient_temp) / (dat$ambient_temp + 240.97))
dat$ea <- dat$es * dat$ambient_percent_rh / 100
dat$VPD <- dat$es - dat$ea


datCEWLWall <- subset(cewl.dat, Habitat == "wall")
datCEWLWallM <- subset(datCEWLWall, Sex == "M")
datCEWLWallF <- subset(datCEWLWall, Sex == "F")
datCEWLNonWall <- subset(cewl.dat, Habitat == "nonwall")
datCEWLNonWallM <- subset(datCEWLNonWall, Sex == "M")
datCEWLNonWallF <- subset(datCEWLNonWall, Sex == "F")

##########################################################################################################
#Analyses testing for differences in individual lizard traits using emmeans 
#emmeans
#SVL
ANCOVASVL <- lm(SVL ~ Habitat + Sex, data = dat)
summary(ANCOVASVL)
SVL_emm <- emmeans(ANCOVASVL, ~ Habitat)
summary(SVL_emm)
SVL_pairwise <- pairs(SVL_emm)
summary(SVL_pairwise)
SVL_emm <- emmeans(ANCOVASVL, ~ Habitat)
SVL_emm_df <- as.data.frame(summary(SVL_emm))

SVL_emm_df$Habitat <- factor(SVL_emm_df$Habitat, levels = c("wall", "nonwall"))

SVL_box <- ggplot(dat %>%
                    mutate(Habitat = fct_relevel(Habitat, "wall")), 
                  aes(x = as.factor(Habitat), y = SVL, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = SVL_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = SVL_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "SVL Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "SVL (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(SVL_box)

#Mass
ANCOVAMass <- lm(Mass ~ Habitat + Sex, data = dat)
summary(ANCOVAMass)
Mass_emm <- emmeans(ANCOVAMass, ~ Habitat)
summary(Mass_emm)
Mass_pairwise <- pairs(Mass_emm)
summary(Mass_pairwise)
Mass_emm <- emmeans(ANCOVAMass, ~ Habitat)
Mass_emm_df <- as.data.frame(summary(Mass_emm))

Mass_emm_df$Habitat <- factor(Mass_emm_df$Habitat, levels = c("wall", "nonwall"))

Mass_box <- ggplot(dat %>%
                    mutate(Habitat = fct_relevel(Habitat, "wall")), 
                  aes(x = as.factor(Habitat), y = Mass, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = Mass_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = Mass_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "Mass Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "Mass (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(Mass_box)

#HeadWidth
ANCOVAHW <- lm(HeadWidth ~ Habitat + Sex, data = dat)
summary(ANCOVAHW)
HW_emm <- emmeans(ANCOVAHW, ~ Habitat)
summary(HW_emm)
HW_pairwise <- pairs(HW_emm)
summary(HW_pairwise)
HW_emm <- emmeans(ANCOVAHW, ~ Habitat)
HW_emm_df <- as.data.frame(summary(HW_emm))

HW_emm_df$Habitat <- factor(HW_emm_df$Habitat, levels = c("wall", "nonwall"))

HW_box <- ggplot(dat %>%
                    mutate(Habitat = fct_relevel(Habitat, "wall")), 
                  aes(x = as.factor(Habitat), y = HeadWidth, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = HW_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = HW_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "HW Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "HW (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(HW_box)

#HeadDepth
ANCOVAHD <- lm(HeadDepth ~ Habitat + Sex, data = dat)
summary(ANCOVAHD)
HD_emm <- emmeans(ANCOVAHD, ~ Habitat)
summary(HD_emm)
HD_pairwise <- pairs(HD_emm)
summary(HD_pairwise)
HD_emm <- emmeans(ANCOVAHD, ~ Habitat)
HD_emm_df <- as.data.frame(summary(HD_emm))

HD_emm_df$Habitat <- factor(HD_emm_df$Habitat, levels = c("wall", "nonwall"))

HD_box <- ggplot(dat %>%
                   mutate(Habitat = fct_relevel(Habitat, "wall")), 
                 aes(x = as.factor(Habitat), y = HeadDepth, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = HD_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = HD_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "HD Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "HD (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(HD_box)

#HeadLength
ANCOVAHL <- lm(HeadLength ~ Habitat + Sex, data = dat)
summary(ANCOVAHL)
HL_emm <- emmeans(ANCOVAHL, ~ Habitat)
summary(HL_emm)
HL_pairwise <- pairs(HL_emm)
summary(HL_pairwise)
HL_emm <- emmeans(ANCOVAHL, ~ Habitat)
HL_emm_df <- as.data.frame(summary(HL_emm))

HL_emm_df$Habitat <- factor(HL_emm_df$Habitat, levels = c("wall", "nonwall"))

HL_box <- ggplot(dat %>%
                   mutate(Habitat = fct_relevel(Habitat, "wall")), 
                 aes(x = as.factor(Habitat), y = HeadLength, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = HL_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = HL_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "HL Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "HL (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(HL_box)

#FemurLength
ANCOVAFL <- lm(FemurLength ~ Habitat + Sex, data = dat)
summary(ANCOVAFL)
FL_emm <- emmeans(ANCOVAFL, ~ Habitat)
summary(FL_emm)
FL_pairwise <- pairs(FL_emm)
summary(FL_pairwise)
FL_emm <- emmeans(ANCOVAFL, ~ Habitat)
FL_emm_df <- as.data.frame(summary(FL_emm))

FL_emm_df$Habitat <- factor(FL_emm_df$Habitat, levels = c("wall", "nonwall"))

FL_box <- ggplot(dat %>%
                   mutate(Habitat = fct_relevel(Habitat, "wall")), 
                 aes(x = as.factor(Habitat), y = FemurLength, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = FL_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = FL_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "FL Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "FL (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(FL_box)


#TibiaLength
ANCOVATL <- lm(TibiaLength ~ Habitat + Sex, data = dat)
summary(ANCOVATL)
TL_emm <- emmeans(ANCOVATL, ~ Habitat)
summary(TL_emm)
TL_pairwise <- pairs(TL_emm)
summary(TL_pairwise)
TL_emm <- emmeans(ANCOVATL, ~ Habitat)
TL_emm_df <- as.data.frame(summary(TL_emm))

TL_emm_df$Habitat <- factor(TL_emm_df$Habitat, levels = c("wall", "nonwall"))

TL_box <- ggplot(dat %>%
                   mutate(Habitat = fct_relevel(Habitat, "wall")), 
                 aes(x = as.factor(Habitat), y = TibiaLength, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = TL_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = TL_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "TL Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "TL (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(TL_box)

#ForearmLength
ANCOVAFAL <- lm(ForearmLength ~ Habitat + Sex, data = dat)
summary(ANCOVAFAL)
FAL_emm <- emmeans(ANCOVAFAL, ~ Habitat)
summary(FAL_emm)
FAL_pairwise <- pairs(FAL_emm)
summary(FAL_pairwise)
FAL_emm <- emmeans(ANCOVAFAL, ~ Habitat)
FAL_emm_df <- as.data.frame(summary(FAL_emm))

FAL_emm_df$Habitat <- factor(FAL_emm_df$Habitat, levels = c("wall", "nonwall"))

FAL_box <- ggplot(dat %>%
                   mutate(Habitat = fct_relevel(Habitat, "wall")), 
                 aes(x = as.factor(Habitat), y = ForearmLength, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = FAL_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = FAL_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "FAL Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "FAL (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(FAL_box)

#BicepLength
ANCOVABL <- lm(BicepLength ~ Habitat + Sex, data = dat)
summary(ANCOVABL)
BL_emm <- emmeans(ANCOVABL, ~ Habitat)
summary(BL_emm)
BL_pairwise <- pairs(BL_emm)
summary(BL_pairwise)
BL_emm <- emmeans(ANCOVABL, ~ Habitat)
BL_emm_df <- as.data.frame(summary(BL_emm))

BL_emm_df$Habitat <- factor(BL_emm_df$Habitat, levels = c("wall", "nonwall"))

BL_box <- ggplot(dat %>%
                    mutate(Habitat = fct_relevel(Habitat, "wall")), 
                  aes(x = as.factor(Habitat), y = BicepLength, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = BL_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = BL_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "BL Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "BL (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(BL_box)

#LimbRatio
ANCOVALR <- lm(LimbRatio ~ Habitat + Sex, data = dat)
summary(ANCOVALR)
LR_emm <- emmeans(ANCOVALR, ~ Habitat)
summary(LR_emm)
LR_pairwise <- pairs(LR_emm)
summary(LR_pairwise)
LR_emm <- emmeans(ANCOVALR, ~ Habitat)
LR_emm_df <- as.data.frame(summary(LR_emm))

LR_emm_df$Habitat <- factor(LR_emm_df$Habitat, levels = c("wall", "nonwall"))

LR_box <- ggplot(dat %>%
                   mutate(Habitat = fct_relevel(Habitat, "wall")), 
                 aes(x = as.factor(Habitat), y = LimbRatio, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = LR_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = LR_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "LR Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "LR (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(LR_box)


Morph_box <- SVL_box + Mass_box + HW_box + HD_box + HL_box + FL_box + TL_box + FAL_box + BL_box + LR_box + plot_layout(nrow = 2)
#ggsave("Morph_box.pdf", plot = Morph_box, width = 30, height = 20, units = "in")

#CEWL
ANCOVACEWL <- lm(CEWL ~ Habitat + Sex, data = cewl.dat)
summary(ANCOVACEWL)
CEWL_emm <- emmeans(ANCOVACEWL, ~ Habitat)
summary(CEWL_emm)
CEWL_pairwise <- pairs(CEWL_emm)
summary(CEWL_pairwise)
CEWL_emm <- emmeans(ANCOVACEWL, ~ Habitat)
CEWL_emm_df <- as.data.frame(summary(CEWL_emm))

CEWL_emm_df$Habitat <- factor(CEWL_emm_df$Habitat, levels = c("wall", "nonwall"))

CEWL_box <- ggplot(cewl.dat %>%
                   mutate(Habitat = fct_relevel(Habitat, "wall")), 
                 aes(x = as.factor(Habitat), y = CEWL, fill = as.factor(Habitat))) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.4, stroke = 0.5, color = "black") +  
  
  geom_errorbar(data = CEWL_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +  
  
  geom_point(data = CEWL_emm_df, aes(x = Habitat, y = emmean, fill = Habitat), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +  
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  
  scale_shape_manual(values = c(24, 21)) +  
  
  labs(
    title = "CEWL Variation Across Habitats with emmeans",
    x = "Habitat",
    y = "CEWL (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")

print(CEWL_box)


##########################################################################################################
#Model Selection Process for Individual Lizards
#CEWL Model Selection
#Global Model of CEWL:
lme.cewl.global2 <- lme(CEWL ~ Habitat + msmt_temp_C + msmt_VPD_kPa + percent_veg_cover + LogSVL + Sex, random = ~ 1 | Site, data = cewl.dat)
summary(lme.cewl.global2)

# Model selection:
cewl.dredge2 <- dredge(lme.cewl.global2)
subset(cewl.dredge2, delta<4)
top.model2 <- get.models(cewl.dredge2, subset = 1)[[1]]
summary(top.model2) 
anova(top.model2)
ms_error <- top.model2$sigma^2

# ANOVA table
aov_table <- anova(top.model2)

# Compute Sum of Squares for each effect
ss_values <- aov_table[, "F-value"] * ms_error * aov_table[, "numDF"]

# Add to table
aov_table$SumSq <- ss_values
aov_table


# Model averaging:
mod1.avg2 <- model.avg(cewl.dredge2, subset = delta<5)
summary(mod1.avg2, get.models(lme.cewl.global2, subset = TRUE))
confint(mod1.avg2)


# Sum of weights of each variable (similar to relative importance):
sw(mod1.avg2)


#PLot with Separate Habitat Regression Lines
CEWL_VPD_PLot <- ggplot(cewl.dat, aes(x = VPD, y = CEWL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(aes(fill = as.factor(Habitat)), size = 2.5, stroke = 0.5, color = "black", alpha = 0.6) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = Habitat), color = "black", se = FALSE, size = 2) +
  geom_smooth(method = "lm", aes(group = Habitat, color = as.factor(Habitat)), se = FALSE, size = 1.5) + 
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Match colors of trend lines
  labs(
    title = "CEWL is higher in urban habitats with lower VPD",
    x = "VPD (kPa)", 
    y = expression(CEWL ~ (g/m^2*h)) 
  ) +
  scale_shape_manual(values = c(24, 21)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1), # Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 14)
    ) +
  guides(color = "none", fill = "none", shape = "none")

#CEWL vs VPD and bar plot
scatter <- ggplot(cewl.dat, aes(x = VPD, y = CEWL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(aes(fill = as.factor(Habitat)), size = 2.5, stroke = 0.4, color = "black", alpha = 0.5) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = Habitat), color = "black", se = FALSE, size = 2) +
  geom_smooth(method = "lm", aes(group = Habitat, color = as.factor(Habitat)), se = FALSE, size = 1.5) + 
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Match colors of trend lines
  labs(
    title = "CEWL is higher in urban habitats with lower VPD",
    x = "VPD (kPa)", 
    y = expression(CEWL ~ (g/m^2*h)) 
  ) +
  scale_shape_manual(values = c(24, 21)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1), # Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 4),
    legend.key.size= unit(0.5, "cm"),
    legend.text = element_text(size = 3)
  )


box <- ggplot(cewl.dat %>%
                mutate(Habitat = fct_relevel(Habitat, "wall")), aes(x = as.factor(Habitat), y = CEWL, fill = as.factor(Habitat))) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.5, size = 1) +  # One box per habitat
  geom_jitter(aes(shape = as.factor(Sex), fill = as.factor(Habitat)), 
              width = 0.2, size = 2.5, alpha = 0.7, stroke = 0.5, color = "black") +  # Points with different shapes for sex
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c(24, 21)) +  # Custom shapes: Triangle for one sex, Circle for another
  labs(
    title = "CEWL Variation Across Habitats",
    x = "Habitat",
    y = expression(CEWL ~ (g/m^2*h)),
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14)
  ) +
  guides(fill = "none", shape = "none")  # Remove the legend for fill and shape from box plot



combined_CEWL_Plot <- ( CEWL_box | scatter)
#ggsave("combined_CEWL_Plot.pdf", width = 12, height = 8, units = "in", dpi = 300)



#Microclimatic variables and their impact on CEWL in PCA
Habitat_PCA <- cewl.dat %>%
  dplyr::select(ambient_temp, percent_veg_cover, VPD)
scaled_Habitat_PCA <- scale(Habitat_PCA)

habitat_pca_result <- prcomp(scaled_Habitat_PCA, center = TRUE, scale. = TRUE)
summary(habitat_pca_result)

cewl.dat.with_pcs <- cewl.dat %>%
  mutate(PC1 = habitat_pca_result$x[, 1],
         PC2 = habitat_pca_result$x[, 2],
         PC3 = habitat_pca_result$x[, 3],
        )

CEWL_PCA <- lm(CEWL ~ PC1 + PC2 + PC3, data = cewl.dat.with_pcs)
summary(CEWL_PCA)

habitat_pca_result$rotation

pc_coefficients <- summary(CEWL_PCA)$coefficients[2:4, 1] 
variable_importance <- habitat_pca_result$rotation %*% pc_coefficients
colnames(variable_importance) <- "Estimated_Impact_on_CEWL"
print(variable_importance)

#Plot CEWL PCA
pca_scores <- as_tibble(habitat_pca_result$x) %>%
  bind_cols(cewl.dat %>% dplyr::select(Habitat))  # Optional: add grouping var like Habitat

# Loadings for PC1 and PC3
loadings <- as.data.frame(habitat_pca_result$rotation[, c(1, 3)])
loadings$Variable <- rownames(loadings)

# Plot PC1 vs PC3
ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC3, color = Habitat), size = 3, alpha = 0.8) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * 2, yend = PC3 * 2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black") +
  geom_text_repel(data = loadings,
                  aes(x = PC1 * 2.2, y = PC3 * 2.2, label = Variable),
                  size = 4) +
  labs(title = "PCA of Habitat Variables: PC1 vs PC3",
       x = paste0("PC1 (", round(summary(habitat_pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC3 (", round(summary(habitat_pca_result)$importance[2, 3] * 100, 1), "% variance)")) +
  theme_minimal() +
  coord_equal()



#Seperate CEWL by habitat and include environmental variables instead of habitat type
#CEWL~ Urban
Urban.CEWL.VPD_dat <- subset(cewl.dat, Habitat == "wall")

lm.Urban.cewl <- lm(CEWL ~ VPD * ambient_temp + percent_veg_cover + LogSVL + Mass + Sex + msmt_temp_C + msmt_VPD_kPa, data = Urban.CEWL.VPD_dat, na.action = "na.fail")
summary(lm.Urban.cewl)
anova(lm.Urban.cewl)

#VPD
Urban.CEWL.VPD <- lm(CEWL ~ VPD, data = Urban.CEWL.VPD_dat)
summary(Urban.CEWL.VPD)

# Model selection:
cewl.dredge2 <- dredge(lm.Urban.cewl)
subset(cewl.dredge2, delta<4)
top.model2 <- get.models(cewl.dredge2, subset = 1)[[1]]
summary(top.model2) 
Anova(top.model2, type = "2")

# Model averaging:
mod1.avg2 <- model.avg(cewl.dredge2, subset = delta<5)
summary(mod1.avg2, get.models(lm.Urban.cewl, subset = TRUE))
confint(mod1.avg2)


# Sum of weights of each variable (similar to relative importance):
sw(mod1.avg2)



#CEWL~ Nonurban
Nonurban.CEWL.VPD_dat <- subset(cewl.dat, Habitat == "nonwall")

lm.Nonurban.cewl <- lm(CEWL ~ VPD * ambient_temp + percent_veg_cover + LogSVL + Mass + Sex + msmt_temp_C + msmt_VPD_kPa, data = Nonurban.CEWL.VPD_dat, na.action = "na.fail")
summary(lm.Nonurban.cewl)
anova(lm.Nonurban.cewl)

#VPD
Nonurban.CEWL.VPD <- lm(CEWL ~ VPD, data = Nonurban.CEWL.VPD_dat)
summary(Nonurban.CEWL.VPD)

# Model selection:
cewl.dredge2 <- dredge(lm.Nonurban.cewl)
subset(cewl.dredge2, delta<4)
top.model2 <- get.models(cewl.dredge2, subset = 1)[[1]]
summary(top.model2) 
Anova(top.model2, type = "2")

# Model averaging:
mod1.avg2 <- model.avg(cewl.dredge2, subset = delta<5)
summary(mod1.avg2, get.models(lm.Nonurban.cewl, subset = TRUE))
confint(mod1.avg2)


# Sum of weights of each variable (similar to relative importance):
sw(mod1.avg2)



##########################################################################################################
#Analyses to test for differences in ind. morphology in PCA space- scaled to SVL and unscaled
#Standardize morphology to body size
morph_dat <- dat

morph_dat$ScaleMass <- residuals(lm(log(morph_dat$Mass) ~ log(morph_dat$SVL), na.action = na.exclude))

morph_dat$ScaleHW <- residuals(lm(log(morph_dat$HeadWidth) ~ log(morph_dat$SVL), na.action = na.exclude))

morph_dat$ScaleHD <- residuals(lm(log(morph_dat$HeadDepth) ~ log(morph_dat$SVL), na.action = na.exclude))

morph_dat$ScaleHL <- residuals(lm(log(morph_dat$HeadLength) ~ log(morph_dat$SVL), na.action = na.exclude))

morph_dat$ScaleFL <- residuals(lm(log(morph_dat$FemurLength) ~ log(morph_dat$SVL), na.action = na.exclude))

morph_dat$ScaleTL <- residuals(lm(log(morph_dat$TibiaLength) ~ log(morph_dat$SVL), na.action = na.exclude))

morph_dat$ScaleBL <- residuals(lm(log(morph_dat$BicepLength) ~ log(morph_dat$SVL), na.action = na.exclude))

morph_dat$ScaleFAL <- residuals(lm(log(morph_dat$ForearmLength) ~ log(morph_dat$SVL), na.action = na.exclude))

morph_dat$ScaleLR <- residuals(lm(log(morph_dat$LimbRatio) ~ log(morph_dat$SVL), na.action = na.exclude))


morph_dat <- subset(morph_dat, select = -c(BagID, Gravid, Morph, Scars, Ticks, Mites, Mass, SVL, HeadWidth, HeadDepth, HeadLength, FemurLength, TibiaLength, BicepLength, ForearmLength, CEWLbt, CEWL, msmt_temp_C, msmt_RH_percent, msmt_e_a_kPa, msmt_VPD_kPa, ambient_percent_rh, ambient_temp, percent_veg_cover, LogMass, LogSVL, LogHW, LogHD, LogHL, LogFL, LogTL, LogBL, LogFAL, LimbRatio, LogLR, SVLZERO))


#Mass
KikiANCOVAMass <- lm(ScaleMass ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVAMass)
anova_Mass <-anova(KikiANCOVAMass)

#Head Width
KikiANCOVAHW <- lm(ScaleHW ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVAHW)
anova_HW <-anova(KikiANCOVAHW)

#Head Depth
KikiANCOVAHD <- lm(ScaleHD ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVAHD)
anova_HD <-anova(KikiANCOVAHD)

#Head Length
KikiANCOVAHL <- lm(ScaleHL ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVAHL) 
anova_HL <-anova(KikiANCOVAHL)

#Femur Length
KikiANCOVAFL <- lm(ScaleFL ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVAFL) 
anova_FL <- anova(KikiANCOVAFL)

#Tibia Length
KikiANCOVATL <- lm(ScaleTL ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVATL) 
anova_TL <- anova(KikiANCOVATL)

#Forearm Length
KikiANCOVAFAL <- lm(ScaleFAL ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVAFAL) 
anova_FAL <- anova(KikiANCOVAFAL)

#Bicep Length
KikiANCOVABL <- lm(ScaleBL ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVABL) 
anova_BL <- anova(KikiANCOVABL)

#Limb Ratio
KikiANCOVALR <- lm(ScaleLR ~ Habitat * Density + Sex, data = morph_dat)
summary(KikiANCOVALR) 
anova_LR <- anova(KikiANCOVALR)


#PCA
#Scaled
#morph_dat <- subset(morph_dat, select = -c(BagID, Gravid, Morph, Scars, Ticks, Mites, Mass, SVL, HeadWidth, HeadDepth, HeadLength, FemurLength, TibiaLength, BicepLength, ForearmLength, CEWLbt, CEWL, msmt_temp_C, msmt_RH_percent, msmt_e_a_kPa, msmt_VPD_kPa, ambient_percent_rh, ambient_temp, percent_veg_cover, LogMass, LogSVL, LogHW, LogHD, LogHL, LogFL, LogTL, LogBL, LogFAL, LimbRatio, LogLR, SVLZERO))

scaled_morph_dat <- subset(morph_dat, select = -c(Island, Site, Habitat, Density, Date, Sex, es, ea, VPD))
scaled_morph_dat <- na.omit(scaled_morph_dat)

scaled_pca <- prcomp(scaled_morph_dat, scale = TRUE)

cleaned_dat <- dat[rownames(scaled_morph_dat), ]

# Create PCA data frame
scaled_pca_data <- data.frame(scaled_pca$x, Habitat = cleaned_dat$Habitat, Sex = cleaned_dat$Sex)

PCA1 <- ggplot(scaled_pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Habitat, color = "black", shape = Sex), size = 3, stroke = 0.5, alpha = 0.6) +  # Black outline, fill by Habitat
  stat_ellipse(aes(color = Habitat), level = 0.95, linewidth = 1) +  # Add ellipses colored by Habitat
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Ellipse colors
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Fill colors for points
  scale_shape_manual(values = c(24, 21)) +  # Shape by Sex (triangle for females, circle for males)
  theme_minimal() +
  labs(title = "PCA of Morphological Traits Scaled to SVL",  # Add title
       x = "PC1", 
       y = "PC2") +
  annotate("text", x = min(scaled_pca_data$PC1), y = max(scaled_pca_data$PC2),  # Add text annotation
           label = paste0("PERMANOVA p = 0.836"), 
           hjust = 0, size = 5) +
  theme(legend.position = "top", 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

# Perform PERMANOVA
pca_scores <- as.data.frame(scaled_pca$x)  # Get PCA scores
permanova_result <- adonis2(pca_scores[, 1:2] ~ cleaned_dat$Habitat, method = "euclidean")
print(permanova_result)

#ggsave("scaled_pca_plot.svg", width = 8, height = 6, units = "in", dpi = 300)



#Not Scaled
unscaled_morph_dat <- subset(dat, select = c(LogMass, LogSVL, LogHW, LogHD, LogHL, LogFL, LogTL, LogBL, LogFAL, LogLR))
unscaled_morph_dat <- na.omit(unscaled_morph_dat)

unscaled_pca <- prcomp(unscaled_morph_dat, scale = TRUE)

cleaned_dat <- dat[rownames(unscaled_morph_dat), ]

# Create PCA data frame
unscaled_pca_data <- data.frame(unscaled_pca$x, Habitat = cleaned_dat$Habitat, Sex = cleaned_dat$Sex)

PCA2 <- ggplot(unscaled_pca_data, aes(x = PC1, y = PC2, fill = Habitat)) +
  geom_point(aes(fill = Habitat, color = "black", shape = Sex), size = 3, stroke = 0.5, alpha = 0.6) +  # Black outline, fill by Habitat
  stat_ellipse(aes(color = Habitat), level = 0.95, linewidth = 1) +  # Color ellipses by Habitat
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Ellipse colors
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Fill colors for points
  scale_shape_manual(values = c(24, 21)) +  # Shape by Sex (triangle for females, circle for males)
  theme_minimal() +
  labs(title = "PCA of Unscaled Morphological Traits",
       x = "PC1",
       y = "PC2") +
  annotate("text", x = min(unscaled_pca_data$PC1), y = max(unscaled_pca_data$PC2),
           label = paste0("PERMANOVA p = 0.033"), 
           hjust = 0, size = 5) +
  theme(legend.position = "top", 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

#ggsave("unscaled_pca_plot.svg", width = 8, height = 6, units = "in", dpi = 300)


pca_scores <- as.data.frame(unscaled_pca$x)  # Get PCA scores
permanova_result <- adonis2(pca_scores[, 1:2] ~ cleaned_dat$Habitat, method = "euclidean")
print(permanova_result)

unscaled_pca$rotation[, 1:2]

summary(unscaled_pca)


combined_PCA_Plot <- (PCA1 | PCA2)
#ggsave("combined_PCA_Plot.svg", width = 8, height = 6, units = "in", dpi = 300)


##################################################################################
#Morphological Model Selection
#Global Model of Morphology:
#SVL
lm.SVL.global <- lm(LogSVL ~ Habitat * Sex + Density, data = dat, na.action = "na.fail")
summary(lm.SVL.global)

# Model selection:
options(na.action = "na.fail")
SVL.dredge <- dredge(lm.SVL.global)
subset(SVL.dredge, delta<4)
SVLtop.model <- get.models(SVL.dredge, subset = 1)[[1]]
summary(SVLtop.model) #Top model SVL ~ | Site
Anova(SVLtop.model, type = 2)

# Model averaging:
SVLmod.avg <- model.avg(SVL.dredge, subset = delta<5)
summary(SVLmod.avg, get.models(lm.SVL.global, subset = TRUE))
confint(SVLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(SVLmod.avg)

SVL_vs_Density.svg <- ggplot(dat, aes(x = Density, y = SVL, fill = as.factor(Habitat), color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2.5, stroke = 0.5, color = "black", alpha = 0.6) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = Habitat), se = FALSE, size = 2, color = "black") +  # Black outline for trend lines
  geom_smooth(method = "lm", aes(group = Habitat, color = as.factor(Habitat)), se = FALSE, linewidth = 1.5) +  # Trend line with Habitat color
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Fill colors by Habitat
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Line colors by Habitat
  scale_shape_manual(values = c(24, 21)) +  # Shapes by Sex
  labs(title = "Body size increases with population density",
       x = "Density (N lizards/100 m transect)", y = "SVL (mm)") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 1), # Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 14)
  )
#ggsave("~/Desktop/SVL_vs_Density.svg", width = 8, height = 6, dpi = 300)

SVL_vs_Density2.svg <- ggplot(dat, aes(x = Density, y = SVL, fill = as.factor(Habitat), 
                                       color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2.5, stroke = 0.5, color = "black", alpha = 0.6) +  # Scatter plot points
  # Black outline trend lines
  geom_smooth(method = "lm", aes(group = interaction(Habitat, Sex)), 
              se = FALSE, size = 2, color = "black") +  # Black outline for trend lines
  # Colored trend lines with different linetypes for Sex
  geom_smooth(method = "lm", aes(group = interaction(Habitat, Sex), linetype = as.factor(Sex)), 
              se = FALSE, size = 1.5) +  # Trend lines by Habitat & Sex with different linetypes
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Fill colors by Habitat
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Line colors by Habitat
  scale_shape_manual(values = c(24, 21)) +  # Shapes by Sex
  scale_linetype_manual(values = c("M" = "solid", "F" = "dotted")) +  # Line types for males and females
  labs(title = "Body size increases with population density",
       x = "Density (N lizards/100 m transect)", y = "SVL (mm)") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),  # Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 1),  # Adds black tick marks
    axis.ticks.length = unit(0.1, "cm")  # Adjusts length of tick marks
  )
#ggsave("~/Desktop/SVL_vs_Density2.svg", width = 8, height = 6, dpi = 300)
  
  
#Head Width
lm.HW.global <- lm(LogHW ~ Habitat + Density + LogSVL + Sex, data = dat, na.action = "na.fail")
summary(lm.HW.global)

# Model selection:
HW.dredge <- dredge(lm.HW.global)
subset(HW.dredge, delta<4)
HWtop.model <- get.models(HW.dredge, subset = 1)[[1]]
summary(HWtop.model) #Top model LogHW ~ SVL | Site
Anova(HWtop.model, type = 2)

# Model averaging:
HWmod.avg <- model.avg(HW.dredge, subset = delta<5)
summary(HWmod.avg, get.models(lm.HW.global, subset = TRUE))
confint(HWmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(HWmod.avg)


ggplot(dat, aes(x = Density, y = SVL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +  # Add regression line
  labs(title = "Scatter Plot of SVL vs Population Density",
       x = "Population Density", y = "SVL") +
  theme_minimal()


#Head Depth
lm.HD.global <- lm(LogHD ~ Habitat + Density + LogSVL + Sex, data = dat, na.action = "na.fail")
summary(lm.HD.global)

# Model selection:
HD.dredge <- dredge(lm.HD.global)
subset(HD.dredge, delta<4)
HDtop.model <- get.models(HD.dredge, subset = 1)[[1]]
summary(HDtop.model) #Top model LogHD ~ SVL | Site
Anova(HDtop.model, type = 2)

# Model averaging:
HDmod.avg <- model.avg(HD.dredge, subset = delta<5)
summary(HDmod.avg, get.models(lm.HD.global, subset = TRUE))
confint(HDmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(HDmod.avg)

ggplot(dat, aes(x = HeadDepth, y = SVL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +  # Add regression line
  labs(title = "Scatter Plot of SVL vs Head Depth",
       x = "HeadDepth", y = "SVL") +
  theme_minimal()


#Head Length
lm.HL.global <- lm(LogHL ~ Habitat + Density + LogSVL + Sex, data = dat, na.action = "na.fail")
summary(lm.HL.global)

# Model selection:
HL.dredge <- dredge(lm.HL.global)
subset(HL.dredge, delta<4)
HLtop.model <- get.models(HL.dredge, subset = 1)[[1]]
summary(HLtop.model) #Top model LogHL ~ SVL | Site
Anova(HLtop.model, type = 2)

# Model averaging:
HLmod.avg <- model.avg(HL.dredge, subset = delta<5)
summary(HLmod.avg, get.models(lm.HL.global, subset = TRUE))
confint(HLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(HLmod.avg)

ggplot(dat, aes(x = HeadLength, y = SVL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +  # Add regression line
  labs(title = "Scatter Plot of SVL vs Head Length",
       x = "HeadLength", y = "SVL") +
  theme_minimal()



#Femur Length
lm.FL.global <- lm(LogFL ~ Habitat + Density + LogSVL + Sex, data = dat, na.action = "na.fail")
summary(lm.FL.global)

# Model selection:
FL.dredge <- dredge(lm.FL.global)
subset(FL.dredge, delta<4)
FLtop.model <- get.models(FL.dredge, subset = 1)[[1]]
summary(FLtop.model) #Top model LogFL ~ SVL | Site
Anova(FLtop.model, type = 2)

# Model averaging:
FLmod.avg <- model.avg(FL.dredge, subset = delta<5)
summary(FLmod.avg, get.models(lm.FL.global, subset = TRUE))
confint(FLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(FLmod.avg)

ggplot(dat, aes(x = FemurLength, y = SVL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +  # Add regression line
  labs(title = "Scatter Plot of SVL vs Femur Length",
       x = "FemurLength", y = "SVL") +
  theme_minimal()



#Tibia Length
lm.TL.global <- lm(LogTL ~ Habitat + Density + LogSVL + Sex, data = dat, na.action = "na.fail")
summary(lm.TL.global)

# Model selection:
TL.dredge <- dredge(lm.TL.global)
subset(TL.dredge, delta<4)
TLtop.model <- get.models(TL.dredge, subset = 1)[[1]]
summary(TLtop.model) #Top model LogTL ~ SVL | Site
Anova(TLtop.model, type = 2)

# Model averaging:
TLmod.avg <- model.avg(TL.dredge, subset = delta<5)
summary(TLmod.avg, get.models(lm.TL.global, subset = TRUE))
confint(TLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(TLmod.avg)

ggplot(dat, aes(x = TibiaLength, y = SVL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +  # Add regression line
  labs(title = "Scatter Plot of SVL vs Tibia Length",
       x = "TibiaLength", y = "SVL") +
  theme_minimal()



#Bicep Length
lm.BL.global <- lm(LogBL ~ Habitat + Density + LogSVL + Sex, data = dat, na.action = "na.fail")
summary(lm.BL.global)

# Model selection:
BL.dredge <- dredge(lm.BL.global)
subset(BL.dredge, delta<4)
BLtop.model <- get.models(BL.dredge, subset = 1)[[1]]
summary(BLtop.model) #Top model LogBL ~ SVL | Site
Anova(BLtop.model, type = 2)

# Model averaging:
BLmod.avg <- model.avg(BL.dredge, subset = delta<5)
summary(BLmod.avg, get.models(lm.BL.global, subset = TRUE))
confint(BLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(BLmod.avg)

ggplot(dat, aes(x = BicepLength, y = SVL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +  # Add regression line
  labs(title = "Scatter Plot of SVL vs Bicep Length",
       x = "BicepLength", y = "SVL") +
  theme_minimal()



#Forearm Length
lm.FAL.global <- lm(LogFAL ~ Habitat + Density + LogSVL + Sex, data = dat, na.action = "na.fail")
summary(lm.FAL.global)

# Model selection:
FAL.dredge <- dredge(lm.FAL.global)
subset(FAL.dredge, delta<4)
FALtop.model <- get.models(FAL.dredge, subset = 1)[[1]]
summary(FALtop.model) #Top model LogFAL ~ SVL | Site
Anova(FALtop.model, type = 2)

# Model averaging:
FALmod.avg <- model.avg(FAL.dredge, subset = delta<5)
summary(FALmod.avg, get.models(lm.FAL.global, subset = TRUE))
confint(FALmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(FALmod.avg)

ggplot(dat, aes(x = ForearmLength, y = SVL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +  # Add regression line
  labs(title = "Scatter Plot of SVL vs Forearm Length",
       x = "ForearmLength", y = "SVL") +
  theme_minimal()


#Limb Ratio
lm.LR.global <- lm(LogLR ~ Habitat + Density + LogSVL + Sex, data = dat, na.action = "na.fail")
summary(lm.LR.global)

# Model selection:
LR.dredge <- dredge(lm.LR.global)
subset(LR.dredge, delta<4)
LRtop.model <- get.models(LR.dredge, subset = 1)[[1]]
summary(LRtop.model) #Top model LogLR ~ SVL | Site
Anova(LRtop.model, type = 2)

# Model averaging:
LRmod.avg <- model.avg(LR.dredge, subset = delta<5)
summary(LRmod.avg, get.models(lm.LR.global, subset = TRUE))
confint(LRmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(LRmod.avg)

ggplot(dat, aes(x = LimbRatio, y = SVL, color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2) +   # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +  # Add regression line
  labs(title = "Scatter Plot of SVL vs Limb Ratio",
       x = "LimbRatio", y = "SVL") +
  theme_minimal()


#Correlation Values
variables <- c("LogHW", "LogHD", "LogHL", "LogFL", "LogTL", "LogBL", "LogFAL", "LogLR")

# Calculate correlations and store in a data frame
correlations <- data.frame(
  Variable = variables,
  Correlation = sapply(variables, function(var) cor(dat$LogSVL, dat[[var]]))
)

# Display the results
print(correlations)



##########################################################################################################
#Testing for collinearity between habitat type and microclimatic variables

#With Habitat
#Ambient Temperature
temp_Habitat_Model <- glm(ambient_temp ~ Habitat, data = dat)
summary(temp_Habitat_Model)

p1 <- ggplot(dat, aes(x = ambient_temp, fill = Habitat)) +
  geom_density(alpha = 0.4) + 
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Adjust transparency for overlap
  labs(title = "Ambient Temperature by Habitat",
       x = "Ambient Temperature", y = "Density") +
  theme_classic()+
  theme(legend.position = "none")
#ggsave("~/Desktop/temp_vs_Habitat.svg", width = 8, height = 6, dpi = 300)

#Ambient Percent Humidity
humid_Habitat_Model <- glm(ambient_percent_rh ~ Habitat, data = dat)
summary(humid_Habitat_Model)

p2 <- ggplot(dat, aes(x = ambient_percent_rh, fill = Habitat)) +
  geom_density(alpha = 0.4) + # Adjust transparency for overlap
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  labs(title = "Ambient Humidity by Habitat",
       x = "Ambient Relative Humidity", y = "Density") +
  theme_classic() +
theme(legend.position = "none")
#ggsave("~/Desktop/humid_vs_Habitat.svg", width = 8, height = 6, dpi = 300)

#VPD
VPD_Habitat_Model <- glm(VPD ~ Habitat, data = dat)
summary(VPD_Habitat_Model)

p3 <- ggplot(cewl.dat, aes(x = VPD, fill = Habitat)) +
  geom_density(alpha = 0.4) +  # Adjust transparency for overlap
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  labs(title = "Vapor Pressue Deficit (VPD) by Habitat",
       x = "VPD", y = "Density") +
  theme_classic() +
theme(legend.position = "none")
#ggsave("~/Desktop/VPD_vs_Habitat.svg", width = 8, height = 6, dpi = 300)

#Percent Vegetation Cover
Veg_Habitat_Model <- glm(percent_veg_cover ~ Habitat, data = dat)
summary(Veg_Habitat_Model)

p4 <- ggplot(dat, aes(x = percent_veg_cover, fill = Habitat)) +
  geom_density(alpha = 0.4) +  # Adjust transparency for overlap
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  labs(title = "Vegetation Cover by Habitat",
       x = "Percent Vegetation Cover", y = "Density") +
  theme_classic() +
theme(legend.position = "none")
#ggsave("~/Desktop/Veg_vs_Habitat.svg", width = 8, height = 6, dpi = 300)


combined_Enviro_Plot <- (p1 | p3 | p4)
#ggsave("~/Desktop/combined_Enviro_Plot.svg", plot = combined_Enviro_Plot, width = 10, height = 8, dpi = 300)

#Vegetation Cover Collinearity with Other Climate Variables
#Ambient Temp
cor(dat$percent_veg_cover, dat$ambient_temp, method ="pearson")
ggplot(dat, aes(x = percent_veg_cover, y = ambient_temp, color = as.factor(Habitat))) +
  geom_point(aes(fill = as.factor(Habitat)), size = 2.5, stroke = 0.5, shape = 21, color = "black", alpha = 0.6) +  # Remove color = "black"
  geom_smooth(method = "lm", aes(group = Habitat), se = FALSE, size = 2, color = "black") +  # Black outline for trend lines
  geom_smooth(method = "lm", aes(group = Habitat), se = FALSE, size = 1.5) +   # Main trend line
  labs(title = "Temperature increases with vegetation cover",
       x = "Vegetation Cover (%)", y = "Temperature (°C)") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 14)
  )
#ggsave("~/Desktop/Veg_vs_Temp.svg", width = 8, height = 6, dpi = 300)

#Ambient relative Humidity
cor(dat$percent_veg_cover, dat$ambient_percent_rh, method ="pearson")
ggplot(dat, aes(x = percent_veg_cover, y = ambient_percent_rh, color = as.factor(Habitat))) +
  geom_point(aes(fill = as.factor(Habitat)), size = 2.5, stroke = 0.5, shape = 21, color = "black", alpha = 0.6) +  # Remove color = "black"
  geom_smooth(method = "lm", aes(group = Habitat), se = FALSE, size = 2, color = "black") +  # Black outline for trend lines
  geom_smooth(method = "lm", aes(group = Habitat), se = FALSE, size = 1.5) +   # Main trend line
  labs(title = "Humidity decreases with vegetation cover",
       x = "Vegetation Cover (%)", y = "Ambient Humidity (%)") +
  theme_minimal()+
  theme(
    axis.line = element_line(color = "black", size = 1),# Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 14)
  )
#ggsave("~/Desktop/Veg_vs_RH.svg", width = 8, height = 6, dpi = 300)

#VPD
cor(cewl.dat$percent_veg_cover, cewl.dat$VPD, method ="pearson")
ggplot(cewl.dat, aes(x = percent_veg_cover, y = VPD, color = as.factor(Habitat))) +
  geom_point(aes(fill = as.factor(Habitat)), size = 2.5, stroke = 0.5, shape = 21, color = "black", alpha = 0.6) +  # Remove color = "black"
  geom_smooth(method = "lm", aes(group = Habitat), se = FALSE, size = 2, color = "black") +  # Black outline for trend lines
  geom_smooth(method = "lm", aes(group = Habitat), se = FALSE, size = 1.5) +   # Main trend line
  labs(title = "VPD increases with vegetation cover",
       x = "Vegetation Cover (%)", y = "VPD (kPa)") +
  theme_minimal()+
  theme(
    axis.line = element_line(color = "black", size = 1),# Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 14)
  )
#ggsave("~/Desktop/Veg_vs_VPD.svg", width = 8, height = 6, dpi = 300)

#GLMs
#CEWL~ Habitat
lm.cewl.global.Habitat <- lm(CEWL ~ Habitat + percent_veg_cover + LogSVL + Mass + Sex + msmt_temp_C + msmt_VPD_kPa, data = cewl.dat, na.action = "na.fail")

summary(lm.cewl.global.Habitat)

# Model selection:
cewl.dredge2 <- dredge(lm.cewl.global.Habitat)
subset(cewl.dredge2, delta<4)
top.model2 <- get.models(cewl.dredge2, subset = 1)[[1]]
summary(top.model2) 
Anova(top.model2, type = "2")

# Model averaging:
mod1.avg2 <- model.avg(cewl.dredge2, subset = delta<5)
summary(mod1.avg2, get.models(lme.cewl.global2, subset = TRUE))
confint(mod1.avg2)


# Sum of weights of each variable (similar to relative importance):
sw(mod1.avg2)

#Testing importance of random variable
#lm.cewl.global.Habitat <- lm(CEWL ~ Habitat + ambient_temp + percent_veg_cover + LogSVL + Sex, data = cewl.dat)
#AIC(lm.cewl.global.Habitat, lme.cewl.global.Habitat)


#MuMIn::r.squaredGLMM(lme.cewl.global.Habitat)


##########################################################################################################
#Gathering Summary Measures of Temp/RH/VPD @ time of CEWL measurement 
#Overall
#Temp 
skim(dat$ambient_temp)

#Relative Humidity
skim(dat$ambient_percent_rh)

#VPD
skim(cewl.dat$VPD)

#Temp 
skim(dat$msmt_temp_C)

#Relative Humidity
skim(dat$msmt_RH_percent)

#VPD
skim(cewl.dat$msmt_VPD_kPa)


#Urban
#Temp 
skim(Urban.CEWL.VPD$ambient_temp)

#Relative Humidity
skim(Urban.CEWL.VPD$ambient_percent_rh)

#VPD
skim(Urban.CEWL.VPD$VPD)

#Nonurban
#Temp 
skim(Nonurban.CEWL.VPD$ambient_temp)

#Relative Humidity
skim(Nonurban.CEWL.VPD$ambient_percent_rh)

#VPD
skim(Nonurban.CEWL.VPD$VPD)

#Habitat Differences in ENV variables
temp <- glm(dat$ambient_temp ~ dat$Habitat)
summary(temp)

vpd <- glm(dat$VPD ~ dat$Habitat)
summary(vpd)

humid <- glm(dat$ambient_percent_rh ~ dat$Habitat)
summary(humid)

veg <- glm(dat$percent_veg_cover ~ dat$Habitat)
summary(veg)

veg2 <- glm(dat$percent_veg_cover ~ dat$Habitat + dat$ambient_temp)
summary(veg2)


#####################################################################################################################################################################################################
#For Kinsey- SVL x CEWL
LogSVL_vs_CEWL.svg <- ggplot(cewl.dat, aes(x = LogSVL, y = CEWL, fill = as.factor(Habitat), color = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2.5, stroke = 0.5, color = "black", alpha = 0.6) +    #Scatter plot points
  geom_smooth(method = "lm", aes(group = Habitat), se = FALSE, size = 2, color = "black") +   #Black outline for trend lines
  geom_smooth(method = "lm", aes(group = Habitat, color = as.factor(Habitat)), se = FALSE, size = 1.5) +   #Trend line with Habitat color
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Fill colors by Habitat
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Line colors by Habitat
  scale_shape_manual(values = c(24, 21)) +  # Shapes by Sex
  labs(title = "CEWL vs LogSVL",
       x = "LogSVL", y = "CEWL") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1), # Solid black axis lines
    panel.grid = element_blank(),
   text = element_text(size = 14)
  )

LogSVL_vs_CEWL2 <- ggplot(cewl.dat, aes(x = LogSVL, y = CEWL, fill = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2.5, stroke = 0.5, color = "black", alpha = 0.6) +  # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", size = 2, se = FALSE) +  # Single trend line for all data
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Fill colors by Habitat
  scale_shape_manual(values = c(24, 21)) +  # Shapes by Sex
  labs(title = "CEWL vs LogSVL",
       x = "LogSVL", y = "CEWL") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 14)
  )

#ggsave("LogSVL_vs_CEWL.svg", width = 12, height = 8, units = "in", dpi = 300)


SVL_vs_CEWL2 <- ggplot(cewl.dat, aes(x = SVL, y = CEWL, fill = as.factor(Habitat), shape = as.factor(Sex))) +
  geom_point(size = 2.5, stroke = 0.5, color = "black", alpha = 0.6) +  # Scatter plot points
  geom_smooth(method = "lm", aes(group = 1), color = "black", size = 2, se = FALSE) +  # Single trend line for all data
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +  # Fill colors by Habitat
  scale_shape_manual(values = c(24, 21)) +  # Shapes by Sex
  labs(title = "CEWL vs SVL",
       x = "SVL", y = "CEWL") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Solid black axis lines
    panel.grid = element_blank(),
    text = element_text(size = 14)
  )












################################################################################################################################################################################################################################
#WorldClim Analyses to ensure macroclimate and habitat type are not confounded

location_dat <- read.csv(("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/data-location.csv"))
head(location_dat)
coords <- data.frame(lat=location_dat[,3], long=location_dat[,4])

#Global Human  Modification index
GHM <- raster("~/Dropbox/SDSU/Brock_Lab/Dissertation/Comparitive_Analyses/gHM/gHM.tif")

coords_df <- data.frame(long = location_dat[,4], lat = location_dat[,3])

coords_sf <- st_as_sf(coords_df, coords = c("long", "lat"), crs = 4326)

coords_sf <- st_transform(coords_sf, crs(GHM))

GHM_values <- raster::extract(GHM, coords_sf)
location_dat$GHM <- GHM_values


#WorldClim Variables
#Precipitation
precip_1 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_01.tif")
precip_2 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_02.tif")
precip_3 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_03.tif")
precip_4 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_04.tif")
precip_5 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_05.tif")
precip_6 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_06.tif")
precip_7 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_07.tif")
precip_8 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_08.tif")
precip_9 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_09.tif")
precip_10 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_10.tif")
precip_11 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_11.tif")
precip_12 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_12.tif")

coordinates(coords) <- ~long+lat
crs(coords) <- crs(raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/wc2.1_30s_prec_01.tif"))

precip_files <- list.files("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/prec/", 
                           pattern = "wc2.1_30s_prec_.*\\.tif$", full.names = TRUE)
precip_stack <- stack(precip_files)

precip_values <- extract(precip_stack, coords, method = "bilinear")
annual_precip <- rowSums(precip_values, na.rm = TRUE)
location_dat$annual_precip <- annual_precip
monthly_avg_precip <- rowMeans(precip_values, na.rm = TRUE)
location_dat$monthly_avg_precip <- monthly_avg_precip

#Solar Radiation
srad_1 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_01.tif")
srad_2 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_02.tif")
srad_3 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_03.tif")
srad_4 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_04.tif")
srad_5 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_05.tif")
srad_6 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_06.tif")
srad_7 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_07.tif")
srad_8 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_08.tif")
srad_9 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_09.tif")
srad_10 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_10.tif")
srad_11 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_11.tif")
srad_12 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_12.tif")

crs(coords) <- crs(raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/wc2.1_30s_srad_01.tif"))

srad_files <- list.files("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/srad/", 
                         pattern = "wc2.1_30s_srad_.*\\.tif$", full.names = TRUE)
srad_stack <- stack(srad_files)

srad_values <- extract(srad_stack, coords, method = "bilinear")

location_dat[paste0("srad_", 1:12)] <- srad_values

monthly_avg_srad <- rowMeans(srad_values, na.rm = TRUE)
location_dat$monthly_avg_srad <- monthly_avg_srad


#Water Vapor Pressure
vapr_1 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_01.tif")
vapr_2 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_02.tif")
vapr_3 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_03.tif")
vapr_4 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_04.tif")
vapr_5 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_05.tif")
vapr_6 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_06.tif")
vapr_7 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_07.tif")
vapr_8 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_08.tif")
vapr_9 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_09.tif")
vapr_10 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_10.tif")
vapr_11 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_11.tif")
vapr_12 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_12.tif")

crs(coords) <- crs(raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/wc2.1_30s_vapr_01.tif"))

vapr_files <- list.files("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/vapr/", 
                         pattern = "wc2.1_30s_vapr_.*\\.tif$", full.names = TRUE)
vapr_stack <- stack(vapr_files)

vapr_values <- extract(vapr_stack, coords, method = "bilinear")

location_dat[paste0("vapr_", 1:12)] <- vapr_values

monthly_avg_vapr <- rowMeans(vapr_values, na.rm = TRUE)
location_dat$monthly_avg_vapr <- monthly_avg_vapr

#Maximum Temperature
tmax_1 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_01.tif")
tmax_2 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_02.tif")
tmax_3 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_03.tif")
tmax_4 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_04.tif")
tmax_5 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_05.tif")
tmax_6 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_06.tif")
tmax_7 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_07.tif")
tmax_8 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_08.tif")
tmax_9 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_09.tif")
tmax_10 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_10.tif")
tmax_11 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_11.tif")
tmax_12 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_12.tif")

crs(coords) <- crs(raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/wc2.1_30s_tmax_01.tif"))

tmax_files <- list.files("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmax/", 
                         pattern = "wc2.1_30s_tmax_.*\\.tif$", full.names = TRUE)
tmax_stack <- stack(tmax_files)

tmax_values <- extract(tmax_stack, coords, method = "bilinear")

location_dat[paste0("tmax_", 1:12)] <- tmax_values

monthly_avg_tmax <- rowMeans(tmax_values, na.rm = TRUE)
location_dat$monthly_avg_tmax <- monthly_avg_tmax


#Minumum Temperature
tmin_1 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_01.tif")
tmin_2 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_02.tif")
tmin_3 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_03.tif")
tmin_4 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_04.tif")
tmin_5 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_05.tif")
tmin_6 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_06.tif")
tmin_7 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_07.tif")
tmin_8 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_08.tif")
tmin_9 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_09.tif")
tmin_10 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_10.tif")
tmin_11 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_11.tif")
tmin_12 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_12.tif")

crs(coords) <- crs(raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/wc2.1_30s_tmin_01.tif"))

tmin_files <- list.files("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tmin/", 
                         pattern = "wc2.1_30s_tmin_.*\\.tif$", full.names = TRUE)
tmin_stack <- stack(tmin_files)

tmin_values <- extract(tmin_stack, coords, method = "bilinear")

location_dat[paste0("tmin_", 1:12)] <- tmin_values

monthly_avg_tmin <- rowMeans(tmin_values, na.rm = TRUE)
location_dat$monthly_avg_tmin <- monthly_avg_tmin


#Average Temperature
tavg_1 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_01.tif")
tavg_2 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_02.tif")
tavg_3 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_03.tif")
tavg_4 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_04.tif")
tavg_5 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_05.tif")
tavg_6 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_06.tif")
tavg_7 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_07.tif")
tavg_8 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_08.tif")
tavg_9 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_09.tif")
tavg_10 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_10.tif")
tavg_11 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_11.tif")
tavg_12 <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_12.tif")

crs(coords) <- crs(raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/wc2.1_30s_tavg_01.tif"))

tavg_files <- list.files("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/tavg/", 
                         pattern = "wc2.1_30s_tavg_.*\\.tif$", full.names = TRUE)
tavg_stack <- stack(tavg_files)

tavg_values <- extract(tavg_stack, coords, method = "bilinear")

location_dat[paste0("tavg_", 1:12)] <- tavg_values

monthly_avg_tavg <- rowMeans(tavg_values, na.rm = TRUE)
location_dat$monthly_avg_tavg <- monthly_avg_tavg

#Calculate VPD
#VPD
location_dat$VPD_1 <- 0.611 * exp((17.502 * location_dat$tmax_1) / (location_dat$tmax_1 + 240.97)) - location_dat$vapr_1
location_dat$VPD_2 <- 0.611 * exp((17.502 * location_dat$tmax_2) / (location_dat$tmax_2 + 240.97)) - location_dat$vapr_2
location_dat$VPD_3 <- 0.611 * exp((17.502 * location_dat$tmax_3) / (location_dat$tmax_3 + 240.97)) - location_dat$vapr_3
location_dat$VPD_4 <- 0.611 * exp((17.502 * location_dat$tmax_4) / (location_dat$tmax_4 + 240.97)) - location_dat$vapr_4
location_dat$VPD_5 <- 0.611 * exp((17.502 * location_dat$tmax_5) / (location_dat$tmax_5 + 240.97)) - location_dat$vapr_5
location_dat$VPD_6 <- 0.611 * exp((17.502 * location_dat$tmax_6) / (location_dat$tmax_6 + 240.97)) - location_dat$vapr_6
location_dat$VPD_7 <- 0.611 * exp((17.502 * location_dat$tmax_7) / (location_dat$tmax_7 + 240.97)) - location_dat$vapr_7
location_dat$VPD_8 <- 0.611 * exp((17.502 * location_dat$tmax_8) / (location_dat$tmax_8 + 240.97)) - location_dat$vapr_8
location_dat$VPD_9 <- 0.611 * exp((17.502 * location_dat$tmax_9) / (location_dat$tmax_9 + 240.97)) - location_dat$vapr_9
location_dat$VPD_10 <- 0.611 * exp((17.502 * location_dat$tmax_10) / (location_dat$tmax_10 + 240.97)) - location_dat$vapr_10
location_dat$VPD_11 <- 0.611 * exp((17.502 * location_dat$tmax_11) / (location_dat$tmax_11 + 240.97)) - location_dat$vapr_11
location_dat$VPD_12 <- 0.611 * exp((17.502 * location_dat$tmax_12) / (location_dat$tmax_12 + 240.97)) - location_dat$vapr_12

location_dat$VPD_avg <- rowMeans(location_dat[, paste0("VPD_", 1:12)], na.rm = TRUE)


#Model Testing

#NDVI and NDMI
NDVI_Hab <- glm(NDVI ~ Habitat, data = location_dat)
summary(NDVI_Hab) #p=0.023

NDMI_Hab <- glm(NDMI ~ Habitat, data = location_dat)
summary(NDMI_Hab) #p=0.158


#Annual Precipitation
#Habitat
annual_precip_Hab <- glm(annual_precip ~ Habitat, data = location_dat)
summary(annual_precip_Hab) #p=0.699

#Monthly Average
avg_precip_Hab <- glm(monthly_avg_precip ~ Habitat, data = location_dat)
summary(avg_precip_Hab) #p=0.699

#Solar Radiation
#Habitat
srad_1_Hab <- glm(srad_1 ~ Habitat, data = location_dat)
summary(srad_1_Hab) #p=0.165

srad_2_Hab <- glm(srad_2 ~ Habitat, data = location_dat)
summary(srad_2_Hab) #p=0.347

srad_3_Hab <- glm(srad_3 ~ Habitat, data = location_dat)
summary(srad_3_Hab) #p=0.079

srad_4_Hab <- glm(srad_4 ~ Habitat, data = location_dat)
summary(srad_4_Hab) #p=0.283

srad_5_Hab <- glm(srad_5 ~ Habitat, data = location_dat)
summary(srad_5_Hab) #p=0.967

srad_6_Hab <- glm(srad_6 ~ Habitat, data = location_dat)
summary(srad_6_Hab) #p=0.718

srad_7_Hab <- glm(srad_7 ~ Habitat, data = location_dat)
summary(srad_7_Hab) #p=0.534

srad_8_Hab <- glm(srad_8 ~ Habitat, data = location_dat)
summary(srad_8_Hab) #p=0.345

srad_9_Hab <- glm(srad_9 ~ Habitat, data = location_dat)
summary(srad_9_Hab) #p=0.147

srad_10_Hab <- glm(srad_10 ~ Habitat, data = location_dat)
summary(srad_10_Hab) #p=0.052

srad_11_Hab <- glm(srad_11 ~ Habitat, data = location_dat)
summary(srad_11_Hab) #p=0.215

srad_12_Hab <- glm(srad_12 ~ Habitat, data = location_dat)
summary(srad_12_Hab) #p=0.092


#Monthly Average
avg_srad_Hab <- glm(monthly_avg_srad ~ Habitat, data = location_dat)
summary(avg_srad_Hab) #p=0.100



#Water Vapor Pressure
#Habitat
vapr_1_Hab <- glm(vapr_1 ~ Habitat, data = location_dat)
summary(vapr_1_Hab) #p=0.917

vapr_2_Hab <- glm(vapr_2 ~ Habitat, data = location_dat)
summary(vapr_2_Hab) #p=0.912

vapr_3_Hab <- glm(vapr_3 ~ Habitat, data = location_dat)
summary(vapr_3_Hab) #p=0.859

vapr_4_Hab <- glm(vapr_4 ~ Habitat, data = location_dat)
summary(vapr_4_Hab) #p=0.847

vapr_5_Hab <- glm(vapr_5 ~ Habitat, data = location_dat)
summary(vapr_5_Hab) #p=0.878

vapr_6_Hab <- glm(vapr_6 ~ Habitat, data = location_dat)
summary(vapr_6_Hab) #p=0.878

vapr_7_Hab <- glm(vapr_7 ~ Habitat, data = location_dat)
summary(vapr_7_Hab) #p=0.876

vapr_8_Hab <- glm(vapr_8 ~ Habitat, data = location_dat)
summary(vapr_8_Hab) #p=0.881

vapr_9_Hab <- glm(vapr_9 ~ Habitat, data = location_dat)
summary(vapr_9_Hab) #p=0.891

vapr_10_Hab <- glm(vapr_10 ~ Habitat, data = location_dat)
summary(vapr_10_Hab) #p=0.922

vapr_11_Hab <- glm(vapr_11 ~ Habitat, data = location_dat)
summary(vapr_11_Hab) #p=0.895

vapr_12_Hab <- glm(vapr_12 ~ Habitat, data = location_dat)
summary(vapr_12_Hab) #p=0.870


#Monthly Average
avg_vapr_Hab <- glm(monthly_avg_vapr ~ Habitat, data = location_dat)
summary(avg_vapr_Hab) #p=0.885


#Maximum Temperature
#Habitat
tmax_1_Hab <- glm(tmax_1 ~ Habitat, data = location_dat)
summary(tmax_1_Hab) #p=0.829

tmax_2_Hab <- glm(tmax_2 ~ Habitat, data = location_dat)
summary(tmax_2_Hab) #p=0.755

tmax_3_Hab <- glm(tmax_3 ~ Habitat, data = location_dat)
summary(tmax_3_Hab) #p=0.742

tmax_4_Hab <- glm(tmax_4 ~ Habitat, data = location_dat)
summary(tmax_4_Hab) #p=0.668

tmax_5_Hab <- glm(tmax_5 ~ Habitat, data = location_dat)
summary(tmax_5_Hab) #p=0.745

tmax_6_Hab <- glm(tmax_6 ~ Habitat, data = location_dat)
summary(tmax_6_Hab) #p=0.694

tmax_7_Hab <- glm(tmax_7 ~ Habitat, data = location_dat)
summary(tmax_7_Hab) #p=0.936

tmax_8_Hab <- glm(tmax_8 ~ Habitat, data = location_dat)
summary(tmax_8_Hab) #p=0.969

tmax_9_Hab <- glm(tmax_9 ~ Habitat, data = location_dat)
summary(tmax_9_Hab) #p=0.718

tmax_10_Hab <- glm(tmax_10 ~ Habitat, data = location_dat)
summary(tmax_10_Hab) #p=0.759

tmax_11_Hab <- glm(tmax_11 ~ Habitat, data = location_dat)
summary(tmax_11_Hab) #p=0.785

tmax_12_Hab <- glm(tmax_12 ~ Habitat, data = location_dat)
summary(tmax_12_Hab) #p=0.816


#Monthly Average
avg_tmax_Hab <- glm(monthly_avg_tmax ~ Habitat, data = location_dat)
summary(avg_tmax_Hab) #p=0.767


#Minimum Temperature
#Habitat
tmin_1_Hab <- glm(tmin_1 ~ Habitat, data = location_dat)
summary(tmin_1_Hab) #p=0.818

tmin_2_Hab <- glm(tmin_2 ~ Habitat, data = location_dat)
summary(tmin_2_Hab) #p=0.824

tmin_3_Hab <- glm(tmin_3 ~ Habitat, data = location_dat)
summary(tmin_3_Hab) #p=0.824

tmin_4_Hab <- glm(tmin_4 ~ Habitat, data = location_dat)
summary(tmin_4_Hab) #p=0.836

tmin_5_Hab <- glm(tmin_5 ~ Habitat, data = location_dat)
summary(tmin_5_Hab) #p=0.863

tmin_6_Hab <- glm(tmin_6 ~ Habitat, data = location_dat)
summary(tmin_6_Hab) #p=0.878

tmin_7_Hab <- glm(tmin_7 ~ Habitat, data = location_dat)
summary(tmin_7_Hab) #p=0.906

tmin_8_Hab <- glm(tmin_8 ~ Habitat, data = location_dat)
summary(tmin_8_Hab) #p=0.895

tmin_9_Hab <- glm(tmin_9 ~ Habitat, data = location_dat)
summary(tmin_9_Hab) #p=0.828

tmin_10_Hab <- glm(tmin_10 ~ Habitat, data = location_dat)
summary(tmin_10_Hab) #p=0.846

tmin_11_Hab <- glm(tmin_11 ~ Habitat, data = location_dat)
summary(tmin_11_Hab) #p=0.806

tmin_12_Hab <- glm(tmin_12 ~ Habitat, data = location_dat)
summary(tmin_12_Hab) #p=0.836



#Monthly Average
avg_tmin_Hab <- glm(monthly_avg_tmin ~ Habitat, data = location_dat)
summary(avg_tmin_Hab) #p=0.843


#Average Temperature
#Habitat
tavg_1_Hab <- glm(tavg_1 ~ Habitat, data = location_dat)
summary(tavg_1_Hab) #p=0.814

tavg_2_Hab <- glm(tavg_2 ~ Habitat, data = location_dat)
summary(tavg_2_Hab) #p=0.796

tavg_3_Hab <- glm(tavg_3 ~ Habitat, data = location_dat)
summary(tavg_3_Hab) #p=0.785

tavg_4_Hab <- glm(tavg_4 ~ Habitat, data = location_dat)
summary(tavg_4_Hab) #p=0.749

tavg_5_Hab <- glm(tavg_5 ~ Habitat, data = location_dat)
summary(tavg_5_Hab) #p=0.817

tavg_6_Hab <- glm(tavg_6 ~ Habitat, data = location_dat)
summary(tavg_6_Hab) #p=0.777

tavg_7_Hab <- glm(tavg_7 ~ Habitat, data = location_dat)
summary(tavg_7_Hab) #p=0.849

tavg_8_Hab <- glm(tavg_8 ~ Habitat, data = location_dat)
summary(tavg_8_Hab) #p=0.880

tavg_9_Hab <- glm(tavg_9 ~ Habitat, data = location_dat)
summary(tavg_9_Hab) #p=0.763

tavg_10_Hab <- glm(tavg_10 ~ Habitat, data = location_dat)
summary(tavg_10_Hab) #p=0.779

tavg_11_Hab <- glm(tavg_11 ~ Habitat, data = location_dat)
summary(tavg_11_Hab) #p=0.798

tavg_12_Hab <- glm(tavg_12 ~ Habitat, data = location_dat)
summary(tavg_12_Hab) #p=0.799



#Monthly Average
avg_tavg_Hab <- glm(monthly_avg_tavg ~ Habitat, data = location_dat)
summary(avg_tavg_Hab) #p=0.795

#VPD
#Habitat
VPD_1_Hab <- glm(VPD_1 ~ Habitat, data = location_dat)
summary(VPD_1_Hab) #p=0.826

VPD_2_Hab <- glm(VPD_2 ~ Habitat, data = location_dat)
summary(VPD_2_Hab) #p=0.650

VPD_3_Hab <- glm(VPD_3 ~ Habitat, data = location_dat)
summary(VPD_3_Hab) #p=0.628

VPD_4_Hab <- glm(VPD_4 ~ Habitat, data = location_dat)
summary(VPD_4_Hab) #p=0.372

VPD_5_Hab <- glm(VPD_5 ~ Habitat, data = location_dat)
summary(VPD_5_Hab) #p=0.714

VPD_6_Hab <- glm(VPD_6 ~ Habitat, data = location_dat)
summary(VPD_6_Hab) #p=0.881

VPD_7_Hab <- glm(VPD_7 ~ Habitat, data = location_dat)
summary(VPD_7_Hab) #p=0.866

VPD_8_Hab <- glm(VPD_8 ~ Habitat, data = location_dat)
summary(VPD_8_Hab) #p=0.852

VPD_9_Hab <- glm(VPD_9 ~ Habitat, data = location_dat)
summary(VPD_9_Hab) #p=0.783

VPD_10_Hab <- glm(VPD_10 ~ Habitat, data = location_dat)
summary(VPD_10_Hab) #p=0.263

VPD_11_Hab <- glm(VPD_11 ~ Habitat, data = location_dat)
summary(VPD_11_Hab) #p=0.674

VPD_12_Hab <- glm(VPD_12 ~ Habitat, data = location_dat)
summary(VPD_12_Hab) #p=0.845

#Monthly Average
avg_VPD_Hab <- glm(VPD_avg ~ Habitat, data = location_dat)
summary(avg_VPD_Hab) #p=0.694

model_list <- list(
  NDMI_Hab,
  annual_precip_Hab, avg_precip_Hab,
  srad_1_Hab, srad_2_Hab, srad_3_Hab, srad_4_Hab, srad_5_Hab, srad_6_Hab, srad_7_Hab,
  srad_8_Hab, srad_9_Hab, srad_10_Hab, srad_11_Hab, srad_12_Hab, avg_srad_Hab,
  vapr_1_Hab, vapr_2_Hab, vapr_3_Hab, vapr_4_Hab, vapr_5_Hab, vapr_6_Hab,
  vapr_7_Hab, vapr_8_Hab, vapr_9_Hab, vapr_10_Hab, vapr_11_Hab, vapr_12_Hab, avg_vapr_Hab,
  tmax_1_Hab, tmax_2_Hab, tmax_3_Hab, tmax_4_Hab, tmax_5_Hab, tmax_6_Hab,
  tmax_7_Hab, tmax_8_Hab, tmax_9_Hab, tmax_10_Hab, tmax_11_Hab, tmax_12_Hab, avg_tmax_Hab,
  tmin_1_Hab, tmin_2_Hab, tmin_3_Hab, tmin_4_Hab, tmin_5_Hab, tmin_6_Hab,
  tmin_7_Hab, tmin_8_Hab, tmin_9_Hab, tmin_10_Hab, tmin_11_Hab, tmin_12_Hab, avg_tmin_Hab,
  tavg_1_Hab, tavg_2_Hab, tavg_3_Hab, tavg_4_Hab, tavg_5_Hab, tavg_6_Hab,
  tavg_7_Hab, tavg_8_Hab, tavg_9_Hab, tavg_10_Hab, tavg_11_Hab, tavg_12_Hab, avg_tavg_Hab,
  VPD_1_Hab, VPD_2_Hab, VPD_3_Hab, VPD_4_Hab, VPD_5_Hab, VPD_6_Hab,
  VPD_7_Hab, VPD_8_Hab, VPD_9_Hab, VPD_10_Hab, VPD_11_Hab, VPD_12_Hab, avg_VPD_Hab
)

max_t <- max(sapply(model_list, function(mod) {
  coefs <- summary(mod)$coefficients
  if ("Habitatwall" %in% rownames(coefs)) {
    abs(coefs["Habitatwall", "t value"])
  } else {
    NA
  }
}), na.rm = TRUE)

max_t


#Does Macroclimate predict microclimate
#Macro vs Micro Climate
combined_dat <- merge(dat, location_dat, by = "Site", all.x = TRUE)
values <- c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")

#annual precip
Macro_V_Micro_1 <- lm(ambient_temp ~ annual_precip, data = combined_dat)
summary(Macro_V_Micro_1)
plot(combined_dat$annual_precip, combined_dat$ambient_temp, 
     col = values,
     pch = 19,
     xlab = "Annual Precipitation (mm)", 
     ylab = "Ambient Temperature (°C)", 
     main = "Ambient Temperature vs. Annual Precipitation")
abline(Macro_V_Micro_1, col = "red", lwd = 2)

Macro_V_Micro_2 <- lm(VPD ~ annual_precip, data = combined_dat)
summary(Macro_V_Micro_2)
plot(combined_dat$annual_precip, combined_dat$VPD, 
     col = values,
     pch = 19,
     xlab = "Annual Precipitation (mm)", 
     ylab = "VPD", 
     main = "VPD vs. Annual Precipitation")
abline(Macro_V_Micro_2, col = "red", lwd = 2)

Macro_V_Micro_3 <- lm(percent_veg_cover ~ annual_precip, data = combined_dat)
summary(Macro_V_Micro_3)
plot(combined_dat$annual_precip, combined_dat$percent_veg_cover, 
     col = values,
     pch = 19,
     xlab = "Annual Precipitation (mm)", 
     ylab = "percent_veg_cover", 
     main = "percent_veg_cover vs. Annual Precipitation")
abline(Macro_V_Micro_3, col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Extract data for each model
ambient_info <- extract_model_info(Macro_V_Micro_1)
vpd_info <- extract_model_info(Macro_V_Micro_2)
veg_info <- extract_model_info(Macro_V_Micro_3)

# Combine into a data frame
annual_precip_results_table <- data.frame(
  Response = c("Ambient Temperature (°C)", "VPD", "Percent Vegetation Cover (%)"),
  Estimate = c(ambient_info[1], vpd_info[1], veg_info[1]),
  Std_Error = c(ambient_info[2], vpd_info[2], veg_info[2]),
  t_value = c(ambient_info[3], vpd_info[3], veg_info[3]),
  p_value = c(ambient_info[4], vpd_info[4], veg_info[4]),
  Adj_R2 = c(ambient_info[5], vpd_info[5], veg_info[5])
)

print(annual_precip_results_table)



#solar radiation
#ambient_temp
srad_ambient_temp_models <- list()
for (i in 1:12) {
  srad_col <- paste0("srad_", i)
  formula <- as.formula(paste("ambient_temp ~", srad_col))
  srad_ambient_temp_models[[srad_col]] <- lm(formula, data = combined_dat)
}

summary(srad_ambient_temp_models[["srad_1"]])
summary(srad_ambient_temp_models[["srad_2"]])
summary(srad_ambient_temp_models[["srad_3"]])
summary(srad_ambient_temp_models[["srad_4"]])
summary(srad_ambient_temp_models[["srad_5"]])
summary(srad_ambient_temp_models[["srad_6"]])
summary(srad_ambient_temp_models[["srad_7"]])
summary(srad_ambient_temp_models[["srad_8"]])
summary(srad_ambient_temp_models[["srad_9"]])
summary(srad_ambient_temp_models[["srad_10"]])
summary(srad_ambient_temp_models[["srad_11"]])
summary(srad_ambient_temp_models[["srad_12"]])

plot(combined_dat$srad_12, combined_dat$ambient_temp,
     pch = 19, col = values,
     xlab = "srad_12", ylab = "Ambient Temperature (°C)",
     main = "Ambient Temp vs srad_12")
abline(srad_ambient_temp_models[["srad_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Loop through each model and extract summary stats
srad_results <- data.frame()
for (month in 1:12) {
  srad_col <- paste0("srad_", month)
  model <- srad_ambient_temp_models[[srad_col]]
  info <- extract_model_info(model)
  srad_results <- rbind(srad_results,
                        data.frame(
                          Month = month,
                          Estimate = info[1],
                          Std_Error = info[2],
                          t_value = info[3],
                          p_value = info[4],
                          Adj_R2 = info[5]
                        ))
}

# View table
print(srad_results)


#VPD
srad_VPD_models <- list()
for (i in 1:12) {
  srad_col <- paste0("srad_", i)
  formula <- as.formula(paste("VPD ~", srad_col))
  srad_VPD_models[[srad_col]] <- lm(formula, data = combined_dat)
}


summary(srad_VPD_models[["srad_1"]])
summary(srad_VPD_models[["srad_2"]])
summary(srad_VPD_models[["srad_3"]])
summary(srad_VPD_models[["srad_4"]])
summary(srad_VPD_models[["srad_5"]])
summary(srad_VPD_models[["srad_6"]])
summary(srad_VPD_models[["srad_7"]])
summary(srad_VPD_models[["srad_8"]])
summary(srad_VPD_models[["srad_9"]])
summary(srad_VPD_models[["srad_10"]])
summary(srad_VPD_models[["srad_11"]])
summary(srad_VPD_models[["srad_12"]])

plot(combined_dat$srad_12, combined_dat$VPD,
     pch = 19, col = values,
     xlab = "srad_12", ylab = "VPD",
     main = "VPD vs srad_12")
abline(srad_VPD_models[["srad_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Loop through each VPD model and extract summary stats
srad_VPD_results <- data.frame()
for (month in 1:12) {
  srad_col <- paste0("srad_", month)
  model <- srad_VPD_models[[srad_col]]
  info <- extract_model_info(model)
  srad_VPD_results <- rbind(srad_VPD_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View table
print(srad_VPD_results)


#Percent Veg Cover
srad_Veg_models <- list()
for (i in 1:12) {
  srad_col <- paste0("srad_", i)
  formula <- as.formula(paste("percent_veg_cover ~", srad_col))
  srad_Veg_models[[srad_col]] <- lm(formula, data = combined_dat)
}

summary(srad_Veg_models[["srad_1"]])
summary(srad_Veg_models[["srad_2"]])
summary(srad_Veg_models[["srad_3"]])
summary(srad_Veg_models[["srad_4"]])
summary(srad_Veg_models[["srad_5"]])
summary(srad_Veg_models[["srad_6"]])
summary(srad_Veg_models[["srad_7"]])
summary(srad_Veg_models[["srad_8"]])
summary(srad_Veg_models[["srad_9"]])
summary(srad_Veg_models[["srad_10"]])
summary(srad_Veg_models[["srad_11"]])
summary(srad_Veg_models[["srad_12"]])

plot(combined_dat$srad_12, combined_dat$percent_veg_cover,
     pch = 19, col = values,
     xlab = "srad_12", ylab = "percent_veg_cover",
     main = "percent_veg_cover vs srad_12")
abline(srad_Veg_models[["srad_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for percent vegetation cover ~ srad models
srad_Veg_results <- data.frame()
for (month in 1:12) {
  srad_col <- paste0("srad_", month)
  model <- srad_Veg_models[[srad_col]]
  info <- extract_model_info(model)
  srad_Veg_results <- rbind(srad_Veg_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View or export the table
print(srad_Veg_results)



#water vapor pressure
#ambient_temp
vapr_ambient_temp_models <- list()
for (i in 1:12) {
  vapr_col <- paste0("vapr_", i)
  formula <- as.formula(paste("ambient_temp ~", vapr_col))
  vapr_ambient_temp_models[[vapr_col]] <- lm(formula, data = combined_dat)
}

summary(vapr_ambient_temp_models[["vapr_1"]])
summary(vapr_ambient_temp_models[["vapr_2"]])
summary(vapr_ambient_temp_models[["vapr_3"]])
summary(vapr_ambient_temp_models[["vapr_4"]])
summary(vapr_ambient_temp_models[["vapr_5"]])
summary(vapr_ambient_temp_models[["vapr_6"]])
summary(vapr_ambient_temp_models[["vapr_7"]])
summary(vapr_ambient_temp_models[["vapr_8"]])
summary(vapr_ambient_temp_models[["vapr_9"]])
summary(vapr_ambient_temp_models[["vapr_10"]])
summary(vapr_ambient_temp_models[["vapr_11"]])
summary(vapr_ambient_temp_models[["vapr_12"]])

plot(combined_dat$vapr_12, combined_dat$ambient_temp,
     pch = 19, col = values,
     xlab = "vapr_12", ylab = "Ambient Temperature (°C)",
     main = "Ambient Temp vs vapr_12")
abline(vapr_ambient_temp_models[["vapr_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create summary table for ambient temp ~ vapr models
vapr_ambient_temp_results <- data.frame()
for (month in 1:12) {
  vapr_col <- paste0("vapr_", month)
  model <- vapr_ambient_temp_models[[vapr_col]]
  info <- extract_model_info(model)
  vapr_ambient_temp_results <- rbind(vapr_ambient_temp_results,
                                     data.frame(
                                       Month = month,
                                       Estimate = info[1],
                                       Std_Error = info[2],
                                       t_value = info[3],
                                       p_value = info[4],
                                       Adj_R2 = info[5]
                                     ))
}

# View or export the table
print(vapr_ambient_temp_results)


#VPD
vapr_VPD_models <- list()
for (i in 1:12) {
  vapr_col <- paste0("vapr_", i)
  formula <- as.formula(paste("VPD ~", vapr_col))
  vapr_VPD_models[[vapr_col]] <- lm(formula, data = combined_dat)
}


summary(vapr_VPD_models[["vapr_1"]])
summary(vapr_VPD_models[["vapr_2"]])
summary(vapr_VPD_models[["vapr_3"]])
summary(vapr_VPD_models[["vapr_4"]])
summary(vapr_VPD_models[["vapr_5"]])
summary(vapr_VPD_models[["vapr_6"]])
summary(vapr_VPD_models[["vapr_7"]])
summary(vapr_VPD_models[["vapr_8"]])
summary(vapr_VPD_models[["vapr_9"]])
summary(vapr_VPD_models[["vapr_10"]])
summary(vapr_VPD_models[["vapr_11"]])
summary(vapr_VPD_models[["vapr_12"]])

plot(combined_dat$vapr_12, combined_dat$VPD,
     pch = 19, col = values,
     xlab = "vapr_12", ylab = "VPD",
     main = "VPD vs vapr_12")
abline(vapr_VPD_models[["vapr_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Generate summary table for VPD ~ vapr models
vapr_VPD_results <- data.frame()
for (month in 1:12) {
  vapr_col <- paste0("vapr_", month)
  model <- vapr_VPD_models[[vapr_col]]
  info <- extract_model_info(model)
  vapr_VPD_results <- rbind(vapr_VPD_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View or export
print(vapr_VPD_results)



#Percent Veg Cover
vapr_Veg_models <- list()
for (i in 1:12) {
  vapr_col <- paste0("vapr_", i)
  formula <- as.formula(paste("percent_veg_cover ~", vapr_col))
  vapr_Veg_models[[vapr_col]] <- lm(formula, data = combined_dat)
}

summary(vapr_Veg_models[["vapr_1"]])
summary(vapr_Veg_models[["vapr_2"]])
summary(vapr_Veg_models[["vapr_3"]])
summary(vapr_Veg_models[["vapr_4"]])
summary(vapr_Veg_models[["vapr_5"]])
summary(vapr_Veg_models[["vapr_6"]])
summary(vapr_Veg_models[["vapr_7"]])
summary(vapr_Veg_models[["vapr_8"]])
summary(vapr_Veg_models[["vapr_9"]])
summary(vapr_Veg_models[["vapr_10"]])
summary(vapr_Veg_models[["vapr_11"]])
summary(vapr_Veg_models[["vapr_12"]])

plot(combined_dat$vapr_12, combined_dat$percent_veg_cover,
     pch = 19, col = values,
     xlab = "vapr_12", ylab = "percent_veg_cover",
     main = "percent_veg_cover vs vapr_12")
abline(vapr_Veg_models[["vapr_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for vapr ~ veg cover
vapr_Veg_results <- data.frame()
for (month in 1:12) {
  vapr_col <- paste0("vapr_", month)
  model <- vapr_Veg_models[[vapr_col]]
  info <- extract_model_info(model)
  vapr_Veg_results <- rbind(vapr_Veg_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View or export
print(vapr_Veg_results)

#monthly max temp
#ambient_temp
tmax_ambient_temp_models <- list()
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  formula <- as.formula(paste("ambient_temp ~", tmax_col))
  tmax_ambient_temp_models[[tmax_col]] <- lm(formula, data = combined_dat)
}


summary(tmax_ambient_temp_models[["tmax_1"]])
summary(tmax_ambient_temp_models[["tmax_2"]])
summary(tmax_ambient_temp_models[["tmax_3"]])
summary(tmax_ambient_temp_models[["tmax_4"]])
summary(tmax_ambient_temp_models[["tmax_5"]])
summary(tmax_ambient_temp_models[["tmax_6"]])
summary(tmax_ambient_temp_models[["tmax_7"]])
summary(tmax_ambient_temp_models[["tmax_8"]])
summary(tmax_ambient_temp_models[["tmax_9"]])
summary(tmax_ambient_temp_models[["tmax_10"]])
summary(tmax_ambient_temp_models[["tmax_11"]])
summary(tmax_ambient_temp_models[["tmax_12"]])

plot(combined_dat$tmax_12, combined_dat$ambient_temp,
     pch = 19, col = values,
     xlab = "tmax_12", ylab = "Ambient Temperature (°C)",
     main = "Ambient Temp vs tmax_12")
abline(tmax_ambient_temp_models[["tmax_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for tmax ~ ambient temperature
tmax_ambient_temp_results <- data.frame()
for (month in 1:12) {
  tmax_col <- paste0("tmax_", month)
  model <- tmax_ambient_temp_models[[tmax_col]]
  info <- extract_model_info(model)
  tmax_ambient_temp_results <- rbind(tmax_ambient_temp_results,
                                     data.frame(
                                       Month = month,
                                       Estimate = info[1],
                                       Std_Error = info[2],
                                       t_value = info[3],
                                       p_value = info[4],
                                       Adj_R2 = info[5]
                                     ))
}

# View or export
print(tmax_ambient_temp_results)


#VPD
tmax_VPD_models <- list()
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  formula <- as.formula(paste("VPD ~", tmax_col))
  tmax_VPD_models[[tmax_col]] <- lm(formula, data = combined_dat)
}

summary(tmax_VPD_models[["tmax_1"]])
summary(tmax_VPD_models[["tmax_2"]])
summary(tmax_VPD_models[["tmax_3"]])
summary(tmax_VPD_models[["tmax_4"]])
summary(tmax_VPD_models[["tmax_5"]])
summary(tmax_VPD_models[["tmax_6"]])
summary(tmax_VPD_models[["tmax_7"]])
summary(tmax_VPD_models[["tmax_8"]])
summary(tmax_VPD_models[["tmax_9"]])
summary(tmax_VPD_models[["tmax_10"]])
summary(tmax_VPD_models[["tmax_11"]])
summary(tmax_VPD_models[["tmax_12"]])

plot(combined_dat$tmax_12, combined_dat$VPD,
     pch = 19, col = values,
     xlab = "tmax_12", ylab = "VPD",
     main = "VPD vs tmax_12")
abline(tmax_VPD_models[["tmax_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for tmax ~ VPD
tmax_VPD_results <- data.frame()
for (month in 1:12) {
  tmax_col <- paste0("tmax_", month)
  model <- tmax_VPD_models[[tmax_col]]
  info <- extract_model_info(model)
  tmax_VPD_results <- rbind(tmax_VPD_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View or export
print(tmax_VPD_results)

#Percent Veg Cover
tmax_Veg_models <- list()
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  formula <- as.formula(paste("percent_veg_cover ~", tmax_col))
  tmax_Veg_models[[tmax_col]] <- lm(formula, data = combined_dat)
}


summary(tmax_Veg_models[["tmax_1"]])
summary(tmax_Veg_models[["tmax_2"]])
summary(tmax_Veg_models[["tmax_3"]])
summary(tmax_Veg_models[["tmax_4"]])
summary(tmax_Veg_models[["tmax_5"]])
summary(tmax_Veg_models[["tmax_6"]])
summary(tmax_Veg_models[["tmax_7"]])
summary(tmax_Veg_models[["tmax_8"]])
summary(tmax_Veg_models[["tmax_9"]])
summary(tmax_Veg_models[["tmax_10"]])
summary(tmax_Veg_models[["tmax_11"]])
summary(tmax_Veg_models[["tmax_12"]])

plot(combined_dat$tmax_12, combined_dat$percent_veg_cover,
     pch = 19, col = values,
     xlab = "tmax_12", ylab = "percent_veg_cover",
     main = "percent_veg_cover vs tmax_12")
abline(tmax_Veg_models[["tmax_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for tmax ~ percent_veg_cover
tmax_Veg_results <- data.frame()
for (month in 1:12) {
  tmax_col <- paste0("tmax_", month)
  model <- tmax_Veg_models[[tmax_col]]
  info <- extract_model_info(model)
  tmax_Veg_results <- rbind(tmax_Veg_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View or export
print(tmax_Veg_results)


#monthly average min temp
#ambient_temp
tmin_ambient_temp_models <- list()
for (i in 1:12) {
  tmin_col <- paste0("tmin_", i)
  formula <- as.formula(paste("ambient_temp ~", tmin_col))
  tmin_ambient_temp_models[[tmin_col]] <- lm(formula, data = combined_dat)
}


summary(tmin_ambient_temp_models[["tmin_1"]])
summary(tmin_ambient_temp_models[["tmin_2"]])
summary(tmin_ambient_temp_models[["tmin_3"]])
summary(tmin_ambient_temp_models[["tmin_4"]])
summary(tmin_ambient_temp_models[["tmin_5"]])
summary(tmin_ambient_temp_models[["tmin_6"]])
summary(tmin_ambient_temp_models[["tmin_7"]])
summary(tmin_ambient_temp_models[["tmin_8"]])
summary(tmin_ambient_temp_models[["tmin_9"]])
summary(tmin_ambient_temp_models[["tmin_10"]])
summary(tmin_ambient_temp_models[["tmin_11"]])
summary(tmin_ambient_temp_models[["tmin_12"]])


plot(combined_dat$tmin_12, combined_dat$ambient_temp,
     pch = 19, col = values,
     xlab = "tmin_12", ylab = "Ambient Temperature (°C)",
     main = "Ambient Temp vs tmin_12")
abline(tmin_ambient_temp_models[["tmin_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for tmin ~ ambient_temp
tmin_ambient_temp_results <- data.frame()
for (month in 1:12) {
  tmin_col <- paste0("tmin_", month)
  model <- tmin_ambient_temp_models[[tmin_col]]
  info <- extract_model_info(model)
  tmin_ambient_temp_results <- rbind(tmin_ambient_temp_results,
                                     data.frame(
                                       Month = month,
                                       Estimate = info[1],
                                       Std_Error = info[2],
                                       t_value = info[3],
                                       p_value = info[4],
                                       Adj_R2 = info[5]
                                     ))
}

# View or export the table
print(tmin_ambient_temp_results)


#VPD
tmin_VPD_models <- list()
for (i in 1:12) {
  tmin_col <- paste0("tmin_", i)
  formula <- as.formula(paste("VPD ~", tmin_col))
  tmin_VPD_models[[tmin_col]] <- lm(formula, data = combined_dat)
}


summary(tmin_VPD_models[["tmin_1"]])
summary(tmin_VPD_models[["tmin_2"]])
summary(tmin_VPD_models[["tmin_3"]])
summary(tmin_VPD_models[["tmin_4"]])
summary(tmin_VPD_models[["tmin_5"]])
summary(tmin_VPD_models[["tmin_6"]])
summary(tmin_VPD_models[["tmin_7"]])
summary(tmin_VPD_models[["tmin_8"]])
summary(tmin_VPD_models[["tmin_9"]])
summary(tmin_VPD_models[["tmin_10"]])
summary(tmin_VPD_models[["tmin_11"]])
summary(tmin_VPD_models[["tmin_12"]])

plot(combined_dat$tmin_12, combined_dat$VPD,
     pch = 19, col = values,
     xlab = "tmin_12", ylab = "VPD",
     main = "VPD vs tmin_12")
abline(tmin_VPD_models[["tmin_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for tmin ~ VPD
tmin_VPD_results <- data.frame()
for (month in 1:12) {
  tmin_col <- paste0("tmin_", month)
  model <- tmin_VPD_models[[tmin_col]]
  info <- extract_model_info(model)
  tmin_VPD_results <- rbind(tmin_VPD_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View or export the table
print(tmin_VPD_results)

#Percent Veg Cover
tmin_Veg_models <- list()
for (i in 1:12) {
  tmin_col <- paste0("tmin_", i)
  formula <- as.formula(paste("percent_veg_cover ~", tmin_col))
  tmin_Veg_models[[tmin_col]] <- lm(formula, data = combined_dat)
}


summary(tmin_Veg_models[["tmin_1"]])
summary(tmin_Veg_models[["tmin_2"]])
summary(tmin_Veg_models[["tmin_3"]])
summary(tmin_Veg_models[["tmin_4"]])
summary(tmin_Veg_models[["tmin_5"]])
summary(tmin_Veg_models[["tmin_6"]])
summary(tmin_Veg_models[["tmin_7"]])
summary(tmin_Veg_models[["tmin_8"]])
summary(tmin_Veg_models[["tmin_9"]])
summary(tmin_Veg_models[["tmin_10"]])
summary(tmin_Veg_models[["tmin_11"]])
summary(tmin_Veg_models[["tmin_12"]])

plot(combined_dat$tmin_12, combined_dat$percent_veg_cover,
     pch = 19, col = values,
     xlab = "tmin_12", ylab = "percent_veg_cover",
     main = "percent_veg_cover vs tmin_12")
abline(tmin_Veg_models[["tmin_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for tmin ~ percent_veg_cover
tmin_Veg_results <- data.frame()
for (month in 1:12) {
  tmin_col <- paste0("tmin_", month)
  model <- tmin_Veg_models[[tmin_col]]
  info <- extract_model_info(model)
  tmin_Veg_results <- rbind(tmin_Veg_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View or export the table
print(tmin_Veg_results)


#monthly avg temp
#ambient_temp
tavg_ambient_temp_models <- list()
for (i in 1:12) {
  tavg_col <- paste0("tavg_", i)
  formula <- as.formula(paste("ambient_temp ~", tavg_col))
  tavg_ambient_temp_models[[tavg_col]] <- lm(formula, data = combined_dat)
}

summary(tavg_ambient_temp_models[["tavg_1"]])
summary(tavg_ambient_temp_models[["tavg_2"]])
summary(tavg_ambient_temp_models[["tavg_3"]])
summary(tavg_ambient_temp_models[["tavg_4"]])
summary(tavg_ambient_temp_models[["tavg_5"]])
summary(tavg_ambient_temp_models[["tavg_6"]])
summary(tavg_ambient_temp_models[["tavg_7"]])
summary(tavg_ambient_temp_models[["tavg_8"]])
summary(tavg_ambient_temp_models[["tavg_9"]])
summary(tavg_ambient_temp_models[["tavg_10"]])
summary(tavg_ambient_temp_models[["tavg_11"]])
summary(tavg_ambient_temp_models[["tavg_12"]])

plot(combined_dat$tavg_12, combined_dat$ambient_temp,
     pch = 19, col = values,
     xlab = "tavg_12", ylab = "Ambient Temperature (°C)",
     main = "Ambient Temp vs tavg_12")
abline(tavg_ambient_temp_models[["tavg_12"]], col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Create table for tavg ~ ambient_temp
tavg_ambient_temp_results <- data.frame()
for (month in 1:12) {
  tavg_col <- paste0("tavg_", month)
  model <- tavg_ambient_temp_models[[tavg_col]]
  info <- extract_model_info(model)
  tavg_ambient_temp_results <- rbind(tavg_ambient_temp_results,
                                     data.frame(
                                       Month = month,
                                       Estimate = info[1],
                                       Std_Error = info[2],
                                       t_value = info[3],
                                       p_value = info[4],
                                       Adj_R2 = info[5]
                                     ))
}

# View or export the table
print(tavg_ambient_temp_results)

#VPD
tavg_VPD_models <- list()
for (i in 1:12) {
  tavg_col <- paste0("tavg_", i)
  formula <- as.formula(paste("VPD ~", tavg_col))
  tavg_VPD_models[[tavg_col]] <- lm(formula, data = combined_dat)
}


summary(tavg_VPD_models[["tavg_1"]])
summary(tavg_VPD_models[["tavg_2"]])
summary(tavg_VPD_models[["tavg_3"]])
summary(tavg_VPD_models[["tavg_4"]])
summary(tavg_VPD_models[["tavg_5"]])
summary(tavg_VPD_models[["tavg_6"]])
summary(tavg_VPD_models[["tavg_7"]])
summary(tavg_VPD_models[["tavg_8"]])
summary(tavg_VPD_models[["tavg_9"]])
summary(tavg_VPD_models[["tavg_10"]])
summary(tavg_VPD_models[["tavg_11"]])
summary(tavg_VPD_models[["tavg_12"]])

plot(combined_dat$tavg_12, combined_dat$VPD,
     pch = 19, col = values,
     xlab = "tavg_12", ylab = "VPD",
     main = "VPD vs tavg_12")
abline(tavg_VPD_models[["tavg_12"]], col = "red", lwd = 2)

tavg_VPD_results <- data.frame()
for (month in 1:12) {
  tavg_col <- paste0("tavg_", month)
  model <- tavg_VPD_models[[tavg_col]]
  info <- extract_model_info(model)
  tavg_VPD_results <- rbind(tavg_VPD_results,
                            data.frame(
                              Month = month,
                              Estimate = info[1],
                              Std_Error = info[2],
                              t_value = info[3],
                              p_value = info[4],
                              Adj_R2 = info[5]
                            ))
}

# View or export the table
print(tavg_VPD_results)


#Percent Veg Cover
tavg_Veg_models <- list()
for (i in 1:12) {
  tavg_col <- paste0("tavg_", i)
  formula <- as.formula(paste("percent_veg_cover ~", tavg_col))
  tavg_Veg_models[[tavg_col]] <- lm(formula, data = combined_dat)
}

summary(tavg_Veg_models[["tavg_1"]])
summary(tavg_Veg_models[["tavg_2"]])
summary(tavg_Veg_models[["tavg_3"]])
summary(tavg_Veg_models[["tavg_4"]])
summary(tavg_Veg_models[["tavg_5"]])
summary(tavg_Veg_models[["tavg_6"]])
summary(tavg_Veg_models[["tavg_7"]])
summary(tavg_Veg_models[["tavg_8"]])
summary(tavg_Veg_models[["tavg_9"]])
summary(tavg_Veg_models[["tavg_10"]])
summary(tavg_Veg_models[["tavg_11"]])
summary(tavg_Veg_models[["tavg_12"]])

plot(combined_dat$tavg_12, combined_dat$percent_veg_cover,
     pch = 19, col = values,
     xlab = "tavg_12", ylab = "percent_veg_cover",
     main = "percent_veg_cover vs tavg_12")
abline(tavg_Veg_models[["tavg_12"]], col = "red", lwd = 2)

tavg_Veg_summary <- data.frame(Model = character(), Estimate = numeric(), Std_Error = numeric(), 
                               t_value = numeric(), p_value = numeric(), Adj_R2 = numeric())

# Loop through each month to summarize the models
for (i in 1:12) {
  # Extract model results
  model <- tavg_Veg_models[[paste0("tavg_", i)]]
  
  # Summary of model
  model_summary <- summary(model)
  
  # Store the relevant values in the results table
  tavg_Veg_summary <- rbind(tavg_Veg_summary, data.frame(
    Model = paste0("tavg_", i),
    Estimate = model_summary$coefficients[2, 1], # Estimate for the predictor
    Std_Error = model_summary$coefficients[2, 2], # Standard error of the predictor
    t_value = model_summary$coefficients[2, 3],  # t-value of the predictor
    p_value = model_summary$coefficients[2, 4],  # p-value of the predictor
    Adj_R2 = model_summary$adj.r.squared         # Adjusted R-squared
  ))
}

# Print the summary table
print(tavg_Veg_summary)

#VPD
#ambient_temp
VPD_ambient_temp_models <- list()
for (i in 1:12) {
  VPD_col <- paste0("VPD_", i)
  formula <- as.formula(paste("ambient_temp ~", VPD_col))
  VPD_ambient_temp_models[[VPD_col]] <- lm(formula, data = combined_dat)
}


summary(VPD_ambient_temp_models[["VPD_1"]])
summary(VPD_ambient_temp_models[["VPD_2"]])
summary(VPD_ambient_temp_models[["VPD_3"]])
summary(VPD_ambient_temp_models[["VPD_4"]])
summary(VPD_ambient_temp_models[["VPD_5"]])
summary(VPD_ambient_temp_models[["VPD_6"]])
summary(VPD_ambient_temp_models[["VPD_7"]])
summary(VPD_ambient_temp_models[["VPD_8"]])
summary(VPD_ambient_temp_models[["VPD_9"]])
summary(VPD_ambient_temp_models[["VPD_10"]])
summary(VPD_ambient_temp_models[["VPD_11"]])
summary(VPD_ambient_temp_models[["VPD_12"]])


plot(combined_dat$VPD_12, combined_dat$ambient_temp,
     pch = 19, col = values,
     xlab = "VPD_12", ylab = "Ambient Temperature (°C)",
     main = "Ambient Temp vs VPD_12")
abline(VPD_ambient_temp_models[["VPD_12"]], col = "red", lwd = 2)

VPD_ambient_temp_summary <- data.frame(Model = character(), Estimate = numeric(), Std_Error = numeric(), 
                                       t_value = numeric(), p_value = numeric(), Adj_R2 = numeric())

# Loop through each month to summarize the models
for (i in 1:12) {
  # Extract model results
  model <- VPD_ambient_temp_models[[paste0("VPD_", i)]]
  
  # Summary of model
  model_summary <- summary(model)
  
  # Store the relevant values in the results table
  VPD_ambient_temp_summary <- rbind(VPD_ambient_temp_summary, data.frame(
    Model = paste0("VPD_", i),
    Estimate = model_summary$coefficients[2, 1], # Estimate for the predictor
    Std_Error = model_summary$coefficients[2, 2], # Standard error of the predictor
    t_value = model_summary$coefficients[2, 3],  # t-value of the predictor
    p_value = model_summary$coefficients[2, 4],  # p-value of the predictor
    Adj_R2 = model_summary$adj.r.squared         # Adjusted R-squared
  ))
}

# Print the summary table
print(VPD_ambient_temp_summary)


#VPD
VPD_VPD_models <- list()
for (i in 1:12) {
  VPD_col <- paste0("VPD_", i)
  formula <- as.formula(paste("VPD ~", VPD_col))
  VPD_VPD_models[[VPD_col]] <- lm(formula, data = combined_dat)
}


summary(VPD_VPD_models[["VPD_1"]])
summary(VPD_VPD_models[["VPD_2"]])
summary(VPD_VPD_models[["VPD_3"]])
summary(VPD_VPD_models[["VPD_4"]])
summary(VPD_VPD_models[["VPD_5"]])
summary(VPD_VPD_models[["VPD_6"]])
summary(VPD_VPD_models[["VPD_7"]])
summary(VPD_VPD_models[["VPD_8"]])
summary(VPD_VPD_models[["VPD_9"]])
summary(VPD_VPD_models[["VPD_10"]])
summary(VPD_VPD_models[["VPD_11"]])
summary(VPD_VPD_models[["VPD_12"]])

plot(combined_dat$VPD_12, combined_dat$VPD,
     pch = 19, col = values,
     xlab = "VPD_12", ylab = "VPD",
     main = "VPD vs VPD_12")
abline(VPD_VPD_models[["VPD_12"]], col = "red", lwd = 2)

VPD_VPD_summary <- data.frame(Model = character(), Estimate = numeric(), Std_Error = numeric(), 
                              t_value = numeric(), p_value = numeric(), Adj_R2 = numeric())

# Loop through each month to summarize the models
for (i in 1:12) {
  # Extract model results
  model <- VPD_VPD_models[[paste0("VPD_", i)]]
  
  # Summary of model
  model_summary <- summary(model)
  
  # Store the relevant values in the results table
  VPD_VPD_summary <- rbind(VPD_VPD_summary, data.frame(
    Model = paste0("VPD_", i),
    Estimate = model_summary$coefficients[2, 1], # Estimate for the predictor
    Std_Error = model_summary$coefficients[2, 2], # Standard error of the predictor
    t_value = model_summary$coefficients[2, 3],  # t-value of the predictor
    p_value = model_summary$coefficients[2, 4],  # p-value of the predictor
    Adj_R2 = model_summary$adj.r.squared         # Adjusted R-squared
  ))
}

# Print the summary table
print(VPD_VPD_summary)


#Percent Veg Cover
VPD_Veg_models <- list()
for (i in 1:12) {
  VPD_col <- paste0("VPD_", i)
  formula <- as.formula(paste("percent_veg_cover ~", VPD_col))
  VPD_Veg_models[[VPD_col]] <- lm(formula, data = combined_dat)
}

summary(VPD_Veg_models[["VPD_1"]])
summary(VPD_Veg_models[["VPD_2"]])
summary(VPD_Veg_models[["VPD_3"]])
summary(VPD_Veg_models[["VPD_4"]])
summary(VPD_Veg_models[["VPD_5"]])
summary(VPD_Veg_models[["VPD_6"]])
summary(VPD_Veg_models[["VPD_7"]])
summary(VPD_Veg_models[["VPD_8"]])
summary(VPD_Veg_models[["VPD_9"]])
summary(VPD_Veg_models[["VPD_10"]])
summary(VPD_Veg_models[["VPD_11"]])
summary(VPD_Veg_models[["VPD_12"]])

plot(combined_dat$VPD_12, combined_dat$percent_veg_cover,
     pch = 19, col = values,
     xlab = "VPD_12", ylab = "percent_veg_cover",
     main = "percent_veg_cover vs VPD_12")
abline(VPD_Veg_models[["VPD_12"]], col = "red", lwd = 2)

VPD_Veg_summary <- data.frame(Model = character(), Estimate = numeric(), Std_Error = numeric(), 
                              t_value = numeric(), p_value = numeric(), Adj_R2 = numeric())

# Loop through each month to summarize the models
for (i in 1:12) {
  # Extract model results
  model <- VPD_Veg_models[[paste0("VPD_", i)]]
  
  # Summary of model
  model_summary <- summary(model)
  
  # Store the relevant values in the results table
  VPD_Veg_summary <- rbind(VPD_Veg_summary, data.frame(
    Model = paste0("VPD_", i),
    Estimate = model_summary$coefficients[2, 1], # Estimate for the predictor
    Std_Error = model_summary$coefficients[2, 2], # Standard error of the predictor
    t_value = model_summary$coefficients[2, 3],  # t-value of the predictor
    p_value = model_summary$coefficients[2, 4],  # p-value of the predictor
    Adj_R2 = model_summary$adj.r.squared         # Adjusted R-squared
  ))
}

# Print the summary table
print(VPD_Veg_summary)

#Combined results table
print(names(annual_precip_results_table))
print(names(srad_results))
print(names(srad_VPD_results))
print(names(srad_Veg_results))
print(names(vapr_ambient_temp_results))
print(names(vapr_VPD_results))
print(names(vapr_Veg_results))
print(names(tmax_ambient_temp_results))
print(names(tmax_VPD_results))
print(names(tmax_Veg_results))
print(names(tmin_ambient_temp_results))
print(names(tmin_VPD_results))
print(names(tmin_Veg_results))
print(names(tavg_ambient_temp_results))
print(names(tavg_VPD_results))
print(names(tavg_Veg_summary))
print(names(VPD_ambient_temp_summary))
print(names(VPD_VPD_summary))
print(names(VPD_Veg_summary))

colnames(srad_results)[1] <- "Response"
colnames(srad_VPD_results)[1] <- "Response"
colnames(srad_Veg_results)[1] <- "Response"
colnames(vapr_ambient_temp_results)[1] <- "Response"
colnames(vapr_VPD_results)[1] <- "Response"
colnames(vapr_Veg_results)[1] <- "Response"
colnames(tmax_ambient_temp_results)[1] <- "Response"
colnames(tmax_VPD_results)[1] <- "Response"
colnames(tmax_Veg_results)[1] <- "Response"
colnames(tmin_ambient_temp_results)[1] <- "Response"
colnames(tmin_VPD_results)[1] <- "Response"
colnames(tmin_Veg_results)[1] <- "Response"
colnames(tavg_ambient_temp_results)[1] <- "Response"
colnames(tavg_VPD_results)[1] <- "Response"
colnames(tavg_Veg_summary)[1] <- "Response"
colnames(VPD_ambient_temp_summary)[1] <- "Response"
colnames(VPD_VPD_summary)[1] <- "Response"
colnames(VPD_Veg_summary)[1] <- "Response"

annual_precip_results_table$model <- "Annual Precipitation"
srad_results$model <- "Solar Radiation"
srad_VPD_results$model <- "Solar Radiation and VPD"
srad_Veg_results$model <- "Solar Radiation and Vegetation"
vapr_ambient_temp_results$model <- "Vapor Pressure and Ambient Temperature"
vapr_VPD_results$model <- "Vapor Pressure and VPD"
vapr_Veg_results$model <- "Vapor Pressure and Vegetation"
tmax_ambient_temp_results$model <- "Max Temp and Ambient Temperature"
tmax_VPD_results$model <- "Max Temp and VPD"
tmax_Veg_results$model <- "Max Temp and Vegetation"
tmin_ambient_temp_results$model <- "Min Temp and Ambient Temperature"
tmin_VPD_results$model <- "Min Temp and VPD"
tmin_Veg_results$model <- "Min Temp and Vegetation"
tavg_ambient_temp_results$model <- "Average Temp and Ambient Temperature"
tavg_VPD_results$model <- "Average Temp and VPD"
tavg_Veg_summary$model <- "Average Temp and Vegetation"
VPD_ambient_temp_summary$model <- "VPD and Ambient Temperature"
VPD_VPD_summary$model <- "VPD and VPD"
VPD_Veg_summary$model <- "VPD and Vegetation"

combined_microVmacro_summary <- rbind(annual_precip_results_table,
                                      srad_results,
                                      srad_VPD_results,
                                      srad_Veg_results,
                                      vapr_ambient_temp_results,
                                      vapr_VPD_results,
                                      vapr_Veg_results,
                                      tmax_ambient_temp_results,
                                      tmax_VPD_results,
                                      tmax_Veg_results,
                                      tmin_ambient_temp_results,
                                      tmin_VPD_results,
                                      tmin_Veg_results,
                                      tavg_ambient_temp_results,
                                      tavg_VPD_results,
                                      tavg_Veg_summary,
                                      VPD_ambient_temp_summary,
                                      VPD_VPD_summary,
                                      VPD_Veg_summary)

write.csv(combined_microVmacro_summary, "combined_microVmacro_summary.csv", row.names = FALSE)

#Macroclimate PCA by habitat
MacroClim_PCA <- location_dat %>%
  dplyr::select(NDVI, NDMI, annual_precip, monthly_avg_srad, 
                monthly_avg_precip, monthly_avg_vapr, 
                monthly_avg_tmax, monthly_avg_tmin, monthly_avg_tavg, VPD_avg)
scaled_MacroClim_PCA <- scale(MacroClim_PCA)

MacroClim_PCA_result <- prcomp(scaled_MacroClim_PCA)

pca_scores <- as.data.frame(MacroClim_PCA_result$x)
pca_scores$Site <- location_dat$Site  # if you have a site column
pca_scores$Habitat <- location_dat$Habitat

habitat_colors <- c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")


ggplot(pca_scores, aes(x = PC1, y = PC2, color = Habitat)) +
  geom_point(size = 3) +
  stat_ellipse(aes(color = Habitat), level = 0.95) +  # Add ellipses (95% confidence)
  scale_color_manual(values = habitat_colors) +
  theme_minimal() +
  labs(title = "PCA of Macroclimate Variables",
       x = paste0("PC1 (", round(summary(MacroClim_PCA_result)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(MacroClim_PCA_result)$importance[2,2] * 100, 1), "%)"),
       color = "Habitat")

#Box Plot of macroclimate variable averages by habitat
long_data <- location_dat %>%
  dplyr::select(NDVI, NDMI, annual_precip, monthly_avg_srad, 
                monthly_avg_precip, monthly_avg_vapr, 
                monthly_avg_tmax, monthly_avg_tmin, monthly_avg_tavg, Habitat, VPD_avg) %>%
  pivot_longer(cols = -Habitat, names_to = "Variable", values_to = "Value")

ggplot(long_data, aes(x = Habitat, y = Value, fill = Habitat)) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "free_y") +  
  scale_fill_manual(values = habitat_colors) +  
  theme_minimal() +
  labs(title = "Boxplots of Macroclimate Variables by Habitat",
       x = "Habitat",
       y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Comparing Microclimate Model fit to macroclimate model fit for ind. lizards
#Microclimates
CEWL_ambient_temp_model <- lm(CEWL ~ ambient_temp, data = na.omit(combined_dat))
summary(CEWL_ambient_temp_model)

CEWL_VPD_model <- lm(CEWL ~ VPD, data = na.omit(combined_dat))
summary(CEWL_VPD_model)

CEWL_Veg_model <- lm(CEWL ~ percent_veg_cover, data = na.omit(combined_dat))
summary(CEWL_Veg_model)

#Macroclimates
#Solar Radiation
srad_models_CEWL <- list()
# Loop through srad_1 to srad_12 and store each model
for (i in 1:12) {
  srad_col <- paste0("srad_", i)
  formula <- as.formula(paste("CEWL ~", srad_col))  # Change ambient_temp to CEWL
  model_data <- combined_dat[complete.cases(combined_dat[, c("CEWL", srad_col)]), ]
  srad_models_CEWL[[srad_col]] <- lm(formula, data = model_data)
}

summary(srad_models_CEWL[["srad_1"]])
summary(srad_models_CEWL[["srad_2"]])
summary(srad_models_CEWL[["srad_3"]])
summary(srad_models_CEWL[["srad_4"]])
summary(srad_models_CEWL[["srad_5"]])
summary(srad_models_CEWL[["srad_6"]])
summary(srad_models_CEWL[["srad_7"]])
summary(srad_models_CEWL[["srad_8"]])
summary(srad_models_CEWL[["srad_9"]])
summary(srad_models_CEWL[["srad_10"]])
summary(srad_models_CEWL[["srad_11"]])
summary(srad_models_CEWL[["srad_12"]])

#Vapor Pressure
vapr_models_CEWL <- list()
# Loop through vapr_1 to vapr_12 and store each model
for (i in 1:12) {
  vapr_col <- paste0("vapr_", i)
  formula <- as.formula(paste("CEWL ~", vapr_col))  # Change ambient_temp to CEWL
  model_data <- combined_dat[complete.cases(combined_dat[, c("CEWL", vapr_col)]), ]
  vapr_models_CEWL[[vapr_col]] <- lm(formula, data = model_data)
}

summary(vapr_models_CEWL[["vapr_1"]])
summary(vapr_models_CEWL[["vapr_2"]])
summary(vapr_models_CEWL[["vapr_3"]])
summary(vapr_models_CEWL[["vapr_4"]])
summary(vapr_models_CEWL[["vapr_5"]])
summary(vapr_models_CEWL[["vapr_6"]])
summary(vapr_models_CEWL[["vapr_7"]])
summary(vapr_models_CEWL[["vapr_8"]])
summary(vapr_models_CEWL[["vapr_9"]])
summary(vapr_models_CEWL[["vapr_10"]])
summary(vapr_models_CEWL[["vapr_11"]])
summary(vapr_models_CEWL[["vapr_12"]])

#Maximum Temp
tmax_models_CEWL <- list()
# Loop through tmax_1 to tmax_12 and store each model
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  formula <- as.formula(paste("CEWL ~", tmax_col))  # Change ambient_temp to CEWL
  model_data <- combined_dat[complete.cases(combined_dat[, c("CEWL", tmax_col)]), ]
  tmax_models_CEWL[[tmax_col]] <- lm(formula, data = model_data)
}

summary(tmax_models_CEWL[["tmax_1"]])
summary(tmax_models_CEWL[["tmax_2"]])
summary(tmax_models_CEWL[["tmax_3"]])
summary(tmax_models_CEWL[["tmax_4"]])
summary(tmax_models_CEWL[["tmax_5"]])
summary(tmax_models_CEWL[["tmax_6"]])
summary(tmax_models_CEWL[["tmax_7"]])
summary(tmax_models_CEWL[["tmax_8"]])
summary(tmax_models_CEWL[["tmax_9"]])
summary(tmax_models_CEWL[["tmax_10"]])
summary(tmax_models_CEWL[["tmax_11"]])
summary(tmax_models_CEWL[["tmax_12"]])

#Minimum Temp
tmin_models_CEWL <- list()
# Loop through tmin_1 to tmin_12 and store each model
for (i in 1:12) {
  tmin_col <- paste0("tmin_", i)
  formula <- as.formula(paste("CEWL ~", tmin_col))  # Change ambient_temp to CEWL
  model_data <- combined_dat[complete.cases(combined_dat[, c("CEWL", tmin_col)]), ]
  tmin_models_CEWL[[tmin_col]] <- lm(formula, data = model_data)
}

summary(tmin_models_CEWL[["tmin_1"]])
summary(tmin_models_CEWL[["tmin_2"]])
summary(tmin_models_CEWL[["tmin_3"]])
summary(tmin_models_CEWL[["tmin_4"]])
summary(tmin_models_CEWL[["tmin_5"]])
summary(tmin_models_CEWL[["tmin_6"]])
summary(tmin_models_CEWL[["tmin_7"]])
summary(tmin_models_CEWL[["tmin_8"]])
summary(tmin_models_CEWL[["tmin_9"]])
summary(tmin_models_CEWL[["tmin_10"]])
summary(tmin_models_CEWL[["tmin_11"]])
summary(tmin_models_CEWL[["tmin_12"]])

#Average Temp
tavg_models_CEWL <- list()
# Loop through tavg_1 to tavg_12 and store each model
for (i in 1:12) {
  tavg_col <- paste0("tavg_", i)
  formula <- as.formula(paste("CEWL ~", tavg_col))  # Change ambient_temp to CEWL
  model_data <- combined_dat[complete.cases(combined_dat[, c("CEWL", tavg_col)]), ]
  tavg_models_CEWL[[tavg_col]] <- lm(formula, data = model_data)
}

summary(tavg_models_CEWL[["tavg_1"]])
summary(tavg_models_CEWL[["tavg_2"]])
summary(tavg_models_CEWL[["tavg_3"]])
summary(tavg_models_CEWL[["tavg_4"]])
summary(tavg_models_CEWL[["tavg_5"]])
summary(tavg_models_CEWL[["tavg_6"]])
summary(tavg_models_CEWL[["tavg_7"]])
summary(tavg_models_CEWL[["tavg_8"]])
summary(tavg_models_CEWL[["tavg_9"]])
summary(tavg_models_CEWL[["tavg_10"]])
summary(tavg_models_CEWL[["tavg_11"]])
summary(tavg_models_CEWL[["tavg_12"]])

#Macro_VPD
VPD_models_CEWL <- list()
# Loop through VPD_1 to VPD_12 and store each model
for (i in 1:12) {
  VPD_col <- paste0("VPD_", i)
  formula <- as.formula(paste("CEWL ~", VPD_col))  # Change ambient_temp to CEWL
  model_data <- combined_dat[complete.cases(combined_dat[, c("CEWL", VPD_col)]), ]
  VPD_models_CEWL[[VPD_col]] <- lm(formula, data = model_data)
}

summary(VPD_models_CEWL[["VPD_1"]])
summary(VPD_models_CEWL[["VPD_2"]])
summary(VPD_models_CEWL[["VPD_3"]])
summary(VPD_models_CEWL[["VPD_4"]])
summary(VPD_models_CEWL[["VPD_5"]])
summary(VPD_models_CEWL[["VPD_6"]])
summary(VPD_models_CEWL[["VPD_7"]])
summary(VPD_models_CEWL[["VPD_8"]])
summary(VPD_models_CEWL[["VPD_9"]])
summary(VPD_models_CEWL[["VPD_10"]])
summary(VPD_models_CEWL[["VPD_11"]])
summary(VPD_models_CEWL[["VPD_12"]])


#Extract and order by R2
extract_r2 <- function(models_list, prefix) {
  data.frame(
    model = names(models_list),
    R2 = sapply(models_list, function(m) summary(m)$r.squared),
    type = prefix
  )
}

micro_models <- list(
  "ambient_temp" = CEWL_ambient_temp_model,
  "VPD" = CEWL_VPD_model,
  "percent_veg_cover" = CEWL_Veg_model
)
micro_r2 <- extract_r2(micro_models, "micro")

srad_r2 <- extract_r2(srad_models_CEWL, "srad")
vapr_r2 <- extract_r2(vapr_models_CEWL, "vapr")
tmax_r2 <- extract_r2(tmax_models_CEWL, "tmax")
tmin_r2 <- extract_r2(tmin_models_CEWL, "tmin")
tavg_r2 <- extract_r2(tavg_models_CEWL, "tavg")
macroVPD_r2 <- extract_r2(VPD_models_CEWL, "macroVPD")

all_r2 <- rbind(micro_r2, srad_r2, vapr_r2, tmax_r2, tmin_r2, tavg_r2, macroVPD_r2)
all_r2 <- all_r2[order(-all_r2$R2), ]

#write.csv(all_r2, "all_r2.csv", row.names = FALSE)



#Microclimate and CEWL ~ Site
#CEWL
combined_dat2 <- combined_dat %>%
  filter(Site != "Zas") %>%
  filter(!is.na(CEWL) & !is.na(Site))

site_levels <- combined_dat2 %>%
  dplyr::select(Site) %>%
  distinct() %>%
  pull(Site)

combined_dat2$Site <- factor(combined_dat2$Site, levels = site_levels)

CEWL_site_Model <- glm(CEWL ~ Site, data =combined_dat2)
summary(CEWL_site_Model)
pairs_CEWL <- pairs(emmeans(CEWL_site_Model, ~ Site)) %>%
  as.data.frame()
pairs_CEWL <- pairs_CEWL %>%
  mutate(
    site1 = sub(" -.*", "", contrast),
    site2 = sub(".*- ", "", contrast)
  )
site_habitat_lookup <- dplyr::select(combined_dat, Site, Habitat.x) %>%
  distinct()
pairs_CEWL <- pairs_CEWL %>%
  left_join(site_habitat_lookup, by = c("site1" = "Site")) %>%
  rename(Habitat1 = Habitat.x) %>%
  left_join(site_habitat_lookup, by = c("site2" = "Site")) %>%
  rename(Habitat2 = Habitat.x)
pairs_CEWL <- pairs_CEWL %>%
  mutate(
    ComparisonType = case_when(
      Habitat1 == "wall" & Habitat2 == "wall" ~ "wall vs wall",
      Habitat1 == "nonwall" & Habitat2 == "nonwall" ~ "nonwall vs nonwall",
      TRUE ~ "wall vs nonwall"
    )
  )
pairs_CEWL
pairs_CEWL %>%
  filter(p.value < 0.05) %>%
  count(ComparisonType)

CEWL_emm_df_site <- emmeans(CEWL_site_Model, ~ Site) %>%
  as.data.frame()

site_habitat_lookup <- combined_dat2 %>%
  dplyr::select(Site, Habitat.x) %>%
  distinct()

site_levels_by_habitat <- site_habitat_lookup %>%
  mutate(Habitat.x = factor(Habitat.x, levels = c("wall", "nonwall"))) %>%
  arrange(Habitat.x) %>%
  pull(Site)

combined_dat2$Site <- factor(combined_dat2$Site, levels = site_levels_by_habitat)
CEWL_emm_df_site$Site <- factor(CEWL_emm_df_site$Site, levels = site_levels_by_habitat)

# Add habitat to CEWL_emm_df_site if not already done
if (!"Habitat.x" %in% colnames(CEWL_emm_df_site)) {
  CEWL_emm_df_site <- CEWL_emm_df_site %>%
    left_join(site_habitat_lookup, by = "Site")
}

# Ensure that the levels of 'Habitat.x' in combined_dat match the manual scales
combined_dat2$Habitat.x <- factor(combined_dat2$Habitat.x, levels = c("nonwall", "wall"))
# Add significance to pairs_CEWL based on p-value < 0.05
pairs_CEWL <- pairs_CEWL %>%
  mutate(Significant = ifelse(p.value < 0.05, "*", ""))  # Only mark significant comparisons

# Remove NA values from pairs_CEWL for the plot
pairs_CEWL <- pairs_CEWL %>% filter(Significant != "")


CEWL_emm_df_site <- CEWL_emm_df_site %>%
  filter(Site %in% levels(combined_dat2$Site)) %>%
  droplevels()

CEWL_emm_df_site$Site <- factor(CEWL_emm_df_site$Site, 
                                levels = levels(combined_dat2$Site))

# Plot
CEWL_box_site <- ggplot(combined_dat2, 
                        aes(x = Site, y = CEWL, fill = Habitat.x)) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = Habitat.x), 
              width = 0.2, size = 2.5, alpha = 0.4, stroke = 0.5, color = "black") +
  
  geom_errorbar(data = CEWL_emm_df_site, 
                aes(x = Site, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +
  
  geom_point(data = CEWL_emm_df_site, 
             aes(x = Site, y = emmean, fill = Habitat.x), 
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +
  
  # Add significance labels (only for significant comparisons)
  #geom_text(data = pairs_CEWL, 
            #aes(x = site1, y = 0.9 * max(combined_dat$CEWL), label = Significant),  # Adjust y position as needed
            #inherit.aes = FALSE, size = 5, color = "red", fontface = "bold") +
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c(24, 21)) +
  
  labs(
    title = "CEWL Variation Across Sites with emmeans",
    x = "Site",
    y = "CEWL (mm)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = "none", shape = "none")

# Show plot
CEWL_box_site



#CEWL vs VPD avg by site
site_summary <- cewl.dat %>%
  group_by(Site, Habitat) %>%  # Keep grouping vars you want
  summarise(
    VPD_mean = mean(VPD, na.rm = TRUE),
    CEWL_mean = mean(CEWL, na.rm = TRUE),
    CEWL_se = sd(CEWL, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


scatter_avg <- ggplot(cewl.dat, 
                      aes(x = VPD, y = CEWL, 
                          color = as.factor(Habitat), 
                          shape = as.factor(Sex))) +
  
  # Raw points
  #geom_point(aes(fill = as.factor(Habitat)), 
             #size = 2.5, stroke = 0.4, color = "black", alpha = 0.5) +
  
  # Site averages (8 points)
  geom_point(data = site_summary, 
             aes(x = VPD_mean, y = CEWL_mean, fill = as.factor(Habitat)), 
             size = 6, shape = 21, color = "black", inherit.aes = FALSE) +
  
  # Error bars (±SE)
  geom_errorbar(data = site_summary, 
                aes(x = VPD_mean, ymin = CEWL_mean - CEWL_se, ymax = CEWL_mean + CEWL_se), 
                width = 0.05, color = "black", size = 0.8, inherit.aes = FALSE) +
  
  # Single overall trend line (no habitat grouping)
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.5, inherit.aes = FALSE,
              aes(x = VPD, y = CEWL)) +
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  
  labs(
    title = "CEWL vs VPD with Site Means and Overall Trend",
    x = "VPD (kPa)", 
    y = expression(CEWL ~ (g/m^2*h))
  ) +
  scale_shape_manual(values = c(24, 21)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10)
  )

scatter_avg


combined_CEWL_Plot <- CEWL_box_site +
  CEWL_box +
  scatter_avg +
  CEWL_VPD_PLot +
  plot_layout(ncol = 2)
  

#Temp
temp_site_Model <- glm(ambient_temp ~ Site, data = combined_dat)
summary(temp_site_Model)
pairs_ambient_temp <- pairs(emmeans(temp_site_Model, ~ Site)) %>%
  as.data.frame()
pairs_ambient_temp <- pairs_ambient_temp %>%
  mutate(
    site1 = sub(" -.*", "", contrast),
    site2 = sub(".*- ", "", contrast)
  )
site_habitat_lookup <- dplyr::select(combined_dat, Site, Habitat.x) %>%
  distinct()
pairs_ambient_temp <- pairs_ambient_temp %>%
  left_join(site_habitat_lookup, by = c("site1" = "Site")) %>%
  rename(Habitat1 = Habitat.x) %>%
  left_join(site_habitat_lookup, by = c("site2" = "Site")) %>%
  rename(Habitat2 = Habitat.x)
pairs_ambient_temp <- pairs_ambient_temp %>%
  mutate(
    ComparisonType = case_when(
      Habitat1 == "wall" & Habitat2 == "wall" ~ "wall vs wall",
      Habitat1 == "nonwall" & Habitat2 == "nonwall" ~ "nonwall vs nonwall",
      TRUE ~ "wall vs nonwall"
    )
  )
pairs_ambient_temp
pairs_ambient_temp %>%
  filter(p.value < 0.05) %>%
  count(ComparisonType)

# Get the emm data for the model
ambient_temp_emm_df_site <- emmeans(temp_site_Model, ~ Site) %>%
  as.data.frame()

# Ensure the Site order is based on Habitat type (wall/nonwall)
site_levels_by_habitat <- site_habitat_lookup %>%
  arrange(desc(Habitat.x)) %>%   # Change to arrange(Habitat.x) if you want nonwall first
  pull(Site)

combined_dat$Site <- factor(combined_dat$Site, levels = site_levels_by_habitat)
ambient_temp_emm_df_site$Site <- factor(ambient_temp_emm_df_site$Site, levels = site_levels_by_habitat)

# Add Habitat to the emm data
if (!"Habitat.x" %in% colnames(ambient_temp_emm_df_site)) {
  ambient_temp_emm_df_site <- ambient_temp_emm_df_site %>%
    left_join(site_habitat_lookup, by = "Site")
}

# Add significance stars for p-value < 0.05
pairs_ambient_temp <- pairs_ambient_temp %>%
  mutate(Significant = ifelse(p.value < 0.05, "*", ""))  # Only mark significant comparisons

# Remove non-significant comparisons
pairs_ambient_temp <- pairs_ambient_temp %>% filter(Significant != "")

# Plot
ambient_temp_box_site <- ggplot(combined_dat, 
                                aes(x = Site, y = ambient_temp, fill = Habitat.x)) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = Habitat.x), 
              width = 0.2, size = 5, alpha = 0.4, stroke = 0.5, color = "black") +
  
  geom_errorbar(data = ambient_temp_emm_df_site, 
                aes(x = Site, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +
  
  geom_point(data = ambient_temp_emm_df_site, 
             aes(x = Site, y = emmean, fill = Habitat.x), 
             size = 10, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +
  
  # Add significance labels (only for significant comparisons)
  #geom_text(data = pairs_ambient_temp, 
           # aes(x = site1, y = 0.9 * max(combined_dat$ambient_temp), label = Significant),  # Adjust y position as needed
            #inherit.aes = FALSE, size = 5, color = "red", fontface = "bold") +
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c(24, 21)) +
  
  labs(
    title = "Ambient Temperature Variation Across Sites with emmeans",
    x = "Site",
    y = "Ambient Temperature (°C)",
    shape = "Sex"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = "none", shape = "none")

# Show plot
ambient_temp_box_site


#Veg
Veg_site_Model <- glm(percent_veg_cover ~ Site, data = combined_dat)
summary(Veg_site_Model)
pairs_percent_veg_cover <- pairs(emmeans(Veg_site_Model, ~ Site)) %>%
  as.data.frame()
pairs_percent_veg_cover <- pairs_percent_veg_cover %>%
  mutate(
    site1 = sub(" -.*", "", contrast),
    site2 = sub(".*- ", "", contrast)
  )
site_habitat_lookup <- dplyr::select(combined_dat, Site, Habitat.x) %>%
  distinct()
pairs_percent_veg_cover <- pairs_percent_veg_cover %>%
  left_join(site_habitat_lookup, by = c("site1" = "Site")) %>%
  rename(Habitat1 = Habitat.x) %>%
  left_join(site_habitat_lookup, by = c("site2" = "Site")) %>%
  rename(Habitat2 = Habitat.x)
pairs_percent_veg_cover <- pairs_percent_veg_cover %>%
  mutate(
    ComparisonType = case_when(
      Habitat1 == "wall" & Habitat2 == "wall" ~ "wall vs wall",
      Habitat1 == "nonwall" & Habitat2 == "nonwall" ~ "nonwall vs nonwall",
      TRUE ~ "wall vs nonwall"
    )
  )
pairs_percent_veg_cover
pairs_percent_veg_cover %>%
  filter(p.value < 0.05) %>%
  count(ComparisonType)

percent_veg_cover_emm_df_site <- emmeans(Veg_site_Model, ~ Site) %>%
  as.data.frame()

# Ensure the Site order is based on Habitat type (wall/nonwall)
site_levels_by_habitat <- site_habitat_lookup %>%
  arrange(desc(Habitat.x)) %>%   # Change to arrange(Habitat.x) if you want nonwall first
  pull(Site)

combined_dat$Site <- factor(combined_dat$Site, levels = site_levels_by_habitat)
percent_veg_cover_emm_df_site$Site <- factor(percent_veg_cover_emm_df_site$Site, levels = site_levels_by_habitat)

# Add Habitat to the emm data
if (!"Habitat.x" %in% colnames(percent_veg_cover_emm_df_site)) {
  percent_veg_cover_emm_df_site <- percent_veg_cover_emm_df_site %>%
    left_join(site_habitat_lookup, by = "Site")
}

# Plot with means and error bars
percent_veg_cover_box_site <- ggplot(combined_dat, 
                                     aes(x = Site, y = percent_veg_cover, fill = Habitat.x)) +
  
  geom_jitter(aes(shape = as.factor(Sex), fill = Habitat.x), 
              width = 0.2, size = 5, alpha = 0.4, stroke = 0.5, color = "black") +  # Jitter plot
  
  # Add error bars (confidence intervals from emm)
  geom_errorbar(data = percent_veg_cover_emm_df_site, 
                aes(x = Site, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +
  
  # Add mean points
  geom_point(data = percent_veg_cover_emm_df_site, 
             aes(x = Site, y = emmean, fill = Habitat.x), 
             size = 10, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +
  
  # Add significance labels (only for significant comparisons)
  #geom_text(data = pairs_percent_veg_cover, 
           # aes(x = site1, y = 0.9 * max(combined_dat$percent_veg_cover), label = "*"),  # Adjust y position as needed
            #inherit.aes = FALSE, size = 5, color = "red", fontface = "bold") +
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c(24, 21)) +
  
  labs(
    title = "Percent Vegetation Cover Variation Across Sites with emmeans",
    x = "Site",
    y = "Percent Vegetation Cover"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = "none", shape = "none")

# Show plot
percent_veg_cover_box_site


#VPD
VPD_site_Model <- glm(VPD ~ Site, data = combined_dat)
summary(VPD_site_Model)
pairs_VPD <- pairs(emmeans(VPD_site_Model, ~ Site)) %>%
  as.data.frame()
pairs_VPD <- pairs_VPD %>%
  mutate(
    site1 = sub(" -.*", "", contrast),
    site2 = sub(".*- ", "", contrast)
  )
site_habitat_lookup <- dplyr::select(combined_dat, Site, Habitat.x) %>%
  distinct()
pairs_VPD <- pairs_VPD %>%
  left_join(site_habitat_lookup, by = c("site1" = "Site")) %>%
  rename(Habitat1 = Habitat.x) %>%
  left_join(site_habitat_lookup, by = c("site2" = "Site")) %>%
  rename(Habitat2 = Habitat.x)
pairs_VPD <- pairs_VPD %>%
  mutate(
    ComparisonType = case_when(
      Habitat1 == "wall" & Habitat2 == "wall" ~ "wall vs wall",
      Habitat1 == "nonwall" & Habitat2 == "nonwall" ~ "nonwall vs nonwall",
      TRUE ~ "wall vs nonwall"
    )
  )
pairs_VPD
pairs_VPD %>%
  filter(p.value < 0.05) %>%
  count(ComparisonType)

VPD_emm_df_site <- emmeans(VPD_site_Model, ~ Site) %>%
  as.data.frame()

# Set site factor levels based on habitat type (to control plotting order)
site_levels_by_habitat <- site_habitat_lookup %>%
  arrange(desc(Habitat.x)) %>%  # Use `arrange(Habitat.x)` if you want nonwalls first
  pull(Site)

combined_dat$Site <- factor(combined_dat$Site, levels = site_levels_by_habitat)
VPD_emm_df_site$Site <- factor(VPD_emm_df_site$Site, levels = site_levels_by_habitat)

# Add habitat type to emmeans output
if (!"Habitat.x" %in% colnames(VPD_emm_df_site)) {
  VPD_emm_df_site <- VPD_emm_df_site %>%
    left_join(site_habitat_lookup, by = "Site")
}

# Plot
VPD_box_site <- ggplot(combined_dat, 
                       aes(x = Site, y = VPD, fill = Habitat.x)) +
  
  # Jittered raw data points
  geom_jitter(aes(shape = as.factor(Sex), fill = Habitat.x),
              width = 0.2, size = 5, alpha = 0.4, stroke = 0.5, color = "black") +
  
  # Error bars from emmeans
  geom_errorbar(data = VPD_emm_df_site, 
                aes(x = Site, ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +
  
  # Mean points
  geom_point(data = VPD_emm_df_site, 
             aes(x = Site, y = emmean, fill = Habitat.x), 
             size = 10, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +
  
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c(24, 21)) +
  
  labs(
    title = "VPD Variation Across Sites with emmeans",
    x = "Site",
    y = "Vapor Pressure Deficit (VPD)"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = "none", shape = "none")

# Show the plot
VPD_box_site


Site_combined_plot <- ambient_temp_box_site + 
  VPD_box_site + 
  percent_veg_cover_box_site + 
  plot_layout(ncol = 3)

Site_combined_plot


#Combine site and habitat microclim plots
Final_Combined_EnviroPlot <- p1 + 
  ambient_temp_box_site + 
  p3 + 
  VPD_box_site + 
  p4 + 
  percent_veg_cover_box_site + 
  plot_layout(ncol = 2)

###############################
#Rasters for Maps
#Annual Precip
#total_precip <- sum(precip_1, precip_2, precip_3, precip_4, precip_5, precip_6, precip_7, precip_8, precip_9, precip_10, precip_11, precip_12)
#writeRaster(total_precip, "total_precipitation.tif", overwrite = TRUE)
total_precip <- raster("~/Dropbox/SDSU/Brock_Lab/UrbanVSNonUrban_Lizards/Data/World_Clim/total_precipitation.tif")

climcol <- colorRampPalette(c("purple", "blue", "skyblue", "green", "lightgreen", "yellow", "orange", "red", "darkred"))
climcol_rev <- colorRampPalette(rev(c("purple", "blue", "skyblue", "green", "lightgreen", "yellow", "orange", "red", "darkred")))
habitat_colors <- c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")


plot(total_precip, xlim=c(25.33, 25.6), ylim=c(36.9, 37.21),
     col=climcol_rev(100), axes=FALSE, ann=FALSE, box=FALSE)
#box()
# Add nonwall (rural) sites - circle
points(location_dat$Long[location_dat$Habitat == "nonwall"],
       location_dat$Lat[location_dat$Habitat == "nonwall"],
       pch = 21, bg = habitat_colors["nonwall"], col = "black", cex = 1.5)

# Add wall (urban) sites - triangle
points(location_dat$Long[location_dat$Habitat == "wall"],
       location_dat$Lat[location_dat$Habitat == "wall"],
       pch = 24, bg = habitat_colors["wall"], col = "black", cex = 1.5)

#text(x = 25.4, y = 37.19, labels = "F. Annual Precipitation", font = 2)



#Tmax
#plot(tmax_5)
plot(tmax_5, xlim=c(25.33, 25.6), ylim=c(36.9, 37.21), zlim=c(19, 28),
     col=climcol(100), axes=FALSE, ann=FALSE, box=FALSE)
#box()  # add a simple box outline

# Add nonwall (rural) sites - circle
points(location_dat$Long[location_dat$Habitat == "nonwall"],
       location_dat$Lat[location_dat$Habitat == "nonwall"],
       pch = 21, bg = habitat_colors["nonwall"], col = "black", cex = 1.5)

# Add wall (urban) sites - triangle
points(location_dat$Long[location_dat$Habitat == "wall"],
       location_dat$Lat[location_dat$Habitat == "wall"],
       pch = 24, bg = habitat_colors["wall"], col = "black", cex = 1.5)

# Add internal title label (adjust x/y as needed for position)
#text(x = 25.375, y = 37.19, labels = "B. May Temperature", font = 2)



#tmax_6
plot(tmax_6, xlim=c(25.33, 25.6), ylim=c(36.9, 37.21), zlim=c(19, 28),
     col=climcol(100), axes=FALSE, ann=FALSE, box=FALSE)
#box()  # add a simple box outline

# Add nonwall (rural) sites - circle
points(location_dat$Long[location_dat$Habitat == "nonwall"],
       location_dat$Lat[location_dat$Habitat == "nonwall"],
       pch = 21, bg = habitat_colors["nonwall"], col = "black", cex = 1.5)

# Add wall (urban) sites - triangle
points(location_dat$Long[location_dat$Habitat == "wall"],
       location_dat$Lat[location_dat$Habitat == "wall"],
       pch = 24, bg = habitat_colors["wall"], col = "black", cex = 1.5)

# Add internal title label (adjust x/y as needed for position)
#text(x = 25.379, y = 37.19, labels = "C. June Temperature", font = 2)


#VPD
vpd_extent <- extent(25.31, 25.6, 36.9, 37.21)


#May
tmax5_crop <- crop(tmax_5, vpd_extent)
vapr5_crop <- crop(vapr_5, vpd_extent)

VPD_May <- (0.611 * exp((17.502 * tmax5_crop) / (tmax5_crop + 240.97))) - vapr5_crop

plot(VPD_May, xlim=c(25.33, 25.6), ylim=c(36.9, 37.21), zlim=c(1, 2), col=climcol(100), axes=FALSE, ann=FALSE, box=FALSE)
#box()  # Add a simple box around the map

# Add nonwall (rural) sites - circle
points(location_dat$Long[location_dat$Habitat == "nonwall"],
       location_dat$Lat[location_dat$Habitat == "nonwall"],
       pch = 21, bg = habitat_colors["nonwall"], col = "black", cex = 1.5)

# Add wall (urban) sites - triangle
points(location_dat$Long[location_dat$Habitat == "wall"],
       location_dat$Lat[location_dat$Habitat == "wall"],
       pch = 24, bg = habitat_colors["wall"], col = "black", cex = 1.5)

# Add internal title label (adjust x and y to your map layout)
#text(x = 25.33, y = 37.19, labels = "D. May VPD", font = 2)


#June
tmax6_crop <- crop(tmax_6, vpd_extent)
vapr6_crop <- crop(vapr_6, vpd_extent)

# Calculate VPD for June
VPD_June <- (0.611 * exp((17.502 * tmax6_crop) / (tmax6_crop + 240.97))) - vapr6_crop

# Plot
plot(VPD_June, xlim=c(25.33, 25.6), ylim=c(36.9, 37.21), zlim=c(1, 2), col=climcol(100), axes=FALSE, ann=FALSE, box=FALSE)
#box()  # Add a simple box around the map

# Add nonwall (rural) sites - circle
points(location_dat$Long[location_dat$Habitat == "nonwall"],
       location_dat$Lat[location_dat$Habitat == "nonwall"],
       pch = 21, bg = habitat_colors["nonwall"], col = "black", cex = 1.5)

# Add wall (urban) sites - triangle
points(location_dat$Long[location_dat$Habitat == "wall"],
       location_dat$Lat[location_dat$Habitat == "wall"],
       pch = 24, bg = habitat_colors["wall"], col = "black", cex = 1.5)

# Add internal title label (adjust x and y to your map layout)
#text(x = 25.34, y = 37.19, labels = "D. June VPD", font = 2)



##############################################################################################################################################################################################################################################################################################################################################################################################################################################
#Population Level Analyses- Taking average values of traits and climates per site and analyzing differences between populations rather than between individual lizards
#Sample size per site
site_counts <- dat %>%
  count(Site, Habitat)  # count by both Site *and* Habitat
print(site_counts)


#Create Site avg dataset
avg_dat <- dat %>%
  group_by(Site) %>%
  summarise(
    Habitat = first(Habitat),  # Assuming one habitat type per site
    
    Density_mean = mean(Density, na.rm = TRUE),
    Density_sd = sd(Density, na.rm = TRUE),
    
    Mass_mean = mean(Mass, na.rm = TRUE),
    Mass_sd = sd(Mass, na.rm = TRUE),
    
    SVL_mean = mean(SVL, na.rm = TRUE),
    SVL_sd = sd(SVL, na.rm = TRUE),
    
    head_width_mean = mean(HeadWidth, na.rm = TRUE),
    head_width_sd = sd(HeadWidth, na.rm = TRUE),
    
    head_depth_mean = mean(HeadDepth, na.rm = TRUE),
    head_depth_sd = sd(HeadDepth, na.rm = TRUE),
    
    head_length_mean = mean(HeadLength, na.rm = TRUE),
    head_length_sd = sd(HeadLength, na.rm = TRUE),
    
    femur_length_mean = mean(FemurLength, na.rm = TRUE),
    femur_length_sd = sd(FemurLength, na.rm = TRUE),
    
    tibia_length_mean = mean(TibiaLength, na.rm = TRUE),
    tibia_length_sd = sd(TibiaLength, na.rm = TRUE),
    
    bicep_length_mean = mean(BicepLength, na.rm = TRUE),
    bicep_length_sd = sd(BicepLength, na.rm = TRUE),
    
    forearm_length_mean = mean(ForearmLength, na.rm = TRUE),
    forearm_length_sd = sd(ForearmLength, na.rm = TRUE),
    
    CEWL_mean = mean(CEWL, na.rm = TRUE),
    CEWL_sd = sd(CEWL, na.rm = TRUE),
    
    msmt_temp_C_mean = mean(msmt_temp_C, na.rm = TRUE),
    msmt_temp_C_sd = sd(msmt_temp_C, na.rm = TRUE),
    
    msmt_RH_percent_mean = mean(msmt_RH_percent, na.rm = TRUE),
    msmt_RH_percent_sd = sd(msmt_RH_percent, na.rm = TRUE),
    
    msmt_VPD_kPa_mean = mean(msmt_VPD_kPa, na.rm = TRUE),
    msmt_VPD_kPa_sd = sd(msmt_VPD_kPa, na.rm = TRUE),
    
    ambient_percent_rh_mean = mean(ambient_percent_rh, na.rm = TRUE),
    ambient_percent_rh_sd = sd(ambient_percent_rh, na.rm = TRUE),
    
    ambient_temp_mean = mean(ambient_temp, na.rm = TRUE),
    ambient_temp_sd = sd(ambient_temp, na.rm = TRUE),
    
    percent_veg_cover_mean = mean(percent_veg_cover, na.rm = TRUE),
    percent_veg_cover_sd = sd(percent_veg_cover, na.rm = TRUE),
    
    LimbRatio_mean = mean(LimbRatio, na.rm = TRUE),
    LimbRatio_sd = sd(LimbRatio, na.rm = TRUE),
    
    VPD_mean = mean(VPD, na.rm = TRUE),
    VPD_sd = sd(VPD, na.rm = TRUE),
    
    .groups = "drop"
  )

head(avg_dat)

avg_dat_urban <- avg_dat[avg_dat$Habitat == "wall", ]
avg_dat_nonurban <- avg_dat[avg_dat$Habitat == "nonwall", ]

#Remove Zas from CEWL analyses. Were only able to sample one individual
avg_dat_no_zas <- avg_dat %>% 
  filter(Site != "Zas")

avg_dat_urban_no_zas <- avg_dat_no_zas[avg_dat_no_zas$Habitat == "wall", ]
avg_dat_nonurban_no_zas <- avg_dat_no_zas[avg_dat_no_zas$Habitat == "nonwall", ]


#Plotting morph averages
traits <- c("SVL", "Mass", "head_depth", "head_width", "head_length",
            "bicep_length", "forearm_length", "femur_length", "tibia_length", "LimbRatio")
avg_dat$Site <- factor(avg_dat$Site, levels = unique(avg_dat$Site[order(avg_dat$Habitat)]))

# Create the plots for each trait, with reordered Site variable
plot_list <- lapply(traits, function(trait) {
  ggplot(avg_dat, aes(x = Site, y = !!sym(paste0(trait, "_mean")), fill = Habitat)) +
    # Add error bars for SD
    geom_errorbar(aes(ymin = !!sym(paste0(trait, "_mean")) - !!sym(paste0(trait, "_sd")),
                      ymax = !!sym(paste0(trait, "_mean")) + !!sym(paste0(trait, "_sd"))),
                  width = 0.2, size = 1, color = "black") +
    # Add points for the means
    geom_point(shape = 21, size = 5, color = "black") +
    # Set custom color scale
    scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
    labs(title = paste(trait, "Mean ± SD"), x = "Site", y = paste(trait, "Mean")) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black", size = 1),
      panel.grid = element_blank(),
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})

# Combine the plots into a grid (e.g., 2 columns)
grid.arrange(grobs = plot_list, ncol = 2)

#Plotting CEWL average
ggplot(avg_dat_no_zas, aes(x = Habitat, y = CEWL_mean, fill = Habitat)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(y = "Average CEWL", x = "Habitat Type") +
  scale_fill_manual(values = c("nonwall" = "#DFBA8A", "wall" = "#1C85F6")) +
  theme(legend.position = "none")

#Model Selection
#Physiology
#Global Model of CEWL ~ Habitat
Avg_CEWL_Hab_Mod <- lm(CEWL_mean ~ Habitat + SVL_mean + msmt_temp_C_mean + msmt_VPD_kPa_mean, data = avg_dat_no_zas)
summary(Avg_CEWL_Hab_Mod)

# Model selection:
avg.cewl.dredge <- dredge(Avg_CEWL_Hab_Mod)
subset(avg.cewl.dredge, delta<4)
top.model <- get.models(avg.cewl.dredge, subset = 1)[[1]]
summary(top.model) 
Anova(top.model, type = 2)

# Model averaging:
mod1.avg <- model.avg(avg.cewl.dredge, subset = delta<5)
summary(mod1.avg, get.models(avg.cewl.dredge, subset = TRUE))
confint(mod1.avg)

# Sum of weights of each variable (similar to relative importance):
sw(mod1.avg)



#Global Model of CEWL ~ Environment
#Urban
Avg_CEWL_Urban_Mod <- lm(CEWL_mean ~ VPD_mean + ambient_temp_mean, data = avg_dat_urban_no_zas)
summary(Avg_CEWL_Urban_Mod)

# Model selection:
avg.cewl.dredge <- dredge(Avg_CEWL_Urban_Mod)
subset(avg.cewl.dredge, delta<4)
top.model <- get.models(avg.cewl.dredge, subset = 1)[[1]]
summary(top.model) 
Anova(top.model, type = 2)

# Model averaging:
mod1.avg <- model.avg(avg.cewl.dredge, subset = delta<5)
summary(mod1.avg, get.models(avg.cewl.dredge, subset = TRUE))
confint(mod1.avg)

# Sum of weights of each variable (similar to relative importance):
sw(mod1.avg)

#Nonurban
Avg_CEWL_Nonurban_Mod <- lm(CEWL_mean ~ VPD_mean, data = avg_dat_nonurban_no_zas)
summary(Avg_CEWL_Nonurban_Mod)

# Model selection:
avg.cewl.dredge <- dredge(Avg_CEWL_Nonurban_Mod)
subset(avg.cewl.dredge, delta<4)
top.model <- get.models(avg.cewl.dredge, subset = 1)[[1]]
summary(top.model) 
Anova(top.model, type = 2)

# Model averaging:
mod1.avg <- model.avg(avg.cewl.dredge, subset = delta<5)
summary(mod1.avg, get.models(avg.cewl.dredge, subset = TRUE))
confint(mod1.avg)

# Sum of weights of each variable (similar to relative importance):
sw(mod1.avg)


#CEWL vs VPD plot
ggplot(avg_dat_no_zas, aes(x = VPD_mean, y = CEWL_mean)) +
  geom_errorbar(aes(ymin = CEWL_mean - CEWL_sd, ymax = CEWL_mean + CEWL_sd), 
                width = 0.05, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    title = "CEWL vs VPD",
    x = "VPD (kPa)",
    y = expression(CEWL ~ (g/m^2*h)),
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )



#Morphology
#Density
Dens_mod <- lm(Density_mean ~ Habitat, data = avg_dat, na.action = "na.fail")
summary(Dens_mod)


#SVL
Avg.SVL.global <- lm(SVL_mean ~ Habitat + Density_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.SVL.global)

# Model selection:
options(na.action = "na.fail")
SVL.dredge <- dredge(Avg.SVL.global)
subset(SVL.dredge, delta<4)
SVLtop.model <- get.models(SVL.dredge, subset = 1)[[1]]
summary(SVLtop.model) 
Anova(SVLtop.model, type = 2)

# Model averaging:
SVLmod.avg <- model.avg(SVL.dredge, subset = delta<5)
summary(SVLmod.avg, get.models(Avg.SVL.global, subset = TRUE))
confint(SVLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(SVLmod.avg)

SVL_dens_lm <- lm(SVL_mean ~ Density_mean, data = avg_dat, na.action = "na.fail")
summary(SVL_dens_lm)

ggplot(avg_dat, aes(x = Density_mean, y = SVL_mean)) +
  geom_errorbar(aes(ymin = SVL_mean - SVL_sd, ymax = SVL_mean + SVL_sd), 
                width = 1, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black", alpha = 0.6) +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    title = "Body size increases with population density",
    x = "Denisty (N lizards/100m transect)",
    y = "SVL (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )



#Head Width
Avg.HW.global <- lm(head_width_mean ~ Habitat + Density_mean + SVL_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.HW.global)

# Model selection:
options(na.action = "na.fail")
HW.dredge <- dredge(Avg.HW.global)
subset(HW.dredge, delta<4)
HWtop.model <- get.models(HW.dredge, subset = 1)[[1]]
summary(HWtop.model) 
Anova(HWtop.model, type = 2)

# Model averaging:
HWmod.avg <- model.avg(HW.dredge, subset = delta<5)
summary(HWmod.avg, get.models(Avg.HW.global, subset = TRUE))
confint(HWmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(HWmod.avg)

ggplot(avg_dat, aes(x = SVL_mean, y = head_width_mean)) +
  geom_errorbar(aes(ymin = head_width_mean - head_width_sd, ymax = head_width_mean + head_width_sd), 
                width = 0.5, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    title = "Head Width vs SVL",
    x = "SVL (mm)",
    y = "Head Width (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )


#Head Depth
Avg.HD.global <- lm(head_depth_mean ~ Habitat + Density_mean + SVL_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.HD.global)

# Model selection:
options(na.action = "na.fail")
HD.dredge <- dredge(Avg.HD.global)
subset(HD.dredge, delta<4)
HDtop.model <- get.models(HD.dredge, subset = 1)[[1]]
summary(HDtop.model) 
Anova(HDtop.model, type = 2)

# Model averaging:
HDmod.avg <- model.avg(HD.dredge, subset = delta<5)
summary(HDmod.avg, get.models(Avg.HD.global, subset = TRUE))
confint(HDmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(HDmod.avg)

ggplot(avg_dat, aes(x = SVL_mean, y = head_depth_mean)) +
  geom_errorbar(aes(ymin = head_depth_mean - head_depth_sd, ymax = head_depth_mean + head_depth_sd), 
                width = 0.5, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    title = "Head depth vs SVL",
    x = "SVL (mm)",
    y = "Head depth (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )



#Head Length
Avg.HL.global <- lm(head_length_mean ~ Habitat + Density_mean + SVL_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.HL.global)

# Model selection:
options(na.action = "na.fail")
HL.dredge <- dredge(Avg.HL.global)
subset(HL.dredge, delta<4)
HLtop.model <- get.models(HL.dredge, subset = 1)[[1]]
summary(HLtop.model) 
Anova(HLtop.model, type = 2)

# Model averaging:
HLmod.avg <- model.avg(HL.dredge, subset = delta<5)
summary(HLmod.avg, get.models(Avg.HL.global, subset = TRUE))
confint(HLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(HLmod.avg)

ggplot(avg_dat, aes(x = SVL_mean, y = head_length_mean)) +
  geom_errorbar(aes(ymin = head_length_mean - head_length_sd, ymax = head_length_mean + head_length_sd), 
                width = 0.5, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    title = "Head length vs SVL",
    x = "SVL (mm)",
    y = "Head length (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )


#Femur Length
Avg.FL.global <- lm(femur_length_mean ~ Habitat + Density_mean + SVL_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.FL.global)

# Model selection:
options(na.action = "na.fail")
FL.dredge <- dredge(Avg.FL.global)
subset(FL.dredge, delta<4)
FLtop.model <- get.models(FL.dredge, subset = 1)[[1]]
summary(FLtop.model) 
Anova(FLtop.model, type = 2)

# Model averaging:
FLmod.avg <- model.avg(FL.dredge, subset = delta<5)
summary(FLmod.avg, get.models(Avg.FL.global, subset = TRUE))
confint(FLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(FLmod.avg)

ggplot(avg_dat, aes(x = SVL_mean, y = femur_length_mean)) +
  geom_errorbar(aes(ymin = femur_length_mean - femur_length_sd, ymax = femur_length_mean + femur_length_sd), 
                width = 0.5, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    title = "Femur length vs SVL",
    x = "SVL (mm)",
    y = "Femur length (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )



#Tibia Length
Avg.TL.global <- lm(tibia_length_mean ~ Habitat + Density_mean + SVL_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.TL.global)

# Model selection:
options(na.action = "na.fail")
TL.dredge <- dredge(Avg.TL.global)
subset(TL.dredge, delta<4)
TLtop.model <- get.models(TL.dredge, subset = 1)[[1]]
summary(TLtop.model) 
Anova(TLtop.model, type = 2)

# Model averaging:
TLmod.avg <- model.avg(TL.dredge, subset = delta<5)
summary(TLmod.avg, get.models(Avg.TL.global, subset = TRUE))
confint(TLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(TLmod.avg)

ggplot(avg_dat, aes(x = SVL_mean, y = tibia_length_mean)) +
  geom_errorbar(aes(ymin = tibia_length_mean - tibia_length_sd, ymax = tibia_length_mean + tibia_length_sd), 
                width = 0.5, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    title = "Tibia length vs SVL",
    x = "SVL (mm)",
    y = "Tibia length (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )


#Bicep Length
Avg.BL.global <- lm(bicep_length_mean ~ Habitat + Density_mean + SVL_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.BL.global)

# Model selection:
options(na.action = "na.fail")
BL.dredge <- dredge(Avg.BL.global)
subset(BL.dredge, delta<4)
BLtop.model <- get.models(BL.dredge, subset = 1)[[1]]
summary(BLtop.model) 
Anova(BLtop.model, type = 2)

# Model averaging:
BLmod.avg <- model.avg(BL.dredge, subset = delta<5)
summary(BLmod.avg, get.models(Avg.BL.global, subset = TRUE))
confint(BLmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(BLmod.avg)

ggplot(avg_dat, aes(x = SVL_mean, y = bicep_length_mean)) +
  geom_errorbar(aes(ymin = bicep_length_mean - bicep_length_sd, ymax = bicep_length_mean + bicep_length_sd), 
                width = 0.5, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    tiBLe = "Bicep length vs SVL",
    x = "SVL (mm)",
    y = "Bicep length (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )


#Forearm Length
Avg.FAL.global <- lm(forearm_length_mean ~ Habitat + Density_mean + SVL_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.FAL.global)

# Model selection:
options(na.action = "na.fail")
FAL.dredge <- dredge(Avg.FAL.global)
subset(FAL.dredge, delta<4)
FALtop.model <- get.models(FAL.dredge, subset = 1)[[1]]
summary(FALtop.model) 
Anova(FALtop.model, type = 2)

# Model averaging:
FALmod.avg <- model.avg(FAL.dredge, subset = delta<5)
summary(FALmod.avg, get.models(Avg.FAL.global, subset = TRUE))
confint(FALmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(FALmod.avg)

ggplot(avg_dat, aes(x = SVL_mean, y = forearm_length_mean)) +
  geom_errorbar(aes(ymin = forearm_length_mean - forearm_length_sd, ymax = forearm_length_mean + forearm_length_sd), 
                width = 0.5, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    tiFALe = "Forearm length vs SVL",
    x = "SVL (mm)",
    y = "Forearm length (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )


#Limb Ratio
Avg.LR.global <- lm(LimbRatio_mean ~ Habitat + Density_mean + SVL_mean, data = avg_dat, na.action = "na.fail")
summary(Avg.LR.global)

# Model selection:
options(na.action = "na.fail")
LR.dredge <- dredge(Avg.LR.global)
subset(LR.dredge, delta<4)
LRtop.model <- get.models(LR.dredge, subset = 1)[[1]]
summary(LRtop.model) 
Anova(LRtop.model, type = 2)

# Model averaging:
LRmod.avg <- model.avg(LR.dredge, subset = delta<5)
summary(LRmod.avg, get.models(Avg.LR.global, subset = TRUE))
confint(LRmod.avg)

# Sum of weights of each variable (similar to relative importance):
sw(LRmod.avg)

ggplot(avg_dat, aes(x = SVL_mean, y = LimbRatio_mean)) +
  geom_errorbar(aes(ymin = LimbRatio_mean - LimbRatio_sd, ymax = LimbRatio_mean + LimbRatio_sd), 
                width = 0.5, color = "black") +  # Plot first (bottom layer)
  geom_smooth(method = "lm", color = "black", size = 1.5, se = FALSE) +  # Second (behind points)
  geom_point(aes(fill = Habitat, shape = Habitat), 
             size = 3, stroke = 1, color = "black") +  # Top layer (points)
  #geom_text(aes(label = Site), hjust = -0.1, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  scale_shape_manual(values = c("wall" = 21, "nonwall" = 21)) +  # Same filled circle
  labs(
    tiLRe = "Limb Ratio vs SVL",
    x = "SVL (mm)",
    y = "Limb Ratio (mm)",
    fill = "Habitat",
    shape = "Habitat"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "right"
  )


###########################################
#Analyses for differences between traits
#emmeans
morph_vars <- c("SVL_mean", "Mass_mean", "head_length_mean", "head_width_mean", "head_depth_mean",
                "femur_length_mean", "tibia_length_mean", "forearm_length_mean", "bicep_length_mean", "LimbRatio_mean")

plot_list <- lapply(morph_vars, function(var) {
  model <- lm(as.formula(paste(var, "~ Habitat + Density_mean")), data = avg_dat)
  emm <- emmeans(model, ~ Habitat)
  emm_df <- as.data.frame(summary(emm))
  pairwise <- summary(pairs(emm))
  p_val <- signif(pairwise$p.value, 3)
  
  # Extract t-value for Habitat (assumes 2-level factor)
  model_summary <- summary(model)
  t_val <- signif(coef(model_summary)[grep("Habitat", rownames(coef(model_summary))), "t value"], 3)
  
  # Reorder habitat for consistency
  emm_df$Habitat <- factor(emm_df$Habitat, levels = c("wall", "nonwall"))
  
  ggplot(avg_dat %>% mutate(Habitat = fct_relevel(Habitat, "wall")),
         aes_string(x = "Habitat", y = var, fill = "Habitat")) +
    geom_jitter(aes(fill = Habitat), width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, shape = 21, color = "black") +
    geom_errorbar(data = emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL),
                  width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +
    geom_point(data = emm_df, aes(x = Habitat, y = emmean, fill = Habitat),
               size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +
    scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
    annotate("text", x = 1.5, y = max(avg_dat[[var]], na.rm = TRUE) * 1.05,
             label = paste0("t = ", t_val, ", p = ", p_val),
             size = 4, hjust = 0.5) +
    labs(
      title = var,
      x = "Habitat",
      y = var
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black", size = 1),
      panel.grid = element_blank(),
      text = element_text(size = 12),
      legend.position = "none"
    )
})

# Combine all plots into a grid
combined_plot <- wrap_plots(plot_list, ncol = 5)
print(combined_plot)


#CEWL
CEWL_model <- lm(CEWL_mean ~ Habitat, data = avg_dat_no_zas)
summary(CEWL_model)

# Get estimated marginal means and pairwise comparisons
CEWL_emm <- emmeans(CEWL_model, ~ Habitat)
CEWL_emm_df <- as.data.frame(summary(CEWL_emm))
CEWL_pairwise <- summary(pairs(CEWL_emm))
CEWL_pval <- signif(CEWL_pairwise$p.value, 3)

# Reorder Habitat levels
CEWL_emm_df$Habitat <- factor(CEWL_emm_df$Habitat, levels = c("wall", "nonwall"))

# Plot
CEWL_plot <- ggplot(avg_dat %>% mutate(Habitat = fct_relevel(Habitat, "wall")),
                    aes(x = Habitat, y = CEWL_mean, fill = Habitat)) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.3, stroke = 0.5, color = "black") +
  geom_errorbar(data = CEWL_emm_df, aes(x = Habitat, ymin = lower.CL, ymax = upper.CL),
                width = 0.2, size = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = CEWL_emm_df, aes(x = Habitat, y = emmean, fill = Habitat),
             size = 6, shape = 21, color = "black", alpha = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  labs(
    title = paste("CEWL_mean by Habitat (p =", CEWL_pval, ")"),
    x = "Habitat",
    y = "CEWL_mean"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.grid = element_blank(),
    text = element_text(size = 14),
    legend.position = "none"
  )

# Print the plot
print(CEWL_plot)




#########################################################################
#Site Microclimate
#Habitat Differences in ENV variables
env_vars <- c("ambient_temp_mean", "VPD_mean", "ambient_percent_rh_mean", "percent_veg_cover_mean")
env_labels <- c("Temperature (°C)", "VPD (kPa)", "Relative Humidity (%)", "Vegetation Cover (%)")

env_plots <- lapply(seq_along(env_vars), function(i) {
  var <- env_vars[i]
  label <- env_labels[i]
  
  model <- glm(as.formula(paste(var, "~ Habitat")), data = avg_dat)
  summary_model <- summary(model)
  p_val <- signif(summary_model$coefficients[2, 4], 3)
  
  ggplot(avg_dat %>% mutate(Habitat = fct_relevel(Habitat, "wall")), 
         aes_string(x = "Habitat", y = var, fill = "Habitat")) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.4, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 5, color = "black") +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, size = 1) +
    scale_fill_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
    labs(
      title = paste(label, "\n(p =", p_val, ")"),
      x = "Habitat",
      y = label
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      text = element_text(size = 12),
      legend.position = "none"
    )
})

# Combine into 2x2 layout
combined_env_plot <- (env_plots[[1]] | env_plots[[2]]) / (env_plots[[3]] | env_plots[[4]])
print(combined_env_plot)





#site Microclimate ~ macro regression
combined_site_avg_dat <- merge(avg_dat, location_dat, by = c("Site", "Habitat"), all.x = TRUE)

#Annual precip
#temp
Avg_Temp.Macro_V_Micro <- lm(ambient_temp_mean ~ annual_precip, data = combined_site_avg_dat)
summary(Avg_Temp.Macro_V_Micro)
plot(combined_site_avg_dat$annual_precip, combined_site_avg_dat$ambient_temp_mean, 
     col = values,
     pch = 19,
     xlab = "Annual Precipitation (mm)", 
     ylab = "Ambient Temperature (°C)", 
     main = "Ambient Temperature vs. Annual Precipitation")
abline(Avg_Temp.Macro_V_Micro, col = "red", lwd = 2)

#VPD
Avg_VPD.Macro_V_Micro <- lm(VPD_mean ~ annual_precip, data = combined_site_avg_dat)
summary(Avg_VPD.Macro_V_Micro)
plot(combined_site_avg_dat$annual_precip, combined_site_avg_dat$VPD_mean, 
     col = values,
     pch = 19,
     xlab = "Annual Precipitation (mm)", 
     ylab = "VPD", 
     main = "VPD vs. Annual Precipitation")
abline(Avg_VPD.Macro_V_Micro, col = "red", lwd = 2)

#veg
Avg_Veg.Macro_V_Micro <- lm(percent_veg_cover_mean ~ annual_precip, data = combined_site_avg_dat)
summary(Avg_Veg.Macro_V_Micro)
plot(combined_site_avg_dat$annual_precip, combined_site_avg_dat$percent_veg_cover_mean, 
     col = values,
     pch = 19,
     xlab = "Annual Precipitation (mm)", 
     ylab = "percent_veg_cover", 
     main = "percent_veg_cover vs. Annual Precipitation")
abline(Avg_Veg.Macro_V_Micro, col = "red", lwd = 2)

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}
# Extract data for each model
ambient_info <- extract_model_info(Avg_Temp.Macro_V_Micro)
vpd_info <- extract_model_info(Avg_VPD.Macro_V_Micro)
veg_info <- extract_model_info(Avg_Veg.Macro_V_Micro)

# Combine into a data frame
AVG.annual_precip_results_table <- data.frame(
  Response = c("Ambient Temperature (°C)", "VPD", "Percent Vegetation Cover (%)"),
  Estimate = c(ambient_info[1], vpd_info[1], veg_info[1]),
  Std_Error = c(ambient_info[2], vpd_info[2], veg_info[2]),
  t_value = c(ambient_info[3], vpd_info[3], veg_info[3]),
  p_value = c(ambient_info[4], vpd_info[4], veg_info[4]),
  Adj_R2 = c(ambient_info[5], vpd_info[5], veg_info[5])
)

print(AVG.annual_precip_results_table)

#solar radiation
#ambient_temp
AVG.srad_ambient_temp_models <- list()
for (i in 1:12) {
  srad_col <- paste0("srad_", i)
  formula <- as.formula(paste("ambient_temp_mean ~", srad_col))
  AVG.srad_ambient_temp_models[[srad_col]] <- lm(formula, data = combined_site_avg_dat)
}

summary(AVG.srad_ambient_temp_models[["srad_1"]])
summary(AVG.srad_ambient_temp_models[["srad_2"]])
summary(AVG.srad_ambient_temp_models[["srad_3"]])
summary(AVG.srad_ambient_temp_models[["srad_4"]])
summary(AVG.srad_ambient_temp_models[["srad_5"]])
summary(AVG.srad_ambient_temp_models[["srad_6"]])
summary(AVG.srad_ambient_temp_models[["srad_7"]])
summary(AVG.srad_ambient_temp_models[["srad_8"]])
summary(AVG.srad_ambient_temp_models[["srad_9"]])
summary(AVG.srad_ambient_temp_models[["srad_10"]])
summary(AVG.srad_ambient_temp_models[["srad_11"]])
summary(AVG.srad_ambient_temp_models[["srad_12"]])

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Loop through each model and extract summary stats
AVG.srad_temp_results <- data.frame()
for (month in 1:12) {
  srad_col <- paste0("srad_", month)
  model <- AVG.srad_ambient_temp_models[[srad_col]]
  info <- extract_model_info(model)
  AVG.srad_temp_results <- rbind(AVG.srad_temp_results,
                        data.frame(
                          Month = month,
                          Estimate = info[1],
                          Std_Error = info[2],
                          t_value = info[3],
                          p_value = info[4],
                          Adj_R2 = info[5]
                        ))
}

# View table
print(AVG.srad_temp_results)


#VPD
AVG.srad_VPD_models <- list()
for (i in 1:12) {
  srad_col <- paste0("srad_", i)
  formula <- as.formula(paste("VPD ~", srad_col))
  AVG.srad_VPD_models[[srad_col]] <- lm(formula, data = combined_dat)
}

# Extract model summaries
extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Summarize all VPD models
AVG.srad_VPD_results <- data.frame()
for (month in 1:12) {
  srad_col <- paste0("srad_", month)
  model <- AVG.srad_VPD_models[[srad_col]]
  info <- extract_model_info(model)
  AVG.srad_VPD_results <- rbind(AVG.srad_VPD_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.srad_VPD_results)

#Veg
AVG.srad_Veg_models <- list()
for (i in 1:12) {
  srad_col <- paste0("srad_", i)
  formula <- as.formula(paste("percent_veg_cover ~", srad_col))
  AVG.srad_Veg_models[[srad_col]] <- lm(formula, data = combined_dat)
}

# Summarize all veg cover models
AVG.srad_Veg_results <- data.frame()
for (month in 1:12) {
  srad_col <- paste0("srad_", month)
  model <- AVG.srad_Veg_models[[srad_col]]
  info <- extract_model_info(model)
  AVG.srad_Veg_results <- rbind(AVG.srad_Veg_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.srad_Veg_results)


#water vapor pressure
#ambient_temp
AVG.vapr_ambient_temp_models <- list()
for (i in 1:12) {
  vapr_col <- paste0("vapr_", i)
  formula <- as.formula(paste("ambient_temp_mean ~", vapr_col))
  AVG.vapr_ambient_temp_models[[vapr_col]] <- lm(formula, data = combined_site_avg_dat)
}

summary(AVG.vapr_ambient_temp_models[["vapr_1"]])
summary(AVG.vapr_ambient_temp_models[["vapr_2"]])
summary(AVG.vapr_ambient_temp_models[["vapr_3"]])
summary(AVG.vapr_ambient_temp_models[["vapr_4"]])
summary(AVG.vapr_ambient_temp_models[["vapr_5"]])
summary(AVG.vapr_ambient_temp_models[["vapr_6"]])
summary(AVG.vapr_ambient_temp_models[["vapr_7"]])
summary(AVG.vapr_ambient_temp_models[["vapr_8"]])
summary(AVG.vapr_ambient_temp_models[["vapr_9"]])
summary(AVG.vapr_ambient_temp_models[["vapr_10"]])
summary(AVG.vapr_ambient_temp_models[["vapr_11"]])
summary(AVG.vapr_ambient_temp_models[["vapr_12"]])

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Loop through each model and extract summary stats
AVG.vapr_temp_results <- data.frame()
for (month in 1:12) {
  vapr_col <- paste0("vapr_", month)
  model <- AVG.vapr_ambient_temp_models[[vapr_col]]
  info <- extract_model_info(model)
  AVG.vapr_temp_results <- rbind(AVG.vapr_temp_results,
                                 data.frame(
                                   Month = month,
                                   Estimate = info[1],
                                   Std_Error = info[2],
                                   t_value = info[3],
                                   p_value = info[4],
                                   Adj_R2 = info[5]
                                 ))
}

# View table
print(AVG.vapr_temp_results)


#VPD
AVG.vapr_VPD_models <- list()
for (i in 1:12) {
  vapr_col <- paste0("vapr_", i)
  formula <- as.formula(paste("VPD ~", vapr_col))
  AVG.vapr_VPD_models[[vapr_col]] <- lm(formula, data = combined_dat)
}

# Extract model summaries
extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Summarize all VPD models
AVG.vapr_VPD_results <- data.frame()
for (month in 1:12) {
  vapr_col <- paste0("vapr_", month)
  model <- AVG.vapr_VPD_models[[vapr_col]]
  info <- extract_model_info(model)
  AVG.vapr_VPD_results <- rbind(AVG.vapr_VPD_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.vapr_VPD_results)

#Veg
AVG.vapr_Veg_models <- list()
for (i in 1:12) {
  vapr_col <- paste0("vapr_", i)
  formula <- as.formula(paste("percent_veg_cover ~", vapr_col))
  AVG.vapr_Veg_models[[vapr_col]] <- lm(formula, data = combined_dat)
}

# Summarize all veg cover models
AVG.vapr_Veg_results <- data.frame()
for (month in 1:12) {
  vapr_col <- paste0("vapr_", month)
  model <- AVG.vapr_Veg_models[[vapr_col]]
  info <- extract_model_info(model)
  AVG.vapr_Veg_results <- rbind(AVG.vapr_Veg_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.vapr_Veg_results)


#monthly max temp
#ambient_temp
AVG.tmax_ambient_temp_models <- list()
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  formula <- as.formula(paste("ambient_temp_mean ~", tmax_col))
  AVG.tmax_ambient_temp_models[[tmax_col]] <- lm(formula, data = combined_site_avg_dat)
}

summary(AVG.tmax_ambient_temp_models[["tmax_1"]])
summary(AVG.tmax_ambient_temp_models[["tmax_2"]])
summary(AVG.tmax_ambient_temp_models[["tmax_3"]])
summary(AVG.tmax_ambient_temp_models[["tmax_4"]])
summary(AVG.tmax_ambient_temp_models[["tmax_5"]])
summary(AVG.tmax_ambient_temp_models[["tmax_6"]])
summary(AVG.tmax_ambient_temp_models[["tmax_7"]])
summary(AVG.tmax_ambient_temp_models[["tmax_8"]])
summary(AVG.tmax_ambient_temp_models[["tmax_9"]])
summary(AVG.tmax_ambient_temp_models[["tmax_10"]])
summary(AVG.tmax_ambient_temp_models[["tmax_11"]])
summary(AVG.tmax_ambient_temp_models[["tmax_12"]])

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Loop through each model and extract summary stats
AVG.tmax_temp_results <- data.frame()
for (month in 1:12) {
  tmax_col <- paste0("tmax_", month)
  model <- AVG.tmax_ambient_temp_models[[tmax_col]]
  info <- extract_model_info(model)
  AVG.tmax_temp_results <- rbind(AVG.tmax_temp_results,
                                 data.frame(
                                   Month = month,
                                   Estimate = info[1],
                                   Std_Error = info[2],
                                   t_value = info[3],
                                   p_value = info[4],
                                   Adj_R2 = info[5]
                                 ))
}

# View table
print(AVG.tmax_temp_results)


#VPD
AVG.tmax_VPD_models <- list()
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  formula <- as.formula(paste("VPD ~", tmax_col))
  AVG.tmax_VPD_models[[tmax_col]] <- lm(formula, data = combined_dat)
}

# Extract model summaries
extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Summarize all VPD models
AVG.tmax_VPD_results <- data.frame()
for (month in 1:12) {
  tmax_col <- paste0("tmax_", month)
  model <- AVG.tmax_VPD_models[[tmax_col]]
  info <- extract_model_info(model)
  AVG.tmax_VPD_results <- rbind(AVG.tmax_VPD_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.tmax_VPD_results)

#Veg
AVG.tmax_Veg_models <- list()
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  formula <- as.formula(paste("percent_veg_cover ~", tmax_col))
  AVG.tmax_Veg_models[[tmax_col]] <- lm(formula, data = combined_dat)
}

# Summarize all veg cover models
AVG.tmax_Veg_results <- data.frame()
for (month in 1:12) {
  tmax_col <- paste0("tmax_", month)
  model <- AVG.tmax_Veg_models[[tmax_col]]
  info <- extract_model_info(model)
  AVG.tmax_Veg_results <- rbind(AVG.tmax_Veg_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.tmax_Veg_results)


#Monthly min temp
#ambient_temp
AVG.tmin_ambient_temp_models <- list()
for (i in 1:12) {
  tmin_col <- paste0("tmin_", i)
  formula <- as.formula(paste("ambient_temp_mean ~", tmin_col))
  AVG.tmin_ambient_temp_models[[tmin_col]] <- lm(formula, data = combined_site_avg_dat)
}

summary(AVG.tmin_ambient_temp_models[["tmin_1"]])
summary(AVG.tmin_ambient_temp_models[["tmin_2"]])
summary(AVG.tmin_ambient_temp_models[["tmin_3"]])
summary(AVG.tmin_ambient_temp_models[["tmin_4"]])
summary(AVG.tmin_ambient_temp_models[["tmin_5"]])
summary(AVG.tmin_ambient_temp_models[["tmin_6"]])
summary(AVG.tmin_ambient_temp_models[["tmin_7"]])
summary(AVG.tmin_ambient_temp_models[["tmin_8"]])
summary(AVG.tmin_ambient_temp_models[["tmin_9"]])
summary(AVG.tmin_ambient_temp_models[["tmin_10"]])
summary(AVG.tmin_ambient_temp_models[["tmin_11"]])
summary(AVG.tmin_ambient_temp_models[["tmin_12"]])

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Loop through each model and extract summary stats
AVG.tmin_temp_results <- data.frame()
for (month in 1:12) {
  tmin_col <- paste0("tmin_", month)
  model <- AVG.tmin_ambient_temp_models[[tmin_col]]
  info <- extract_model_info(model)
  AVG.tmin_temp_results <- rbind(AVG.tmin_temp_results,
                                 data.frame(
                                   Month = month,
                                   Estimate = info[1],
                                   Std_Error = info[2],
                                   t_value = info[3],
                                   p_value = info[4],
                                   Adj_R2 = info[5]
                                 ))
}

# View table
print(AVG.tmin_temp_results)


#VPD
AVG.tmin_VPD_models <- list()
for (i in 1:12) {
  tmin_col <- paste0("tmin_", i)
  formula <- as.formula(paste("VPD ~", tmin_col))
  AVG.tmin_VPD_models[[tmin_col]] <- lm(formula, data = combined_dat)
}

# Extract model summaries
extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Summarize all VPD models
AVG.tmin_VPD_results <- data.frame()
for (month in 1:12) {
  tmin_col <- paste0("tmin_", month)
  model <- AVG.tmin_VPD_models[[tmin_col]]
  info <- extract_model_info(model)
  AVG.tmin_VPD_results <- rbind(AVG.tmin_VPD_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.tmin_VPD_results)

#Veg
AVG.tmin_Veg_models <- list()
for (i in 1:12) {
  tmin_col <- paste0("tmin_", i)
  formula <- as.formula(paste("percent_veg_cover ~", tmin_col))
  AVG.tmin_Veg_models[[tmin_col]] <- lm(formula, data = combined_dat)
}

# Summarize all veg cover models
AVG.tmin_Veg_results <- data.frame()
for (month in 1:12) {
  tmin_col <- paste0("tmin_", month)
  model <- AVG.tmin_Veg_models[[tmin_col]]
  info <- extract_model_info(model)
  AVG.tmin_Veg_results <- rbind(AVG.tmin_Veg_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.tmin_Veg_results)


#Monthly avg temp
#ambient_temp
AVG.tavg_ambient_temp_models <- list()
for (i in 1:12) {
  tavg_col <- paste0("tavg_", i)
  formula <- as.formula(paste("ambient_temp_mean ~", tavg_col))
  AVG.tavg_ambient_temp_models[[tavg_col]] <- lm(formula, data = combined_site_avg_dat)
}

summary(AVG.tavg_ambient_temp_models[["tavg_1"]])
summary(AVG.tavg_ambient_temp_models[["tavg_2"]])
summary(AVG.tavg_ambient_temp_models[["tavg_3"]])
summary(AVG.tavg_ambient_temp_models[["tavg_4"]])
summary(AVG.tavg_ambient_temp_models[["tavg_5"]])
summary(AVG.tavg_ambient_temp_models[["tavg_6"]])
summary(AVG.tavg_ambient_temp_models[["tavg_7"]])
summary(AVG.tavg_ambient_temp_models[["tavg_8"]])
summary(AVG.tavg_ambient_temp_models[["tavg_9"]])
summary(AVG.tavg_ambient_temp_models[["tavg_10"]])
summary(AVG.tavg_ambient_temp_models[["tavg_11"]])
summary(AVG.tavg_ambient_temp_models[["tavg_12"]])

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Loop through each model and extract summary stats
AVG.tavg_temp_results <- data.frame()
for (month in 1:12) {
  tavg_col <- paste0("tavg_", month)
  model <- AVG.tavg_ambient_temp_models[[tavg_col]]
  info <- extract_model_info(model)
  AVG.tavg_temp_results <- rbind(AVG.tavg_temp_results,
                                 data.frame(
                                   Month = month,
                                   Estimate = info[1],
                                   Std_Error = info[2],
                                   t_value = info[3],
                                   p_value = info[4],
                                   Adj_R2 = info[5]
                                 ))
}

# View table
print(AVG.tavg_temp_results)


#VPD
AVG.tavg_VPD_models <- list()
for (i in 1:12) {
  tavg_col <- paste0("tavg_", i)
  formula <- as.formula(paste("VPD ~", tavg_col))
  AVG.tavg_VPD_models[[tavg_col]] <- lm(formula, data = combined_dat)
}

# Extract model summaries
extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Summarize all VPD models
AVG.tavg_VPD_results <- data.frame()
for (month in 1:12) {
  tavg_col <- paste0("tavg_", month)
  model <- AVG.tavg_VPD_models[[tavg_col]]
  info <- extract_model_info(model)
  AVG.tavg_VPD_results <- rbind(AVG.tavg_VPD_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.tavg_VPD_results)

#Veg
AVG.tavg_Veg_models <- list()
for (i in 1:12) {
  tavg_col <- paste0("tavg_", i)
  formula <- as.formula(paste("percent_veg_cover ~", tavg_col))
  AVG.tavg_Veg_models[[tavg_col]] <- lm(formula, data = combined_dat)
}

# Summarize all veg cover models
AVG.tavg_Veg_results <- data.frame()
for (month in 1:12) {
  tavg_col <- paste0("tavg_", month)
  model <- AVG.tavg_Veg_models[[tavg_col]]
  info <- extract_model_info(model)
  AVG.tavg_Veg_results <- rbind(AVG.tavg_Veg_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.tavg_Veg_results)


#VPD
#ambient_temp
AVG.VPD_ambient_temp_models <- list()
for (i in 1:12) {
  VPD_col <- paste0("VPD_", i)
  formula <- as.formula(paste("ambient_temp_mean ~", VPD_col))
  AVG.VPD_ambient_temp_models[[VPD_col]] <- lm(formula, data = combined_site_avg_dat)
}

summary(AVG.VPD_ambient_temp_models[["VPD_1"]])
summary(AVG.VPD_ambient_temp_models[["VPD_2"]])
summary(AVG.VPD_ambient_temp_models[["VPD_3"]])
summary(AVG.VPD_ambient_temp_models[["VPD_4"]])
summary(AVG.VPD_ambient_temp_models[["VPD_5"]])
summary(AVG.VPD_ambient_temp_models[["VPD_6"]])
summary(AVG.VPD_ambient_temp_models[["VPD_7"]])
summary(AVG.VPD_ambient_temp_models[["VPD_8"]])
summary(AVG.VPD_ambient_temp_models[["VPD_9"]])
summary(AVG.VPD_ambient_temp_models[["VPD_10"]])
summary(AVG.VPD_ambient_temp_models[["VPD_11"]])
summary(AVG.VPD_ambient_temp_models[["VPD_12"]])

extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Loop through each model and extract summary stats
AVG.VPD_temp_results <- data.frame()
for (month in 1:12) {
  VPD_col <- paste0("VPD_", month)
  model <- AVG.VPD_ambient_temp_models[[VPD_col]]
  info <- extract_model_info(model)
  AVG.VPD_temp_results <- rbind(AVG.VPD_temp_results,
                                 data.frame(
                                   Month = month,
                                   Estimate = info[1],
                                   Std_Error = info[2],
                                   t_value = info[3],
                                   p_value = info[4],
                                   Adj_R2 = info[5]
                                 ))
}

# View table
print(AVG.VPD_temp_results)


#VPD
AVG.VPD_VPD_models <- list()
for (i in 1:12) {
  VPD_col <- paste0("VPD_", i)
  formula <- as.formula(paste("VPD ~", VPD_col))
  AVG.VPD_VPD_models[[VPD_col]] <- lm(formula, data = combined_dat)
}

# Extract model summaries
extract_model_info <- function(model) {
  sm <- summary(model)
  estimate <- round(sm$coefficients[2, 1], 2)
  se <- round(sm$coefficients[2, 2], 2)
  t_val <- round(sm$coefficients[2, 3], 2)
  p_val <- round(sm$coefficients[2, 4], 3)
  adj_r2 <- round(sm$adj.r.squared, 3)
  return(c(estimate, se, t_val, p_val, adj_r2))
}

# Summarize all VPD models
AVG.VPD_VPD_results <- data.frame()
for (month in 1:12) {
  VPD_col <- paste0("VPD_", month)
  model <- AVG.VPD_VPD_models[[VPD_col]]
  info <- extract_model_info(model)
  AVG.VPD_VPD_results <- rbind(AVG.VPD_VPD_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.VPD_VPD_results)

#Veg
AVG.VPD_Veg_models <- list()
for (i in 1:12) {
  VPD_col <- paste0("VPD_", i)
  formula <- as.formula(paste("percent_veg_cover ~", VPD_col))
  AVG.VPD_Veg_models[[VPD_col]] <- lm(formula, data = combined_dat)
}

# Summarize all veg cover models
AVG.VPD_Veg_results <- data.frame()
for (month in 1:12) {
  VPD_col <- paste0("VPD_", month)
  model <- AVG.VPD_Veg_models[[VPD_col]]
  info <- extract_model_info(model)
  AVG.VPD_Veg_results <- rbind(AVG.VPD_Veg_results,
                                data.frame(
                                  Month = month,
                                  Estimate = info[1],
                                  Std_Error = info[2],
                                  t_value = info[3],
                                  p_value = info[4],
                                  Adj_R2 = info[5]
                                ))
}

# View results
print(AVG.VPD_Veg_results)

#Combined results table
print(names(AVG.annual_precip_results_table))
print(names(AVG.srad_temp_results))
print(names(AVG.srad_VPD_results))
print(names(AVG.srad_Veg_results))
print(names(AVG.vapr_temp_results))
print(names(AVG.vapr_VPD_results))
print(names(AVG.vapr_Veg_results))
print(names(AVG.tmax_temp_results))
print(names(AVG.tmax_VPD_results))
print(names(AVG.tmax_Veg_results))
print(names(AVG.tmin_temp_results))
print(names(AVG.tmin_VPD_results))
print(names(AVG.tmin_Veg_results))
print(names(AVG.tavg_temp_results))
print(names(AVG.tavg_VPD_results))
print(names(AVG.tavg_Veg_results))
print(names(AVG.VPD_temp_results))
print(names(AVG.VPD_VPD_results))
print(names(AVG.VPD_Veg_results))

colnames(AVG.srad_temp_results)[1] <- "Response"
colnames(AVG.srad_VPD_results)[1] <- "Response"
colnames(AVG.srad_Veg_results)[1] <- "Response"
colnames(AVG.vapr_temp_results)[1] <- "Response"
colnames(AVG.vapr_VPD_results)[1] <- "Response"
colnames(AVG.vapr_Veg_results)[1] <- "Response"
colnames(AVG.tmax_temp_results)[1] <- "Response"
colnames(AVG.tmax_VPD_results)[1] <- "Response"
colnames(AVG.tmax_Veg_results)[1] <- "Response"
colnames(AVG.tmin_temp_results)[1] <- "Response"
colnames(AVG.tmin_VPD_results)[1] <- "Response"
colnames(AVG.tmin_Veg_results)[1] <- "Response"
colnames(AVG.tavg_temp_results)[1] <- "Response"
colnames(AVG.tavg_VPD_results)[1] <- "Response"
colnames(AVG.tavg_Veg_results)[1] <- "Response"
colnames(AVG.VPD_temp_results)[1] <- "Response"
colnames(AVG.VPD_VPD_results)[1] <- "Response"
colnames(AVG.VPD_Veg_results)[1] <- "Response"

AVG.annual_precip_results_table$model <- "Annual Precipitation"
AVG.srad_temp_results$model <- "Solar Radiation and Ambient Temperature"
AVG.srad_VPD_results$model <- "Solar Radiation and VPD"
AVG.srad_Veg_results$model <- "Solar Radiation and Vegetation"
AVG.vapr_temp_results$model <- "Vapor Pressure and Ambient Temperature"
AVG.vapr_VPD_results$model <- "Vapor Pressure and VPD"
AVG.vapr_Veg_results$model <- "Vapor Pressure and Vegetation"
AVG.tmax_temp_results$model <- "Max Temp and Ambient Temperature"
AVG.tmax_VPD_results$model <- "Max Temp and VPD"
AVG.tmax_Veg_results$model <- "Max Temp and Vegetation"
AVG.tmin_temp_results$model <- "Min Temp and Ambient Temperature"
AVG.tmin_VPD_results$model <- "Min Temp and VPD"
AVG.tmin_Veg_results$model <- "Min Temp and Vegetation"
AVG.tavg_temp_results$model <- "Average Temp and Ambient Temperature"
AVG.tavg_VPD_results$model <- "Average Temp and VPD"
AVG.tavg_Veg_results$model <- "Average Temp and Vegetation"
AVG.VPD_temp_results$model <- "VPD and Ambient Temperature"
AVG.VPD_VPD_results$model <- "VPD and VPD"
AVG.VPD_Veg_results$model <- "VPD and Vegetation"

combined_AVGmicroVmacro_summary <- rbind(AVG.annual_precip_results_table,
                                         AVG.srad_temp_results,
                                         AVG.srad_VPD_results,
                                         AVG.srad_Veg_results,
                                         AVG.vapr_temp_results,
                                         AVG.vapr_VPD_results,
                                         AVG.vapr_Veg_results,
                                         AVG.tmax_temp_results,
                                         AVG.tmax_VPD_results,
                                         AVG.tmax_Veg_results,
                                         AVG.tmin_temp_results,
                                         AVG.tmin_VPD_results,
                                         AVG.tmin_Veg_results,
                                         AVG.tavg_temp_results,
                                         AVG.tavg_VPD_results,
                                         AVG.tavg_Veg_results,
                                         AVG.VPD_temp_results,
                                         AVG.VPD_VPD_results,
                                         AVG.VPD_Veg_results)

#write.csv(combined_AVGmicroVmacro_summary, "combined_AVGmicroVmacro_summary.csv", row.names = FALSE)



#################################################
#Micro vs Macro Model Comparisons
#Remove Zas
combined_site_avg_dat_no_zas <- subset(combined_site_avg_dat, Site != "Zas")

#Microclimates
AVG.CEWL_ambient_temp_model <- lm(CEWL_mean ~ ambient_temp_mean, data = combined_site_avg_dat_no_zas)
summary(AVG.CEWL_ambient_temp_model)

AVG.CEWL_VPD_model <- lm(CEWL_mean ~ VPD_mean, data = combined_site_avg_dat_no_zas)
summary(AVG.CEWL_VPD_model)

AVG.CEWL_Veg_model <- lm(CEWL_mean ~ percent_veg_cover_mean, data = combined_site_avg_dat_no_zas)
summary(AVG.CEWL_Veg_model)

#Macroclimates
#Solar Radiation
AVG.srad_models_CEWL <- list()
# Loop through srad_1 to srad_12 and store each model
for (i in 1:12) {
  srad_col <- paste0("srad_", i)
  formula <- as.formula(paste("CEWL_mean ~", srad_col))  # Change ambient_temp to CEWL
  AVG.srad_models_CEWL[[srad_col]] <- lm(formula, data = combined_site_avg_dat_no_zas)
}

summary(AVG.srad_models_CEWL[["srad_1"]])
summary(AVG.srad_models_CEWL[["srad_2"]])
summary(AVG.srad_models_CEWL[["srad_3"]])
summary(AVG.srad_models_CEWL[["srad_4"]])
summary(AVG.srad_models_CEWL[["srad_5"]])
summary(AVG.srad_models_CEWL[["srad_6"]])
summary(AVG.srad_models_CEWL[["srad_7"]])
summary(AVG.srad_models_CEWL[["srad_8"]])
summary(AVG.srad_models_CEWL[["srad_9"]])
summary(AVG.srad_models_CEWL[["srad_10"]])
summary(AVG.srad_models_CEWL[["srad_11"]])
summary(AVG.srad_models_CEWL[["srad_12"]])


#Vapor Pressure
AVG.vapr_models_CEWL <- list()
# Loop through vapr_1 to vapr_12 and store each model
for (i in 1:12) {
  vapr_col <- paste0("vapr_", i)
  formula <- as.formula(paste("CEWL_mean ~", vapr_col))  # Change ambient_temp to CEWL
  AVG.vapr_models_CEWL[[vapr_col]] <- lm(formula, data = combined_site_avg_dat_no_zas)
}

summary(AVG.vapr_models_CEWL[["vapr_1"]])
summary(AVG.vapr_models_CEWL[["vapr_2"]])
summary(AVG.vapr_models_CEWL[["vapr_3"]])
summary(AVG.vapr_models_CEWL[["vapr_4"]])
summary(AVG.vapr_models_CEWL[["vapr_5"]])
summary(AVG.vapr_models_CEWL[["vapr_6"]])
summary(AVG.vapr_models_CEWL[["vapr_7"]])
summary(AVG.vapr_models_CEWL[["vapr_8"]])
summary(AVG.vapr_models_CEWL[["vapr_9"]])
summary(AVG.vapr_models_CEWL[["vapr_10"]])
summary(AVG.vapr_models_CEWL[["vapr_11"]])
summary(AVG.vapr_models_CEWL[["vapr_12"]])



#Maximum Temp
AVG.tmax_models_CEWL <- list()
# Loop through tmax_1 to tmax_12 and store each model
for (i in 1:12) {
  tmax_col <- paste0("tmax_", i)
  formula <- as.formula(paste("CEWL_mean ~", tmax_col))  # Change ambient_temp to CEWL
  AVG.tmax_models_CEWL[[tmax_col]] <- lm(formula, data = combined_site_avg_dat_no_zas)
}

summary(AVG.tmax_models_CEWL[["tmax_1"]])
summary(AVG.tmax_models_CEWL[["tmax_2"]])
summary(AVG.tmax_models_CEWL[["tmax_3"]])
summary(AVG.tmax_models_CEWL[["tmax_4"]])
summary(AVG.tmax_models_CEWL[["tmax_5"]])
summary(AVG.tmax_models_CEWL[["tmax_6"]])
summary(AVG.tmax_models_CEWL[["tmax_7"]])
summary(AVG.tmax_models_CEWL[["tmax_8"]])
summary(AVG.tmax_models_CEWL[["tmax_9"]])
summary(AVG.tmax_models_CEWL[["tmax_10"]])
summary(AVG.tmax_models_CEWL[["tmax_11"]])
summary(AVG.tmax_models_CEWL[["tmax_12"]])


#Minimum Temp
AVG.tmin_models_CEWL <- list()
# Loop through tmin_1 to tmin_12 and store each model
for (i in 1:12) {
  tmin_col <- paste0("tmin_", i)
  formula <- as.formula(paste("CEWL_mean ~", tmin_col))  # Change ambient_temp to CEWL
  AVG.tmin_models_CEWL[[tmin_col]] <- lm(formula, data = combined_site_avg_dat_no_zas)
}

summary(AVG.tmin_models_CEWL[["tmin_1"]])
summary(AVG.tmin_models_CEWL[["tmin_2"]])
summary(AVG.tmin_models_CEWL[["tmin_3"]])
summary(AVG.tmin_models_CEWL[["tmin_4"]])
summary(AVG.tmin_models_CEWL[["tmin_5"]])
summary(AVG.tmin_models_CEWL[["tmin_6"]])
summary(AVG.tmin_models_CEWL[["tmin_7"]])
summary(AVG.tmin_models_CEWL[["tmin_8"]])
summary(AVG.tmin_models_CEWL[["tmin_9"]])
summary(AVG.tmin_models_CEWL[["tmin_10"]])
summary(AVG.tmin_models_CEWL[["tmin_11"]])
summary(AVG.tmin_models_CEWL[["tmin_12"]])


#Average Temp
AVG.tavg_models_CEWL <- list()
# Loop through tavg_1 to tavg_12 and store each model
for (i in 1:12) {
  tavg_col <- paste0("tavg_", i)
  formula <- as.formula(paste("CEWL_mean ~", tavg_col))  # Change ambient_temp to CEWL
  AVG.tavg_models_CEWL[[tavg_col]] <- lm(formula, data = combined_site_avg_dat_no_zas)
}

summary(AVG.tavg_models_CEWL[["tavg_1"]])
summary(AVG.tavg_models_CEWL[["tavg_2"]])
summary(AVG.tavg_models_CEWL[["tavg_3"]])
summary(AVG.tavg_models_CEWL[["tavg_4"]])
summary(AVG.tavg_models_CEWL[["tavg_5"]])
summary(AVG.tavg_models_CEWL[["tavg_6"]])
summary(AVG.tavg_models_CEWL[["tavg_7"]])
summary(AVG.tavg_models_CEWL[["tavg_8"]])
summary(AVG.tavg_models_CEWL[["tavg_9"]])
summary(AVG.tavg_models_CEWL[["tavg_10"]])
summary(AVG.tavg_models_CEWL[["tavg_11"]])
summary(AVG.tavg_models_CEWL[["tavg_12"]])


#Macro VPD
AVG.VPD_models_CEWL <- list()
# Loop through VPD_1 to VPD_12 and store each model
for (i in 1:12) {
  VPD_col <- paste0("VPD_", i)
  formula <- as.formula(paste("CEWL_mean ~", VPD_col))  # Change ambient_temp to CEWL
  AVG.VPD_models_CEWL[[VPD_col]] <- lm(formula, data = combined_site_avg_dat_no_zas)
}

summary(AVG.VPD_models_CEWL[["VPD_1"]])
summary(AVG.VPD_models_CEWL[["VPD_2"]])
summary(AVG.VPD_models_CEWL[["VPD_3"]])
summary(AVG.VPD_models_CEWL[["VPD_4"]])
summary(AVG.VPD_models_CEWL[["VPD_5"]])
summary(AVG.VPD_models_CEWL[["VPD_6"]])
summary(AVG.VPD_models_CEWL[["VPD_7"]])
summary(AVG.VPD_models_CEWL[["VPD_8"]])
summary(AVG.VPD_models_CEWL[["VPD_9"]])
summary(AVG.VPD_models_CEWL[["VPD_10"]])
summary(AVG.VPD_models_CEWL[["VPD_11"]])
summary(AVG.VPD_models_CEWL[["VPD_12"]])

#Extract and order by R2
extract_r2 <- function(models_list, prefix) {
  data.frame(
    model = names(models_list),
    R2 = sapply(models_list, function(m) summary(m)$r.squared),
    type = prefix
  )
}

micro_models <- list(
  "ambient_temp" = AVG.CEWL_ambient_temp_model,
  "VPD" = AVG.CEWL_VPD_model,
  "percent_veg_cover" = AVG.CEWL_Veg_model
)
micro_r2 <- extract_r2(micro_models, "micro")

srad_r2 <- extract_r2(AVG.srad_models_CEWL, "srad")
vapr_r2 <- extract_r2(AVG.vapr_models_CEWL, "vapr")
tmax_r2 <- extract_r2(AVG.tmax_models_CEWL, "tmax")
tmin_r2 <- extract_r2(AVG.tmin_models_CEWL, "tmin")
tavg_r2 <- extract_r2(AVG.tavg_models_CEWL, "tavg")
macroVPD_r2 <- extract_r2(AVG.VPD_models_CEWL, "macroVPD")

all_r2 <- rbind(micro_r2, srad_r2, vapr_r2, tmax_r2, tmin_r2, tavg_r2, macroVPD_r2)
all_r2 <- all_r2[order(-all_r2$R2), ]


###########################################################
#Site Average Morphology PCA
#Unscaled
morpho_vars <- avg_dat %>%
  dplyr::select(ends_with("mean")) %>%
  dplyr::select(Mass_mean, SVL_mean, head_width_mean, head_depth_mean,
                head_length_mean, femur_length_mean, tibia_length_mean,
                bicep_length_mean, forearm_length_mean, LimbRatio_mean)

morpho_pca <- prcomp(morpho_vars, scale. = TRUE)
pca_scores <- as.data.frame(morpho_pca$x[, 1:2])  # Only keep PC1 and PC2
pca_scores$Site <- avg_dat$Site
pca_scores$Habitat <- avg_dat$Habitat

adonis_result <- adonis2(morpho_vars ~ avg_dat$Habitat, method = "euclidean")
adonis_result

# Plot PCA with coloring for Habitat (Urban vs Nonurban) and add p-value from Adonis test
avg_unscaled_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Habitat)) +
  geom_point(size = 4, alpha = 0.6) +
  stat_ellipse(aes(group = Habitat), linetype = "solid", linewidth = 1) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  labs(
    title = paste("PCA of Unscaled Morphological Traits (Adonis p = 0.669)"),
    x = "PC1", y = "PC2"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    text = element_text(size = 12),
    legend.position = "none"
  )

# Perform ANOVA for PC1 and PC2 separately to see if there's a significant difference between Habitat groups
pc1_aov <- summary(aov(PC1 ~ Habitat, data = pca_scores))
pc2_aov <- summary(aov(PC2 ~ Habitat, data = pca_scores))

# Print ANOVA results
print(pc1_aov)
print(pc2_aov)



#Unscaled
morpho_scaled <- morpho_vars %>%
  mutate(across(-SVL_mean, ~ .x / SVL_mean))
morpho_scaled_for_pca <- morpho_scaled %>% dplyr::select(-SVL_mean)

morpho_scaled_pca <- prcomp(morpho_scaled_for_pca, scale. = TRUE)
pca_scores <- as.data.frame(morpho_scaled_pca$x[, 1:2])  # Only keep PC1 and PC2
pca_scores$Site <- avg_dat$Site
pca_scores$Habitat <- avg_dat$Habitat

adonis_result2 <- adonis2(morpho_scaled_for_pca ~ avg_dat$Habitat, method = "euclidean")
adonis_result2

avg_scaled_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Habitat)) +
  geom_point(size = 4, alpha = 0.6) +
  stat_ellipse(aes(group = Habitat), linetype = "solid", linewidth = 1) +
  scale_color_manual(values = c("wall" = "#1C85F6", "nonwall" = "#DFBA8A")) +
  labs(
    title = paste("PCA of Morphological Traits Scaled to SVL (Adonis p = 0.962)"),
    x = "PC1", y = "PC2"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    text = element_text(size = 12),
    legend.position = "none"
  )

# Perform ANOVA for PC1 and PC2 separately to see if there's a significant difference between Habitat groups
pc1_aov <- summary(aov(PC1 ~ Habitat, data = pca_scores))
pc2_aov <- summary(aov(PC2 ~ Habitat, data = pca_scores))

# Print ANOVA results
print(pc1_aov)
print(pc2_aov)

combined_Avg_PCA_Plot <- (avg_scaled_pca | avg_unscaled_pca) + 
  plot_layout(nrow = 1)
































