#### Mugwort Proj ----
# Code by N.R. Sommer
# Updated 16 July 2025

#### Packages ----
library(lme4)
library(parameters)
library(MuMIn)
library(tidyverse)
library(ggthemes)
library(survival)
library(survminer)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(patchwork)

#### Elemental content ----

cn_data <- read.csv("Data/mugwort_cn.csv")

cn_clean <- cn_data %>%
  mutate(
    Vegetation = factor(Vegetation, levels = c("Goldenrod", "Mugwort", "Goldenrod and Mugwort")),
    Treatment = factor(Treatment, levels = c("Control", "Herbivore", "Predator")),
    Block = factor(Block)
  )


cn_long <- cn_clean %>%
  filter(Vegetation %in% c("Goldenrod", "Mugwort")) %>%  # Only single species treatments
  select(Cage, Vegetation, Treatment, Block, 
         Gold.N, Gold.C, Gold.C.N, 
         Mug.N, Mug.C, Mug.C.N,
         Soil.N, Soil.C) %>%
  pivot_longer(
    cols = c(Gold.N, Gold.C, Gold.C.N, Mug.N, Mug.C, Mug.C.N),
    names_to = c("Species", "Variable"),
    names_pattern = "(Gold|Mug)\\.(N|C|C\\.N)",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  pivot_wider(
    names_from = Variable,
    values_from = Value,
    values_fn = mean  # Take mean of any duplicates
  ) %>%
  mutate(
    Species = factor(Species, levels = c("Gold", "Mug"), 
                    labels = c("Goldenrod", "Mugwort"))
  )

c_model <- lm(C ~ Species + Soil.C, data = cn_long)
summary(c_model)

n_model <- lm(N ~ Species + Soil.N, data = cn_long)
summary(n_model)




##### Figures ----

p_carbon <- ggplot(cn_long, aes(x = Species, y = C, fill = Species)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(aes(fill = Species), position = position_jitter(width = 0.2), alpha = 0.6, shape = 21) +
  labs(x = "Species", y = "Carbon Content (%)", fill = "Species") +
  scale_fill_manual(values = c("Goldenrod" = "#F4D03F", "Mugwort" = "#2E5A88")) +
  theme_few() +
  theme(legend.position = "none") +
  coord_flip()

p_nitrogen <- ggplot(cn_long, aes(x = Species, y = N, fill = Species)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(aes(fill = Species), position = position_jitter(width = 0.2), alpha = 0.6, shape = 21) +
  labs(x = "Species", y = "Nitrogen Content (%)", fill = "Species") +
  scale_fill_manual(values = c("Goldenrod" = "#F4D03F", "Mugwort" = "#2E5A88")) +
  theme_few() +
  theme(legend.position = "none") +
  coord_flip()

p_carbon
p_nitrogen

#ggsave("Figures/carbon_content.png", p_carbon, dpi = 300, width = 8, height = 3)
#ggsave("Figures/nitrogen_content.png", p_nitrogen, dpi = 300, width = 8, height = 3)


#### Plant structure ----

leafscan <- read.csv("Data/LeafScan.csv")
stemdata <- read.csv("Data/StemData.csv")

# Total leaf area / total leaf numbers = "per leaf area"
# (Per leaf area * total number of intersections) / stem length

stemdata$aggregate <- ((stemdata$Total_Leaf.Area/stemdata$Total_Leaf.Num)*stemdata$Total_Intersections)/stemdata$Stem_Length

stemdata_clean <- stemdata %>%
  filter(!is.na(aggregate), !is.na(Total_Leaf.Area), !is.na(Total_Intersections)) %>%
  mutate(Species = factor(Species))


# 1. Stem Complexity (aggregate metric) ----
complexity_model <- lm(aggregate ~ Species, data = stemdata_clean)
summary(complexity_model)
anova(complexity_model)

complexity_sim <- simulateResiduals(complexity_model, n = 1000)
plot(complexity_sim)
testResiduals(complexity_sim)

# 2. Total Leaf Area ----
leaf_area_model <- lm(Total_Leaf.Area ~ Species, data = stemdata_clean)
summary(leaf_area_model)
anova(leaf_area_model)

area_sim <- simulateResiduals(leaf_area_model, n = 1000)
plot(area_sim)
testResiduals(area_sim)

# 3. Total Intersections ----
intersections_model <- lm(Total_Intersections ~ Species, data = stemdata_clean)
summary(intersections_model)
anova(intersections_model)

intersections_sim <- simulateResiduals(intersections_model, n = 1000)
plot(intersections_sim)
testResiduals(intersections_sim)

# Emmeans
complexity_emm <- emmeans(complexity_model, ~ Species)
pairs(complexity_emm)

leaf_area_emm <- emmeans(leaf_area_model, ~ Species)
pairs(leaf_area_emm)

intersections_emm <- emmeans(intersections_model, ~ Species)
pairs(intersections_emm)

# 4. Figures ----
p1 <- stemdata_clean %>%
  ggplot() + 
  geom_boxplot(aes(y = aggregate, x = Species, fill = Species)) +
  geom_jitter(aes(y = aggregate, x = Species, fill = Species), width = 0.2, alpha = 0.6, shape = 21) +
  labs(x = "", y = "Complexity Index") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F")) +
  theme_few() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0, size = 16, face = "bold")) +
  ggtitle("A")

p2 <- stemdata_clean %>%
  ggplot() + 
  geom_boxplot(aes(y = Total_Leaf.Area, x = Species, fill = Species)) +
  geom_jitter(aes(y = Total_Leaf.Area, x = Species, fill = Species), width = 0.2, alpha = 0.6, shape = 21) +
  labs(x = "", y = "Total Leaf Area (cm²)") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F")) +
  theme_few() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0, size = 16, face = "bold")) +
  ggtitle("B")

p3 <- stemdata_clean %>%
  ggplot() + 
  geom_boxplot(aes(y = Total_Intersections, x = Species, fill = Species)) +
  geom_jitter(aes(y = Total_Intersections, x = Species, fill = Species), width = 0.2, alpha = 0.6, shape = 21) +
  labs(x = "", y = "Total Stem Intersections") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F")) +
  theme_few() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0, size = 16, face = "bold")) +
  ggtitle("C")


p1
p2
p3

# Create combined structural complexity figure
# Combine panels vertically with equal heights and proper spacing
combined_structural_figure <- (p1 / p2 / p3) + 
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.margin = margin(10, 10, 10, 10))

# Save the combined structural figure
ggsave("Figures/combined_structural_figure.png", combined_structural_figure, 
       dpi = 300, width = 6, height = 9, units = "in")

#ggsave("Figures/stem_complexity.png", p1, dpi = 300, width = 4, height = 6)
#ggsave("Figures/leaf_area.png", p2, dpi = 300, width = 4, height = 6)
#ggsave("Figures/stem_intersections.png", p3, dpi = 300, width = 4, height = 6)


#### Mass-Gain Lab Experiment ----

massgain <- read.csv("Data/MassGain.csv")

# Cleaning
massgain_clean <- massgain %>%
  mutate(
    Date = as.Date(Date, format = "%d-%b-%y"),
    Days_From_Start = as.numeric(Date - min(Date)),
    Measurement = factor(Measurement, levels = 1:6),
    Molted = grepl("Molt", Notes, ignore.case = TRUE),
    Status = case_when(
      grepl("Dead|Escaped", Notes, ignore.case = TRUE) ~ "Event",
      TRUE ~ "Censored"
    )
  )

# 1. Overall Growth rate  ----

# overall growth rate for each individual
growth_rates <- massgain_clean %>%
  group_by(ID, Type) %>%
  arrange(Days_From_Start) %>%
  summarise(
    Initial_Mass = first(Mass),
    Final_Mass = last(Mass),
    Days_Survived = max(Days_From_Start),
    Total_Growth = Final_Mass - Initial_Mass,
    Growth_Rate = Total_Growth / Days_Survived,
    Total_Molts = sum(Molted, na.rm = TRUE),
    Molting_Rate = Total_Molts / Days_Survived,
    .groups = 'drop'
  )

growth_rates %>% group_by(Type) %>% 
      summarise(mean_growth = mean(Growth_Rate, na.rm = TRUE),
                sd_growth = sd(Growth_Rate, na.rm = TRUE),
                n = n())

# Differences in growth rate between vegetation diets
growth_model <- lm(Growth_Rate ~ Type, data = growth_rates)
anova(growth_model)
summary(growth_model)

# 2. Nonlinear growth----
model_linear <- lmer(Mass ~ Type * Days_From_Start + (1|ID), 
                    data = massgain_clean)

model_quad <- lmer(Mass ~ Type * poly(Days_From_Start, 2) + (1|ID), 
                   data = massgain_clean)

model_cubic <- lmer(Mass ~ Type * poly(Days_From_Start, 3) + (1|ID), 
                    data = massgain_clean)

# Compare 
model_comparison <- anova(model_linear, model_quad, model_cubic)
model_comparison

best_model <- ifelse(which.min(c(AIC(model_linear), AIC(model_quad), AIC(model_cubic))) == 1,
                     "linear", 
                     ifelse(which.min(c(AIC(model_linear), AIC(model_quad), AIC(model_cubic))) == 2,
                            "quadratic", "cubic"))

final_model <- if(best_model == "linear") model_linear else if(best_model == "quadratic") model_quad else model_cubic

summary(final_model)


# diagnostics for non-linear model
sim_final <- simulateResiduals(final_model, n = 1000)
plot(sim_final)
testResiduals(sim_final)

# diagnostics for growth rate model
sim_growth <- simulateResiduals(growth_model, n = 1000)
plot(sim_growth)
testResiduals(sim_growth)

# Parametric assessment of vegetation type differences
anova(final_model)
emm <- emmeans(final_model, ~ Type | Days_From_Start, 
               at = list(Days_From_Start = c(0, 5, 11)))  # Start, middle, end
pairs(emm)



# Interaction between vegetation type and time
emm_interaction <- emmeans(final_model, ~ Type * Days_From_Start, 
                          at = list(Days_From_Start = c(0, 5, 11)))
pairs(emm_interaction, by = "Days_From_Start")

# 3. Molting ----

molting_model <- lm(Molting_Rate ~ Type, data = growth_rates)
anova(molting_model)
summary(molting_model)

# Diagnostics
sim_molting <- simulateResiduals(molting_model, n = 1000)
plot(sim_molting)
testResiduals(sim_molting)

molting_emm <- emmeans(molting_model, ~ Type)
pairs(molting_emm)

# 4. Survival ---

survival_data <- massgain_clean %>%
  filter(!grepl("Escaped", Notes, ignore.case = TRUE)) %>%
  group_by(ID, Type) %>%
  summarise(
    Time = max(Days_From_Start),
    Status = if(any(grepl("Dead", Notes, ignore.case = TRUE))) 1 else 0,
    .groups = 'drop'
  ) %>%
  filter(!is.na(Time))

surv_model <- survfit(Surv(Time, Status) ~ Type, data = survival_data)
summary(surv_model)

surv_diff <- survdiff(Surv(Time, Status) ~ Type, data = survival_data)
surv_diff

# Cox proportional hazards model
cox_model <- coxph(Surv(Time, Status) ~ Type, data = survival_data)
summary(cox_model)


# 4. Figures ----

# Growth trajectories
p1 <- ggplot(massgain_clean, aes(x = Days_From_Start, y = Mass, color = Type, group = ID)) +
  geom_line(alpha = 0.3, size = 0.5) +
  geom_point(alpha = 0.6, size = 1) +
  stat_summary(aes(group = Type), geom = "line", fun = mean, size = 2) +
  stat_summary(aes(group = Type), geom = "point", fun = mean, size = 3) +
  facet_wrap(~Type) +
  labs(x = "Days from Start", y = "Mass (g)") +
  scale_color_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F", "Grass" = "#72874EFF")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 14)) +
  ggtitle("A")

# Molting rate boxplot
p2 <- ggplot(growth_rates, aes(x = Type, y = Molting_Rate, fill = Type)) +
  geom_boxplot() +
  geom_point(aes(fill = Type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.6, shape = 21) +
  labs(x = "Vegetation Type", y = "Molting Rate (molts/day)") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F", "Grass" = "#72874EFF")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  ggtitle("B") +
  coord_flip()

# Survival curves
p3_survplot <- ggsurvplot(surv_model, data = survival_data,
                          pval = FALSE, conf.int = TRUE,
                          xlab = "Days from Start", ylab = "Survival Probability", 
                          palette = c("#2E5A88", "#F4D03F", "#453947FF"),
                          font.title = 18, font.x = 16, font.y = 16, font.tickslab = 14)

# Extract the plot and add panel label
p3 <- p3_survplot$plot + 
  theme(plot.title = element_text(hjust = 0, size = 18, face = "bold")) +
  ggtitle("C")


p1
p2
p3

# Create combined diet figure
# Combine panels vertically with equal widths and proper spacing
combined_diet_figure <- (p1 / p2 / p3) + 
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.margin = margin(10, 10, 10, 10))

# Save the combined diet figure
ggsave("Figures/combined_diet_figure.png", combined_diet_figure, 
       dpi = 300, width = 10, height = 12, units = "in")

ggsave("Figures/growth_trajectories.png", p1, dpi = 300, width = 8, height = 4)
ggsave("Figures/molting_rates.png", p2, dpi = 300, width = 8, height = 4)
ggsave("Figures/survival_curves.png", p3, dpi = 300, width = 8, height = 4)


#### Mesocosm Survival ----
surv <- read.csv("Data/Survival.csv")

# Clean
surv_clean <- surv %>%
  mutate(
    Date = as.Date(Julian_Date, origin = as.Date("2020-01-01")),
    Days_From_Start = as.numeric(Date - min(Date))
  )

# Interval survival rates (accounts for restocking)
interim_survival <- surv_clean %>%
  group_by(Cage_ID) %>%
  arrange(Days_From_Start) %>%
  mutate(
    Grasshoppers_End = Grasshopper_Count,
    Grasshoppers_Start = ifelse(Grasshoppers_End > 3, 
                               Grasshoppers_End,
                               ifelse(lag(Grasshopper_Count, default = 3) < 3, 
                                     3, 
                                     lag(Grasshopper_Count, default = 3))), 
    # interval survival rate (as end/start)
    Interim_Survival_Rate = ifelse(Grasshoppers_Start == 0, 
                                  ifelse(Grasshoppers_End == 0, 1, 0),
                                  Grasshoppers_End / Grasshoppers_Start)
  ) %>%
  filter(!is.na(lag(Grasshopper_Count))) %>%   # remove first observation (no previous data for interim calculation)
  ungroup()

# Transform for beta distribution
interim_survival <- interim_survival %>%
  mutate(
    Interim_Survival_Transformed = ifelse(Interim_Survival_Rate == 0, 0.001,
                                         ifelse(Interim_Survival_Rate == 1, 0.999, 
                                                Interim_Survival_Rate))
  )

# Beta GLMM for interim survival rates
interim_glmm_beta <- glmmTMB(Interim_Survival_Transformed ~ Veg_Type * Pred_Type + (1|Cage_ID), 
                            data = interim_survival, 
                            family = beta_family(link = "logit"))

# Model summary
summary(interim_glmm_beta)

print(interim_glmm_beta$sdr$pdHess) # convergence

interim_sim <- simulateResiduals(interim_glmm_beta, n = 1000)
plot(interim_sim)
testResiduals(interim_sim) # not great diagnostics; however there is a lack of significant effects which is unlikely to change interpretation

# Post-hoc comparisons for interim survival
interim_emm <- emmeans(interim_glmm_beta, ~ Veg_Type * Pred_Type)
pairs(interim_emm)

##### Figures ----

survival_trajectories <- surv_clean %>%
  group_by(Veg_Type, Pred_Type, Days_From_Start) %>%
  summarise(
    Mean_Count = mean(Grasshopper_Count, na.rm = TRUE),
    SE_Count = sd(Grasshopper_Count, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  filter(!is.na(Mean_Count))

p_survival_trajectories <- ggplot(survival_trajectories, aes(x = Days_From_Start, y = Mean_Count, color = Veg_Type, linetype = Pred_Type)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Mean_Count - SE_Count, ymax = Mean_Count + SE_Count, fill = Veg_Type), alpha = 0.2, color = NA) +
  facet_wrap(~Veg_Type) +
  labs(x = "Days from Start", y = "Mean Grasshopper Count", 
       color = "Vegetation Type", linetype = "Predator Treatment", fill = "Vegetation Type") +
  scale_color_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F", "Mixed" = "#453947FF")) +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F", "Mixed" = "#453947FF")) +
  theme_few() +
  theme(legend.position = "bottom") +
  coord_flip()


p_survival_trajectories

#ggsave("Figures/survival_trajectories.png", p_survival_trajectories, dpi = 300, width = 8, height = 6)



#### Biomass effects ----

initial_cover <- read.csv("Data/InitialCover.csv")
plant_harvest <- read.csv("Data/PlantHarvest.csv")

# Cleaning
final_biomass <- plant_harvest %>%
  group_by(Cage_ID, Plant_Type) %>%
  summarise(
    Final_Total_Biomass = sum(Mass, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = Plant_Type,
    values_from = Final_Total_Biomass,
    names_prefix = "Final_"
  ) %>%
  mutate(
    Final_Solidago = coalesce(Final_Solidago, 0),
    Final_Artemisia = coalesce(Final_Artemisia, 0),
    Final_Grass = coalesce(Final_Grass, 0),
    Total_Final_Biomass = Final_Solidago + Final_Artemisia + Final_Grass
  )

initial_cover_clean <- initial_cover %>%
  mutate(
    Artemisia_Initial_Cover = Artesima_PercCover,
    Solidago_Initial_Cover = Solidago_PercCover,
    Grass_Initial_Cover = Grass_PercCover,
    Total_Initial_Cover = Artesima_PercCover + Solidago_PercCover + Grass_PercCover,
    Veg_Type = factor(Veg_Type, levels = c("Artemisia", "Solidago", "Mixed")),
    Pred_Type = factor(Pred_Type, levels = c("Herbivore", "Predator"))
  )

biomass_analysis <- final_biomass %>%
  left_join(initial_cover_clean, by = "Cage_ID") %>%
  select(Cage_ID, Veg_Type, Pred_Type, Total_Final_Biomass, 
         Total_Initial_Cover, Artemisia_Initial_Cover, Solidago_Initial_Cover, Grass_Initial_Cover,
         Final_Solidago, Final_Artemisia, Final_Grass)

# Species-specific responses --

# Grass 
grass_model <- lm(Final_Grass ~ Pred_Type * Veg_Type + Grass_Initial_Cover, data = biomass_analysis)
summary(grass_model)

# Solidago
solidago_model <- lm(Final_Solidago ~ Pred_Type * Veg_Type + Solidago_Initial_Cover, data = biomass_analysis)
summary(solidago_model)

# Artemisia
artemisia_model <- lm(Final_Artemisia ~ Pred_Type * Veg_Type + Artemisia_Initial_Cover, data = biomass_analysis)
summary(artemisia_model)

# Total biomass (net effect)
total_model <- lm(Total_Final_Biomass ~ Pred_Type * Veg_Type + Total_Initial_Cover, data = biomass_analysis)
summary(total_model)

# Model diagnostics
grass_sim <- simulateResiduals(grass_model, n = 1000)
plot(grass_sim)
solidago_sim <- simulateResiduals(solidago_model, n = 1000)
plot(solidago_sim)
artemisia_sim <- simulateResiduals(artemisia_model, n = 1000)
plot(artemisia_sim)
total_sim <- simulateResiduals(total_model, n = 1000)
plot(total_sim)


# Emmeans
grass_emm <- emmeans(grass_model, ~ Pred_Type * Veg_Type)
pairs(grass_emm)
solidago_emm <- emmeans(solidago_model, ~ Pred_Type * Veg_Type)
pairs(solidago_emm)
artemisia_emm <- emmeans(artemisia_model, ~ Pred_Type * Veg_Type)
pairs(artemisia_emm)
total_emm <- emmeans(total_model, ~ Pred_Type * Veg_Type)
pairs(total_emm)


# Grasshopper feeding preferences (herbivore-only) -


# Artemisia vs Grass
artemisia_treatment_long <- biomass_analysis %>%
  filter(Veg_Type == "Artemisia") %>%
  select(Cage_ID, Pred_Type, Final_Artemisia, Final_Grass) %>%
  pivot_longer(cols = c(Final_Artemisia, Final_Grass),
               names_to = "Species", values_to = "Final_Biomass") %>%
  mutate(Species = gsub("Final_", "", Species))

artemisia_preference_model <- lm(Final_Biomass ~ Species * Pred_Type, data = artemisia_treatment_long)
summary(artemisia_preference_model)

# Solidago vs Grass
solidago_treatment_long <- biomass_analysis %>%
  filter(Veg_Type == "Solidago") %>%
  select(Cage_ID, Pred_Type, Final_Solidago, Final_Grass) %>%
  pivot_longer(cols = c(Final_Solidago, Final_Grass),
               names_to = "Species", values_to = "Final_Biomass") %>%
  mutate(Species = gsub("Final_", "", Species))

solidago_preference_model <- lm(Final_Biomass ~ Species * Pred_Type, data = solidago_treatment_long)
summary(solidago_preference_model)

# mixed treatment
mixed_biomass_long <- biomass_analysis %>%
  filter(Veg_Type == "Mixed") %>%
  select(Cage_ID, Pred_Type, Final_Artemisia, Final_Solidago, Final_Grass) %>%
  pivot_longer(cols = c(Final_Artemisia, Final_Solidago, Final_Grass),
               names_to = "Species", values_to = "Final_Biomass") %>%
  mutate(Species = gsub("Final_", "", Species))

mixed_biomass_model <- lm(Final_Biomass ~ Species * Pred_Type, data = mixed_biomass_long)
summary(mixed_biomass_model)


##### Figures ----

p_grass <- ggplot(biomass_analysis, aes(x = Pred_Type, y = Final_Grass, fill = Veg_Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(aes(fill = Veg_Type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.6, shape = 21) +
  labs(x = "", y = "Grass Biomass (g)", fill = "Plant Community") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F", "Mixed" = "#453947FF")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  ggtitle("A") +
  coord_flip()


p_solidago <- biomass_analysis %>%
  filter(Veg_Type %in% c("Solidago", "Mixed")) %>%
  ggplot(aes(x = Pred_Type, y = Final_Solidago, fill = Veg_Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(aes(fill = Veg_Type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.6, shape = 21) +
  labs(x = "", y = "Solidago Biomass (g)", fill = "Plant Community") +
  scale_fill_manual(values = c("Solidago" = "#F4D03F", "Mixed" = "#453947FF")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  ggtitle("B") +
  coord_flip()

p_artemisia <- biomass_analysis %>%
  filter(Veg_Type %in% c("Artemisia", "Mixed")) %>%
  ggplot(aes(x = Pred_Type, y = Final_Artemisia, fill = Veg_Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(aes(fill = Veg_Type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.6, shape = 21) +
  labs(x = "", y = "Artemisia Biomass (g)", fill = "Plant Community") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Mixed" = "#453947FF")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  ggtitle("C") +
  coord_flip()

p_total <- ggplot(biomass_analysis, aes(x = Pred_Type, y = Total_Final_Biomass, fill = Veg_Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(aes(fill = Veg_Type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.6, shape = 21) +
  labs(x = "", y = "Total Biomass (g)", fill = "Plant Community") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F", "Mixed" = "#453947FF")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  ggtitle("D") +
  coord_flip()

p_grass
p_solidago
p_artemisia
p_total

# Create combined biomass figure
# Add legends to the right of each individual plot
p_grass_with_legend <- p_grass + 
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Plant Community"))

p_solidago_with_legend <- p_solidago + 
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Plant Community"))

p_artemisia_with_legend <- p_artemisia + 
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Plant Community"))

p_total_with_legend <- p_total + 
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Plant Community"))

# Combine panels vertically with equal widths and proper spacing
combined_biomass_figure <- (p_grass_with_legend / p_solidago_with_legend / p_artemisia_with_legend / p_total_with_legend) + 
  plot_layout(heights = c(1, 1, 1, 1)) &
  theme(plot.margin = margin(10, 10, 10, 10))

# Save the combined biomass figure
# Wider figure to accommodate legends on the right side
ggsave("Figures/combined_biomass_figure.png", combined_biomass_figure, 
       dpi = 300, width = 10, height = 16, units = "in")

#ggsave("Figures/grass_biomass_response.png", p_grass, dpi = 300, width = 10, height = 4)
#ggsave("Figures/solidago_biomass_response.png", p_solidago, dpi = 300, width = 10, height = 4)
#ggsave("Figures/artemisia_biomass_response.png", p_artemisia, dpi = 300, width = 10, height = 4)
#ggsave("Figures/total_biomass_response.png", p_total, dpi = 300, width = 10, height = 4)

#### Habitat domains ----

habitat_data <- read.csv("Data/HabitatDomains.csv")


habitat_clean <- habitat_data %>%
  mutate(
    Veg_Treatment = factor(Veg_Treatment, levels = c("Artemisia", "Solidago", "Mixed")),
    Pred_Treatment = factor(Pred_Treatment, levels = c("Herbivore", "Predator")),
    Animal = factor(Animal, levels = c("MEFE", "PIMI")),
    Perch = factor(Perch, levels = c("Artemisia", "Solidago", "Grass", "Mesh", "Soil")),
    Beh = factor(Beh, levels = c("Sitting", "Foraging", "Moving", "Mating", "Molting")),
    Obs_ID = paste(Cage, Animal, Time, sep = "_")
  ) %>%
  filter(!is.na(Height), !is.na(Width))

# 1. Are there differences in grasshopper height between vegetation and predator treatments?

grasshopper_data <- habitat_clean %>% 
  filter(Animal %in% c("MEFE", "PIMI"))

height_model <- lmer(Height ~ Pred_Treatment * Veg_Treatment + (1|Cage), data = grasshopper_data)

summary(height_model)
anova(height_model)

# Diagnostics
height_sim <- simulateResiduals(height_model, n = 1000)
plot(height_sim)
testResiduals(height_sim)

# Emmeans: Do vegetation communities affect grasshopper responses to predators?
height_emm <- emmeans(height_model, ~ Pred_Treatment | Veg_Treatment)
pairs(height_emm)

# 2. Are there differences in the vegetation grasshoppers use between predator/non-predator treatments?

grasshopper_perch <- grasshopper_data %>%
  filter(Perch %in% c("Artemisia", "Solidago", "Grass")) %>%
  group_by(Cage, Pred_Treatment, Veg_Treatment) %>%
  summarise(
    total_obs = n(),
    artemisia_use = sum(Perch == "Artemisia", na.rm = TRUE),
    solidago_use = sum(Perch == "Solidago", na.rm = TRUE),
    grass_use = sum(Perch == "Grass", na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    artemisia_prop = artemisia_use / total_obs,
    solidago_prop = solidago_use / total_obs,
    grass_prop = grass_use / total_obs
  )

artemisia_use_model <- lm(artemisia_prop ~ Pred_Treatment * Veg_Treatment, data = grasshopper_perch)
solidago_use_model <- lm(solidago_prop ~ Pred_Treatment * Veg_Treatment, data = grasshopper_perch)
grass_use_model <- lm(grass_prop ~ Pred_Treatment * Veg_Treatment, data = grasshopper_perch)

summary(artemisia_use_model)
summary(solidago_use_model)
summary(grass_use_model)

# Diagnostics
artemisia_use_sim <- simulateResiduals(artemisia_use_model, n = 1000)
plot(artemisia_use_sim)
testResiduals(artemisia_use_sim)

solidago_use_sim <- simulateResiduals(solidago_use_model, n = 1000)
plot(solidago_use_sim)
testResiduals(solidago_use_sim)

grass_use_sim <- simulateResiduals(grass_use_model, n = 1000)
plot(grass_use_sim)
testResiduals(grass_use_sim)


# Emmeans: Do vegetation communities affect grasshopper habitat use responses to predators?
artemisia_emm <- emmeans(artemisia_use_model, ~ Pred_Treatment | Veg_Treatment)
pairs(artemisia_emm)

solidago_emm <- emmeans(solidago_use_model, ~ Pred_Treatment | Veg_Treatment)
pairs(solidago_emm)

grass_emm <- emmeans(grass_use_model, ~ Pred_Treatment | Veg_Treatment)
pairs(grass_emm)



# 3. When grasshoppers are foraging, what substrate are they foraging on?
# Is this different between predator and vegetation treatments?

grasshopper_foraging <- grasshopper_data %>%
  filter(Beh == "Foraging") %>%
  group_by(Cage, Pred_Treatment, Veg_Treatment) %>%
  summarise(
    foraging_obs = n(),
    artemisia_forage = sum(Perch == "Artemisia", na.rm = TRUE),
    solidago_forage = sum(Perch == "Solidago", na.rm = TRUE),
    grass_forage = sum(Perch == "Grass", na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    artemisia_forage_prop = artemisia_forage / foraging_obs,
    solidago_forage_prop = solidago_forage / foraging_obs,
    grass_forage_prop = grass_forage / foraging_obs
  )

artemisia_forage_model <- lm(artemisia_forage_prop ~ Pred_Treatment * Veg_Treatment, data = grasshopper_foraging)
solidago_forage_model <- lm(solidago_forage_prop ~ Pred_Treatment * Veg_Treatment, data = grasshopper_foraging)
grass_forage_model <- lm(grass_forage_prop ~ Pred_Treatment * Veg_Treatment, data = grasshopper_foraging)

summary(artemisia_forage_model)
summary(solidago_forage_model)
summary(grass_forage_model)

# Emmeans: Do vegetation communities affect grasshopper foraging substrate responses to predators?
artemisia_forage_emm <- emmeans(artemisia_forage_model, ~ Pred_Treatment | Veg_Treatment)
pairs(artemisia_forage_emm)

solidago_forage_emm <- emmeans(solidago_forage_model, ~ Pred_Treatment | Veg_Treatment)
pairs(solidago_forage_emm)

grass_forage_emm <- emmeans(grass_forage_model, ~ Pred_Treatment | Veg_Treatment)
pairs(grass_forage_emm)


##### Figures ----

p_height <- ggplot(grasshopper_data, aes(x = Pred_Treatment, y = Height, fill = Veg_Treatment)) +
  geom_violin(alpha = 0.7) +
  geom_point(aes(fill = Veg_Treatment), position = position_jitter(width = 0.2), alpha = 0.3, size = 0.5, shape = 21) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black", shape = 16) +
  facet_wrap(~Veg_Treatment) +
  labs(x = "", y = "Height (cm)", fill = "Vegetation Type") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F", "Mixed" = "#453947FF")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(margin = margin(r = 5)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  ggtitle("A")

perch_long <- grasshopper_perch %>%
  pivot_longer(cols = c(artemisia_prop, solidago_prop, grass_prop),
               names_to = "Plant_Type", values_to = "Proportion") %>%
  mutate(Plant_Type = gsub("_prop", "", Plant_Type),
         Plant_Type = factor(Plant_Type, levels = c("artemisia", "solidago", "grass")))

p_habitat_use <- ggplot(perch_long, aes(x = Pred_Treatment, y = Proportion, fill = Plant_Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(aes(fill = Plant_Type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.3, size = 0.5, shape = 21) +
  facet_wrap(~Veg_Treatment) +
  labs(x = "", y = "Proportion of Time", fill = "Plant Type") +
  scale_fill_manual(values = c("artemisia" = "#2E5A88", "solidago" = "#F4D03F", "grass" = "#72874EFF"),
                    labels = c("artemisia" = "Artemisia", "solidago" = "Solidago", "grass" = "Grass")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 14)) +
  ggtitle("B") +
  coord_flip()

# Use data from mixed treatments only
mixed_habitat_raw <- grasshopper_perch %>%
  filter(Veg_Treatment == "Mixed") %>%
  pivot_longer(cols = c(artemisia_prop, solidago_prop, grass_prop),
               names_to = "Plant_Type", values_to = "Proportion") %>%
  mutate(Plant_Type = gsub("_prop", "", Plant_Type),
         Plant_Type = factor(Plant_Type, levels = c("artemisia", "solidago", "grass")),
         Treatment = "Habitat Use")

mixed_forage_raw <- grasshopper_foraging %>%
  filter(Veg_Treatment == "Mixed") %>%
  pivot_longer(cols = c(artemisia_forage_prop, solidago_forage_prop, grass_forage_prop),
               names_to = "Plant_Type", values_to = "Proportion") %>%
  mutate(Plant_Type = gsub("_forage_prop", "", Plant_Type),
         Plant_Type = factor(Plant_Type, levels = c("artemisia", "solidago", "grass")),
         Treatment = "Foraging")

mixed_functional_data <- bind_rows(mixed_habitat_raw, mixed_forage_raw) %>%
  mutate(Plant = factor(Plant_Type, levels = c("artemisia", "solidago", "grass"),
                       labels = c("Artemisia", "Solidago", "Grass")))

p_functional_separation <- ggplot(mixed_functional_data, aes(x = Treatment, y = Proportion, fill = Plant)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge", alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.9), width = 0.25) +
  facet_wrap(~Pred_Treatment) +
  labs(x = "", y = "Proportion", fill = "Plant Type") +
  scale_fill_manual(values = c("Artemisia" = "#2E5A88", "Solidago" = "#F4D03F", "Grass" = "#72874EFF")) +
  theme_few() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 14)) +
  ggtitle("C") +
  coord_flip()

p_height
p_habitat_use
p_functional_separation

# Create combined paneled figure
# Add legends to the right of each individual plot
p_height_with_legend <- p_height + 
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Vegetation Type"))

p_habitat_use_with_legend <- p_habitat_use + 
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Plant Type"))

p_functional_separation_with_legend <- p_functional_separation + 
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Plant Type"))

# Combine panels vertically with equal widths and proper spacing
final_figure <- (p_height_with_legend / p_habitat_use_with_legend / p_functional_separation_with_legend) + 
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.margin = margin(10, 10, 10, 10))

# Display the final figure
final_figure

# Save the combined figure for publication
# Wider figure to accommodate legends on the right side
ggsave("Figures/combined_habitat_figure.png", final_figure, 
       dpi = 300, width = 10, height = 12, units = "in")

#ggsave("Figures/grasshopper_height.png", p_height, dpi = 300, width = 10, height = 6)
#ggsave("Figures/grasshopper_habitat_use.png", p_habitat_use, dpi = 300, width = 10, height = 6)
#ggsave("Figures/functional_separation.png", p_functional_separation, dpi = 300, width = 10, height = 6)
