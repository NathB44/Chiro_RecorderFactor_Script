setwd("C:/Users/nbesn/Desktop/BOULOBOULO/Stages/Master1/CESCO/Data et analyses/génération jdd")
source("fun_for_bat_gen.R")

library(ggplot2)
library(tidygam)
library(mgcv)
library(dplyr)
library(MASS)
library(DHARMa)
library(glmmTMB)

windowsFonts(Avenir = windowsFont("Avenir Next LT Pro"))
design_tcs <- read.csv2("design_tcs.csv")
mycol <- c("adm"= "#E69F00","bcd"="#56B4E9","blg"="#009E73","sm4"="#D55E00")

design_tcs$sensi <- as.factor(design_tcs$sensi)
design_tcs$recorder <- as.factor(design_tcs$recorder)
design_tcs$ID <- as.factor(design_tcs$ID)
summary(design_tcs)

#######################################
design_tcs <- design_tcs %>%
  group_by(recorder) %>%
  mutate(relative_contacts=(contacts-min(contacts))/(max(contacts)-min(contacts)))

gam_model <- gam(relative_contacts~s(sensi_val, by=recorder, k=11), data=design_tcs, fx=T)

model_p <- predict_gam(gam_model)
model_p

model_p <- model_p %>%
  group_by(recorder) %>%
  mutate(sensi_scaled=(sensi_val-min(sensi_val))/(max(sensi_val)-min(sensi_val))) %>%
  ungroup()

model_p <- model_p %>%
  rename(fit = relative_contacts,
    lower=lower_ci,
    upper=upper_ci)

p_gam <- ggplot(model_p, aes(x = sensi_val, y = fit, color = recorder)) +
  geom_line(linewidth = 1.25) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = recorder), alpha = 0.15, color = NA) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  labs(x = "Sensibilité",
    y = "Contacts relatifs",
    color = "Détecteurs",
    fill = "Détecteurs") +
  theme_minimal(base_family = "Avenir") +
  theme(axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(vjust = 2),
    panel.grid.major.y = element_line(linewidth = 0.8))
p_gam

ggsave("courbes_sensi_gam.jpeg", plot = p_gam, width = 20, height = 20, units = "cm")

###########EXRACTION DES EQUATIONS : gam sort des droites donc équivalent à des glm :

glm_addit_poiss <- glm(relative_contacts~recorder+sensi_val, data = design_tcs, family = poisson)
glm_inter_poiss <- glm(relative_contacts~recorder*sensi_val, data = design_tcs, family = poisson)
glm_addit_nb <- glm.nb(relative_contacts~recorder+sensi_val, data = design_tcs)
glm_inter_nb <- glm.nb(relative_contacts~recorder*sensi_val, data = design_tcs)

AIC(glm_addit_nb, glm_inter_nb)
AIC(glm_addit_poiss, glm_inter_poiss)

#le modèle additif en negative binomial a la meilleure valeur d'AIC

summary(glm_addit_nb)
simulateResiduals(glm_addit_nb, plot = T)

#on recupere les coeff :

coef(summary(glm_addit_nb))
coefs <- coef(glm_addit_nb)
recorder <- c("adm", "bcd", "blg", "sm4")
sapply(recorder, get_eq) #fonction permettant de donner les équations automatiquement

#################### graphique : ###########################

a <- coefs["sensi_val"]
# tableau avec intercepts par détecteur
equations <- data.frame(recorder=recorder, intercept=sapply(recorder, eq.df), slope=a)

# grille de valeurs de sensis
sensi_vals <- seq(0, 1.5, length.out = 100)

# construction des courbes
curve_data <- equations %>%
  group_by(recorder) %>%
  do({data.frame(sensi_val = sensi_vals,
      contacts = exp(.$slope * sensi_vals + .$intercept))})

p_eq <- ggplot(curve_data, aes(x = sensi_val, y = contacts, color = recorder)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = mycol) +
  labs(title = "Courbes de sensibilité", x = "Sensibilité", y = "Nombre de contacts relatifs", color = "Détecteurs :") +
  theme_minimal(base_family = "Avenir") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))
p_eq

ggsave("courbes_sensi_relativ_simple.jpeg", plot = p_eq, width = 20, height = 20, units = "cm")

#######
########################   AJUSTEMENT DES COURBES
#######

design <- read.csv2("design.csv")
design$recorder <- as.factor(design$recorder)

#### Passage en relatif au sm4 :
ref_value <- mean(design$contacts[design$ID == "sm4_high"], na.rm = TRUE)
design$relative_contacts <- design$contacts / ref_value

obs_means <- design %>%
  filter(grepl("low", ID)) %>%
  group_by(recorder) %>%
  summarise(mean_obs = mean(relative_contacts, na.rm = TRUE)) %>%
  ungroup()

# Intercept du modèle glm_addit_nb pour chaque détecteur
model_means <- equations %>%
  mutate(mean_model = exp(intercept + slope * mean(sensi_vals))) %>%
  dplyr::select(recorder, mean_model)

# Fusion avec les moyennes observées
adjustments <- left_join(model_means, obs_means, by = "recorder") %>%
  mutate(adj_factor = mean_obs / mean_model)

# Mise à jour des équations avec les facteurs de correction
equations_adj <- left_join(equations, adjustments, by = "recorder")

# Reconstruction des courbes ajustées
curve_data_adj <- equations_adj %>%
  group_by(recorder) %>%
  do({
    data.frame(
      sensi_val = sensi_vals,
      contacts = exp(.$slope * sensi_vals + .$intercept) * .$adj_factor
    )
  }) %>%
  ungroup()

p_eq_adj <- ggplot(curve_data_adj, aes(x = sensi_val, y = contacts, color = recorder)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = mycol) +
  labs(title = "Courbes de sensibilité ajustées", x = "Sensibilité", y = "Nombre de contacts relatifs", color = "Détecteurs :") +
  theme_minimal(base_family = "Avenir") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))
p_eq_adj

# Sauvegarde
ggsave("courbes_sensi_relativ_ajustees.jpeg", plot = p_eq_adj, width = 20, height = 20, units = "cm")


