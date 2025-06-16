library(glmmTMB)
library(emmeans)
library(ggplot2)
library(viridis)
library(DHARMa)
library(multcompView)
library(multcomp)

design <- read.csv2("design.csv")

#### GLMM ####

glmm_contacts <- glmmTMB(contacts~dist+ID+(1|site), family=nbinom2, data=(design))
summary(glmm_contacts)

#### Chart representation (ggplot) ####

mycol <- c("adm"= "#E69F00","bcd"="#56B4E9","blg"="#009E73","sm4"="#D55E00") 

windowsFonts(Avenir = windowsFont("Avenir Next LT Pro"))

phoc <- emmeans(glmm_contacts, pairwise ~ ID, adjust = "tukey") #very adapted in our situation
phoc
plot(phoc$emmeans) + labs(title = "Estimated marginal means (GLMM)",
                          x = "Log-scale estimated mean", y = "Recorder + Sensitivity (ID)") + theme_minimal(base_family = "Avenir")

####PLOT :

#associating labels to significant pairs (a=0.05) :
phoc_labels <- cld(phoc$emmeans, Letters = letters, alpha = 0.05, adjust = "tukey")
summary(phoc_labels)
phoc_labels$ID <- as.factor(phoc_labels$ID)
phoc_labels$.group <- gsub(" ", "", phoc_labels$.group)

#graphique :
p3 <- ggplot(design, aes(x = ID, y = contacts, fill = recorder)) + 
  geom_boxplot(outlier.shape = NA, color = "grey25") + 
  geom_jitter(width = 0.2, shape = 21, color = "white", alpha = 0.4, stroke = 0.4) +
  scale_fill_manual(values = mycol, 
                    labels = c("adm" = "Audiomoth", "bcd" = "Batcorder",
                               "blg" = "Batlogger", "sm4" = "SM4BAT")) +
  theme_minimal(base_family = "Avenir") + 
  labs(x = "Associations Détecteurs/Réglages",
       y = "Nombre de passage de chauves-souris par nuit", fill = "Détecteurs :") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t=15)),
    axis.title.y = element_text(vjust = 2),
    panel.grid.major.y = element_line(linewidth = 0.80)) +
  scale_x_discrete(labels = c(
    "adm_low" = "Faible",
    "adm_high" = "Fort",
    "bcd_low" = "Faible",
    "bcd_high" = "Fort",
    "blg_low" = "Faible",
    "blg_high" = "Fort",
    "sm4_low" = "Faible",
    "sm4_high" = "Fort"))

p4 <- p3 + geom_text(data = phoc_labels, inherit.aes = F,
                    aes(x = ID, y = max(design$contacts)-50, label = .group),
                    vjust = 0, size = 5)
p4

ggsave("siginificant_ID_groups_1.jpeg", plot = p4, width = 20, height = 20, units = "cm")

######Tableau des intercepts du GLMM######

# Extraire les coefficients du modèle conditionnel
estim <- summary(glmm_contacts)$coeff$cond

# Transformer en data frame avec uniquement les noms et estimates
tab_estim <- data.frame(
  param = rownames(estim),
  estim_glmm = estim[, "Estimate"])

# Mise en forme du tableau

tab_estim$param <- gsub("(Intercept)", "adm_high", tab_estim$param)
tab_estim$param <- gsub("[()]", "", tab_estim$param)
tab_estim$param <- gsub("^ID", "", tab_estim$param)

#Rétablir les bonnes valeurs des estimates en fonction de l'intercept
# Récupérer la première valeur

ref_value <- tab_estim$estim_glmm[1]

# Faire la soustraction de la valeur de l'intercept - les estimates des autres recorder

tab_estim$custom <- ref_value + tab_estim$estim_glmm
tab_estim$custom[1] <- ref_value #Rajouter la valeur de l'intercepet dans la première cellule

############Post hoc################ 

emm_ID <- emmeans(glmm_contacts, ~ ID) # Obtenir les moyennes marginales estimées pour le facteur ID
summary(emm_ID)

post_hoc <- pairs(emm_ID) # Effectuer des comparaisons par paires entre les niveaux de ID
summary(post_hoc)

#################Tableau_post_h##############

#Création d'un data frame du post hoc afin de pouvoir sélectionner les paires significativement différentes

tab_post_h<-as.data.frame(post_hoc)

# Séparation de la colonne "contrast" en deux colonnes. strsplit() pour découper chaque chaîne. do.call(rbind, ...) pour convertir la liste en matrice, que l'on ajoute ensuite comme deux nouvelles colonnes.
contrasts_split <- do.call(rbind, strsplit(as.character(tab_post_h$contrast), " - "))

# Ajouter les colonnes séparées au data.frame
tab_post_h$ID_A <- contrasts_split[,1]
tab_post_h$ID_B <- contrasts_split[,2]

tab_post_h <- tab_post_h[,-c(2,3,4,5)] #Suppression des colonnes non-utiles

#############Facteur de correction ###################

#Création d'une colonne contenant les facteurs pour chaque recorder sans celui du recorder présumé le plus faible (ici batcorder low) afin de s'y reférer
design$estimate_glmm<-tab_estim$custom[match(design$ID, tab_estim$param)]


# Extraire les lignes où ID == "bcd_low"
bcd_vals <- design[design$ID == "bcd_low", c("site", "estimate_glmm")]

# S'assurer que le vecteur est bien nommé par site
bcd_par_site <- setNames(bcd_vals$estimate_glmm, bcd_vals$site)

# Créer la nouvelle colonne en copiant la valeur de bcd_low par site
design$estimate_bcd_low <- bcd_par_site[as.character(design$site)]

#Division des estimates des recoders par celui du bcd_low
design$factors<-exp(design$estimate_glmm)/exp(design$estimate_bcd_low)

#Mise à l'exponentiel de l'estimate obtenu afin d'avoir la bonne valeur (glmm valeurs en log) et donc obtention du facteur de correction à appliquer au nombre de contatcs

#design$factors<-exp(design$estimate_log)

#Remettre l'estimeste de bcd_low à 1 afin de ne pas appliquer de facteur de correction au recorder de reférence
design$factors[design$ID == "bcd_low"] <- 1

#Application du facteur de correction

design$cts_with_f <-design$contacts/design$factors
design$cts_with_f <-round(design$cts_with_f)

#####GLMM de vérification de la correction du nombre de contact 

glmm_verif <- glmmTMB(cts_with_f~dist+ID+(1|site), family=nbinom2, data=(design))

summary(glmm_verif)

#Posthoc avec facterus de correction

f_post_hoc <- emmeans(glmm_verif, pairwise ~ ID, adjust= "tukey") # Obtenir les moyennes marginales estimées pour le facteur ID

summary(f_post_hoc)

###Représentation graphique (ggplot)

max_y_value <- max(design$cts_with_f, na.rm = TRUE)
p5 <- ggplot(design, aes(x = ID, y = cts_with_f, fill = recorder)) + 
  geom_boxplot(outlier.shape = NA, color = "grey25") + 
  geom_jitter(width = 0.2, shape = 21, color = "white", alpha = 0.4, stroke = 0.4) +
  scale_fill_manual(values = mycol, 
                    labels = c("adm" = "Audiomoth", "bcd" = "Batcorder",
                               "blg" = "Batlogger", "sm4" = "SM4BAT")) +
  theme_minimal(base_family = "Avenir") + 
  labs(x = "Associations Détecteurs/Réglages",
       y = "Nombre de passage de chauves-souris par nuit", fill = "Détecteurs :") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t=15)),
    axis.title.y = element_text(vjust = 2),
    panel.grid.major.y = element_line(linewidth = 0.80)) +
  scale_x_discrete(labels = c(
    "adm_low" = "Faible",
    "adm_high" = "Fort",
    "bcd_low" = "Faible",
    "bcd_high" = "Fort",
    "blg_low" = "Faible",
    "blg_high" = "Fort",
    "sm4_low" = "Faible",
    "sm4_high" = "Fort"))
p5

#associating labels to significant pairs (a=0.05) :
phoc_labels2 <- cld(f_post_hoc$emmeans, Letters = letters, alpha = 0.05, adjust = "tukey")
summary(phoc_labels2)
#clean up:
phoc_labels2$ID <- as.factor(phoc_labels2$ID)
phoc_labels2$.group <- gsub(" ", "", phoc_labels2$.group)

p6 <- p5 + geom_text(data = phoc_labels2, inherit.aes = F,
                     aes(x = ID, y = max(design$cts_with_f)-8, label = .group),
                     vjust = 0, size = 5)
p6

ggsave("non_siginificant_ID_groups_adjusted.jpeg", plot = p6, width = 20, height = 20, units = "cm")
#########Exporter les data frame en csv
dir.create("data_output")
write.csv2(design, "adj_design.csv", row.names = FALSE ) 

125/0.09032073
125/1.956552
