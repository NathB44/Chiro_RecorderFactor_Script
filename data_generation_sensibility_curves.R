set.seed(0) #set a seed to generate the same data

source("fun_for_bats.R")

#### Definition des variables : ####

recorder <- c("blg", "bcd", "sm4", "adm") #Quatre détecteurs étudiés
sensi_tcs <- c("low4","low3","low2","low1","low0","medium","high0","high1","high2","high3","high4") #Sensibilité qui sera utile à prendre en compte pour la génération aléatoire

#### Paramètres pour chaque detecteurs ####

settings <- data.frame(
  recorder = recorder,
  mean = c(120, 24, 150, 80),    
  var = c(5, 2, 4, 4)
)

#### Dataframe avec l'ID pour chaque combinaison (detecteur + reglage) ####

nrep <- 1
recorder_settings_tcs <- expand.grid(
  recorder = recorder,
  sensi = sensi_tcs,
  rep = 1:nrep
)
recorder_settings_tcs$ID <- paste0(recorder_settings_tcs$recorder, "_", recorder_settings_tcs$sensi) #fusion des "rec" et "sensi"

#### Définition des valeurs de sensi

sensi_effect_tcs <- data.frame(
  sensi = sensi_tcs,
  sensi_val = c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5) # À ajuster si besoin
)

#### Structure qui servira à générer :

design_tcs <- merge(sensi_effect_tcs, recorder_settings_tcs, by = "sensi") 
design_tcs$sensi <- factor(design_tcs$sensi, levels = sensi_tcs, ordered = T)

#### Génération des valeurs 

design_tcs$contacts <- mapply(summon_tcs, rec = design_tcs$recorder, sensi_val = design_tcs$sensi_val)

#### On met les contacts en relatif : 

design_tcs$relative_contacts <- design_tcs$contacts/design_tcs$contacts[design_tcs$ID=="sm4_high4"] 

### Le jeu de donnée est prêt à être utilisé !

write.csv2(design_tcs, "design_tcs.csv", row.names = F)

#########################################
