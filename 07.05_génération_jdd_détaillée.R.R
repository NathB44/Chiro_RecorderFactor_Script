source("fun_for_bat_gen.R")

#### Definition des variables : ####

recorder <- c("blg", "bcd", "sm4", "adm") #Quatre détecteurs étudiés
sensi <- c("high", "low") #Sensibilité qui sera utile à prendre en compte pour la génération aléatoire
class <- c("A", "B", "C", "D") #Classes de végétation définies dans la même idée, on verra dessous
sites_par_class <- 4  #On a chaque classe représentée par 4 sites
n_sites <- length(class) * sites_par_class #donc n sites : n = nombre de catégories de classes * le nombre de sites qui les représentent

#### Paramètres pour chaque detecteurs ####

settings <- data.frame(
  recorder = recorder,
  mean = c(120, 24, 150, 80),    
  var = c(5, 2, 4, 4))

#### Dataframe avec l'ID pour chaque combinaison (detecteur + reglage) ####

recorder_settings <- expand.grid(
  recorder = recorder,
  sensi = sensi
)
recorder_settings$ID <- paste0(recorder_settings$recorder, "_", recorder_settings$sensi) #fusion des "rec" et "sensi"

# Structure qui servira à générer :

design <- expand.grid(
  site = 1:n_sites,
  ID = recorder_settings$ID
)
design <- merge(design, recorder_settings, by = "ID") #Jointure pour récupérer detecteur et reglage à partir de l’ID

#### Attribution des distances de vegetation ####

design$class <- rep(class, each = sites_par_class) #nouvelle colonne "classe" dans design

# On a donc : 4 fois un type de classe (a,b,c,d)  par combinaison de détecteurs/sensi = 2 (reglages) * 4 (recorders) * 4 (classes) = 32 combinaisons * 4 (réplicats) = 128 rows

#### Ajout d'une variable continue de distance selon les classes ####

## étape 1 définir les sites dans un nv tableau 

dist_for_sites <- data.frame(site = unique(design$site), class = rep(class, each = sites_par_class)) 
# Une colonne site qui aura le meme nombre de lignes qu'il y a de sites dans "design"
# et une colonne classe où on répète l'opération antérieure de définition des classes
# /!\ à noter que le nom "site" doit être le même que dans le design en prévision de la fusion future

## étape 2 définir les plages de distances pour chq classe :

dist_ranges <- list(A = 0:7, B = 8:14, C = 15:21, D = 22:29) 
# plages en mètres définies par le protocole

## étape 3 génération aléatoire de valeur correspondant aux plages

dist_for_sites$dist <- mapply(dst, dist_for_sites$class) 
# indication de quelle colonne doit etre étiquetée pour se voir attribuer une distance avec fonction dst cf:"fun_for_bat_gen.R"

## étape 4 fusion par site 

design <- merge(design, dist_for_sites, by = "site")
summary(design)

# Renommer la colonne classx et supprimer celle en double 
design <- design[,-6]
names(design)[names(design) == "class.x"] <- "class"

## étape 5 : appication d'une fonction d'ajustement de distance f cf: fun_for_bat_gen.R :

dist_values <- 0:29 #variable continue de 0 à 29 m

#simple visualisation de f  
plot(dist_values, f(dist_values), type = "l", col = "violetred4", lwd = 4,
     ylab = "Effet de la distance", xlab = "Distance (m)",
     main = "Effet de la distance par rapport à la végétation")

design$adjusted_values <- mapply(f, design$dist) #application de f sur les distances (dist=x)

#### Génération des valeurs pour les réglages "low" ####

design$contacts[design$sensi == "low"] <- mapply(summon, #application de la fonction summon, cf: fun_for_bat_gen.R
rec = design$recorder[design$sensi == "low"], #definition de sur quelle colonne et modalité chaque fonction s'applique
adj = design$adjusted_values[design$sensi == "low"] )

#on observe dans design que les valeurs lows sont maintenant remplies

#### Génération des valeurs pour les réglages "high" ####

design <- high_gen(design, recorder, sensi_effect) 
#application de la fonction high_gen cf:"fun_for_gen"

design$contacts <- round(design$contacts) 
#on arrondit au plus proche car on enregistre pas de demi-pipistrelles

### Le jeu de donnée est prêt à être utilisé ! super !!

#### Aperçu ####

head(design)
print(design$contacts)

plot(contacts~recorder, data=design, outline=F)
plot(contacts~sensi, data=design, outline=F)
plot(contacts~ID, data=design, main = "nombre de contacts selon les détecteurs et leurs sensibilités", xlab= "détecteurs et sensibilité (high/low)" , outline=F)

#### Enregistrement du design en csv ####

write.csv2(design, "design.csv", row.names = F)

## /!\ Si vous souhaitez d'autres valeurs dans votre csv, réexecutez le code ; ou utilisez set.seed()

