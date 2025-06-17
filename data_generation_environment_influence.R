source("fun_for_bats.R")

#### Description of variables : ####

recorder <- c("blg", "bcd", "sm4", "adm") # 4 recorders of the study
sensi <- c("high", "low") # 2 senbilities useful for random generation
class <- c("A", "B", "C", "D") # 4 classes for random generation (see later)
sites_par_class <- 4  # each class has 4 sites
n_sites <- length(class) * sites_par_class # 4*4 = 16 sites (n)

#### Recorders' settings ####

settings <- data.frame( # settings useful for random generation
  recorder = recorder,
  mean = c(120, 24, 150, 80),    
  var = c(5, 2, 4, 4))

#### Dataframe with l'ID for each association recorder+setting ####

recorder_settings <- expand.grid(recorder = recorder, sensi = sensi)
recorder_settings$ID <- paste0(recorder_settings$recorder, "_", recorder_settings$sensi) #ID = associations

#### Generation structure creation ####

design <- expand.grid(site = 1:n_sites, ID = recorder_settings$ID)
design <- merge(design, recorder_settings, by = "ID") # settings recovery

#### Giving out vegetation distances ####

design$class <- rep(class, each = sites_par_class)
# we now have : 4 times a class type (a,b,c,d) per association of recorder/sensibility =2 (sensibilites) * 4 (recorders) * 4 (classes) = 32 associations * 4 (replication) = 128 rows

#### Adding continuous variable of distance according to vegetation classes #### 

## STEP 1 definition of sites in a new dataframe 

dist_for_sites <- data.frame(site = unique(design$site), class = rep(class, each = sites_par_class)) 

## STEP 2 ditance ranges for each class

dist_ranges <- list(A = 0:7, B = 8:14, C = 15:21, D = 22:29) 
# ranges in meters

## STEP 3 summoning random values to suit ranges

dist_for_sites$dist <- mapply(dst, dist_for_sites$class) #applyinf "dst" function in fun.for.bats.R

## STEP 4 merger per sites

design <- merge(design, dist_for_sites, by = "site")
summary(design)

#Cleaning 
design <- design[,-6]
names(design)[names(design) == "class.x"] <- "class"

## STEP 5 : application of fun : distance adjusment "f" cf: fun_for_bats.R : 

dist_values <- 0:29 # 0 to 29 meters = continuous variable

# ploting distance effect
plot(dist_values, f(dist_values), type = "l", col = "violetred4", lwd = 4,
     ylab = "Distance effect", xlab = "Distance (m)",
     main = "Distance effect in relation to vegetation")

design$adjusted_values <- mapply(f, design$dist) #application of f on distances (dist=x)

#### Generation of values for "low" settings ####

design$contacts[design$sensi == "low"] <- mapply(summon, #application of "summon" fun, cf: fun_for_bats.R
rec = design$recorder[design$sensi == "low"], 
adj = design$adjusted_values[design$sensi == "low"] )

#low values in design$contacts are now filled

#### Generation of values for "high" settings ####

design <- high_gen(design, recorder, sensi_effect) 
#application of "high_gen" fun, cf: fun_for_bats.R

design$contacts <- round(design$contacts) 
#round to closest (rbinom generate integer numbers, but here, these integer have been multiplied)

### The data set is now ready to use ! great !!

#### Overview ####

head(design)
print(design$contacts)

plot(contacts~recorder, data=design, outline=F)
plot(contacts~sensi, data=design, outline=F)
plot(contacts~ID, data=design, main = "Number of contacts according to recorders and their sensibilities", xlab= "Recorder and Sensibilities (high/low)" , outline=F)

#### design into csv ####

write.csv2(design, "design.csv", row.names = F)

## /!\ If you want other values in your csv, execute the scipt agin ; or use set.seed()

