#Use set.seed() to generate the same data

source("fun_for_bats.R")

#### Variables : ####

recorder <- c("blg", "bcd", "sm4", "adm") # 4 recorders studied
sensi_tcs <- c("low4","low3","low2","low1","low0","medium","high0","high1","high2","high3","high4") # Sensibilities for random generation

#### Same recorder settings as in "design" ####

settings <- data.frame(recorder = recorder,
  mean = c(120, 24, 150, 80),    
  var = c(5, 2, 4, 4))

#### Dataframe with ID for each association (recorder + settings) ####
recorder_settings_tcs <- expand.grid(recorder = recorder, sensi = sensi_tcs)
recorder_settings_tcs$ID <- paste0(recorder_settings_tcs$recorder, "_", recorder_settings_tcs$sensi) #merger of rec and sensi

#### Sensibility values ####

sensi_effect_tcs <- data.frame(sensi = sensi_tcs,
  sensi_val = c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5) # adjust if needed
)

#### Generation structure ####

design_tcs <- merge(sensi_effect_tcs, recorder_settings_tcs, by = "sensi") 
design_tcs$sensi <- factor(design_tcs$sensi, levels = sensi_tcs, ordered = T)

#### Generation of contacts ####

design_tcs$contacts <- mapply(summon_tcs, rec = design_tcs$recorder, sensi_val = design_tcs$sensi_val)
# Application of "summon_tcs" fun, cf: fun_for_bats.R

#### Getting relative to sm4 ####

design_tcs$relative_contacts <- design_tcs$contacts/design_tcs$contacts[design_tcs$ID=="sm4_high4"] 

### DATASET READY TO USE
write.csv2(design_tcs, "design_tcs.csv", row.names = F)
