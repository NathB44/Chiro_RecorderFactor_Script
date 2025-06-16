#### TAG AND RANDOM DISTANCE FUNCTION ####

dst <- function(dst) { sample(dist_ranges[[dst]], 1) } 

# generation of random values according to the ranges with sample() for each label "dst"

####ADJUSTED DISTANCE TO VEGETATION####

f <- function(x) {1/(0.4*log(0.5*x+1.5)+0.8)}

#### ADJUSTED SENSITIVITY ####

sensi_effect <- function(sens) { if (sens == "high") return(1.5) else return(0.6) } 

# association of an arbitrary value to recorder's sensitivity high/low

#### CONTACTS GENERATION FOR LOW SENSORS ####

summon <- function(rec, adj) {
    
    base_mean <- settings$mean[settings$recorder == rec] 
    #calculated mean from settings set in the beginning
    
    mu <- base_mean * adj * sensi_effect("low") 
    #in µ : calculated mean before * effect applied on distances * effect applied on low sensitivities
    #µ is the mean used for the generation
    
    size <- settings$var[settings$recorder == rec] 
    #size = the opposite of dispersion (the higher size is, the closer we get from a Poisson distribution)
    #size = recorder's variance set, used as a proxy of distribution
    
    rnbinom(1, mu = mu, size = size) 
    #contact generation, row per row thanks to mapply() for each value :
    # 1 for randomizing ; mean sets on µ ; and no need for sd here, since our dispersion is set on variance with size
}
  
#### LOOP : ADJUSTING CONTACTS FOR HIGH SENSORS ####

high_gen <- function (design, recorder, sensi_effect) {
  
  for (s in unique(design$site)) { 
  #creation of a loop with "for" that will apply on sites (s label)
    
    for (r in recorder) { 
    #in s loop, creation of another one to work on recorders (r label)
    
    row_low  <- which(design$site == s & design$recorder == r & design$sensi == "low") 
    row_high <- which(design$site == s & design$recorder == r & design$sensi == "high")
    #creation of two objects row_x that include for each row their sensitivity modality
    
    value_low <- design$contacts[row_low] 
    #integrating the low contacts generated sooner in a new object
    
    # Generation of "high" from "low" x factor
    value_high <- value_low * (sensi_effect("high") / sensi_effect("low")) 
    #in this object, we multiply low values with the quotient of the sensitivity factor 
    #for instance : val_low * (1.2 / 0.8) = val_low * 1.5) 
    
    design$contacts[row_high] <- value_high 
    #including high calculated values to its correspondancy in design
    }
  } 
  #The loop applies before on each recorder, then in each site, row after row 
  #(otherwise we would have the same values everywhere)
  
return(design)
}

#### CONTACTS GENERATION FOR SENSIBILITY CURVES' TEST ####

summon_tcs <- function(rec, sensi_val) {
  
  base_mean <- settings$mean[settings$recorder == rec] 
  
  mu_tcs <- base_mean * sensi_val
  #in µ : calculated mean before * sensibilities factor
  #µ is the mean used for the generation
  
  size <- settings$var[settings$recorder == rec] 
  
  rnbinom(1, mu = mu_tcs, size = size) 
  #contact generation, row per row thanks to mapply() for each value :
  # 1 for randomizing ; mean sets on µ ; and no need for sd here, since our dispersion is set on variance with size
}

#### GATHERING COEFF OF GLM ####

get_eq <- function(coef1) {
  b <- coefs["(Intercept)"]
  if (coef1 != "adm") {
    b <- b + coefs[paste0("recorder", coef1)]
  }
  a <- coefs["sensi_val"]
  eq <- paste0("y_", coef1, " = exp(", round(a, 3), "x", round(b, 3), ")")
  return(eq)
}

#### PLOTING EQUATION WITH DATAFRAME ####

eq.df <- function(eq_df) {b <- coefs["(Intercept)"]
  if (eq_df != "adm") {b <- b + coefs[paste0("recorder", eq_df)]}
  return(b)}

#### ADJUSTING COEFF FOR THE LAST GLMM ####

adj_coeff <- function(adj_coeff) {
  b <- coefs["(Intercept)"]
  if (rec != "adm") {
    b <- b + coefs[paste0("recorder", rec)]
  }
  return(as.numeric(b))
}
