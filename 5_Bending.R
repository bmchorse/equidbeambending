# Bending analysis
## B. McHorse, 2016

## Set up the workspace and call relevant libraries
rm(list=ls())
library(dplyr)

bbvars <- read.csv("BBVars.csv", header=T)

####### Bending analysis -------------------------------------------------------------
beamanalysis <- function(data, scaling, type) {   # data is data frame with file info, 
      # genus/species, element, and the beam bending inputs like Fgrf. i.e., bbvars.csv
      # default is to not scale the grf
      
      # Pull variables from input data frame
      Genus <- data$Genus
      Species <- data$Species
      Element <- data$Element
      A <- data$A
      l <- data$l
      IML <- data$IML
      yML <- data$yML
      
      if(scaling=="Body") {
            Fgrf <- data$Fgrf
      }   else if(scaling=="TRI") { 
            Fgrf <- data$Fgrftoedness
      } else { 
            Fgrf <- data$Fgrf
            warning('No scaling specified; defaulting to body weight load for Fgrf.')
      }
      
      R <- data$R
      r <- data$r
      beta <- data$beta
      
      if(type=="Normal") {
            theta <- data$theta_norm
      }   else if(type=="Performance") { 
            theta <- data$theta_perf
      } else {
            theta <- data$theta_norm
            warning('No theta specified; defaulting to 5 degrees (normal locomotion).')
      }
      
      ## Calculate beam bending  
      Fm <- (Fgrf*R)/r
      h = l/2
      sigmac <- -(Fm * cos(beta*pi/180) + 
                        Fgrf * cos(theta*pi/180))/A # Compression stress; 
                        # convert degrees to radians because that's how R wants them
      sigmabAP <- ((Fm * sin(beta*pi/180) - 
                          Fgrf * sin(theta*pi/180)) * h * yML)/IML # Bending stress 
                        # in Anteroposterior direction (about the ML axis)
      
      tensionAP <- sigmac + abs(sigmabAP) # Tension on anterior surface of the bone:
                                    # Total stress reduced by tension counteracting compression
      compressionAP <- sigmac - abs(sigmabAP) # Compression on posterior surface of bone:
                                    # Total stress increased by compression + compression
      
      output <- as.data.frame(cbind(sigmac, sigmabAP, tensionAP, compressionAP))
      colnames(output) <- c("Compression stress", "Bending stress (AP)", 
                            "Anterior surface stress", "Posterior surface stress")
      outmpa <- round(output/100, 2) # These are in MPa
      info <- select(data, File:Element)
      results <- cbind(info, outmpa)
      return(results)
}

# Run beam analysis.
outputlist <- list(NA)
loads <- c("Body", "TRI")
conditions <- c("Normal", "Performance")

k = 1 # set iterator for output list
for(i in 1:length(loads)) {
      for(j in 1:length(conditions)) {
            bendresults <- beamanalysis(bbvars, loads[i], conditions[j]) 
            # that lets us iterate over both loads and both performance conditions
            bendresults$load <- loads[i] # Add a load column
            bendresults$condition <- conditions[j] # Add a condition column

            # this gives us a list of all the object names we'll want to call later
            outputlist[[k]] <- bendresults
            k = k + 1
      }
}

results <- outputlist %>% Reduce(full_join,.) # Merge all the outputs into a data frame.
write.csv(results, "./Results/bendingresults.csv", row.names=FALSE)
