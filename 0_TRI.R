## Calculations of TRI
# B. McHorse, 2016

## Set up the workspace and call relevant libraries
rm(list=ls())
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("TRIdata.csv", header=T)

### Method 1: take genus averages for digit length and calculate TRI
# (bigger error bars, but works when digits aren't associated, e.g., for Parahippus)
averages <-  ddply(df, c("Genus", "WhichDigit"), function(x) {
      n <- length(x$Genus)
      lmean <- mean(x$Max.articular.length, na.rm = TRUE)
      data.frame(n = n, l = lmean)
})
TRI <-  ddply(averages, "Genus", function(x) { # for each genus,
      lratio <- (x[x$WhichDigit == "side",][,4])/(x[x$WhichDigit == "center",][,4]) # divide side digit length (column 4) by center digit length
      data.frame(TRI = lratio)
})

### Method 2: calculate TRI for each specimen, then take a genus average for TRI
# (smaller error bars; works for most species, and results are similar to Method 1)
specimenTRI <- ddply(df, c("Genus", "Specimen.number"), function(x) { # for each specimen,
      meanside <- mean((x[x$WhichDigit == "side",][,8])) # take the average of the side digit length (col 7), then
      meancenter <- mean((x[x$WhichDigit == "center",][,8])) # take the average of the center digits, then
      TRI <- meanside/meancenter # divide avg side digit length by avg center digit length
      data.frame(TRI = TRI)
})

specimenTRI <- specimenTRI[complete.cases(specimenTRI$TRI),] # remove any NaNs

# Now take genus averages for TRI
specTRI <- ddply(specimenTRI, "Genus", function(x) {
      avgTRI <- mean(x$TRI)
      data.frame(TRI = avgTRI)
})

# We'll use the specimen-wise method (second method), except where it is not 
# available (e.g., Parahippus has no associated specimens). In that case,
# we will take the overall genus averages (first method).

saveTRI <- rbind(specTRI, TRI[TRI$Genus == "Parahippus",])

write.csv(saveTRI, "TRI-specimenwise.csv", row.names=FALSE)
