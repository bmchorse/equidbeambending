# Table creation from results
## B. McHorse, 2016

## Set up the workspace and call relevant libraries
rm(list=ls())

library(dplyr)
library(tidyr)
library(htmlTable)

results <- read.csv("./Results/bendingresults.csv")
basics <- read.csv("bbvars.csv") %>% select(Genus, Species, avgTRI, mass) %>% distinct()

## Safety factors - for Table 1
# Get summary of SF (min, max, mean, sd) for POSTERIOR surface stress in REGULAR
# and PERFORMANCE locomotion, for both BODY and TRI loads
sf <- results %>% select(Genus, Species, Element, Posterior.surface.stress, 
                         load, condition) %>%
      filter(Element == "MCIII") %>%
      group_by(load, condition) %>%
      mutate(SF = abs(200/Posterior.surface.stress)) %>%
      summarize(Mean = round(mean(SF),1), SD = round(sd(SF), 1), 
                Min = round(min(SF), 1), Max = round(max(SF), 1))

htmlTable(sf)

## Metacarpal vs. Metatarsal difference - for Table 2
# Mean strain on posterior surface of MCIII minus mean strain on same for MTIII
# Shown for normal and performance trials
# Using ONLY genera for which we have the same specimen's MC and MT 
difftable <- results %>% separate(File, c("Number", "Binomial", "El"), 
                             sep = "_") %>%
      select(Number, Genus, Species, Element, Posterior.surface.stress, 
             load, condition) %>%
      filter(Number != "MCZNoNum") %>% # Remove non-numbered specimens; may not come from same individual
      filter(load == "TRI") %>% group_by(Number, Genus, 
                                              Element, condition) %>% 
      summarise(AvgPostStress = round(mean(Posterior.surface.stress), 1)) %>%
      spread(Element, AvgPostStress) %>% 
      mutate(Difference = round(MCIII - MTIII, 1), 
             PercDiff = round((Difference/MCIII*100), 1)) %>%
      filter(!is.na(Difference)) %>% arrange(desc(condition))

htmlTable(difftable)

## Overall summary - for general reference
oversumm <- results %>% select(Genus, Species, Element, Posterior.surface.stress,
                               load, condition) %>%
      filter(load == "TRI") %>% group_by(Genus, Element, condition) %>% 
      summarise(AvgPostStress = round(mean(Posterior.surface.stress), 1)) %>%
      spread(Element, AvgPostStress)  %>%
      arrange(desc(condition))
htmlTable(oversumm)

## Small summary of posterior stress - for reference when writing
smallsumm <- results %>% select(Genus, Species, Element, Posterior.surface.stress,
                                load, condition) %>%
      filter(Element == "MCIII") %>% group_by(load, condition) %>% 
      summarise(MeanStress = round(mean(Posterior.surface.stress), 1),
                SDStress = round(sd(Posterior.surface.stress), 1),
                MinStress = round(max(Posterior.surface.stress), 1),
                MaxStress = round(min(Posterior.surface.stress), 1)) %>% 
                      # we use min() and max() in opposite order than is intuitive here,
                      # because the highest stress values are in compression.
      arrange(desc(condition))
htmlTable(smallsumm)

## Small summary of anterior stress - for reference when writing
smallsumm2 <- results %>% select(Genus, Species, Element, Anterior.surface.stress,
                                load, condition) %>%
      filter(Element == "MCIII") %>% group_by(load, condition) %>% 
      summarise(MeanStress = round(mean(Anterior.surface.stress), 1),
                SDStress = round(sd(Anterior.surface.stress), 1),
                MinStress = round(max(Anterior.surface.stress), 1),
                MaxStress = round(min(Anterior.surface.stress), 1)) %>% 
      # we use min() and max() in opposite order than is intuitive here,
      # because the highest stress values are in compression.
      arrange(desc(condition))
htmlTable(smallsumm2)



