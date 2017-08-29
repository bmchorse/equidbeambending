# Process bone geometry data from BoneJ output files
## B. McHorse, 2016

## Set up the workspace and call relevant libraries
rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot) # Aplies pleasant default ggplot aesthetics

# Read in the data. 
path <- "./Data/"
# Pulls all .csv files from the directory in, names the object based on filename.
files <- list.files(path=path, pattern="*.csv")
for(file in files)
{
      perpos <- which(strsplit(file, "")[[1]]==".")
      assign(
            gsub(" ","",substr(file, 1, perpos-1)), 
            read.csv(paste(path,file,sep="")))
}

## Matching up metadata with filenames
bonejrecords$filename <- NA # Blank filename column
bonejrecords <- filter(bonejrecords, 
                       !(Species == "niobrarensis" & Element == "MCIII")) 
# Remove the two-part problem child scan

for (i in 1:length(bonejrecords[,1])) {
      checknum <- paste0(as.character(bonejrecords[i, grep("^SpecNo$", colnames(bonejrecords))]), "_") # checks specimen number 
      # and adds underscore to fix the regular vs. B numbering problem
      checkname <- as.character(bonejrecords[i, grep("^Species$", colnames(bonejrecords))]) # checks species column
      checkelement <- as.character(bonejrecords[i, grep("^Element$", colnames(bonejrecords))])
      grepout <- grepl(checknum, files, fixed=TRUE) & 
            grepl(checkname, files, fixed=TRUE) & 
            grepl(checkelement, files, fixed=TRUE)
      bonejrecords[i, grep("^filename$", colnames(bonejrecords))] <- files[grepout] # the grep bit 
      # assigns it to the 'filename' column without needing numerical index
}

bonejrecords$filename <- gsub(".csv", "", bonejrecords$filename, fixed=TRUE)
# Now the 'filename' column can be cross-referenced with object names in the R
# environment, so we can access those files to pull variable information out. 

## Get the variables from the files.
getvars <- function(records) { # the bonejrecords file, which must include a list of the object names (see above code)
      # make a blank df for results to go into
      vars <- data.frame(File=NA, Genus=NA, Species=NA, Element=NA,
                         A=NA, l=NA, IAP=NA, IML=NA, yAP=NA, yML=NA, mass=NA)
      
      for (i in 1:length(records$filename)) {
            slicelength.m <- records$voxDepth.mm[i]*0.1 # Convert mm to cm, to get voxel depth
            
            df <- noquote(records$filename[i]) # get the object name for this row
            x <- get(df) # call the object from the global environment!
            l <- length(x$Slice)*slicelength.m
            x$iap.m <- x$IAP*0.0001 # mm^4 to cm^4
            x$iml.m <- x$IML*0.0001 # Convert mm^4 to cm^4
            
            if(is.na(records$manualslice[i]) == TRUE) {
                  IAP <- x$iap.m[length(x$Slice)/2] # This is at midshaft only; rounds down if between slices
                  IML <- x$iml.m[length(x$Slice)/2] # IML is *about the ML axis*, so in the AP direction.
                  yAP <- x$RAP[length(x$Slice)/2]*0.1 # Convert radius about AP axis to cm
                  yML <- x$RML[length(x$Slice)/2]*0.1 # Same but for ML axis
                  A <- x$CSA[length(x$Slice)/2]*0.01 # Convert cross-sectional area from mm^2 to cm^2
            } else {
                  manualslice <- bonejrecords$manualslice[i]
                  IAP <- x$iap.m[manualslice]
                  IML <- x$iml.m[manualslice] 
                  yAP <- x$RAP[manualslice]*0.1 
                  yML <- x$RML[manualslice]*0.1 
                  A <- x$CSA[manualslice]*0.01 
            }
            Genus <- paste(noquote(records$Genus[i]))
            Species <- paste(noquote(records$Species[i]))
            File <- records$filename[i]
            Element <- paste(noquote(records$Element[i]))
            Mass <- paste(noquote(records$mass[i]))
            
            vars[i,] <- c(File, Genus, Species, Element, A, l, IAP, IML, yAP, yML, Mass)
      }
      return(vars)
      
}

# Get second moment of area along the whole bone, and save separately.
getIG <- function(records) { # the bonejrecords file, which must include a list of the object names
      for (i in 1:length(records$filename)) {
            slicelength.m <- records$voxDepth.mm[i]*0.1 # Convert mm to cm, to get voxel depth
            df <- noquote(records$filename[i]) # get the object name for this row
            x <- get(df) # call the object from the global environment!
            genus <- records$Genus[i]
            l <- length(x$Slice)*slicelength.m
            x$iap.m <- x$IAP*0.0001 # mm^4 to cm^4
            x$iml.m <- x$IML*0.0001
            x$A.m <- x$CSA*0.01 # mm^4 to cm^4
            
            out <- x %>% select(Slice, A = A.m, IAP = iap.m, IML = iml.m) %>% 
                  mutate(PercentOfBone = Slice/length(Slice),
                         normA = A^(1/2),
                         normIAP = IAP^(1/4),
                         normIML = IML^(1/4),
                         Genus = genus) %>% 
                  select(Genus,PercentOfBone,A, normA, IAP,IML,normIAP,normIML)
            write.csv(out, paste("./Data/IG/", gsub("-", "", df), "_I", ".csv", 
                                 sep=""), # note gsub to remove hyphens from filename
                      quote=FALSE,
                      row.names=FALSE)
      }
}

bbvars <- getvars(bonejrecords)
getIG(bonejrecords) # This will save internal geometry along the bone length.

# Because the Equus bone is a composite made with different voxel sizes, 
# AND was done in centimeters, we need to calculate length manually.
EniobrarensisMCIIIproximal_13478 <- read.csv("./SeparateEniobFiles/EniobrarensisMCIIIproximal_13478.csv")
EniobrarensisMCIIIdistal_13478 <- read.csv("./SeparateEniobFiles/EniobrarensisMCIIIdistal_13478.csv")

# Update slice numbers for the two-part scan and save as one single file
EniobrarensisMCIIIdistal_13478$Slice <- EniobrarensisMCIIIproximal_13478[1:340,4]+342
MCZ13478_Eniobrarensis_MCIII <- rbind(EniobrarensisMCIIIproximal_13478, EniobrarensisMCIIIdistal_13478)
write.csv(EniobrarensisMCIII_13478, "./Data/MCZ13478_Eniobrarensis_MCIII.csv")

lengthofEquusbone <- length(EniobrarensisMCIIIproximal_13478$Slice)*0.0098042 + length(EniobrarensisMCIIIdistal_13478$Slice)*0.0168072
lengthofEquusbone.m <- lengthofEquusbone*2.54 # convert inches to cm
getEquusvars <- function(x) {
      l <- lengthofEquusbone.m
      x$iap.m <- x$IAP*41.62 # inches^4 to cm^4
      x$iml.m <- x$IML*41.62 
      x$A.m <- x$CSA*6.452 # in^2 to cm^2
      
      # Save the full second moment of area data
      out <- x %>% select(Slice, A = A.m, IAP = iap.m, IML = iml.m) %>% 
            mutate(PercentOfBone = Slice/length(Slice),
                   normA = A^(1/2),
                   normIAP = IAP^(1/4),
                   normIML = IML^(1/4)) %>% 
            select(PercentOfBone,A, normA, IAP, IML, normIAP, normIML) %>% na.omit
      write.csv(out, "./Data/IG/MCZ13478_Eniobrarensis_MCIII_I.csv", quote=FALSE,
                row.names=FALSE)
      
      IAP <- x$iap.m[length(x$Slice)/2] # This is at midshaft only; rounds down if between slices
      IML <- x$iml.m[length(x$Slice)/2] # Note that IML is *about the ML axis*, so in the AP direction.
      yAP <- x$RAP[length(x$Slice)/2]*2.54 # Radius about the AP axis midshaft (so, mediolateral radius), convert to cm
      yML <- x$RML[length(x$Slice)/2]*2.54 # Same but for ML axis
      A <- x$CSA[length(x$Slice)/2]*6.452 # Convert cross-sectional area from in^2 to cm^2
      
      vars <- as.data.frame(cbind(A, l, IAP, IML, yAP, yML))
      return(vars)
}

Equus <- getEquusvars(MCZ13478_Eniobrarensis_MCIII) # We'll set slicelength to 0 since we'll override bone length anyway
Equus$l <- lengthofEquusbone.m
Equus2 <- data.frame(File = "MCZ13478_Eniobrarensis_MCIII", Genus = "Equus", Species = "niobrarensis",
                     Element = "MCIII")
Equusvars <- cbind(Equus2, Equus)
Equusvars$mass <- 350 # Manually add the mass
totalbbvars <- rbind(bbvars, Equusvars)

## Update TRI. 
TRI <- read.csv("TRI-specimenwise.csv", header=T)
bbvars <- left_join(totalbbvars, TRI)

## Add variables for beam bending analysis (see README for variable descriptions)

bbvars$beta <- 0 # Muscle effectively acts in line with long axis of bone
bbvars$theta_norm <- 5 # For normal locomotion
bbvars$theta_perf <- 20 # For high performance (accelerating, cutting, jumping)

bbvars <- bbvars %>% mutate_at(vars(Genus:Element), funs(factor)) %>% 
      mutate_at(vars(A:mass), funs(as.numeric)) # Set data types appropriately

# r and R will be scaled according to calculated values, using published Equus
# data as a starting point. 
Equusr <- 3 # 3cm is a published r for Equus
EquusR <- 10 # 10cm is a published R for Equus
Equusmass <- 450 # 450kg is the mass of the Equus for which r and R are known
bbvars$r <- Equusr*(bbvars$mass/Equusmass)^0.43 # r proportional to mass^(0.43)

# Now we will calculate EMA, which scales approximately with mass^(0.258),
# beginning with Equus. From this, we will back-calculate R, because EMA = r/R.
EquusEMA <- 3/10 # r/R
bbvars$EMA <- EquusEMA*(bbvars$mass/Equusmass)^0.258
bbvars$R <- bbvars$r/bbvars$EMA

# Add columns for the other beam bending variables.
# Fgrf will be set equivalent to a body-weight load. 
# We will have one regular Fgrf and one Fgrf scaled by TRI.
bbvars$Fgrf <- bbvars$mass*9.8 # Body weight load on a single leg.
bbvars$Fgrftoedness <- bbvars$Fgrf*(1/(2*bbvars$TRI+1)) # Body weight reduced proportionately by TRI.

write.csv(bbvars, "BBVars.csv", quote = FALSE, row.names = FALSE)
