## Analyzing second moment of area and cross-sectional area across the length of the bone
## B. McHorse, 2016

# Setup workspace
rm(list=ls())
library(ggplot2)
library(reshape2)
library(cowplot) # Applies pleasant default ggplot aesthetics
library(dplyr)

bbvars <- read.csv("BBVars.csv", header=T)
tipcolors <- read.csv("./Extra/TipColors-fromTRIContMap.csv", col.names = c("Genus", "Color"))

## Internal Geometry
# Read in cross-sectional area (A) and second moment of area (I) data
path <- "./Data/IG/"
# Pulls all .csv files from the directory in, names the object based on filename.
files <- list.files(path=path, pattern="*.csv")
for(file in files)
{
      perpos <- which(strsplit(file, "")[[1]]==".")
      assign(
            gsub(" ","",substr(file, 1, perpos-1)), 
            read.csv(paste(path,file,sep="")))
}

# Add genus names to equus files
MCZ13478_Eniobrarensis_MCIII_I$Genus <- "Equus"
MCZ13478_Equusniobrarensis_MTIII_I$Genus <- "Equus"

# Get a list of our dataframes
objs <- ls()
# Separate out metacarpals from metatarsals (and drop other misc objects)
MC <- objs[grepl("MCIII", objs)]
MT <- objs[grepl("MTIII", objs)]

## Metacarpal plotting
# Combine data from all taxa into one long frame
plotMC <- data.frame(Genus = NA, PercentOfBone = NA, normA = NA, normIML = NA, normIAP = NA)

for(i in 1:length(MC)) {
            x <- get(MC[i])
            final <- x %>% select(Genus, PercentOfBone, normA, normIML, normIAP)
            plotMC <- rbind(plotMC, final)
}
plotMC <- plotMC[-1,]
plotMC$Genus <- as.factor(plotMC$Genus)

# Read in TRI values
TRIvalues <- bbvars %>% group_by(Genus) %>% summarise(TRI = mean(TRI))

# Sort the data to be plotted (in order of TRI)
plotdata <- left_join(plotMC, TRIvalues)
plotdata$Genus <- factor(plotdata$Genus,
                         levels = unique(plotdata$Genus[order(plotdata$TRI,decreasing=TRUE)]))
plotdata$PercentOfBone <- plotdata$PercentOfBone*100 # Make an actual percent
# Combine with TRI
levels(tipcolors$Color)[2] <- "#000000" # Replace tapir color with black for easier viewing
tipjoin <- left_join(TRIvalues, tipcolors) %>% arrange(desc(TRI))
tips <- tipjoin$Color
names(tips) <- tipjoin$Genus
# Trim to only those present
tipkeep <- names(tips) %in% levels(plotMC$Genus) # Get True/False for colors found in trimmed data
plottip <- as.character(tips[tipkeep]) # Drop the FALSE values

# Plot each internal geometry variable along the length of the bone
mcA <- ggplot(data=plotdata, aes(x=PercentOfBone, y=normA, color=Genus)) +
      geom_smooth(method="lm", formula = y ~ poly(x, 15), se=FALSE) +
      scale_color_manual(values=plottip) +
      labs(x = "Percent of Total Bone Length", y = "Normalized A") +
      theme(legend.text = element_text(face = "italic"))

mcIML <- ggplot(data=plotdata, aes(x=PercentOfBone, y=normIML, color=Genus)) +
      geom_smooth(method="lm", formula = y ~ poly(x, 15), se=FALSE) +
      scale_color_manual(values=plottip) +
      labs(y = expression(Normalized~I[ML]), x = "Percent of Total Bone Length") +
      theme(legend.text = element_text(face = "italic"))

mcIAP <- ggplot(data=plotdata, aes(x=PercentOfBone, y=normIAP, color=Genus)) +
      geom_smooth(method="lm", formula = y ~ poly(x, 15), se=FALSE) +
      scale_color_manual(values=plottip) +
      labs(y = expression(Normalized~I[AP]), x = "Percent of Total Bone Length") +
      theme(legend.text = element_text(face = "italic"))


## Now do all of same, but for metatarsals
plotMT <- data.frame(Genus = NA, PercentOfBone = NA, normA = NA, normIML = NA, normIAP = NA)

for(i in 1:length(MT)) {
      x <- get(MT[i])
      final <- x %>% select(Genus, PercentOfBone, normA, normIML, normIAP)
      plotMT <- rbind(plotMT, final)
}
plotMT <- plotMT[-1,]
plotMT$Genus <- as.factor(plotMT$Genus)

# Sort the data to be plotted (in order of TRI)
plotdataMT <- left_join(plotMT, TRIvalues)
plotdataMT$Genus <- factor(plotdataMT$Genus,
                         levels = unique(plotdataMT$Genus[order(plotdataMT$TRI,decreasing=TRUE)]))
plotdataMT$PercentOfBone <- plotdataMT$PercentOfBone*100

# Combine with TRI
tipkeepMT <- names(tips) %in% levels(plotMT$Genus) # Get True/False for colors found in trimmed data
plottipMT <- as.character(tips[tipkeepMT]) # Drop the FALSE values

# Plot internal geometry
mtA <- ggplot(data=plotdataMT, aes(x=PercentOfBone, y=normA, color=Genus)) +
      geom_smooth(method="lm", formula = y ~ poly(x, 15), se=FALSE) +
      scale_color_manual(values=plottipMT) +
      labs(x = "Percent of Total Bone Length", y = "Normalized A") +
      theme(legend.text = element_text(face = "italic"))
mtIML <- ggplot(data=plotdataMT, aes(x=PercentOfBone, y=normIML, color=Genus)) +
      geom_smooth(method="lm", formula = y ~ poly(x, 15), se=FALSE) +
      scale_color_manual(values=plottipMT) +
      labs(y = expression(Normalized~I[ML]), x = "Percent of Total Bone Length") +
      theme(legend.text = element_text(face = "italic"))
mtIAP <- ggplot(data=plotdataMT, aes(x=PercentOfBone, y=normIAP, color=Genus)) +
      geom_smooth(method="lm", formula = y ~ poly(x, 15), se=FALSE) +
      scale_color_manual(values=plottipMT) +
      labs(y = expression(Normalized~I[AP]), x = "Percent of Total Bone Length") +
      theme(legend.text = element_text(face = "italic"))

# Read in the allometry figures

source("Allometry.R")

allometry = plot_grid(allomA, allomI, labels="AUTO", align='v', ncol=2)
fig2 = plot_grid(allometry, mcIML, labels=c('', 'C'), ncol=1, rel_heights = c(1,1.3))


# Save files
pdf(file = "./Figures/Figure2.pdf", onefile = TRUE, width=8, height=6)
fig2
dev.off()

pdf(file = "./Supplemental/IntGeom_MC.pdf", onefile = TRUE, width=7, height=8)
plot_grid(mcA, mcIAP, ncol = 1, labels='AUTO')
dev.off()

pdf(file = "./Supplemental/IntGeom_MT.pdf", onefile = TRUE, width=7, height=12)
plot_grid(mtA, mtIAP, mtIML, ncol = 1, labels='AUTO')
dev.off()
