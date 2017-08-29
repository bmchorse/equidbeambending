# Plotting of the bending results
## B. McHorse, 2016

## Set up the workspace and call relevant libraries
rm(list=ls())
getwd()
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)

results <- read.csv("./Results/bendingresults.csv")
bbvars <- read.csv("BBVars.csv") %>% select(Genus, Species, TRI)
tipcolors <- read.csv("./Extra/TipColors-fromTRIContMap.csv", col.names = c("Genus", "Color"))

# Order the factors by TRI so it plots nicely
results$Genus <- factor(results$Genus, 
                       levels = unique(results$Genus[order(bbvars$TRI, 
                                                          decreasing=TRUE)]))

# Melt results and get the mean stress
mresults <- melt(results) %>% group_by(Genus, Element, load, condition, variable) %>% 
      summarise(se = sqrt(var(value)/n()), value = mean(value))

mresults$variable <- gsub("\\.", " ", mresults$variable) # Replace periods with spaces
mresults$variable <- tools::toTitleCase(mresults$variable) # Convert to title case

# Create an object so we can color genus names according to TRI values
TRIvalues <- bbvars %>% group_by(Genus) %>% summarise(TRI = mean(TRI))
tipjoin <- left_join(TRIvalues, tipcolors) %>% arrange(desc(TRI))
tips <- tipjoin$Color
names(tips) <- tipjoin$Genus

# Plotting function for results
basicplot <- function(data, cond, el, title, tipcolors = tips) {
      # First filter by condition and element
      stressvariables = c("Anterior Surface Stress", "Posterior Surface Stress")
      
      CE_filter <- data %>% filter(condition == paste(cond), 
                                   Element == paste(el)) 
      # Then select only the rows containing the variable(s) specified
      var_filter <- CE_filter[CE_filter$variable %in% stressvariables,]
      
      # Now trim the tip colors vector to only those found in the trimmed data
      var_filter$Genus <- droplevels(var_filter$Genus)
      tipkeep <- names(tips) %in% levels(var_filter$Genus) # Get True/False for colors found in trimmed data
      plottip <- as.character(tips[tipkeep]) # Drop the FALSE values
      
      # Print out a summary of the surface stresses
      datasumm <- var_filter %>% group_by(variable) %>%
            summarise(min = min(value), max=max(value), mean=mean(value))
      print(datasumm)
      # Adjust the fracture stress line
      if(datasumm[1,4] > 0 & datasumm[2,4] > 0) {
            fracstress <- 180
      } else if(datasumm[1,4] < 0 & datasumm[2,4] < 0) {
            fracstress <- 200
      } else {
            fracstress <- c(180, 200)
      }
      print(fracstress)
      SF <- c(fracstress/4, fracstress/2)
      myplot <- ggplot(var_filter, aes(Genus, abs(value), fill=load)) +
            geom_bar(position = "dodge", stat="identity", colour="black") +
            geom_errorbar(aes(ymin=abs(value)-abs(se), ymax=abs(value)+abs(se)),
                          width=0.25, position=position_dodge(.9), color='grey50') +
            # Add stress fracture for tension AND compression; later delete the irrelevant one.
            geom_hline(yintercept=fracstress, color="firebrick", linetype="longdash", size=1) +
            geom_hline(yintercept=SF, color='grey50', size=1) +
            scale_fill_manual(name = "Load", values = c("black", "white")) +
            scale_y_continuous(expand = c(0,0)) +
            facet_wrap(~ variable, nrow = length(stressvariables)) +
            ylab("Stress (MPa)") +
            ggtitle(paste(title)) +
            #Commented-out code below prints the stress value above each bar.
            # geom_text(aes(label=round(value, digits=0)),
            #           position = position_dodge(width = 0.9),
            #           color="gray50", vjust=-.5, size=2.8) +
            theme(plot.title = element_text(size = 14, face = "bold"),
                  text = element_text(size = 12),
                  axis.text.x=element_text(colour = plottip, face = "italic", 
                                           angle = 45, vjust = 0.5),
                  strip.text.x = element_text(size = 12),
                  strip.background = element_blank(),
                  legend.background = element_rect(size=.5),
                  legend.key = element_rect(size = 5),
                  legend.key.size = unit(1, 'lines'))
      return(myplot)
}

## Supplemental Materials plots
# Plot 1: Metacarpals, normal locomotion
plot1 <- basicplot(mresults, cond = "Normal", el = "MCIII", 
          title = "Stress in MCIII during normal locomotion")

# Plot 2: Metatarsals, normal locomotion
plot2 <- basicplot(mresults, cond = "Normal", el = "MTIII", 
                   title = "Stress in MTIII during normal locomotion")

# Plot 3: Metacarpals, performance
plot3 <- basicplot(mresults, cond = "Performance", el = "MCIII", 
                   title = "Stress in MCIII during performance locomotion")

# Plot 4: Metatarsals, performance
plot4 <- basicplot(mresults, cond = "Performance", el = "MTIII", 
                   title = "Stress in MTIII during performance locomotion")

# Save PDFs of figures.
# Figure 3
pdf(file = "./Figures/Figure3.pdf", onefile = TRUE)
plot3
dev.off()

# Supplemental figures
pdf(file = "./Supplemental/Metacarpal-stresses.pdf", onefile = TRUE)
plot_grid(plot1, labels='A')
plot_grid(plot3, labels='B')
dev.off()

pdf(file = "./Supplemental/Metatarsal-stresses.pdf", onefile = TRUE)
plot_grid(plot2, labels='A')
plot_grid(plot4, labels='B')
dev.off()