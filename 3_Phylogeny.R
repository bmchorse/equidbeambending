## Phylogeny data for paper, showing TRI.
# B. McHorse, 2016

## Set up the workspace and call relevant libraries
rm(list=ls())
library(ape)
library(phytools)
library(dplyr)

# Get tree and TRI data.
tree <- read.nexus("./Trees/treeforpaper.nex")
bbvars <- read.csv("bbvars.csv")
avgTRI <- bbvars %>% select(Genus, TRI) %>% group_by(Genus) %>%
      summarise(TRI = mean(TRI))
bb <- data.frame(avgTRI[,2], row.names = avgTRI$Genus)

# Trim the tree to just taxa that we have TRI values for.
taxainboth<-intersect(tree$tip.label,rownames(bb)) # find the taxa that are found in both the tree and the data
dropset<-setdiff(tree$tip.label,taxainboth)# taxa in tree but not data
prunedtree<-drop.tip(tree,tip=dropset)
bb2 <- bb[row.names(bb) == taxainboth, , FALSE] # Selects all columns, sets drop rownames to False
plot(prunedtree) # Make sure all is well
axisPhylo() # Add the time axis

bb <- bb2$TRI
names(bb) <- rownames(bb2)

# Save the modified tree to file if needed in the future.
writeNexus(prunedtree, "./Trees/trimmedtreeforpaper.nex")

# Plot TRI on the phylogeny and save to file.
pdf(file = "./Figures/Figure1.pdf", onefile = TRUE, width=8, height=8)
TRI <- contMap(prunedtree, bb, res=200,fsize=c(1.25,1.1), 
               lims=c(0, 1), lwd=7, leg.txt="Toe Reduction Index")
dev.off()

# Get exact tip colors to use for reference later.
# From http://blog.phytools.org/2015/05/how-to-get-colors-of-tips-in-plotted.html
gettipcolors <-function(obj, tip){
      jj <- which(obj$tree$tip.label == tip)
      kk <- which(obj$tree$edge[,2] == jj)
      setNames(obj$cols[names(obj$tree$maps[[kk]])[length(obj$tree$maps[[kk]])]],
               NULL)
}
colors <- sapply(TRI$tree$tip.label, gettipcolors, obj = TRI)
colors
write.csv(colors, "./Extra/TipColors-fromTRIContMap.csv")

