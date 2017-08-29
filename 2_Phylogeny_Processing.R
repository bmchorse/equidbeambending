## Processing phylogeny data for paper.
# B. McHorse, 2016

# We use two composite trees to construct the tree for Figure 1 in this paper.

# Tree 1 is from: Fraser, D., Gorelick, R. & Rybczynski, N. 2015.
# "Macroevolution and climate change influence phylogenetic community assembly 
# of North American hoofed mammals." Biol J Linn Soc Lond 114, 485–494. 
# (doi:10.1111/bij.12457)

# Tree 2 is from: Jones, K. E. 2016 "New insights on equid locomotor evolution 
# from the lumbar region of fossil horses." Proc. R. Soc. B, 20152947. 
# (doi:10.1098/rspb.2015.2947)

# Both are time-scaled. Here we drop species, combine the trees, 
# do some trimming and splicing to get the tapirs where we want them,
# adjust, and save a genus tree to be imported into 
# 3_Phylogeny.R for plotting.

rm(list=ls())

require(phytools)
require(geiger)
require(paleotree)

frasertree <- read.nexus("./Trees/Fraser_Perissodactyl_timescaled_aba_edited.nex")
jonestree <- read.nexus("./Trees/Jones_Equid_timescaled.nex")

getgenustree <- function(x) {
      # Remove species names
      treegenus <- gsub('[_][A-z ]*', '' , x$tip.label)
      
      # Remove the periods that got left behind from sp. entries
      treegenus <- gsub("[.]", "", treegenus) 
      x$tip.label <- treegenus 
      
      # Get the indices (which()) that say TRUE for a duplicated() tip label
      dup.indices <- which(duplicated(x$tip.label)) 
      
      # Get rid of the duplicates.
      treetrim <- drop.tip(x, dup.indices) 
      return(treetrim)
}

# Trim both trees to genera only.
ftreetrim <- getgenustree(frasertree) 
jtreetrim <- getgenustree(jonestree)

# Add Miohippus
jtreetrim <- bind.tip(jtreetrim, "Miohippus", edge.length=10, where=15, position=2)
# Drop some duplicated taxa
jtreedrop <- drop.tip(jtreetrim, tips(jtreetrim, node=19)) 

# Splice the trees together
combotrim <- bind.tree(jtreedrop, ftreetrim, where=1)
plot(combotrim)
edgelabels(round(combotrim$edge.length,2),adj=c(0.5,-0.2),frame="none") # Show branch lengths

# Now the tapir node is in a weird place. Let's trim it and put it back where it goes.
tapirnode <- extract.clade(combotrim, node=65); tapirnode$root.edge <- 3 
combodrop <- drop.tip(combotrim, tips(combotrim, node=65)) # Drop tapir node
plot(combodrop) # Check to make sure the tapirs are gone
combodrop$root.edge <- 4 # Add root length
edgelabels(round(combodrop$edge.length,2),adj=c(0.5,-0.2),frame="none") # Show tree 
comboplus <- bind.tree(combodrop, tapirnode, position=3) # Stitch tapirs back on
comboplus$edge.length[5] <- 4 # Adjust the branch length after Mesohippus 

# Drop the perissodactyl genera we don't need
tipswedontwant <- c("Tapiravus", "Miotapirus", "Teleoceras", "Peraceras", "Aphelops", "Moropus")
tree <- drop.tip(comboplus, tipswedontwant) 
tree <- ladderize(tree, right=FALSE) # Reorganize the tree for easy viewing

tree <- drop.tip(tree, "Onohippidium") # Synonym of Hippidion
tree$edge.length[1] <- 65 # Make tapir go to present-day

# Plot the final version.
plot(tree)
edgelabels(round(tree$edge.length,2),adj=c(0.5,-0.2),frame="none") # Show branch lengths
axisPhylo() # Show the time scale

writeNexus(tree, "./Trees/treeforpaper.nex")
