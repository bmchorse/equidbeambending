# Testing allometry in internal geometry variables (A and I).
# Sourced from 2_InternalGeometry.R
## B. McHorse, 2016

# Setup workspace
detach("package:dplyr", unload=TRUE) # Because IntGeom will have dplyr loaded,
# and smatr calls plyr
library(smatr)
library(ape)
library(phytools)
library(dplyr)

colorvals <- as.character(tipcolors$Color)
prunedtree <- read.nexus("./Trees/trimmedtreeforpaper.nex")

internal <- bbvars %>% group_by(Genus) %>% dplyr::summarise(A = log(mean(A)),
                                                            IML = log(mean(IML)),
                                                            m = log(mean(mass)))
# Make separate named objects for PIC
# Note that names(data) and tree$tip.labels should be in the same order
A <- internal$A
names(A) <- internal$Genus
I <- internal$IML
names(I) <- internal$Genus
mass <- internal$m
names(mass) <- internal$Genus

## Phylogenetic independent contrasts -> regression
picA <- pic(A, prunedtree)
picI <- pic(I, prunedtree)
picM <- pic(mass, prunedtree)

#The argument -1 forces intercept through 0, which should be done for PIC regressions
(allomAcontrast <- sma(picA ~ picM -1, method='MA', slope.test=2/3)) 
(allomIcontrast <- sma(picI ~ picM -1, method='MA', slope.test=4/3))

plot(allomAcontrast)
plot(allomIcontrast)

allomplot <- function(ma, yvar, testslope) {
      
      slopetest = ma$groupsummary[[11]]
      
      ggplot(data=ma$data, aes(x = ma$data[[2]], 
                               y = ma$data[[1]])) + 
            geom_point(size=2) +
            geom_abline(slope = ma$groupsummary[[5]], intercept=0, color="Navy",
                        size=1.5) +
            geom_abline(slope=ma$groupsummary[[6]], intercept=0, color="gray") + # lower slope CI
            geom_abline(slope=ma$groupsummary[[7]], intercept=0, color="gray") + # upper slope CI
            xlab("Mass contrasts") + 
            ylab(paste(yvar, " contrasts")) +
            geom_abline(slope=testslope, linetype="longdash", intercept=0, 
                        color="black", size=1, alpha=0.9) #+
            #labs(title = paste("Slope = ", round(ma$groupsummary[[5]], 3),
                               #"\nSlope test p = ", round(ma$groupsummary$Slope_test_p, 3)))
      
}



allomA <- allomplot(allomAcontrast, "A", 2/3)
allomI <- allomplot(allomIcontrast, "I", 4/3)

# Diagnostic plots
plot(allomAcontrast, which='residual')
plot(allomAcontrast, which='qq')      
      
plot(allomIcontrast, which='residual')
plot(allomIcontrast, which='qq')      

# NOTE: these figures are called from BB_InternalGeometry.R.