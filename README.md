# Bending performance and digit reduction in fossil equids
## A beam-bending analysis of metapodials through evolutionary time

This project uses CT scans of fossil equid (and one extant tapir) metapodials to assess how resistance to normal forces of locomotion - i.e., bending and axial compression - changes through evolutionary time as the side toes are reduced and body mass changes.

The paper can be found published in Proceedings of the Royal Society B as, "[Mechanics of evolutionary digit reduction in fossil horses (Equidae)](http://rspb.royalsocietypublishing.org/content/284/1861/20171174)". 
### Required packages (with version used in this analysis)
R version 3.3.2 (2016-10-31)

ape (4.0), cowplot (0.7.0), dplyr (0.5.0), geiger (2.0.6), ggplot2 (2.2.0), htmlTable (1.7), paleotree (2.7), phytools (0.5-38), plyr (1.8.4), reshape2 (1.4.2), tidyr (0.6.0)

### Scripts (in order)

0. Calculating TRI values - **0_TRI.R**
    - Inputs: TRIdata.csv
1. Reading in and cleaning data - **1_Processing.R**
    - Inputs: *./Data/bonejrecords.csv*, *TRI-specimenwise.csv*, and all BoneJ outputs for specimens from ./Data
    - Outputs: *BBVars.csv*, ./Data/IG (internal geometry data files)
2. Reading in and processing trees - **2_Phylogeny_Processing.R**
    - Inputs: *./Trees/Fraser_Perissodactyl_timescaled_aba_edited.nex*, *./Trees/Jones_Equid_timescaled.nex*
    - Outputs: *./Trees/treeforpaper.nex*
3. Producing phylogeny visualization of TRI - **3_Phylogeny.R**
    - Inputs: *BBVars.csv*, *./Trees/treeforpaper.nex* (time-scaled .nex tree file of Equidae)
    - Outputs: *trimmedtreeforpaper.nex* (Equidae tree with only taxa from this study), *./Figures/Figure1.pdf*, *./Extra/TipColors-fromTRIContMap.csv* (TRI colors to show relative digit state in later figures) 
4. Calculating internal geometry - **4_InternalGeometry.R**
    - Inputs: *BBvars.csv*, *./Extra/TipColors-fromTRIContMap.csv*, **Allometry.R**, all csv files in ./Data/IG
    - Outputs: *./Figures/Figure2.pdf*, *./Supplemental/IntGeom_MC.pdf*, *./Supplemental/IntGeom_MT.pdf*
5. Beam analysis - **5_Bending.R**
    - Inputs: *BBVars.csv*
    - Outputs: *./Results/bendingresults.csv*
6. Results graphs - **6_Plotting.R**
    - Inputs: *BBVars.csv*, *./Results/bendingresults.csv*, *./Extra/TipColors-fromTRIContMap.csv*
    - Outputs: *./Results/Figure3.pdf*, *./Supplemental/AllSurfaceStresses.pdf* 
7. Tables - **7_Tables.R**
    - Inputs: *BBVars.csv*, *./Results/bendingresults.csv*
    - Outputs: HTML summary tables of results (not directly saved)

Finally, Allometry.R is sourced by Step 4. 

### Data: ./Data (all files other than bonejrecords.csv), ./Data/IG (internal geometry data), ./SeparateEniobFiles
Spreadsheets of bone slice geometry, output from BoneJ analysis of CT stacks. ./Data/IG contains A and I along the length of each bone, rather than at midshaft (see variables described below). 
Within BoneJ, the anteroposterior and mediolateral orientations have been set for each bone.

**BoneJ output variables used in this study**
- Slice number (of downsampled CT stack)
- CSA: Cross-sectional area (mm<sup>2</sup>). Later called A.
- IAP: Second moment of area about anteroposterior axis (i.e., in mediolateral direction) (mm<sup>4</sup>)
- IML: Second moment of area about mediolateral axis (i.e., in anteroposterior direction) (mm<sup>4</sup>)
- RAP: Radius from center to outer surface of bone in anteroposterior direction (mm). Later called yAP.
- RML: Radius from center to outer surface of bone in mediolateral direction (mm). Later called yML.

### Metadata: ./Data/bonejrecords.csv
Spreadsheet containing information about each specimen scanned. 

**Variables**
- Collection: The collection where the fossil is reposited 
- SpecNo: Specimen number
- Genus
- Species
- voxSize: Voxel depth from CT scan (mm); used with Increment to calculate voxel depth
- Increment: n, where every nth slice from the CT scan is used during downsampling
- voxDepth: Voxel depth (mm); used with Slice number to calculate bone length
- Element: Metacarpal (MCIII) or Metatarsal (MTIII)
- manualslice: If the exact halfway point of the bone is broken or unsuitable for the scan, this provides a slice number to override
- mass: Estimated body mass of the genus or species
- mass.otherspecies: If a body mass is not available for one species in the genus but is for another, we indicate here which replacement species was used
- masssource: Literature reference for body mass estimate
- Side: Which side of the body the bone is from
- Notes: Any important notes about the specimen or scan

### Analysis Files
#### 0_TRI.R
Calculates toe reduction index (TRI) from lengths of proximal phalanges for the genera in this study.

1. Calculate TRI using genus average lengths of center and side digits (less robust method)
2. Calculate TRI for each specimen, then take genus average (more robust but requires toes to be associated with each other)
3. Save TRI values using second method for all genera except *Parahippus*, which uses the first method

#### 1_Processing.R
Processes BoneJ output files. Note that the *Equus* metacarpal is handled separately from the other files in most cases; it was scanned in two parts and units are in inches, thus needs particular treatment when calculating bone length, calling filenames, etc. 

1. Reads in all .csv files in the Data directory
2. Matches filenames with entries in metadata sheet
3. Calculates bone length 
4. Converts all units and retrieves variables from BoneJ output file for beam bending
5. Saves data files A and I along the length of the bone for all specimens
6. Handles two-part *Equus* metacarpal scan
7. Adds genus-average TRI to variable spreadsheet
8. Adds angle variables for beam bending
9. Calculates moment arm variables for beam bending based on scaling relationships
10. Saves variable spreadsheet BBvars.csv

#### BBvars.csv
This file is generated by the **BB_Processing.R** script.

##### Variables
In addition to those carried over from *bonejrecords.csv* and the BoneJ output files, both described above, the variables included in this file are:
- File: The filename of the datasheet for cross-referencing
- beta: The angle from the extensor muscle to the long axis of the bone. Here set to 0.
- theta_norm: The angle from the ground reaction force to the long axis of the bone, for normal locomotion.
- theta_perf: Same as above, for high-performance locomotion (e.g., accelerating or jumping).
- r: The moment arm of the extensor muscles (cm)
- EMA: Effective mechanical advantage, or r/R. Unitless.
- R: The moment arm of the ground reaction force (cm)
- Fgrf: The magnitude of the ground reaction force (N)
- Fgrftoedness: Same, but reduced proportional to TRI (higher TRI = larger side toes = reduced load on third digit)

#### 2_Phylogeny_Processing.R
Processes two time-scaled phylogenies to combine them for our purposes.

#### 3_Phylogeny.R
Plots TRI on a phylogeny of the genera used in this beam bending study.

1. Reads in a time-scaled phylogenetic tree of Equidae
2. Reads in TRI data for the genera in this study
3. Prunes phylogenetic tree to just the taxa from this study
4. Plots a continuous trait map of TRI onto the phylogeny
5. Saves the trait map phylogeny figure
6. Saves the exact color (= TRI) at the tip for each genus

#### 4_InternalGeometry.R
Uses internal geometry variables to calculate and plot second moment of area (I) and cross-sectional area (A) along the length of each bone. 

#### 5_Bending.R
Performs the beam bending calculations and saves results.

1. Read in the variable spreadsheet (*BBVars.csv*)
2. Define a bending function that uses engineering equations to calculate stress. See methods in manuscript for details.
  - **Inputs**
    - data: a dataframe that has all relevant metadata (BBvars.csv)
    - scaling: calls either body-weight load or body-weight scaled by TRI
    - type: calls either normal locomotion or high performance locomotion theta
  - **Outputs**
    - 'results' dataframe of compressive stress, bending stress for anteroposterior and mediolateral trials, and summed surface stresses for each trial
3. Set loads and conditions (includes both loads and both conditions)
4. Loop the bending analysis over all conditions and loads; write results to an object
5. Merge all results objects into a single data frame
6. Save results data frame as *bendingresults.csv*

#### 6_Plotting.R
Creates graphs for publication. 

1. Read in results from *bendingresults.csv* 
2. Define a bar graph function
3. Make bar graphs of stress at midshaft for various loading and performance conditions
4. Save graphs as PDF for either Figure 3 or for Supplemental Info

Note: some post-editing was performed in Illustrator to add TRI-colored labels, shift line for fracture stress depending on whether total loading is tensile or compressive, etc. No data values were changed.

#### 7_Tables.R
Creates an html table for Tables 1 and 2 in the paper, as well as other tables that were used to reference while writing. 

1. Read in results from *bendingresults.csv*
2. Create Table 2 (same-specimen difference between metacarpal and metatarsal stress)
3. Create table of average posterior stress for TRI-scaled loads
4. Create table summarizing metacarpal stresses
5. Create Table 1 (summary of safety factors in metacarpals)

Note that these tables are not automatically saved. The HTML output can be copied and pasted into Microsoft Word. 

#### Allometry.R
Tests for allometry in the scaling of cross-sectional area and second moment of area. 

1. Read in internal geometry and phylogeny data
2. Log-scale body mass and internal geometry variables
3. Calculate phylogenetic independent contrasts (PIC) to correct for evolutionary relatedness
4. Perform major axis regression to test for significant difference from isometry

### Results and other folders

#### ./Extra
Contains TRI hex codes for reference.

#### ./Figures
Contains figures output from plotting and phylogeny files, as well as externally created figures.

#### ./Manuscript
Contains manuscript files.

#### ./Results
Contains *results.csv* and a Word version of the html tables generated by **BB_Tables.R**.

#### ./Supplemental
Contains supplemental figures.
