# FourierAnalysis

ImageJ Plugins for extracting dynamics from movies. Compatible with ImageJ 1.49 and java 1.8.0

DDM is composed of:

DICF calculation: DDM_Classical_5.java
DICF Fit: CorrFctFit_PartSwimmerFit_SecDiffLM.java
Fit utility: SecantDiffLevenbergMarquartFitter.java

PhiDM is implemented by: 

PhiDM_Version4_1.java

Local PhiDM (PIV method) is:

PhiDM_LocalVersion_Filter_MultiCore.java
with the utility: VelocityDisplayRoi2.java

Installation:

1. utilities

a. copy the packages/ folder at the root of ImageJ (.../ImageJ/)

b. open ImageJ.cfg file at the root with a text editor

c. modify '-cp ij.jar ij.ImageJ' into '-cp ij.jar;packages ij.ImageJ' and save

2. Main plugin

d. copy the FourierAnalysis/ folder in the .../ImageJ/Plugin/ folder 
