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

I. Installation

1. Utilities

Copy the packages/ folder at the root of ImageJ (.../ImageJ/)

Make it recognizable by ImageJ:

    WINDOWS :
    ---------
	- modify the .cfg file of imageJ in the folling fashion :
	- it should look something like that :
		.
		jre\bin\javaw.exe
		-Xmx2500m -cp ij.jar  ij.ImageJ

	- add ";packages" after "ij.jar"
	- so you get something like that :
		.
		jre\bin\javaw.exe
		-Xmx2500m -cp ij.jar;packages   ij.ImageJ

     MAC OS :
     --------
	- ctrl-click on ImageJ.app and choose "show package contents", 
	- double click on the file info.plist
	- modify the properties of Java/ClassPath : You should add something like : 
		$APP_PACKAGE/../packages/
	- so that the list of pathes looks like : 
		$APP_PACKAGE/../ij.jar:$APP_PACKAGE/../packages:/System/Library/Java/Extensions/QTJava.jar 

     LINUX :
     ------- 
	- edit the "run" file and modify the call to java so that it doesn't execute ij.jar. 
	- add ij.jar and packages in the classpath : 
		./jre/bin/java -Xmx512m -cp ij.jar:./packages ij.ImageJ

2. Main plugin

Copy the FourierAnalysis/ folder in the .../ImageJ/Plugin/ folder 

II. Use

II.0. Data acquisition

The data analysis is much faster if the movie features square images with power of 2 sizes. It is highly recommanded to use this 
configuration. The algorithms are more efficient if the typical motion between 2 frames is (much) smaller than 1 pixel. This
should be the criterion to decide the frame rate at which the movie is recorded. If the exposure time can be several time smaller than the time between 2 frames, it is best.

II.1. DDM 

Refer to Wilson et al. Phys. Rev. Lett. 106 018101 (2011) for theoretical background/

II.1.a. Computing the DICF

Load the movie file in ImageJ. Launch DDM_Classical_5. Parameters are:

Maximum_lag_time: the maximum lag time for which the differential intermediate correlation function (DICF) is computed. Default 
is length_of_the_movie/5.

Number_of_Points_per_decade: The DICF is computed for logarithmically distributed lag times. This is the number of lag times in a decade for which the DICF is computed. Default is 20.

Averaging_every_n_frames: 1 + Number of frames skipped when running the temporal average. SET TO 1.

The Plugin then propose a name for the output file, in which the DICF will be saved for all wave number, and runs. Wait until completion.  

II.1.b. Fitting the DICF

The fitting stage allows to extract diffusion coefficient, fraction of swimmers, mean and std of the swimming speed.

Run the CorrFctFit_PartSwimmerFit_SecDiffLM plugin. Select the file you just created as input. Guesstimate of the parameters are
required to start the fitting procedure. They need to be realistic but not accurate.

 - Guess_D: Diffusion coefficient 

 - Guess_mean v: mean swimming speed

 - Guess_std dev v: std dev of the swimming speed (take it smaller than the mean...)

 - Guess_fractionSwimm: fraction of swimmers

Run and save the output file.

The output file contains the fit of the parameters, for the formula indicated in the first line, for each wave number (q)
separately. Chances are that the fits are innaccurate for very small and very large wave numbers, because the DCIF is not 
entirely decorrelated and noise level is reached, respectively. Descard them. Fitting parameters should be constant in the 
intermediate range of wave numbers, or your suspension does not contain swimming bacteria.

You may want to modify the CorrFctFit_PartSwimmerFit_SecDiffLM file to fit the DCIF to other types of motion. Feel free to 
contact us if you need help to do this. 


II.2. PhiDM 
