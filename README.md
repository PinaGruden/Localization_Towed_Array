# Localization Towed Array
 
 This repo contains code for the localization with ambiguity surfaces for towed array data.
 
 It is developed with Matlab version 2022a. (if using older versions commands clim and FontSize for plotting do not exist)
 
 ## How to use
 
 1. To do *simulations*
 
 Go to folder ./Simulations/ and use *towedarray_LS_Simulations.m* for Simulation scenarios. 

Available simulations:
 
- Simulation 1: Stationary source, moving array, no noise on measurements

- Simulation 2: Stationary source, moving array, noisy measurements 

- Simulation 3: Moving source, moving array, no noise on measurements - Incorporates Dilation of Surfaces

- Simulation 4: Moving source, moving array with more sensors, no noise on measurements

- Simulation 5: Two stationary sources, moving array, no noise
 
 
2. To use with *Real data*

This packages assumes that you have run **TDOA tracking master** package first and that you have extracted TDOA tracks available. Optionally, if you have Pamguard detections available, you should have run **Extract Pamguard detections** package before running this current package.

Before running the package specify the paths and parameters for your application by modifying the following scripts: 

- specify_paths.m - Specify folders where data is located and where results should be saved to. The expected data format are .mat files obtained from **TDOA tracking master** and optionally **Extract Pamguard detections** packages.
- specify_parameters.m - This is where you specify parameters for localization. Change any parameters in the sections labeled “ CHANGABLE:” as needed. 


Then run the package by running:
1) A0_SetUp.m - this reads all relevant data into Matlab workspace, it also selects TDOA tracks and relevant hydrophone positions that will be used to obtain final localizations.
2) A1_Localize.m - this localizes a selected group of animals (TDOA track). It requires user input on which of the track fragments one wants to localize (since one track is often fragmented into separate track fragments). Results are displayed as Cartesian coordinates and latitude/longitude in the command prompt and as plots.


## Output

1. Otuput of simulated script are localization estimates and plots of the results.

2. Otuput of real data processing is an localization estimate for a selected TDOA track (or fragments of the same TDOA track), and plots of the results. 

 
 
 


