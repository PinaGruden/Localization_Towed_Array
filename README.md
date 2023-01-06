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
2) A1_Localize.m - this localizes user selected group(s) of animals (TDOA track(s)). It requires user input on which of the track fragments one wants to localize (since one track is often fragmented into separate track fragments). Results are saved into a table and contain information of localization (Cartesian coordinates and latitude/longitude), and perpendicular distance to boat trackline for localized sources. Localized groups are also plotted.


## Output

1. Otuput of simulated script are localization estimates and plots of the results.

2. Otuput of real data processing is a table containing localization estimate and perpendicular distance to boat trackline for each selected TDOA track (or fragments of the same TDOA track), and plots of the results. Specifically the saved results are:

- AStotal: final ambiguity surface, a 1 x M cell array where M is a 
          number tracked sources/groups. Each cell is an A x B matrix 
          where A is number of y coordinates, and B is number of x coordinates. 
- ASdilatetotal: final ambiguity surface with dilation, same dimensions 
                  as AStotal. 
- AStotal_hyperbolas: final surface with intersecting hyperbolas, same
                      dimensions as AStotal.
- Loc_table: table with localization information with 7 columns:
             ~ 'TrackID'- selected track indices corresponding to entries
                          of input structure 'Tracks'.
             ~ 'Loc_m'- localization estimate in Cartesian coordinates 
                        with respect to start of boat track (in m) 
             ~ 'Loc_LatLong' - localization estimate as latitude &
                              longitude (in decimal degrees)
             ~ 'Loc_m_dilated' - localization estimate from dilated AS   
                              in Cartesian coordinates with respect to  
                              start of boat track (in m)
             ~ 'Loc_LatLong_dilated'- localization estimate from dilated 
                                      AS as latitude & longitude (in 
                                      decimal degrees)
             ~ 'distance_m' - perpendicular distance from localization 
                             estimate ('Loc_m') to boat trackline (in m)
             ~ 'distance_m_dilated' - perpendicular distance from 
                                  localization estimate from dilated AS 
                                  ('Loc_m_dilated') to boat trackline (in m)

 
 
 


