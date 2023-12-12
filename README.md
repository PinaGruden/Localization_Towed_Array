# Localization Towed Array
 
 This package contains code for the model-based localization of marine mammals from towed array data. Ambiguity surfaces (probabilistic indicators of source positions) (Nosal, 2013; Barkley et al., 2021) are used to perform the localization.

 *References*
- Eva-Marie Nosal. Methods for tracking multiple marine mammals with widebaseline
passive acoustic arrays. _The Journal of the Acoustical Society of America_,
134(3):2383–2392, 2013.
- Yvonne M Barkley, Eva-Marie Nosal, and Erin M Oleson. Model-based localization
of deep-diving cetaceans using towed line array acoustic data. _The Journal
of the Acoustical Society of America_, 150(2):1120–1132, 2021.
 

## 1. Required Matlab version and toolboxes

This package was developed with **Matlab version 2022a (9.12)**. 

**IMPORTANT**: Do not use earlier versions of Matlab since certain functions have changed behavior. For example, function `ismembertol` has different behavior between version 2022a and older versions, and the processing will not work correctly if using older versions. Also, earlier versions of Matlab do not have certain functions such as `clim` and `FontSize`, so these will not work. 


The package uses the following Matlab toolboxes:
- *Signal Processing Toolbox*
- *Mapping Toolbox*
- *Image Processing Toolbox*

## 2. Package contents

This package contains functions and scripts to perform model-based localization with ambiguity surfaces. You can either use simulated data or real data. The code for simulations is available in a folder **./Simulations/**. The rest of the package is for real data. 

The package contains the following functions, scripts, and files in the main folder:
- `A0_GetGPSData.m` - main script that reads all relevant GPS data for a given encounter from a larger GPS database, and computes sensor positions for each time step.
- `A0_SetUp.m` - main script that reads all relevant data into the Matlab workspace, selects TDOA tracks and relevant sensor positions that will be used to obtain final localization.
- `A1_Localize.m` - main script that localizes sources using a model-based ambiguity surface approach.
- `README.md` - readme document specifying package usage.
- `specify_parameters.m` - specifies parameters for towed array localization.
- `specify_paths.m`  - specifies paths to folders where data is located and results saved to.


The package contains the following folders in the main folder:

1) **./AS_Localization_code/** - contains code to extract tables from SQLite3 database, and extract detections information from tables. Functions included are:
- `computeAS.m` - computes an ambiguity surface (AS) for a given source based on the measured and modeled time difference of arrivals (TDOAs).
- `findpeaks2D.m` - finds peaks in the ambiguity surface and computes their widths.
- `GetGPSandPhonePositions.m` - creates a table with GPS data and sensor positions for a given time-frame of the encounter.
- `localize_tracks.m` - localizes a selected TDOA track(s) and computes perpendicular distance between localized source and trackline.
- `M2LatLon.m` - converts position in Meters to Lat/Lon.
- `plot_AS.m` - plots the final ambiguity surface, boat and array positions.
- `plot_tracks.m` - plots all and selected TDOA tracks against Pamguard detections (if available) or against cross-correlograms.
- `select_tracks.m` - allows user to select which TDOA track / fragments they want to localize.
- `simulate_array_pos.m` - simulates sensor positions assuming the boat moves in a straight line specified by the slope and interstect.

2) **./Simulations/** - code to perform localization with a model-based ambiguity surface approach using simulated data and scenarios.
3) **./Test_example/** - contains data to use for testing the package works correctly.

 
 ## 3. How to use
 
 ### 3.1. To use with *Simulated data*
 
 Go to folder **./Simulations/** and use `towedarray_LS_Simulations.m` for Simulation scenarios. 

Available simulations:
 
- Simulation 1: Stationary source, moving array, no noise on measurements

- Simulation 2: Stationary source, moving array, noisy measurements 

- Simulation 3: Moving source, moving array, no noise on measurements - Incorporates Dilation of Surfaces

- Simulation 4: Moving source, moving array with more sensors, no noise on measurements

- Simulation 5: Two stationary sources, moving array, no noise
 
 
### 3.2.  To use with *Real data*

This packages assumes that you have run **TDOA tracking master** package first and that you have extracted the available TDOA tracks. Optionally, if you have Pamguard detections available, you should have run **Extract Pamguard detections** package before running this current package.

First, specify paths and parameters you want to use for localization - see Section 3.2.1. After this, you can run the localization - see Section 3.2.2. 

#### 3.2.1. Modify

- `specify_paths.m` - Specify folders where data is located and where results should be saved to. The expected data format are .mat files obtained from **TDOA tracking master** and optionally **Extract Pamguard detections** packages.  If GPS and sensor depth information are available these are expected to be in .csv format.

- `specify_parameters.m` - This is where you specify parameters for localization. Scroll down to change any parameters in the sections labeled “ CHANGABLE:” as needed. 

#### 3.2.2. Run

Then run the package by running:

1) **[Optional]** `A0_GetGPSData.m`: If real GPS data is available, then this reads all relevant GPS data into Matlab workspace, and computes sensor positions. It saves results into a .csv table (at the location specified in specify_paths.m).

2) `A0_SetUp.m` - this reads all relevant data into the Matlab workspace. It also selects TDOA tracks and relevant hydrophone positions that will be used to obtain final localizations.

3) `A1_Localize.m` - this localizes user selected group(s) of animals (TDOA track(s)). It requires user input on which of the track fragments one wants to localize (since one track is often fragmented into separate track fragments). Results are saved into a table and contain localization estimates (Cartesian coordinates and latitude/longitude), and perpendicular distance to the array for localized sources. Localized groups are also plotted.


## 4. Output

1. Otuput of simulated script are localization estimates and plots of the results.

2. Output of Real data processing is a table containing localization estimates (in Cartesian coordinates and latitude/longitude formats) and the perpendicular distances to the array for each selected TDOA track (or fragments of the same TDOA track), and plots of the results. The variables that are saved are:

- Table containing localization information (**Loc_table**)-each source is in a separate row. The table contains the following columns:
  - 'TrackID'- selected track indices corresponding to entries
                          of input structure 'Tracks'.
  - 'Loc_m'- localization estimate in Cartesian coordinates 
                        with respect to start of boat track (in m) 
  - 'Loc_LatLong' - localization estimate as latitude and
                              longitude (in decimal degrees)
  - 'Loc_m_dilated' - localization estimate from dilated AS   
                              in Cartesian coordinates with respect to  
                              start of boat track (in m)
  - 'Loc_LatLong_dilated'- localization estimate from dilated 
                                      AS as latitude & longitude (in 
                                      decimal degrees)
  - 'distance_m' - perpendicular distance from localization 
                             estimate ('Loc_m') to the array (in m)
  - 'distance_m_dilated' - perpendicular distance from 
                                  localization estimate from dilated AS 
                                  ('Loc_m_dilated') to the array (in m)

- Final ambiguity surface (**AStotal**)- a 1 x M cell array where M is a 
          number tracked sources/groups. Each cell is an A x B matrix 
          where A is number of y coordinates, and B is number of x coordinates. 
- Final ambiguity surface with dilation (**ASdilatetotal**)- a 1 x M cell array where M is a number tracked sources/groups. Each cell is an A x B matrix where A is number of y coordinates, and B is number of x coordinates.
- Final surface with intersecting hyperbolas (**AStotal_hyperbolas**)- a 1 x M cell array where M is a number tracked sources/groups. Each cell is an A x B matrix where A is number of y coordinates, and B is number of x coordinates.
  

 
 
 


