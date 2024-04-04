# Localization Towed Array
 
 This package contains code for the model-based localization of marine mammals from towed array data. Ambiguity surfaces (probabilistic indicators of source positions) (Nosal, 2013; Barkley et al., 2021) are used to perform the localization. **Note, it is assumed that the boat moves in a straight line and no turns are present.**

 *References*
- Eva-Marie Nosal. Methods for tracking multiple marine mammals with widebaseline
passive acoustic arrays. _The Journal of the Acoustical Society of America_,
134(3):2383–2392, 2013.
- Yvonne M Barkley, Eva-Marie Nosal, and Erin M Oleson. Model-based localization
of deep-diving cetaceans using towed line array acoustic data. _The Journal
of the Acoustical Society of America_, 150(2):1120–1132, 2021.

Copyright (c) 2024, Pina Gruden 


## 1. Required Matlab version and toolboxes

This package was developed with **Matlab version 2022a (9.12)**. 

---

**IMPORTANT**: **Do not use earlier versions of Matlab** since certain functions have changed behavior. For example, function `ismembertol` has different behavior between version 2022a and older versions, and the processing will not work correctly if using older versions. Also, earlier versions of Matlab do not have certain functions such as `clim` and `FontSize`, so these will not work. Also, **do not use Matlab version 2023b (23.2.0.2485118)**, since function `ismembertol` also changed behavior and the package will not work correctly.

---


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
- `specify_paths.m`  - specifies paths to folders where data is located and results are saved to.


The package contains the following folders in the main folder:

1) **./AS_Localization_code/** - contains code to compute ambiguity surfaces and provide localization based on that. Functions included are:
- `computeAS.m` - computes an ambiguity surface (AS) for a given source based on the measured and modeled time difference of arrivals (TDOAs).
- `findpeaks2D.m` - finds peaks in the ambiguity surface and computes their widths.
- `get_xyrange.m` - limits the search range to a narrower area around the source location.
- `GetGPSandPhonePositions.m` - creates a table with GPS data and sensor positions for a given time-frame of the encounter.
- `localize_track.m` - localizes a selected TDOA track(s) and computes perpendicular distance between localized source and trackline.
- `localizetracks.m` - localizes user selected TDOA tracks/track segments. It gives a visual plot of the localization and offers an option for keeping/discarding localization and localizing mulitple tracks.
- `M2LatLon.m` - converts position in meters to lat/lon.
- `make\_2Dgrids.m` - creates a 2-D grid and also outputs a structure of grid parameters.
- `obtain\_ambiguitysurface.m` - computes an ambiguity surface for a given source for a specified grid and surface parameters. 
- `plot_AS.m` - plots the final ambiguity surface, boat and array positions.
- `plot_tracks.m` - plots all and selected TDOA tracks against Pamguard detections (if available) or against cross-correlograms.
- `select_tracks.m` - allows user to select which TDOA track / fragments they want to localize.
- `simulate_array_pos.m` - simulates sensor positions assuming the boat moves in a straight line specified by the slope and interstect.

2) **./Simulations/** -  code to perform localization with a model-based ambiguity surface approach using simulated data and scenarios.
3) **./Test_example/** - contains data to use for testing that the package works correctly. **IMPORTANT**: some of the test files were too big to be stored on GitHub. These are data for the folders `Crosscorrelograms` and `Raw_GPSandArray_info`. Each of these folders contains `README.txt` that contains a link where the files can be downloaded. 

 
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
To use with *Test data*- refer to Section 2, point 3).

For all other real data: This packages assumes that you have run **TDOA tracking master** package first and that you have extracted the available TDOA tracks. Optionally, if you have Pamguard detections available, you should have run **Extract Pamguard detections** package before running this current package.

This package expects a certain folder structure, that is outlined in Section 3.2.1.

First, specify paths and parameters you want to use for localization - see Section 3.2.2. After this, you can run the localization - see Section 3.2.3. 

#### 3.2.1. Expected folder structure
It is expected that each encounter will be in a separate folder and have its results (and intermediate files) saved in a separate folder from other encounters. So for **each encounter** one should have:

- a separate folder for Extracted TDOA tracks - i.e., where you stored results from ''TDOA tracking master'' package
- (Optional) a separate folder for Pamguard detections -  i.e., where you stored results from ''Extract Pamguard detections'' package 
- a separate folder for raw GPS data  
- a separate folder for extracted table from GPS data 
- a separate folder for Final localization results 

Note, most of these folders (apart from essential folders holding your data- i.e. extracted TDOA tracks) will be automatically created on a path you specify if they do not already exist. 

#### 3.2.2. Modify

Before running the package, specify the paths and parameters for your application by modifying the following scripts: 

- `specify_paths.m` - Specify folders where data is located and where results should be saved. The expected folder structure is outlined in 3.2.1. The expected data format are `.mat` files obtained from ''TDOA tracking master''  (e.g., `<Encounter>_Results.mat` and `<Encounter>_<whistles/clicks>_rawCrossCorrelogram_ALL.mat`) and optionally ''Extract Pamguard detections'' \\ 
    (e.g., `<Encounter>_Extracted_Annotated<Clicks/Whistles>.mat`) packages. If GPS and sensor depth information are available, these are expected to be in `.csv` format.

- `specify_parameters.m` - This is where you specify parameters for localization. Scroll down to change any parameters in the sections labeled ''CHANGABLE:'' as needed. The parameters are documented in the function.

#### 3.2.3. Run

Before running scripts below, your ''Current Folder'' must be navigated to the main folder of this package (then paths will get added automatically) OR add the main folder and subfolders to path manually. Either works!

Then run the package by running:

1) **[Optional]** `A0_GetGPSData.m`: If real GPS data is available, this reads all relevant GPS data into the Matlab workspace and computes sensor positions. It saves results into a `csv` table (at the location specified in `specify_paths.m`).

2) `A0_SetUp.m`: This reads all relevant data into the Matlab workspace. It also selects TDOA tracks and relevant hydrophone positions that will be used to obtain final localizations. `A0_SetUp.m` also generates figures that are used in the next step - `A1_Localize.m`. *Do not close these figures.*

3) `A1_Localize.m`: This localizes user selected group(s) of animals (TDOA track(s)). It requires user input on which track or which of the track fragments (since one track is often fragmented into separate track fragments) one wants to localize. This is interactive part where user is asked questions and is expected to interact through the commandline. The following questions are asked:

     - a. Which track/tracks you want to localize? Look at Figure 1 and 2 for track numbers.
     - b. Do you want to keep this localization? (1 for yes, 0 for no)
     - c. Localize another group? (1 for yes, 0 for no)
 
    Note, if there are not sufficient bearing changes in the selected tracks to achieve the localization (in question a.), then the user is prompted to try again.
    If the user wants to stop localizing, they need to choose 0 at question c. If the user wants to stop localizing at question a., then they should select something they localized already, then at the questions b. and c. select 0.
   
Results are saved into a table and contain localization estimates (Cartesian coordinates and latitude/longitude) and the perpendicular distance to the array for localized sources. Localized groups are also plotted.
    In case the process exists unexpectedly (e.g. computer crash, force quit) the temporary `TempResults.mat` is saved in the specified Results folder. Note, no other output is saved in this case and if the user wishes to keep this result, they must rename this file, otherwise it will be re-written next time `A1_Localize.m` is run.

## 4. Output

The package outputs:

1. Output of Simulated scenarios are localization estimates and plots of the results.

2. Output of Real data processing is a table containing localization estimates (in Cartesian coordinates and latitude/longitude formats) and the perpendicular distances to the array for each selected TDOA track (or fragments of the same TDOA track), and plots of the results. The variables that are saved in a `.mat` file (`Results_<Encounter>_StartTime_<First recording Date>_<First recording Time>_AnalyzeTime_<Today's Date>_<Time when saved>`) are:

- Table containing localization information (`Loc_table`) - each source is in a separate row. The table contains the following columns:
  - 'TrackID'- selected track indices;
  - 'Loc_m'- localization estimate in Cartesian coordinates with respect to start of boat track (in m);
  - 'Loc_LatLong' - localization estimate as latitude and longitude (in decimal degrees);
  - 'Loc_m_dilated' - localization estimate from dilated AS in Cartesian coordinates with respect to start of boat track (in m);
  - 'Loc_LatLong_dilated'- localization estimate from dilated AS as latitude and longitude (in decimal degrees);
  - 'distance_m' - perpendicular distance from localization estimate (\texttt{Loc\_m}) to the array (in m);
  - 'distance_m_dilated' - perpendicular distance from localization estimate from dilated AS ('Loc_m_dilated') to the array (in m).

- Final ambiguity surface (`AStotal`) - a 1 x M cell array where M is a number tracked sources/groups. Each cell is an A x B matrix where A is number of y coordinates, and B is number of x coordinates.
- Final ambiguity surface with dilation (`ASdilatetotal`) -  a 1 x M cell array where M is a number tracked sources/groups. Each cell is an A x B matrix where A is number of y coordinates, and B is number of x coordinates.
- Final surface with intersecting hyperbolas (`AStotal_hyperbolas`) - a 1 x M cell array where M is a number tracked sources/groups. Each cell is an A x B matrix where A is number of y coordinates, and B is number of x coordinates.


Example plot that is obtained as the output of this package is shown in Fig. 1.

![Results_Lasker_AC109_StartTime_20170912_165826_AnalyzeTime_20240223_1608](https://github.com/PinaGruden/Localization_Towed_Array/assets/62533526/2da87f33-5e40-4ec9-89da-8468b6fda7b0)
Fig.1 Example output plot for ''Localization towed array'' package.

 
 
 


