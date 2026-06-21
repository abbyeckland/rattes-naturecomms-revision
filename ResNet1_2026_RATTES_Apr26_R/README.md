# RATTES mlr model component
<u>Overview:</u> this code is for multiple linear regression model development and predictions for the second RATTES model component. This code reads in input files derived from ResNet (Hurst et al., 2025), dam parameter datasets (flow-conditioned parameter grids (FCPG) from Barnhart et al. (2020) and LakeCat from Hill et al. (2018)), and the sediment-contributing drainage area (A_sedMLR) timeseries, parameterized using the initial trap efficiency and design capacity for all ResNet reservoirs. The code pre-processes the data to ensure all requirements for multiple linear regression are met. Then, several MLR models are built using different combinations of parameter data. Then, the final MLR model to make predictions is selected and parameterized to estimate sedimentation rates for all unsurveyed RATTES reservoirs for each year between the first year after the prediction year (*yrp*), which is the year with the first available capacity data, and the end of the time loop (year 2050).

<u>Scripts:</u> R scripts and Python Jupyter Notebooks are utilized in this analysis. The R script **1_MLR_modeltraining.R** is the script that builds the MLR model. The R script **2_MLR_predictions.R** is the script where predictions are made using the model generated in script 1. 
- For **1_MLR_modeltraining.R**, in line 24, the user must specify the version out of two options: "rattes_v1p2_silt" or "rattes_v1p2_3param_silt". These are the two MLR models that can be accessed without changing any additional code in the script(s). Version "rattes_v1p2_silt" sets up the RATTES model with all four parameters (A_sedMLR, T, TE, and S), while "rattes_v1p2_3param_silt" excludes the T (time) parameter.
- For **2_MLR_predictions.R**, in line 24, the user must specify the version out of four options: "rattes_v1p2_silt", "rattes_v1p2_sand", "rattes_v1p2_clay", and "rattes_v1p2_3param_silt". Each option calls different input files, which were generated using different kappa values for the various grain sizes.
The two additional R scripts (**readin_lakecat_data.R** and **compute_wSAatdam.R**) are called within the main two scripts. The only Python code used in this analysis is contained within a Python Notebook (fcpg_data/obtain_FCPG_data.ipynb) which was utilized to first obtain the FCPG data, which are used as potential input data to the MLR model. This notebook does not need to be rerun, but is provided for reference.

<u>Input files:</u> The main input files are located in the folder "SYfiles", including: "MLR_Input_GRAINSIZE.csv", which is derived from ResNet, and "MLR_SAatDamInitial_GRAINSIZE.csv", which contains the sediment-contributing drainage area for each reservoir in ResNet calculated using design trap efficiencies and capacities. 

<u>Recommended steps for running code:</u>
1. Run **1_MLR_modeltraining.R** with version "rattes_v1p2_silt".
2. Run **2_MLR_predictions.R** with version "rattes_v1p2_silt". Run again with version "rattes_v1p2_sand". Run again with version "rattes_v1p2_clay".
3. Now that the code for the main 4-parameter model and predictions has been run, rerun code for the reduced model:
    - 3a. Run **1_MLR_modeltraining.R** with version "rattes_v1p2_3param_silt".
    - 3b. Run **2_MLR_predictions.R** with version "rattes_v1p2_3param_silt". 

<u>Other notes:</u>
- Upon running **MLR_toShare.R**, additional folders with figures, results dataframes, model summaries, etc will be automatically created. See the comments in the code for more details on these derivations.
- Zip files in folder "lakecat_data" must be unziped prior to running any code. The resulting csv files must be moved into the "lakecat_data" folder.
