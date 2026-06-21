RATTES3_readme - this code completes the modeling for the endmember scenario assuming all clay sediment.  The ModelComboAnalysis_Aug2025clay.m code includes 5 functions: CombineMLRSY.m, CompareSY_MLRAug25.m, DeltaTerm.m, RiverSedDAinput.m, SDArivers.m

RATTES 3 combine the sediment yield and multiple linear regression model components, back projects the multiple linear regression model results for instances where the only available survey data post-date dam closure, and completes model analysis. 
The model calls output from RATTES1, assuming file structure is unchanged. 
You will need to add MLR input files... adding the prefix mlr from RATTES2 output. The prefix is added to differentiate with final model results:
 (example: mlr_cap_rattes_v1p1_clay.csv) 
You will need to add a copy of DeltaIDs.csv from the RATTES_final/ModelCombo_RATTES3 folder. 
-creates .mat file outputs and csv output
