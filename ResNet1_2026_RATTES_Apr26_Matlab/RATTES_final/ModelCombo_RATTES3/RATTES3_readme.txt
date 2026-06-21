RATTES3_readme
RATTES 3 combine the sediment yield and multiple linear regression model components, back projects the multiple linear regression model results for instances where the only available survey data post-date dam closure, and completes model analysis. 
-Figures 3-5 in the paper are created using this code
-Table 1 in the paper is created using this code
-Supplemental delta data table is created using this code
-Supplemental figures are created using this code. 
-The RATTESqueries function prints data queries which exist in the paper
TO RUN:
-Using output from RATTES1 (MATLAB code) and RATTES2 (R code) run RATTES3 
-The RATTES1 input is pulled into the code assuming the folder "ModelCombo_RATTES3" is stored in the same parent folder as "SYModel_MLRinput_RATTES1"
-The RATTES2 input files are generated in R. Add the prefix "mlr_" from RATTES2 output and use to run RATTES3. 
	-The prefix is added to differentiate with final model results:
 	(example: mlr_cap_rattes_v1p1.csv) 
-creates .mat file outputs and csv output

Empty folders exist to be populated with data in the code. 
