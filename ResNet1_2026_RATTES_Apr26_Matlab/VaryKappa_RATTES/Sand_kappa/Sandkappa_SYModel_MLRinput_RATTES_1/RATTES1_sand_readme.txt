# SedYield_MLRinput_Jan2026_sand is the first sedimentation model component- this one assumes all reservoir sediment is sand and adjusts the kappa value in the Brown equation for clay.
This contains functions: AddInputsUpdates, RankTagApr25, SetCalcMatrix, FirstYieldFunctionApr25,FinalYieldFunctionApr25,MLRinputApr25
-Creates files MLR_Input%Date%.csv and MLR_SAatDamInitial%Date%.csv used to run MLR model for clay
-Creates mat file loaded into final MATLAB model sand component to combine MLR model and SY model results

You will need to add the ResNetInput_SitesCanada_052225.csv from the RATTES_final/SYModel_MLRInput_RATTES1 folder. 