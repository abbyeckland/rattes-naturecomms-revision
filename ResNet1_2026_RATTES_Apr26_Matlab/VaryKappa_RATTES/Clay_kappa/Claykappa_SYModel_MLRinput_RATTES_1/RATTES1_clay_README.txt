# SedYield_MLRinput_Jan2026_clay is the first sedimentation model component- this one assumes all reservoir sediment is clay and adjust the kappa value in the Brown equation for clay.

# To run this code, you will need to add a copy of ResNetInput_SitesCanada_052225.csv. This is contained in the RATTES_final/SYModel_MLRinput_RATTES_1 folder.

This contains functions: AddInputsUpdates, RankTagApr25, SetCalcMatrix, FirstYieldFunctionApr25,FinalYieldFunctionApr25,MLRinputApr25
-Creates files MLR_Input%Date%.csv and MLR_SAatDamInitial%Date%.csv used to run MLR model
-Creates "SYModelResults_MLRinput.mat" loaded into final MATLAB model component to combine MLR model and SY model results for clay