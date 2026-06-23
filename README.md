# RATTES v1.2

This repository contains the code required to reproduce the sediment yield (SY) and multiple linear regression (MLR) model components of RATTES. 

## Repository structure
#### ResNet1_2026_RATTES_Apr26_Matlab 
This folder contains the MATLAB code used for the SY component of RATTES. The code was developed using MATLAB version 2024b and subsequently verified using MATLAB 2025b. Testing was performed on a Windows-based computer. Typical run times are approximately 15 minutes for Model 1 (SY model) and 30 minutes for Model 3 (model combination and final analysis).

#### ResNet1_2026_RATTES_Apr26_R 
This folder contains the R code used for the MLR component of RATTES. The code was developed and tested using RStudio version 2025.05.1+513. Testing was performed on an Apple Mac computer. Typical run times are approximately 5 minutes for the model-training script and 30-40 minutes for the model-prediction script for each model configuration.

## Software requirements

To use this repository, clone or download the Github repository and run the code in the appropriate software environment.
- MATLAB: The SY component requires MATLAB 2024b or later. Information on obtaining MATLAB is available from Mathworks: https://www.mathworks.com/campaigns/products/trials.html.
- R and RStudio: The MLR component requires R and RStudio. Installation instructions are available here: https://rstudio-education.github.io/hopr/starting.html.

## Reproducing the results
Detailed instructions for running each model component and reproducing the results presented in the manuscript are provided within the corresponding folders of this repository.
