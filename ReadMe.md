ReadMe file for the supplemental code to Alther, Fronhofer & Altermatt (2021) 'Dispersal rate and riverine network connectivity shape the genetic structure of freshwater amphipod metapopulations' (DOI:)
Date: 2020-12-10
Author: Roman Alther, Eawag, Duebendorf, Switzerland
************************************************************************************

The ZIP folder contains three subfolders ('01_Simulation','02_Data_prep','03_Analysis') and three corresponding R scripts. All required data is provided.

If you only want to produce the figures from the paper, just run step 3 (see below).

In order to reconstruct the whole analysis, you will need to run three R scripts in consecutive order. Attention: If you rearrange files or folders you would need to change the R scripts accordingly (tidy work). We therefore suggest not touching the folder structure as given by the zipped file.

1) First run 'make_sim_package.R' in the '01_Simulation' folder. This will initiate the stochastic simulations, run in C++. The underlying network is in the folder 'indiv_sim/input/'. Depending on the number of available cores, the simulations take a few hours to weeks. We ran the analysis on a Ubuntu system on 18 cores simultaneously, resulting in approx. 3 days runtime. The simulations will generate about 60 GB of data, organized hierarchically in folders.
2) Run 'PopGenNet_Prep_20201210.R' in the '02_Data_prep' folder. This script relies on the output from the stochastic simulations and will prepare the data for subsequent analysis. The required microsatellite data and spatial information is provided as 'Gammarus_all_microsat_data_Rhein_2018.txt' and 'Gfos_Network.RData'. The preparation wil result in an 'output' folder with 18 Rdata files with summarized statistics from the stochastic simulations. Also it will create an 'PopGenNet_data.RData' file with additional data. Both, folder and files are needed for subsequent analyses.
3) Finally run 'PopGenNet20201210.R' in the '03_Analysis' folder. The script relies on two custom-made packages (OpenSwissRiverPlot and MultiPanel). They will be automatically installed from the web or can be installed from the included .tgz files (lines 31 and 37). If you did not run the previous two steps, you may rely on the existing 'PopGenNet_data.RData' in '02_Data_prep' and 'PopGenData.RData' in '03_Analysis'. Make sure to set the 'existing_data' parameter to TRUE (default). Otherwise the script will throw an error on line 278. The analysis script will produce all the figures from the paper, organized in a folder 'Analysis_DATE' and a subfolder 'SuppFigs'. Figures are prepared as vector graphics (PDF), but can be changed to pixelized figures by setting the parameter 'PDF' to FALSE (line 60).