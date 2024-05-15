# synthLethal

A repository for analysis of synthetic lethal data in the BHK lab.

All the datasets that were used are too big to store them in Github, so please contact us at amadeu.carceller.pardina@gmail.com and we would be happy to share it. Without it, no networks can be created, but existing ones can be used.

In order to evaluate our networks:

1. Download the zip file of the repository and extract it. 
2. Open RStudio
3. Source the following files in the R folder: load_SL_pairs.R, enlight_validation.R, calculate_score.R, and auxiliary_functions.R
4. Install the corresponding packages if necessary.
5. Set the working directory to the synthLethal folder
   ```
   setwd('C:/Users/you/synthLethal')
   ```
6. Run the following code, which will directly plot the results of the clinical evaluations in your plot window:
   ```
    SL_pairs <- load_SL_pairs()
    results <- validate_enlight(SL_pairs)
   ```
The rest of the networks are also stored in the data folder. The create_network.R function is used to create the synthetic lethal networks. With all the data available, it is possible to run:
 ```
    SL_pairs <- create_network()
    results <- validate_enlight(SL_pairs)
   ```
Scripts used for network analysis are available at scripts/scripts_by_Amadeu
