# Ukraine_AirRaids

The code in this repository can be used to replicate the results in the publication "Civil Response to Government Alerts Declines During Russian Invasion of Ukraine" (Forthcoming in PNAS), conditional on obtaining the required mobile device data from mobility data vendor [Veraset](https://veraset.com). All datasets considered confidential or sensitive have been removed from the replication package. To access the remaining datasets used in the publication, please refer to the accompanying OpenIPCSR repository. 


The code is made up of three folders:

1. `spark`: Code to parse the raw mobility data, should be run on a Spark cluster (e.g. AWS EMR) as the data is too large for local processing. 
2. `cleaning`: Code to clean and merge the various datasets used in the analysis. 
3. `analysis`: Code to perform the analysis and generate the figures in the paper. 


------ 
`spark` 
* `grab_veraset_v2.py`: This code loads the raw Veraset data from S3, filters on Ukraine, and optionally calculates "dwells" for each device.
* `runEventStudies.ipynb`: This code estimates the main event studies reported in the paper in parallel on the Spark cluster and saves them to S3. The resulting files should be downloaded locally (to the data/final folder) and the rest of the analysis can then be done on one's device, using the code in the other folders.
* `homelocs_vs.py`: Code to calculate the home locations of the devices, which are used for the "two-wall rule" and the representativeness robustness tests in the paper. Requires the "dwells" to be calculated. 
* `bootstrap_solved.sh`: Bash script that can be run as a bootstrap step on AWS EMR to make sure all the necessary Python libraries are installed on the Spark cluster. Provided for convenience. This code has some potentially superfluous infrastructure management that was required because of a bug in AWS, which may or may not have been resolved since. 

-----
`cleaning`
* `clean_air_raids.py`: Cleans and prepares the raw air alerts data downloaded from Volodymyr Agafonkin's repository (see link in the paper). 

-----
`analysis`
* `eventstudy.R`: This code produces the majority of the figures in the paper, including all variations of the main event study plot. Open it from the `analysis.Rproj` directory to make sure the paths are set correctly.
* `frontline.ipynb`: This code downloads the territorial control data from the VIINA database (see paper) and backs out an estimated frontline from it. 
* `representative.py`: This code produces the figures for the representativeness exercise in the SI, comparing the geographic density of the device data using the estimated home locations to WorldPop and Census data.
* `counterfactuals.do`: This code produces the estimated counterfactual results and figures in the paper, as well as the bootstrap estimates. 
* `alert_systems_plot.ipynb`: This code produces the worldwide alert systems figure in the paper. 




