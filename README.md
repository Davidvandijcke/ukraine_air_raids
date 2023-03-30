# Ukraine_AirRaids

The code in this repository can be used to replicate the results in the publication "Civil Response to Government Alerts Declines During Russian Invasion of Ukraine" (Forthcoming in PNAS), conditional on obtaining the required mobile device data from mobility data vendor [Veraset](https://veraset.com). All datasets considered confidential or sensitive have been removed from the replication package. To access the remaining datasets used in the publication, please refer to the accompanying OpenIPCSR repository. 


The code is made up of three folders:

1. `spark`: Code to parse the raw mobility data, should be run on a Spark cluster (e.g. AWS EMR) as the data is too large for local processing. 
2. `cleaning`: Code to clean and merge the various datasets used in the analysis. 
3. `analysis`: Code to perform the analysis and generate the figures in the paper. 





