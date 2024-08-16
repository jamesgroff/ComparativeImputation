# Comparative Imputation of Missing Values (MVs) Codes
This is the repository of the codes, results and data for the Master project and thesis paper entitled "Comparative evaluation of imputaion methods in SILAC-based phosphoproteomics expereiments".

This is an accompanied repository for purposes of code sharing, reproducing and using the scripts to follow the general workflow of this project. 

![image](https://github.com/user-attachments/assets/41c7a72a-6760-4e88-87d6-88452882776a)

This will be updated if fixes or changes need to be made. Feel free to announce issues if things do not work as intended.

# How To Use the scripts

In the data folder, you have example data files of "proteinGroups.txt", "Phospho (STY)Sites.txt" and "Phospho (STY)Sites_normed.txt", as well as "overviewConditions.txt" all coming from the data of the following paper after automatic MaxQuant pipeline analysis based on default parameters in "overviewExample.docx": https://doi.org/10.1016/j.celrep.2021.109762

For those who work under Dengjel's laboratory in the University of Fribourg, there is also a direct access of those files that contains more details than this repository on the servers under "data/JMG".
For those who don't, contact the laboratory to have access of the specific files with the code "ZH_01" for this purpose.

Make sure that the session or working directory is under the Data folder before running the R codes.

You should then run the scripts in this order on either version to get the full steps and intermediary results you would find in the corresponding folder "Results": 
"DescriptionMV.R" -> "BenchmarkStatistics.R" -> "GenerateMatrixMV.R" -> "ImputationMV.R" (or any other imputation scripts) -> "VisualizeMetrics.R"

The major differences of the different imputation scripts are as follows: 
- "ImputationMV.R": Full data considered, no filters
- "ImputationWithoutFullRowMV.R": All the observations that only has MVs are discarded from the full data
- "ImputationWithoutFullOrNearMV.R": All the observations that has at most one complete value are discarded from the full data
- "ImputationWithFoldChangeMV.R": Same as the last one but with an added 1.5-fold change filter on the test statistics classification

An additional script called "CheckHistogramMetrics.R" is available as an option for a visual quality check of your data distribution.

There are also codes for the first two steps with customizable parameters that you can modify on top of the scripts, they are entitled "DescriptionMVwithParameters.R" and "BenchmarkStatisticsWithParameters.R".
