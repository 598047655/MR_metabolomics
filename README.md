**README**

**Overview**
This repository contains R scripts and a CSV file for generating various plots and analyses related to metabolite data in cerebrospinal fluid (CSF). The scripts include Mendelian Randomization analyses, forest plots, and circular heatmaps to visualize the results.

**Files**

**CODE.R**   **Description**: This script performs Mendelian Randomization analysis specifically for 338 CSF metabolites. It extracts exposure and outcome data, harmonizes the data, runs the MR analysis, and generates various plots, including scatter plots, funnel plots, forest plots, and leave-one-out plots. It also includes heterogeneity and pleiotropy tests. The analysis for plasma metabolites is similar and hence not repeated here.
**Usage**: Ensure all necessary libraries are installed and load the required data files. Execute the script to perform the full analysis and generate the results and plots.

**CSF_Forest_Plot_CODE.R**
**Description**: This script generates a forest plot for CSF metabolite data. It uses the forestploter library to create a visual representation of the odds ratios and confidence intervals for various metabolites.
**Usage**: Load the me.csv data file and run the script to generate and save the forest plot as a PNG file.

**me.csv**
**Description**: This CSV file contains the data required for generating the forest plot in the CSF_Forest_Plot_CODE.R script. It includes information on metabolites, their mean values, and confidence intervals.
**Usage**: Ensure this file is in the same directory as the CSF_Forest_Plot_CODE.R script or provide the correct path within the script.

**Circular_Heatmap_CODE.R**
**Description**: This script generates a circular heatmap to visualize consistent significance across different Mendelian Randomization methods. The heatmap highlights the metabolites with significant and consistent beta values across various MR methods.
**Usage**: Load the necessary data and run the script to generate the circular heatmap. Customize the paths and parameters as needed within the script.

****How to Use**

**CODE.R:****
Install required libraries if not already installed.
Modify the script to point to the correct data files and parameters.
Execute the script to perform the analysis and generate the plots.

**CSF_Forest_Plot_CODE.R:**
Make sure the me.csv file is in the correct location.
Install the forestploter library if not already installed.
Run the script to create and save the forest plot.

**me.csv:**
No action needed other than ensuring it is in the correct directory for the CSF_Forest_Plot_CODE.R script.

**Circular_Heatmap_CODE.R:**
Ensure the necessary libraries are installed.
Adjust file paths and parameters as needed.
Run the script to generate the circular heatmap.

**Notes**
Customize file paths and parameters within the scripts to fit your data and analysis needs.
Ensure all necessary libraries are installed before running the scripts.
Check for any additional dependencies that might be required for specific functionalities in the scripts.

**Contact**
For any questions or issues, please contact the author at [22218425@zju.edu.cn].
