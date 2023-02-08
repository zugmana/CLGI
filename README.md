# CLGI

This repository contain data and code related to the work:
"Country-level gender inequality is associated with structural differences in the brains of women and men"
DOI: TBD

A full list with information on where to find the rawdata used in this manuscript can be found in the supplementary material. We cannot share the rawdata for individual participants included, however we provide a file that allows to run the meta-regression.

The provided code will load the result from escalc using the "Least Squares" option and additional info needed for the meta-regression. 
The steps to organize the rawdata (not available in this repository), filter out outliers and run the escalc function are in the script named "Data_prepration.R". This scripts output the file dataformeta.Rdata (available in /tables).
This file is then loaded for the meta-regression in "Thickness-GII.R"
This sample code is for the Cortical thickness analysis



  
