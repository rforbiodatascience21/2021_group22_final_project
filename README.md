# Proteomics of SARS-Cov-2-infected host cells
Exam project of group 22 in the course R for Bio Data Science (22100).

**FLOW CHART**  
### Repository overview  
Input the repository figure when done

## Description
This project analyses proteomics data from a study of the effects of SARS-Cov-2 infection on host cells. The aim is to perform a full, tidy data analysis of the proteomics data.

### Data
The LC-MS/MS proteomics data was obtained from the article [Proteomics of SARS-CoV-2-infected host cells reveals therapy targets](https://www.nature.com/articles/s41586-020-2332-7?fbclid=IwAR3HEcdWjX3-4zTxGjXoiOtb2ol6iBMM6zt4uZ-ycECLEuu31KNJT_5uqaQ) by Bojkova et al. 

### Usage
This project includes a data analysis, a presentation and a shiny app. To generate the full data analysis and presentation, run the script 00_doit.R found in '/R'. This will install all required packages, before running all scripts in the appropriate order. Any data generated can then be found in '/data' and results and figures in '/results'.  

#### Shiny app
To run the shiny app follow the [link](https://signe-hjort.shinyapps.io/shiny/).

### Dependencies
All project code is written in R, and requires the following packages:  
* tidyverse  
* shiny  
* shinythemes  
* patchwork  
* readxl  
* broom  
* ggrepel  

### Authors
This project was made by:  
* Alban Laus Obel Slabowska (s200347)  
* Astrid Rasmussen (s193254)  
* Jens Waaben (s205972)  
* Signe Katrine Hjort (s174587)  
* Signe Fischer Ravn (s175410)
