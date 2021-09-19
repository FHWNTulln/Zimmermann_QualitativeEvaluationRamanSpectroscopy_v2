# Spectra Data Analyser

This tool is an extension of the tool developed by Alexander Amberger.   
(https://github.com/FHWNTulln/Amberger_QualitativeEvaluationRamanSpectroscopy)

The script **Spectra_Data_Analyser_PLOT_ONLY.Rmd** delivers a plot of the original spectra only.   
It works very fast and the user can get a first overview how the spectra looks like.  
That´s helpful to decide which wavenumber areas should be reduced/excluded.

The script **Spectra_Data_Analyser.R** is an interactive tool for stepwise analysis of the spectra data. 
Here, the user can adjust parameters for data pretreatment as well as for the PCA. Outliers can be selected during the execution of the script and should be written down for later analysis. Finally, a PCA-LDA model can be build and validated using repeated double cross validation (RDCV).

The final step of the analysis can be done with the script **Spectra_Data_Analyser.Rmd** where all selected parameters from the previous steps can be set and beforehand selected outliers can be entered into lists.  
The script delivers a html file which can be stored for further analysis steps and for later comparisons.

All parameters which can be set are described in the scripts directly.

Steps of the analysis are:  
* read spectra data from multiple files  
* reduce wavenumber range  
* exclude wavenumber area  
* normalize data/calculate baseline  
* PCA  
* select outliers  
* PCA
* Cross-validated PCA-LDA

All necessary parameters can be set in the scripts and are included in the output files as well.  
There are 10 data pretreatment methods included in the scripts.  

From the R package *prospectr*:  
* Standard Normal Variate  
* Detrend  

Detailed information can be found at https://www.rdocumentation.org/packages/prospectr/versions/0.2.0/topics/prospectr-package  

From the R package *baseline*:
* Asymmetric Least Squares  
* Fill Peaks  
* Iterative Restricted Least Squares  
* Median Window  
* Modified Polynomial Fitting  
* Simultaneous Peak Detection and Baseline Correction  
* Robust Baseline Estimation  
* Rolling Ball  

Detailed information can be found at https://cran.r-project.org/web/packages/baseline/index.html

The PCA is done with the R function *prcomp*.  
Detailed information can be found at https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp  

For outlier selection, the R function *identify* is used.  
Detailed information can be found at https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/identify  

LDA is performed using the function *lda* from the package *MASS*.  
Detailed information can be found at https://cran.r-project.org/web/packages/MASS/index.html
