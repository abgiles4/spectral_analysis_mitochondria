# spectral_analysis_mitochondria
These R scripts can be used to perform multi-wavelength spectral analysis on absorbance data collected in mitochondria and heart.
To start, create a .csv file as follows: the first column = wavelength, the second column = dark current, the third column = incident light, and all other columns should be transmission data where each column represents one sample.
Use 2024_readcsv.R to convert transmission data to absorbance and then use the appropriate script to analyze absorbance data.
2024_nls_mito.R can be used to analyse data from mitochondria or heart data from Mb deficient mice.
2024_nls_Mb.R can be used to analyse heart data from Wt mice.
2024_nls_time.R can be used to perform time course analyses using data from mitochondria or Mb deficient mice.

To perform these analyses, you must collect a reference library to be used for analyses in the same system used for data collection. In order for the library to be compatible with the existing scripts, the library should be formatted as follows: column 1 = wavelength, column 3 = MbD, column 4 = MbO, column 7 = bH, column 8 = bL, column 10 = reduced cytochrome c1, column 13 = reduced cytochrome c. All unlisted columns can be used for additional reference spectra. The wavelength column of the reference library must match the transmission data.

These scripts were written in R version 4.4.0 and use packages: minpack.lm, ggplot2, zoo.

Nonlinear least-squares regression was conducted using routines for nonlinear least-squares optimisation in the minpack.lm package available in R (Elzhov, T. V. et al R Interface to the Levenberg-Marquardt Nonlinear Least-Squares Algorithm Found in MINPACK, Plus Support for Bounds, version 1.2-4, September 11, 2023).
