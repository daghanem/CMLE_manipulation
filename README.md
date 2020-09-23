# CMLE_manipulation: The files in this repository help practitioners apply the censored maximum likelihood estimation in Ghanem, Shen and Zhang (2020) available free at https://www.journals.uchicago.edu/doi/full/10.1086/709649

1.	The two Stata ado files, gb2lfit_mod.ado and gb2lfit_II_mod.ado, included in the folder are supporting functions for the censored MLE. 
2.	Since we cannot disclose the original data that we use in the paper, we include a subfolder DEMO with a demonstration dataset (demodata.dta), which is one year of alternated Shanghai PM10 concentration data. We altered the original data by adding white noise. 
3.	The file Stata demo.do then uses this demo dataset to reproduce the figures and estimators that we reported in the paper with the real data.
