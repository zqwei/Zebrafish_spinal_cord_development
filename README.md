# README #

###Current Stage###
1. Remove and plot *possible* duplicated units in the dataset
2. Covariance analysis: clustering
	* ZW: YW, Please check whether these results make sense to you
3. Covariance analysis: factor analysis

###Dataset list###
* Data_Dre_E1_BTXinjHuCH2BGCaMP6f_TL_20140818_045650_corrected_signal (2014-10-01)
* (-) Data_Dre_E1_HuCH2BGCaMP6f_0_20141006_041947_corrected (2014-11-01)
* Data_Dre_E1_HuCH2BGCaMP6f_0_20141006_041947_corrected (2015-01-14)
* Data_Dre_E1_HuCGCaMP6f-mnx1TagRFP_0-1_20150410_032910.corrected (2015-07-20)

### What is this repository for? ###
* This repository is for data analysis code, figures, and manuscript of **Yinan**'s dataset.
* version: 1.0

### Analysis list ###
* Covariance analysis: clustering
	* ZW: Re-index neurons based clustering
	* ZW: Covariance without ungroup neurons (ungroup neurons are those with low correlation to the others)
	* ZW: Development of ungroup neurons
* Covariance analysis: factor analysis
	* SD: The variance of neurons explained by the FA change over time (number of factors is fixed)
	* SD: The likelihood of the data at different time points as a function of the model created at a particular time point (number of factors is fixed)
* Time-series analysis: linear dynamical system
* Structure analysis: model prediction of functional connectivity vs anatomical structures
* Comparison between different datasets

### Structure of directories ###
* .gitignore  -- file to be ignored in git update
* TempDat [ignored] -- a temporary storage place for data without duplicated units.
* Plots [ignored] -- figures
* Text -- manuscripts, slides etc.
* Codes -- analysis codes
* Codes/Func -- shared functions across analysis
* Codes/Data_Analysis  -- each file stands for a code for a specific goal of analysis
* Codes/Old_code -- older version of codes, which are depreciated
* Codes/SD_code -- codes previously contributed from Shaul Druckmann
* Data -- directories of datasets

### Contribution/comments guidelines ###
* Text                             -- Leave the comments in *Latex* file beginning with initials like ZW, YW, SD, PK
* Code                           -- Keep a comment in the commit in detail before pushing to git