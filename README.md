# README #

## Redirect the git from bitbucket to github
In command-line terminal:
```bash
git remote rename origin bitbucket
git remote add origin https://github.com/zqwei/Zebrafish_spinal_cord_development.git
git remote remove bitbucket
git pull
```

## Current Stage
Please always put the current request of analysis to [Issue List](https://github.com/zqwei/Zebrafish_spinal_cord_development/issues).
1. Remove and plot *possible* duplicated units in the dataset
2. Covariance analysis: clustering
3. Covariance analysis: factor analysis

## Dataset list
1. 20140818 (no Atlas)
1. 20141006 (no Atlas)
1. 20150410 (**Example data in figures**)
1. 20150417
1. 20160312
1. 20160328
1. 20161004
1. 20161026
1. 20161027
1. 20161028
1. 20170111 (Lineage)
1. 20170112
1. 20170126
1. 20170201
1. 20170202
1. 20170216
1. 20170315
1. 20170323
1. 20170412 (MO)
1. 20170503 (MO)
1. 20170517 (MO)

### What is this repository for? ###
This repository is for data analysis code, figures, and manuscript of zebrafish developping spinal cord imaging data.

### Analysis list ###
* Covariance analysis: clustering
	* Re-index neurons based clustering
	* Covariance without ungroup neurons (ungroup neurons are those with low correlation to the others)
	* Development of ungroup neurons
* Covariance analysis: factor analysis
	* The variance of neurons explained by the FA change over time (number of factors is fixed)
	* The likelihood of the data at different time points as a function of the model created at a particular time point (number of factors is fixed)
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
