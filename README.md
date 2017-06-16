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
* data with functional imaging + atlas annotation  
3\. 20150410 (Example Dataset for Figures)  
4\. 20150417  
5\. 20160312  
7\. 20161004  
8\. 20161026  
12\. 20170112  

* data with functional imaging + dev info (birth time)  
6\. 20160328  
9\. 20161027 (missing birthtime)  
11\. 20170111 (missing birthtime)  
14\. 20170201 (missing birthtime)  

* data with functional imaging + IHC info (islet.mat)  
10\. 20161028  
13\. 20170126  
15\. 20170202  
16\. 20170216  

* data with functional imaging but no atlas information  
1\. 20140818  
2\. 20141006  

* data with MO_slc6a9 injection  
17\. 20170315  
18\. 20170323  
19\. 20170412  
20\. 20170503  
21\. 20170517  

### What is this repository for? ###
This repository is for data analysis code, figures, and manuscript of zebrafish developping spinal cord imaging data.

### Analysis list ###
%% Need to write

### Structure of directories ###
* .gitignore  -- file to be ignored in git update
* TempDat [ignored] -- a temporary storage place for data without duplicated units.
* Plots [ignored] -- figures
* Codes -- analysis codes
* Codes/Func -- shared functions across analysis
* Codes/Data_Analysis  -- each file stands for a code for a specific goal of analysis
* Data -- directories of datasets (raw data can be downloaded from https://www.dropbox.com/sh/g7wahoj4o3f0j04/AAAqjndWq2mjHgcEmqpvoumBa?dl=0)
* Atlas_Analysis -- Code to project neuronal location in recording to the common zebrafish atlas space
