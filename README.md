# proj_TwoPhaseCure

## Introduction

This repository houses the R code used for simulation studies within the project titled "Two-phase designs with failure time processes subject to non-susceptibility." Inside, a folder labeled *utility* contains all utility functions. Additionally, there's a running script named *sim_script.R* designed for executing simulations.

## utility

This folder contains five *.r* files. 

*design.r* and *design_withV.r*: contains utitlity functions for phase II designs including stratified and residual-dependent sampling schemes;

*em_coxph.r*: contains functions for esitmation via EM algorithm;

*avar_louis.r*: contains functions for asymptotic variance computation for EM;

*opt_str_design_gau.r*: contains functions for optimal stratified designs, which are not feasible in practice but implementable in simulation studies. 


## sim_script.R

This script is designed to execute an example simulation study. By modifying the data configurations, it can replicate all the simulation results outlined in the manuscript. The script employs a for-loop to conduct repeated simulations. Within each simulation iteration, three key steps are executed: 1) data generation; 2) implementation of the two-phase design; and 3) estimation and inference.

## sample_dt.phII.csv

A toy data set showing the data structure after a phase II sub-sampling: 

- del: event status;

- time: observed times;

- x: expensive covariate, NA for missing;

- v: inexpensive covariate; 

- id: subject id;

- r: a binary indicator suggesting if the subject was selected into the phase II sub-sample. 
