Name
====
MATLAB implementation for ``Fast Iterative Method for SOAV Minimization Problem with 
Linear Equality and Box Constraints and Its Linear Convergence''

## Overview
The supplementary PDF material and codes used for the numerical experiments 
are goint to be uploaded to this repository. 

## Description
### Supplementary PDF Material
``supplementary_material.pdf'' illustrates some simple numerical examples of 
the proximal operations of the SOAV-type cost 
by utilizing the proposed bisection algorithm and table based algorithms. 
### Codes for Numerical Examples
- Main File:
  1. main_run_LP_ADMM.m
- Sub-Routine Files:
  1. sub_MATaverage.m
  1. sub_instance_make.m
  1. sub_simulate_LP_ADMM.m
- SOAV Minimization Functions
  1. soav_bisec.cpp
  1. soav_table.cpp
  1. soav_conventional.cpp
  1. soav_LP_QP.m
