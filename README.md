# Summary
Companion Matlab files implementing the analyses described in the summited paper:\

Bae SA, Androulakis IP. Mathematical modeling informs the impact of circadian oscillatory characteristics and meal patterns on insulin secretion. 2019.

# Abstract

# Instructions
Most of the folders start with + because they are "package folders." Drop the +hem folder anywhere in MATLAB's search path (go to File -> Set Path to change this) and then you can call the scripts and functions in this package from anywhere.\

+hem/+model: Files that collectively output the model solutions. Start from test_model.m\
+hem/+data: Data files\
+hem/+plot: Plotting functions\
+hem/+util: Any other functions: parameter estimation, post-processing, ...\
+hem/param: Model parameters, read by +model/run.m. 
