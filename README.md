# CCWSIM

A guide for users 
Manuscript: “CCWSIM: fast and efficient CCSIM patch-based algorithm using Discrete Wavelet Transform to simulate categorical variables”
Requirements and important points: 
The present MATLAB scripts were written for the multiple-point geostatistical simulation by CCWSIM algorithm. These scripts perform categorical 2D simulation for binary variables. Furthermore, the scripts handle up to three level of the wavelet decomposition. We tried all the input parameters be sufficiently explained in the scripts for running the algorithm. Before running the script, consider some essential points as follow: 
For running the main CCWSIM script these subroutines are necessary: 
1- hd_resize_2D: It changes the simulation grid size based on wavelet decomposition level and relocates the hard data.
2- mincut_func: It finds the best minimum error boundary in simulation.
3- mincut: It finds the best minimum error boundary in simulation.
4- combine_2D: It merges desire pattern with previous pattern based on best minimum error boundary. 
5- hist_cat: It matches the facies proportion for categorical variable.
These subroutines are available as an open source codes in: https://github.com/SCRFpublic/MS_CCSIM.
For the simulation in the first level of the wavelet decomposition, the size of ti, simulation grid, OL and T should be even. 
For the simulation in the second level of the wavelet decomposition the size of ti, simulation grid, OL and T should be a factor of four. 
For the simulation in the third level of the wavelet decomposition the size of ti, simulation grid, OL and T should be a factor of eight. 
