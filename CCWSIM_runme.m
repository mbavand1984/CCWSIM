% ---------------------------------------------------------------------------------------------------------
% This script can perform 2D multiple-point simulation for categorical
% binary variables through CCWSIM algorithm.
 
% Manuscript: "CCWSIM: fast and efficient CCSIM patch-based algorithm using
% Discrete Wavelet Transform to simulate categorical variables"
  
%-----------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS:

%   ti             : A 2D training image should be imported. Three TIs with size are available:
%                     ti_Channel, ti_StoneWall and ti_Delta
%   hd             : for unconditional simulation according to the size of ti, a NaN matrix should be defined for example NaN(2000),               
%                     for conditional simulation, a simulation grid with conditional data should be imported.                                                               
%                     Three HDs are available according to available TIs,respectively: cnd_Channel, cnd_StoneWall and cnd_Delta.
%                     Also, to test the performance of the algorithm in the case of dense hard data,
%                     a file called cnd_Channel_Dense is prepared. 
%   CT             : determination of the Co-Template area size. normally, [2 2] is a suitable size
%                     for finding hard data in Co-Template area.       
%   OL             : size of the overlap area between the desired pattern and previous simulated 
%                     pattern. This parameter should be determined based on wavelet decomposition level. For detail, information
%                     users are invited to see the README file.                                 
%   cand           : number of candidate patterns that are selected based on minimum similarity distance.
%   fc             : 0 equal to no need for facies matching. Otherwise, the facies proportion should be provided  
%                     for example: for a 2 facies TI, this matrix should be like [0.28 0.72])
%   mrp            : Multiple random path flag (1:on, 0:off). 
%   T_vibration    : Multiple template size flag (1:on, 0:off)
%   prop           : the proportion of sorted candidates that are considered for finding best-matched pattern based on 
%                     hard data inside the template and in the co-template area.
%   real_numb      : number of realizations
%   wavelet_level  : number of wavelet decomposition. This script can handle up to three levels of the DWT.
%--------------------------------------------------------------------------------------------------------------------------------------------                                 
% Output Parameters
% - C: output Simulation grid 
% - error: The number mismatch HD ( zero means no mismatch)
%---------------------------------------------------------------------------------------------------------------------------------------------
% This program needs the following subroutines:
% - hd_resize_2D: It changes the simulation grid size based on wavelet decomposition level and relocates the hard data.
% - mincut_func: It finds the best minimum error boundary in simulation.
% - mincut: It finds the best minimum error boundary in simulation.
% - combine_2D: It merges desire pattern with previous pattern based on best minimum error boundary. 
% - hist_cat: It matches the facies proportion for categorical variable.
% - CCWSIM_main: It simulates the categorical binary variables in simulation grid.
% - CCWSIM_2D: It shows the simulation results and runtime. 
%--------------------------------------------------------------------------
%% 2D  multiple-point simulation through CCWSIM algorithm.

clear;clc
close all
load ti_Channel; 
%  hd = NaN(2000);       
load cnd_Channel
CT = [2 2];             
OL =32;                
T = 200;                
cand = 20;             
fc = 0;              
mrp = 1;                
T_vibration = 0;        
prop = 0.1;             
real_numb =1;        
wavelet_level = 3;      
 
[C , error] = CCWSIM_2D(ti, hd, T, OL, CT, prop, cand, fc, mrp, T_vibration, wavelet_level, real_numb);

