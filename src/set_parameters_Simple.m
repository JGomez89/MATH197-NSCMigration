global modelType timeStep cancer_center chemo_sensitivity dmax dstochastic alpha alpha4chmtx beta4dist FolderName1 FolderName2;

modelNum = 2;

chemo_sensitivity = 10;
has_cancer = 0;
dstochastic = 8; 
alpha4chmtx = 1;
beta4dist = 1;
cancer_radius = 1000;
coh_limit = 0.1;  
chmtx_limit = 0.15;
% timeStep = 588;

FolderName1 = '\Figures\CancerNA\'; %Save plots to this folder (Must pre-exist)
FolderName2 = '\Figures\Cancer\'; %Save plots to this folder (Must pre-exist)