global modelType timeStep cancer_center chemo_sensitivity dmax dstochastic alpha beta4chmtx beta4dist FolderName1 FolderName2;

modelNum = 2;

% Test values
dstochasticArr = [1 2 4 8];
beta4distArr = [1 2 3 4];


chemo_sensitivity = 10;
% dstochastic = 8; %TRY: 1,2,4,8
alpha = 1; 
beta4chmtx = 5;
% beta4dist = 1; %TRY: 1,2,3,4
cancer_radius = 1000;
coh_limit = 0.1;  
% timeStep = 588;

FolderName1 = '/Figures/CancerNA/'; %Save plots to this folder (Must pre-exist)
FolderName2 = '/Figures/Cancer/'; %Save plots to this folder (Must pre-exist)
