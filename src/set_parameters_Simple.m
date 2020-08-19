global has_cancer modelType timeStep cancer_center dstochastic alpha4chmtx beta4dist FolderName1 FolderName2;

%%% Test values
dstochasticArr =    [1 2 5 8];
beta4distArr =      [1 2 3 4];
alpha4chmtxArr =    [1 2 3 4];

%%% Standard values
modelNum =          2;
has_cancer =        0;
dstochastic =       2;
alpha4chmtx =       1;
beta4dist =         1;
cancer_radius =     1000;
coh_limit =         0.1;
chmtx_limit =       0.15;
% timeStep =          588;

FolderName1 = '/Figures/CancerNA/';     %Save plots to this folder (Must pre-exist)
FolderName2 = '/Figures/Cancer/';       %Save plots to this folder (Must pre-exist)