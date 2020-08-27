global has_cancer d_w d_g TIMESTEP  chemo_sensitivity alpha4chmtx beta4dist FolderName1 FolderName2;

%%% Test values
beta4distArr =      [1 2 3 4];
alpha4chmtxArr =    [1 2 3 4];

%%% Standard values
modelNum =          2;
has_cancer =        0;
chemo_sensitivity = 1;
alpha4chmtx =       1;
beta4dist =         1;
coh_limit =         0.1;
chmtx_limit =       0.3;
cancer_radius =     1000;
% TIMESTEP =          588;
d_w =               5;
d_g =               2.5;

FolderName1 = '/Figures/CancerNA/';     %Save plots to this folder (Must pre-exist)
FolderName2 = '/Figures/Cancer/';       %Save plots to this folder (Must pre-exist)