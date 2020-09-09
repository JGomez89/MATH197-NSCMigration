global inj_center cancer_center has_cancer d_w d_g chemo_sensitivity alpha4chmtx beta4dist FolderName1 FolderName2;

%%% Standard values
n_seeds =           10;                 %Total number of paths generated
Finaltimestep =     5000;               %Num of steps for each cell
% TIMESTEP = ???                          %Experimentally calculated realtime to tick
inj_center =        [400, 270, 180];      %In corpus callosum
cancer_center =     [160, 280, 180];    %Near contralateral corpus callosum 

%%% Conditional
has_cancer =        0;
modelNum =          2;

%%% Variables
chemo_sensitivity = 1;                  %chmtx_bias sensitivity
alpha4chmtx =       1;                  %alpha(beta=1)  parameter for WM vs chemo
beta4dist =         1;                  %beta(alpha=1)  parameter for distance
coh_limit =         0.4;                %coherency threshold
chmtx_limit =       0.1;                %chemotaxis threshold
d_w =               5;                  %standard step size on white matter
d_g =               0.05;               %standard step size on grey matter

FolderName1 = '/Figures/CancerNA/';     %Save plots to this folder when there is no cancer (Must pre-exist)
FolderName2 = '/Figures/Cancer/';       %Save plots to this folder when there is cancer (Must pre-exist)