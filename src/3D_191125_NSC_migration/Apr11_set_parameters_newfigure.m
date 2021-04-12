global inj_center cancer_center has_cancer d_w d_g chemo_sensitivity alpha4chmtx beta4dist FolderName1 FolderName2 modelType CONVERT2MICRON cancer_size;

%%% Standard values
n_seeds =           500;                %Total number of paths generated
Finaltimestep =     10000;              %Num of steps for each cell
% TIMESTEP =          ???;                %Experimentally calculated tick per real time
CONVERT2MICRON =    13.5;               %Avg μm per pixel in our figure
coh_limit =         0.4;                %coherency threshold
chmtx_limit =       0.1;                %chemotaxis threshold
inj_center =        [100 319 66];       %[400, 270, 180];    %In corpus callosum
cancer_center =     [300 400 175];      %[160, 280, 180];    %Near contralateral corpus callosum 
d_w =               5;                  %Reference step size on white matter
cancer_size =       [100/CONVERT2MICRON, 400/CONVERT2MICRON, 100/CONVERT2MICRON]; %radius [x y z]

%%% Conditional
has_cancer =        1;
modelNum =          2;

%%% Variables
alpha4chmtx =       1;                  %alpha(beta=1)  parameter for WM vs chemo
beta4dist =         1;                  %beta(alpha=1)  parameter for distance
chemo_sensitivity = 10;                 %chmtx_bias sensitivity
d_g =               0.1;                %reference step size on grey matter

%%% Test Values
alpha4chmtxArr =    [1 2 3 4];
beta4distArr =      [1 2 3 4];
chemo_sensArr =     [1 2 4 8];
d_gArr =            [0.01 0.1 1 2.5];

FolderName1 = '/Figures/CancerNA/';     %Save plots to this folder when there is no cancer (Must pre-exist)
FolderName2 = '/Figures/Cancer/';       %Save plots to this folder when there is cancer (Must pre-exist)