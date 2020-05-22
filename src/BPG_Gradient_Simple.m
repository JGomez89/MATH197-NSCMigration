% BIDIRECTOINAL PATH GENERATION

%%_________________
%Using 5 voxel kernerl blur orientation and coherency map 

%% load data 
% clear all
% directory = pwd;

% WM = imread(strcat(pwd,'/WM_M43.tif'));
% WM = logical(WM.*(WM>10));                                                                              % Threshold image to remove background
% ori_map = 180*imread(strcat(pwd,'/Orientation_5vox_blur.tif'))/pi;   % -90 to 90                                    % Orientation map converted to degrees (-90 to 90 deg)
% temp = zeros(size(ori_map));
% temp(ori_map>0) = ori_map(ori_map>0);                                                                   % positive angles are ascribed to 0-90 degrees
% temp(ori_map<0) = 180+ori_map(ori_map<0);                                                               % Negative angles are ascribed to 90 - 180 deg
% ori_map = temp; clear temp;
% coh_map = imread(strcat(pwd,'/Coherency_5vox_blur.tif')).*WM;

ori_map = 180*imread(strcat(pwd,'/Orientation_5vox_blur.tif'))/pi;   % -90 to 90                                    % Orientation map converted to degrees (-90 to 90 deg)
ori_map = -ori_map; 

%% parameters 
dmax = 10;                                                                                                  % Step size of each simulation step
coh_limit = 0.1;                                                                                        % Coherency limit for terminating path


%% Initialize a matrix of same size as the white matter image and initialize a seed point (inj center) around which the seeds for each simulation are placed
[X,Y] = meshgrid(1:size(WM,2),1:size(WM,1));
inj_center = [2150,1000];
%inj_center = [5100,730];
seed_ROI = sqrt((X-inj_center(2)).^2 + (Y-inj_center(1)).^2)<250;                                       % All points within 500 pixel distance from inj_center
seed_ind = find(seed_ROI);                                                                              % Indices of all these points
n_seeds = 20;                                                                                          % Total number of paths generated
WM_path = double(WM);  
% Storing the coorinate, angle of eigen vector and coherency value for both
% positive and negative direction paths that start from a seed point
p = struct('coord1',{},'angle1',[],'coherency1',[],'coord2',{},'angle2',[],'coherency2',[]);             % the three fields 



% cancer injection
cancer_center = [1500, 3500];
cancer_sd = 1500;
concentration = exp((-(X - cancer_center(2)).^2 - (Y - cancer_center(1)).^2)./cancer_sd^2);

cgradX = zeros(size(ori_map));
cgradX(2:end-1,:) = concentration(3:end,:) - concentration(1:end-2,:);

cgradY = zeros(size(ori_map));
cgradY(:,2:end-1) = concentration(:,3:end) - concentration(:,1:end-2);

% find largest vector magnitude
L = (cgradX.^2 + cgradY.^2).^(.5);
maxL = max(max(L));

% normalize gradient
cgradX = cgradX./maxL;
cgradY = cgradY./maxL;

Tstart = tic; 

%% for each seed 
for seed_loop = 1:n_seeds

    %%%% initialize first two steps 
    seed = seed_ind(round(rand*length(seed_ind)));                                                       % Choose a random point index to start the simulation
    [seed(1,1),seed(1,2)] = ind2sub(size(seed_ROI),seed);                                               % Coordinate of this starting point
    p(seed_loop).coord1(:,1) = seed';                                                                      % Starting location stored in data matrix
    ind = sub2ind(size(WM),seed(1,1),seed(1,2));                                                      % Index of this point in the white matter image
    % If coherency map data is present then use the angle indicated by the
    % data, else use a random angle
%     if coh_map(ind) > coh_limit 
        p(seed_loop).angle1(1) = ori_map(ind);                                                         % Starting orientation
%     else
%         p(seed_loop).angle1(1) = -90 + 2*rand*90;                                                       % Random angle generated if no coherency data
%     end

    %%%% second step 
    i = 2; 
    ind1 = sub2ind(size(WM),round(p(seed_loop).coord1(1,i-1)),round(p(seed_loop).coord1(2,i-1)));                                             % Index of previous coordinate point
            
    angle = ori_map(ind1);
    p(seed_loop).angle1(i) = angle;            
    p(seed_loop).coherency1(i) = coh_map(ind1);
    
    eigen_vect = [sind(angle),cosd(angle)];                                                     % Calculate the eigen vector (EV) direction

    dstep = (rand(1)*2-1)*dmax; %* double(coh_map(ind1)); 
    
    p(seed_loop).coord1(:,i) = p(seed_loop).coord1(:,i-1) + dstep*eigen_vect';                           % Find next coordinate point at step d along EV direction

    
    path = zeros(size(WM));
    
    flag1 = 0; 
    for i = 3:100000                                                                                    % Run the loop a large number of steps            
        if ~flag1  
            seedx = round(p(seed_loop).coord1(1,i-1)); 
            seedy = round(p(seed_loop).coord1(2,i-1)); 
            
            ind1 = sub2ind(size(WM),seedx,seedy);                                             % Index of previous coordinate point
            
            if coh_map(ind1) %> coh_limit                                                                           % Get Eigen vector angle at the previous coordinate point
                angle = ori_map(ind1);
                p(seed_loop).angle1(i) = angle; 
                
            else
                angle_sampled = -90+2*rand*90;
                angle = rem(angle_sampled + p(seed_loop).angle1(i-1),360);                              % Real value of angle if sampled angle + original angle > 360 deg
                p(seed_loop).angle1(i) = angle_sampled;
                
            end 
            
            % move with less magnitude if no WM present 
            if WM(ind1) 
                dmax = 5; 
            else
                dmax = 0.5; 
            end                
            
            p(seed_loop).coherency1(i) = coh_map(ind1);
            eigen_vect = [sind(angle),cosd(angle)];                                                     % Calculate the eigen vector (EV) direction
            
            dstep = rand(1)*dmax ; %%%%%%* double(coh_map(ind1)); 
           
            % Find next coordinate point at step d along EV direction
            p(seed_loop).coord1(:,i) = p(seed_loop).coord1(:,i-1) ... 
                - sign(dot( p(seed_loop).coord1(:,i-2) - p(seed_loop).coord1(:,i-1), eigen_vect ) + eps ) * dstep*eigen_vect';                           % Find next coordinate point at step d along EV direction
                % let it move away from the direction it came from. 
%                 + dstep*eigen_vect';   
            
%             p(seed_loop).coord1(1,i) = p(seed_loop).coord1(1,i) + dstep * cgradX(seednow(1,1),seednow(1,2));
%             p(seed_loop).coord1(2,i) = p(seed_loop).coord1(2,i) + dstep * cgradY(seednow(1,1),seednow(1,2));
            
% % the path is out of bounds or hit max steps then don't proceed and change flag
%             
% % ------------------------------------------                
            if all(p(seed_loop).coord1(:,i)' - size(WM)<=-1)  && all(round(p(seed_loop).coord1(:,i)')>1) % && ~path(round(p(seed_loop).coord1(1,i)),round(p(seed_loop).coord1(2,i)))
% %                 % to interpolate the path to integers 
% %                 xpts = round( linspace(p(seed_loop).coord1(1,i-1),p(seed_loop).coord1(1,i),10) );
% %                 ypts = round( linspace(p(seed_loop).coord1(2,i-1),p(seed_loop).coord1(2,i),10) );                
%                 
            else
                flag1 = 1;
            end
% % ------------------------------------------                
        end
    end
        
end

toc( Tstart ); 

%% Plotting the paths
% figure; imagesc(WM_path); axis equal;
%save 2D_5vox_blur_100k.mat
%edit_indices = [109,93,303,416,133,112,378];

% get the coordinate data from each coord1 and coord2 data and plot them.
X_coords = zeros(1,n_seeds);
Y_coords = zeros(1,n_seeds);
figure; imagesc(coh_map);   colormap gray;  hold on;
for i = 1:n_seeds
  
    coord_data1 = p(i).coord1; %cell2mat(p(i).coord1);
    l1 = size(coord_data1,2);
    if( l1 <= 5000 )
    X_coords = coord_data1(1, :);   
    Y_coords = coord_data1(2, :); 
    else
    X_coords = coord_data1(1,1:ceil(l1/5000):l1)';    
    Y_coords = coord_data1(2,1:ceil(l1/5000):l1)'; 
    end         
    plot(Y_coords,X_coords,'.', 'markersize', 2); %axis equal
        
end
% plot(Y_coords,X_coords,'.', 'markersize', 2); axis equal
% saveas(gcf,strcat('NSCMigration_pStay(',num2str(pStay),')_pCurr(',num2str(pCurr),')_deathParam(',int2str(deathParam),').fig'))

disp( 'done' ); 
