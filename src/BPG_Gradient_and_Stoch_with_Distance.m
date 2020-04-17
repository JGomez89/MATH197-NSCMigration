% BIDIRECTOINAL PATH GENERATION

%%_________________
%Using 5 voxel kernerl blur orientation and coherency map 

clear all
directory = pwd;

WM = imread(strcat(pwd,'/WM_M43.tif'));
WM = logical(WM.*(WM>10));                                                                              % Threshold image to remove background
ori_map = 180*imread(strcat(pwd,'/Orientation_5vox_blur.tif'))/pi;                                      % Orientation map converted to degrees (-90 to 90 deg)
temp = zeros(size(ori_map));
temp(ori_map>0) = ori_map(ori_map>0);                                                                   % positive angles are ascribed to 0-90 degrees
temp(ori_map<0) = 180+ori_map(ori_map<0);                                                               % Negative angles are ascribed to 90 - 180 deg
ori_map = temp; clear temp;

coh_map = imread(strcat(pwd,'/Coherency_5vox_blur.tif')).*WM;
d = 5;                                                                                                  % Step size of each simulation step
coh_limit = 0.1;                                                                                        % Coherency limit for terminating path

%% Initialize a matrix of same size as the white matter image and initialize a seed point (inj center) around which the seeds for each simulation are placed
[X,Y] = meshgrid(1:size(WM,2),1:size(WM,1));
inj_center = [2150,1000];
%inj_center = [5100,730];
seed_ROI = sqrt((X-inj_center(2)).^2 + (Y-inj_center(1)).^2)<500;                                       % All points within 500 pixel distance from inj_center
seed_ind = find(seed_ROI);                                                                              % Indices of all these points
n_seeds = 500;                                                                                          % Total number of paths generated
WM_path = double(WM);  
% Storing the coorinate, angle of eigen vector and coherency value for both
% positive and negative direction paths that start from a seed point
p = struct('coord1',{},'angle1',[],'coherency1',[],'coord2',{},'angle2',[],'coherency2',[]);             % the three fields 


% boxplot of distances
num_intervals = 10;
x_intervals = zeros(2*n_seeds,num_intervals);
y_intervals = x_intervals;
timeIntervals = 30;
num_iterations = 300;

% cancer injection
cancer_center = [3500, 2500];
cancer_param = 1250;
concentration = exp((-(X - cancer_center(2)).^2 - (Y - cancer_center(1)).^2)./cancer_param^2);

cgradY = zeros(size(ori_map));
cgradY(2:end-1,:) = concentration(3:end,:) - concentration(1:end-2,:);

cgradX = zeros(size(ori_map));
cgradX(:,2:end-1) = concentration(:,3:end) - concentration(:,1:end-2);

% find largest vector magnitude
L = (cgradX.^2 + cgradY.^2).^(.5);
maxL = max(max(L));

% normalize gradient
cgradX = cgradX./maxL;
cgradY = cgradY./maxL;


% Parameters for choosing movement
a = .5;     % temperary ratio
pStay = a;
pCurr = (1-pStay)*a;
pN = (1 - pStay - pCurr) / 8;

% Parameters for choosing death
deathParam = 800;


%%
coord = zeros(1,2,n_seeds);                                                                              % coordinates for each point for all paths generated in simulation
for seed_loop = 1:n_seeds
% seed_loop
    seed = seed_ind(round(rand*length(seed_ind)));                                                       % Choose a random point index to start the simulation
    [seed(1,1),seed(1,2)] = ind2sub(size(seed_ROI),seed);                                               % Coordinate of this starting point
    p(seed_loop).coord1{1} = seed;                                                                      % Starting location stored in data matrix
    p(seed_loop).coord2{1} = seed;
    ind = sub2ind(size(WM),seed(1,1),seed(1,2));                                                      % Index of this point in the white matter image
    % If coherency map data is present then use the angle indicated by the
    % data, else use a random angle
    if coh_map(ind)
        p(seed_loop).angle1(1) = -ori_map(ind);                                                         % Starting orientation
        p(seed_loop).angle2(1) = -ori_map(ind);
    else
        p(seed_loop).angle1(1) = -90 + 2*rand*90;                                                       % Random angle generated if no coherency data
        p(seed_loop).angle2(1) = -90 + 2*rand*90;
    end
    
    path = zeros(size(WM));
    flag1 = 0;  flag2 = 0;                                                                              % Flag indicating whether the path should continue 
    
    % Determine how long the cell will live
    stepsMax = poissinv(rand(1),deathParam);

    
    for i = 2:20000                                                                                    % Run the loop a large number of steps

        if ~flag1                                                                                       % Positive direction path
            
            seedx = round(p(seed_loop).coord1{i-1}(1)); 
            seedy = round(p(seed_loop).coord1{i-1}(2)); 
            
            seednow = [round(p(seed_loop).coord1{i-1}(1)), round(p(seed_loop).coord1{i-1}(2))]; 
            ind1 = sub2ind(size(WM),seednow(1),seednow(2));                                             % Index of previous coordinate point
            
            pRand = rand(1);

            if coh_map(ind1)                                                                            % Get Eigen vector angle at the previous coordinate point

                
                % Choose from Intervals
                %[0,p0][p0,1-8pn][1-8pn,1-7pn]
                if pRand <= pStay
                   
                elseif pRand <= 1- 8*pN     % Choose curr
                    % Keep the previous eigenvalue
                    
                elseif pRand <= 1-7*pN      % Choose top-left neighbor 
                    seedx= seedx -1;
                    seedy= seedy +1;

                elseif pRand <= 1-6*pN      % Choose above neighbor
                    seedy= seedy +1;
                
                elseif pRand <= 1-5*pN      % etc...
                    seedx= seedx +1;
                    seedy= seedy +1;
                
                elseif pRand <= 1-4*pN
                    seedx= seedx +1;
                    
                elseif pRand <= 1-3*pN
                    seedx= seedx +1;
                    seedy= seedy -1;
                    
                elseif pRand <= 1-2*pN
                    seedy= seedy -1;
                    
                elseif pRand <= 1-pN
                    seedx= seedx -1;
                    seedy= seedy -1;  
                    
                else % Last option
                    seedx= seedx -1;
                    
                end
                ind1 = sub2ind(size(WM),seedx,seedy);         % Index of previous coordinate point

                angle = -ori_map(ind1);
                p(seed_loop).angle1(i) = angle;                
                
            else
                angle_sampled = -90+2*rand*90;
                angle = rem(angle_sampled + p(seed_loop).angle1(i-1),360);                              % Real value of angle if sampled angle + original angle > 360 deg
                p(seed_loop).angle1(i) = angle_sampled;
            end
            
            p(seed_loop).coherency1(i) = coh_map(ind1);
            eigen_vect = [sind(angle),cosd(angle)];                                                     % Calculate the eigen vector (EV) direction
            
            sensitivityMag = 10;    % Grad variable
                        
            if pRand <= pStay
               p(seed_loop).coord1{i} = p(seed_loop).coord1{i-1};                           % Find next coordinate point at step d along EV direction
            else
               p(seed_loop).coord1{i} = p(seed_loop).coord1{i-1} - d*eigen_vect;                           % Find next coordinate point at step d along EV direction
            end
            
            p(seed_loop).coord1{i}(1) = p(seed_loop).coord1{i}(1) + sensitivityMag * cgradY(seednow(1,1),seednow(1,2));
            p(seed_loop).coord1{i}(2) = p(seed_loop).coord1{i}(2) + sensitivityMag * cgradX(seednow(1,1),seednow(1,2));

% % the path is out of bounds or hit max steps then don't proceed and change flag
%             
% % ------------------------------------------                
            if all(p(seed_loop).coord1{i} - size(WM)<=-1)  && all(round(p(seed_loop).coord1{i})>1)  && ~path(round(p(seed_loop).coord1{i}(1)),round(p(seed_loop).coord1{i}(2)))                  
    
                if i < stepsMax
                xpts = linspace(p(seed_loop).coord1{i-1}(1),p(seed_loop).coord1{i}(1),d);
                ypts = linspace(p(seed_loop).coord1{i-1}(2),p(seed_loop).coord1{i}(2),d);

                else
                    flag1 = 1;
                    
                end
                if mod(i,timeIntervals) == 0
                        x_intervals(seed_loop, i/timeIntervals) = p(seed_loop).coord1{i}(1);
                        y_intervals(seed_loop, i/timeIntervals) = p(seed_loop).coord1{i}(2);
                end   
            else
                flag1 = 1;
                x_intervals(seed_loop, find(x_intervals(seed_loop,:) == 0, 1):end) = p(seed_loop).coord1{i}(1);
                y_intervals(seed_loop, find(y_intervals(seed_loop,:) == 0, 1):end) = p(seed_loop).coord1{i}(2);
            end
% % ------------------------------------------                
        end
        
        
        
        if ~flag2
                        
            seedx = round(p(seed_loop).coord2{i-1}(1)); 
            seedy = round(p(seed_loop).coord2{i-1}(2)); 
            
            seednow = [round(p(seed_loop).coord2{i-1}(1)), round(p(seed_loop).coord2{i-1}(2))]; 
            ind2 = sub2ind(size(WM),round(p(seed_loop).coord2{i-1}(1)),round(p(seed_loop).coord2{i-1}(2)));
            
            pRand = rand(1);

            if coh_map(ind2)                                                                            % Get Eigen vector angle at the previous coordinate point
                
                % Choose from Intervals
                %[0,p0][p0,1-8pn][1-8pn,1-7pn]
                if pRand <= pStay
                   
                elseif pRand <= 1- 8*pN     % Choose curr
                    % Keep the previous eigenvalue
                    
                elseif pRand <= 1-7*pN      % Choose top-left neighbor 
                    seedx= seedx -1;
                    seedy= seedy +1;

                elseif pRand <= 1-6*pN      % Choose above neighbor
                    seedy= seedy +1;
                
                elseif pRand <= 1-5*pN      % etc...
                    seedx= seedx +1;
                    seedy= seedy +1;
                
                elseif pRand <= 1-4*pN
                    seedx= seedx +1;
                    
                elseif pRand <= 1-3*pN
                    seedx= seedx +1;
                    seedy= seedy -1;
                    
                elseif pRand <= 1-2*pN
                    seedy= seedy -1;
                    
                elseif pRand <= 1-pN
                    seedx= seedx -1;
                    seedy= seedy -1;  
                    
                else % Last option
                    seedx= seedx -1;
                    
                end
                ind2 = sub2ind(size(WM),seedx,seedy);         % Index of previous coordinate point

                angle = -ori_map(ind2);
                p(seed_loop).angle2(i) = angle;                
                
            else
                angle_sampled = -90+2*rand*90;
                angle = rem(angle_sampled + p(seed_loop).angle2(i-1),360);                              % Real value of angle if sampled angle + original angle > 360 deg
                p(seed_loop).angle2(i) = angle_sampled;
            end
            
            p(seed_loop).coherency2(i) = coh_map(ind2);
            eigen_vect = [sind(angle),cosd(angle)];                                                     % Calculate the eigen vector (EV) direction
            
            
            if pRand <= pStay
               p(seed_loop).coord2{i} = p(seed_loop).coord2{i-1};                           % Find next coordinate point at step d along EV direction
            else
               p(seed_loop).coord2{i} = p(seed_loop).coord2{i-1} + d*eigen_vect;                           % Find next coordinate point at step d along EV direction
            end
            
            
% ------------------------------------------                
            if all(p(seed_loop).coord2{i}-size(WM)<=-1)  && all(round(p(seed_loop).coord2{i})>1)  && ~path(round(p(seed_loop).coord2{i}(1)),round(p(seed_loop).coord2{i}(2)))                  % the path is out of bounds or path a loop has formed
                if i < stepsMax
                    xpts = linspace(p(seed_loop).coord2{i-1}(1),p(seed_loop).coord2{i}(1),d);
                    ypts = linspace(p(seed_loop).coord2{i-1}(2),p(seed_loop).coord2{i}(2),d);
                    
                else
                    flag2 = 1;
                end
                if mod(i,timeIntervals) == 0
                    x_intervals(seed_loop + n_seeds, i/timeIntervals) = p(seed_loop).coord2{i}(1);
                    y_intervals(seed_loop + n_seeds, i/timeIntervals) = p(seed_loop).coord2{i}(2);
                end
            else
                flag2 = 1;
                x_intervals(seed_loop + n_seeds, find(x_intervals(seed_loop + n_seeds,:) == 0, 1):end) = p(seed_loop).coord2{i}(1);
                y_intervals(seed_loop + n_seeds, find(y_intervals(seed_loop + n_seeds,:) == 0, 1):end) = p(seed_loop).coord2{i}(2);
            end
% ------------------------------------------                
        end    
    end
    
    
%    % Initialize the path indices on the path image with each path having
%    a different index. This is just another way for plotting on same plot.
%    WM_path(logical(path)) = .1*seed_loop;
end

%% Plotting the paths
% figure; imagesc(WM_path); axis equal;
%save 2D_5vox_blur_100k.mat
%edit_indices = [109,93,303,416,133,112,378];


% calculate distances and plot
distance = sqrt((inj_center(1) - x_intervals).^2 + (inj_center(2) - y_intervals).^2);
mean(distance)
figure()
boxplot(distance)

% determine whether NSC made it within a certain radius to cancer center
radius = 1000;
numAtCancer = 0;
numAtCancerGraph = zeros(num_intervals, 1);
xAxis = zeros(num_intervals, 1);

for i = 1:num_intervals
    for j = 1:2*n_seeds
        if (sqrt((cancer_center(1) - x_intervals(j, i)).^2 + (cancer_center(2) - y_intervals(j, i)).^2) < radius)
            numAtCancer = numAtCancer + 1;
        end
    end
    numAtCancerGraph(i, 1) = numAtCancer/(2*n_seeds)*100;
    xAxis(i, 1) = i;
    numAtCancer = 0;
end

figure()
plot(xAxis,numAtCancerGraph,'r')
xlabel('Time Intervals')
ylabel('Percent of NSC that Reach Cancer Site')

% get the coordinate data from each coord1 and coord2 data and plot them.
X_coords = zeros(1,n_seeds);
Y_coords = X_coords;
for i = 1:n_seeds
    coord_data1 = cell2mat(p(i).coord1);
    l1 = length(coord_data1);
    X_coords(1:l1/2,i) = coord_data1(1:2:l1)';    
    Y_coords(1:l1/2,i) = coord_data1(2:2:l1)'; 
    coord_data2 = cell2mat(p(i).coord2);
    l2 = length(coord_data2);
    X_coords(l1/2+1:(l1+l2)/2,i) = coord_data2(1:2:l2)';    
    Y_coords(l1/2+1:(l1+l2)/2,i) = coord_data2(2:2:l2)'; 
end
figure; imagesc(coh_map);   colormap gray;  hold on;
plot(Y_coords,X_coords,'.'); axis equal
for i = 1:num_intervals
    figure; imagesc(coh_map);   colormap gray;  hold on;
    plot(y_intervals(1:end,i),x_intervals(1:end,i),'.'); axis equal
end
% saveas(gcf,strcat('NSCMigration_pStay(',num2str(pStay),')_pCurr(',num2str(pCurr),')_deathParam(',int2str(deathParam),').fig'))


