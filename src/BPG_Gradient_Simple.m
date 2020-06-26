% BIDIRECTOINAL PATH GENERATION

%%_________________
%Using 5 voxel kernerl blur orientation and coherency map 

% Info about 'NSC_atlas_Vikram/WM_M43.tif' 
% conversion 1000 pixel [1444 ?m] radius
% Total size 12.840 mm 7.226 mm [8892  5004 pixel]

% LM-NSC008 cells
% Median distances at 3, 6, and 9 months post-injection were 437 ?m, 1030 ?m, and 902 ?m, respectively. 
% 66, 72, and 59 percent of NSCs were found within the WM (Fig 3D).
% median distances between NSCs and the nearest WM/GM interface were 263 ?m, 118 ?m, and 87 ?m. 


%% load data 
% clear all
% directory = pwd;
if( ~exist('coh_map','var') ) 
    WM = imread(strcat('NSC_atlas_Vikram/WM_M43.tif'));
    WM = logical(WM.*(WM>10));                                                                              % Threshold image to remove background
    ori_map = 180*imread(strcat('NSC_atlas_Vikram/Orientation_5vox_blur.tif'))/pi;   % -90 to 90                                    % Orientation map converted to degrees (-90 to 90 deg)
    temp = zeros(size(ori_map));
    temp(ori_map>0) = ori_map(ori_map>0);                                                                   % positive angles are ascribed to 0-90 degrees
    temp(ori_map<0) = 180+ori_map(ori_map<0);                                                               % Negative angles are ascribed to 90 - 180 deg
    ori_map = temp; clear temp;
    coh_map = imread(strcat('NSC_atlas_Vikram/Coherency_5vox_blur.tif')).*WM;
    
    clear WM 
end

%% set parameters 
set_parameters_Simple()                                                                                      % Coherency limit for terminating path
n_seeds = 1000;                                                                                                % Total number of paths generated
Finaltimestep = 5000; 

%% Initialize a matrix of same size as the white matter image and initialize a seed point (inj center) around which the seeds for each simulation are placed
[inj_center, seed_ind] = set_initial(coh_map); 


% Storing the coorinate, angle of eigen vector and coherency value for both
% positive and negative direction paths that start from a seed point
p = struct('coord',{},'angle',[],'coherency',[]);             % the three fields 


%% cancer injection
cancer_center = [1500, 3500]; 
cancer_sd = 1000; 
[concentration, cgradX, cgradY] = set_cancer( cancer_center, cancer_sd, coh_map ); 



%% 
Tstart = tic; 

%% for each seed 
for seed_loop = 1:n_seeds

    %%%% initialize first two steps 
    seed = seed_ind(round(rand*length(seed_ind)));                                                       % Choose a random point index to start the simulation
    [seed(1,1),seed(1,2)] = ind2sub(size(coh_map),seed);                                               % Coordinate of this starting point
    p(seed_loop).coord(:,1) = seed';                                                                      % Starting location stored in data matrix
    ind = sub2ind(size(coh_map),seed(1,1),seed(1,2));                                                      % Index of this point in the white matter image
    % If coherency map data is present then use the angle indicated by the
    % data, else use a random angle
    p(seed_loop).angle(1) = -ori_map(ind);                                                         % Starting orientation

    %%%% second step 
    i = 2; 
    ind1 = sub2ind(size(coh_map),round(p(seed_loop).coord(1,i-1)),round(p(seed_loop).coord(2,i-1)));                                             % Index of previous coordinate point
            
    angle = -ori_map(ind1);
    p(seed_loop).angle(i) = angle;            
    p(seed_loop).coherency(i) = coh_map(ind1);
    
    eigen_vect = [sind(angle),cosd(angle)];                                                     % Calculate the eigen vector (EV) direction
    dstep = (rand(1)*2-1)*dmax; %* double(coh_map(ind1));     
    p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*eigen_vect';                           % Find next coordinate point at step d along EV direction
    
%     path = zeros(size(coh_map));

%stochasticity in IC of NSC 
    dmax  = 5; 
    dstochastic = 5; 
    beta4dist = 2; 
    dmax = betainv( rand(1), 1, beta4dist )*dmax; % *(1+beta4dist); 
%     dstochastic = dstochastic*rand(1); 
    
    flag1 = 0; 
    for i = 3:Finaltimestep                                                                                    % Run the loop a large number of steps            
        if ~flag1  %% while NSC stays in the 
            seedx = round(p(seed_loop).coord(1,i-1)); 
            seedy = round(p(seed_loop).coord(2,i-1)); 
            
            ind1 = sub2ind(size(coh_map),seedx,seedy);                                             % Index of previous coordinate point
            
            if coh_map(ind1) %> coh_limit                                                                           % Get Eigen vector angle at the previous coordinate point
                angle = -ori_map(ind1);
                d = dmax; 
                
            else
                angle_sampled = -90+2*rand*90;
%                 angle_sampled = -45+2*rand*45;
                angle = rem( angle_sampled + p(seed_loop).angle(i-1) ,360);                              % Real value of angle if sampled angle + original angle > 360 deg
                % move with less magnitude if no WM present 
                
%                 d = dmax/dstochastic; 
                d = dstochastic*rand(1); 
                
            end 
            p(seed_loop).angle(i) = angle; 
            
            p(seed_loop).coherency(i) = coh_map(ind1);
            eigen_vect = [sind(angle),cosd(angle)];                                                     % Calculate the eigen vector (EV) direction
            
            %%%%%% deteministic 
            dstep = d; 
            %%%%%% uniform[0, 1] 
%             dstep = rand(1)*2*d; 
            %%%%%% beta[1, bb], if bb==1, uniform 
%             beta4dist = 5; 
%             dstep = betainv( rand(1), 1, beta4dist )*d *(1+beta4dist) ; 
            
           
            %% model 1 
            % Find next coordinate point at step d along EV direction 
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) ... 
                - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * dstep*eigen_vect';                           % Find next coordinate point at step d along EV direction
            % let it move away from the direction it came from. 
            
            %% model 2 
            % add chemotaxis with (1-w) weight
%             w = betainv( rand(1), beta4chmtx, 1 ) ; 
% 
%             eigen_vect = - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * eigen_vect'; 
%             chmtx_vect = [cgradY(seedx,seedy); cgradX(seedx,seedy)]; 
%             p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) ... 
%                 + dstep*( w*eigen_vect + (1-w)*chmtx_vect ); 
            
            %% model 3 
            %  move along chemotaxis once it is near by cancer 
%             if( concentration(seedx,seedy) <= chmtx_limit )
%                 
%             p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) ... 
%                 - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * dstep*eigen_vect';   
% 
%             else
%             
%             p(seed_loop).coord(1,i) = p(seed_loop).coord(1,i) + dstep * cgradY(seedx,seedy);
%             p(seed_loop).coord(2,i) = p(seed_loop).coord(2,i) + dstep * cgradX(seedx,seedy);
% 
%             end 
            

% % the path is out of bounds or hit max steps then don't proceed and change flag
%             
% % ------------------------------------------                
            if all(p(seed_loop).coord(:,i)' - size(coh_map)<=-1)  && all(round(p(seed_loop).coord(:,i)')>1) % && ~path(round(p(seed_loop).coord(1,i)),round(p(seed_loop).coord(2,i)))
% %                 % to interpolate the path to integers 
% %                 xpts = round( linspace(p(seed_loop).coord(1,i-1),p(seed_loop).coord(1,i),10) );
% %                 ypts = round( linspace(p(seed_loop).coord(2,i-1),p(seed_loop).coord(2,i),10) );                
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
xfinal = zeros(1,n_seeds); 
yfinal = zeros(1,n_seeds);
figure; imagesc(coh_map);   colormap gray;  axis equal; hold on;
for i = 1:n_seeds
    
    coord_data1 = p(i).coord; 
    l1 = size(coord_data1,2);
    if( l1 <= 5000 )
    X_coords = coord_data1(1, :);   
    Y_coords = coord_data1(2, :); 
    else
    X_coords = coord_data1(1,1:ceil(l1/5000):l1)';    
    Y_coords = coord_data1(2,1:ceil(l1/5000):l1)'; 
    end         
    plot(Y_coords,X_coords,'.', 'markersize', 2);  
    
    xfinal(i) = coord_data1(1, end); 
    yfinal(i) = coord_data1(2, end); 
end
% %%% to plot final location 
plot(yfinal, xfinal,'ow', 'markersize', 4, 'markerfacecolor', 'w'); 
saveas( gcf, strcat( '200626_model1_trajectoryAll_d', int2str(dmax), '_', num2str(dmax/dstochastic), '.jpg') ); 

nprint = 1; 
if( nprint ) 

    % %%% to plot distribution each time  
    for i = 1:n_seeds 
        if( size( p(i).coord, 2 ) < Finaltimestep ) % if stopped before Finaltimestep
            ind = size( p(i).coord, 2 ); 
            p(i).coord(1,(ind+1):Finaltimestep) = p(i).coord(1,ind); 
            p(i).coord(2,(ind+1):Finaltimestep) = p(i).coord(2,ind);         
        end
        coord(:,:,i) = p(i).coord';         
    end
    
    Tstep = floor( Finaltimestep/4 ); 
    figure; 
        subplot( 1, 5, 1 ); imagesc(coh_map);   colormap gray;  hold on; title( 'time step 0' )
        plot(squeeze( coord(1,2,:)), squeeze( coord(1,1,:)),'ow', 'markersize', 4, 'markerfacecolor', 'w');     
    for n = 1:4 
        subplot( 1, 5, n+1 ); imagesc(coh_map);   colormap gray;  hold on; title( int2str(Tstep*n) )  
        plot(squeeze( coord(Tstep*n,2,:)), squeeze( coord(Tstep*n,1,:)),'ow', 'markersize', 4, 'markerfacecolor', 'w'); 
    end 
    saveas( gcf, strcat( '200626_model1_trajectory_d', int2str(dmax), '_', num2str(dmax/dstochastic), '.jpg') ); 
    % plot density only when enough n_seeds 
%     if( n_seeds > 100 ) 
%         figure; 
%         for n = 1:5 
%             subplot( 1, 5, n ); 
%             [X,Y] = meshgrid(1:(size(coh_map,2)/150):size(coh_map,2),1:(size(coh_map,1)/150):size(coh_map,1));
%             h = ksdensity( [squeeze( coord(Tstep*n,2,:)), squeeze( coord(Tstep*n,1,:))], [X(:),Y(:)] ); 
%             surfo( X, Y, reshape( h, size(X,1), size(X,2) ) ); 
%             set(gca, 'YDir','reverse')
%         end 
%     end 
    
end 

for n = 1:length(p) 
%%%% distance from the center of injection site 
    distInit(n,:) = sqrt( sum ( (p(n).coord - inj_center'*ones(1,size(p(n).coord,2))).^2, 1 ) );  
%%%% distance from its Initial injection location 
%     distInit(n,:) = sqrt( sum ( (p(n).coord - p(n).coord(:,1)*ones(1,size(p(n).coord,2))).^2, 1 ) );  
end 
figure; boxplot( distInit(:, [1:1000:4001, 5000])*1.444 ); ylabel( 'Distance from injection site (\mu m)' ); 
hold on; plot( [(8+84/(5000)):(84/(5000)):92], median( distInit )*1.444 )
hold on; plot( [(8+84/(5000)):(84/(5000)):92], mean( distInit )*1.444 )
set(gca,'XTick', [(8):(84/(5)):92], 'XTicklabel',{'step 0', '1000',  '2000',  '3000',  '4000',  '5000'});
saveas( gcf, strcat( '200626_model1_distboxplot_d', int2str(dmax), '_', num2str(dmax/dstochastic), '.jpg') ); 

disp( 'done' ); 


function [concentration, cgradX, cgradY] = set_cancer( cancer_center, cancer_sd, coh_map ) 
    [X,Y] = meshgrid(1:size(coh_map,2),1:size(coh_map,1));

    concentration = exp((-(X - cancer_center(2)).^2 - (Y - cancer_center(1)).^2)./cancer_sd^2);

    cgradY = zeros(size(X));
    cgradY(2:end-1,:) = concentration(3:end,:) - concentration(1:end-2,:);

    cgradX = zeros(size(X));
    cgradX(:,2:end-1) = concentration(:,3:end) - concentration(:,1:end-2);

    % find largest vector magnitude
    L = (cgradX.^2 + cgradY.^2).^(.5);
    maxL = max(max(L));

    % normalize gradient
    cgradX = cgradX./maxL;
    cgradY = cgradY./maxL;

end 

function [inj_center, seed_ind] = set_initial(coh_map) 

[X,Y] = meshgrid(1:size(coh_map,2),1:size(coh_map,1));
inj_center = [2150,1000];
%inj_center = [5100,730];

seed_sd = 250; 
seed_ROI = sqrt((X-inj_center(2)).^2 + (Y-inj_center(1)).^2)<seed_sd;                                       % All points within 500 pixel distance from inj_center
seed_ind = find(seed_ROI);  


end 
