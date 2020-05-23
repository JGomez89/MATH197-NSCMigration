% BIDIRECTOINAL PATH GENERATION

%%_________________
%Using 5 voxel kernerl blur orientation and coherency map 

%% load data 
% clear all
% directory = pwd;
load( 'NSC_atlas_Vikram/180slice' );

eigen_map = EV_slice; 
coh_map   = FA_slice; % fractional anisotropy 
ori_map = atand(squeeze( EV_slice(:,:,2) )./squeeze( EV_slice(:,:,1) )); 
clear EV_slice FA_slice 

%% parameters 
set_parameters_Simple()                                                                                      % Coherency limit for terminating path

%% Initialize a matrix of same size as the white matter image and initialize a seed point (inj center) around which the seeds for each simulation are placed
[X,Y] = meshgrid(1:size(coh_map,2),1:size(coh_map,1));
inj_center = [400,270]; % in corpus callosum
seed_ROI = sqrt((X-inj_center(2)).^2 + (Y-inj_center(1)).^2)<10;                                       % All points within 500 pixel distance from inj_center
seed_ind = find(seed_ROI);                                                                              % Indices of all these points
n_seeds = 10;                                                                                          % Total number of paths generated

% Storing the coorinate, angle of eigen vector and coherency value for both
% positive and negative direction paths that start from a seed point
p = struct('coord1',{},'angle1',[],'coherency1',[],'coord2',{},'angle2',[],'coherency2',[]);             % the three fields 



%% cancer injection
cancer_center = [160,280]; % near contralateral corpus callosum 
cancer_sd = 20;
[concentration, cgradX, cgradY] = set_cancer( cancer_center, cancer_sd, X, Y ); 

chmtx_limit = 0.1; 


%% 
Tstart = tic; 
Finaltimestep = 10000; 

%% for each seed 
for seed_loop = 1:n_seeds

    %%%% initialize first two steps 
    seed = seed_ind(round(rand*length(seed_ind)));                                                       % Choose a random point index to start the simulation
    [seed(1,1),seed(1,2)] = ind2sub(size(seed_ROI),seed);                                               % Coordinate of this starting point
    p(seed_loop).coord1(:,1) = seed';                                                                      % Starting location stored in data matrix
    ind = sub2ind(size(coh_map),seed(1,1),seed(1,2));                                                      % Index of this point in the white matter image
    % If coherency map data is present then use the angle indicated by the
    % data, else use a random angle
    p(seed_loop).angle1(1) = ori_map(seed(1,1),seed(1,2)); 

    %%%% second step 
    i = 2; 
%     ind1 = sub2ind(size(WM),round(p(seed_loop).coord1(1,i-1)),round(p(seed_loop).coord1(2,i-1)));                                             % Index of previous coordinate point
            
%     angle = -ori_map(ind1);
%     p(seed_loop).angle1(i) = angle;            
%     p(seed_loop).coherency1(i) = coh_map(ind1);
    
    eigen_vect = squeeze( eigen_map( round(p(seed_loop).coord1(1,i-1)),round(p(seed_loop).coord1(2,i-1)), : )); 
    dstep = (rand(1)*2-1)*dmax; %* double(coh_map(ind1));     
    p(seed_loop).coord1(:,i) = p(seed_loop).coord1(:,i-1) + dstep*eigen_vect;                           % Find next coordinate point at step d along EV direction
    
%     path = zeros(size(WM));
    
    flag1 = 0; 
    for i = 3:Finaltimestep                                                                                    % Run the loop a large number of steps            
        if ~flag1  
            seedx = round(p(seed_loop).coord1(1,i-1)); 
            seedy = round(p(seed_loop).coord1(2,i-1)); 
            
            ind1 = sub2ind(size(coh_map),seedx,seedy);                                             % Index of previous coordinate point
            
            if coh_map(ind1) > coh_limit                                                                           % Get Eigen vector angle at the previous coordinate point
%                 angle = -ori_map(ind1);
%                 p(seed_loop).angle1(i) = angle; 
                eigen_vect = squeeze( eigen_map( seedx, seedy, : ) ); 
                d = dmax; 
                
            else
% %                 angle_sampled = -90+2*rand*90;
                angle = -45+2*rand*45;
                if ~isnan(ori_map(seedx, seedy))
                    angle = rem( angle + ori_map(seedx, seedy) ,360);     
                end
                p(seed_loop).angle1(i) = angle;
                eigen_vect = [cosd(angle), sind(angle)]'; 

                % move with less magnitude if no WM present 
                d = 0; %dmax/100; 
                
            end 
            
            p(seed_loop).coherency1(i) = coh_map(ind1);
            
            
            %%%%%% uniform[0, 1] 
%             dstep = rand(1)*dmax ; 
            %%%%%% beta[1, beta], if beta==1, uniform 
            dstep = betainv( rand(1), 1, beta4dist )*d ; 
            dstep = d; 
            
            
           
            %% model 1 
            % Find next coordinate point at step d along EV direction 
            p(seed_loop).coord1(:,i) = p(seed_loop).coord1(:,i-1) ... 
                - sign(dot( p(seed_loop).coord1(:,i-2) - p(seed_loop).coord1(:,i-1), eigen_vect ) + eps ) * dstep*eigen_vect;                           % Find next coordinate point at step d along EV direction
            % let it move away from the direction it came from. 
            
            p(seed_loop).coord1(:,i) =   p(seed_loop).coord1(:,i) + 0.1*randn(2,1); 
            
            %% model 2 
%             % add chemotaxis with (1-w) weight
%             w = betainv( rand(1), beta4chmtx, 1 ) ; 
% 
%             eigen_vect = - sign(dot( p(seed_loop).coord1(:,i-2) - p(seed_loop).coord1(:,i-1), eigen_vect ) + eps ) * eigen_vect; 
%             chmtx_vect = [cgradY(seedx,seedy); cgradX(seedx,seedy)]; 
%             p(seed_loop).coord1(:,i) = p(seed_loop).coord1(:,i-1) ... 
%                 + dstep*( w*eigen_vect + (1-w)*chmtx_vect );
            
            %% model 3 
            %  move along chemotaxis once it is near by cancer 
%             if( concentration(seedx,seedy) <= chmtx_limit )
%                 
%             p(seed_loop).coord1(:,i) = p(seed_loop).coord1(:,i-1) ... 
%                 - sign(dot( p(seed_loop).coord1(:,i-2) - p(seed_loop).coord1(:,i-1), eigen_vect ) + eps ) * dstep*eigen_vect;   
% 
%             else
%             
%             p(seed_loop).coord1(1,i) = p(seed_loop).coord1(1,i) + dstep * cgradY(seedx,seedy);
%             p(seed_loop).coord1(2,i) = p(seed_loop).coord1(2,i) + dstep * cgradX(seedx,seedy);
% 
%             end 
            
            if( mod( i, 100 ) == 0 )
                i;
            end

% % the path is out of bounds or hit max steps then don't proceed and change flag
%             
% % ------------------------------------------                
            seednow = round(p(seed_loop).coord1(:,i)); 
            if all( seednow' - size(coh_map)<=-1)  && all( seednow' >1) && ~isnan(ori_map(seednow(1),seednow(2)))
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
xfinal = zeros(1,n_seeds); 
yfinal = zeros(1,n_seeds);
figure; imagesc(coh_map);   colormap gray;  axis equal;  hold on;
for i = 1:n_seeds
    
    coord_data1 = p(i).coord1; 
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

nprint = 0; 
if( nprint ) 

    % %%% to plot distribution each time  
    for i = 1:n_seeds 
        if( size( p(i).coord1, 2 ) < Finaltimestep ) 
            ind = size( p(i).coord1, 2 ); 
            p(i).coord1(1,(ind+1):Finaltimestep) = p(i).coord1(1,ind); 
            p(i).coord1(2,(ind+1):Finaltimestep) = p(i).coord1(2,ind);         
        end
        coord(:,:,i) = p(i).coord1';         
    end
    Tstep = floor( Finaltimestep/5 ); 
    figure; 
    for n = 1:5 
        subplot( 1, 5, n ); imagesc(coh_map);   colormap gray;  hold on; 
        plot(squeeze( coord(Tstep*n,2,:)), squeeze( coord(Tstep*n,1,:)),'ow', 'markersize', 4, 'markerfacecolor', 'w'); 
    end 
    
    

end 

disp( 'done' ); 


function [concentration, cgradX, cgradY] = set_cancer( cancer_center, cancer_sd, X, Y ) 

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

