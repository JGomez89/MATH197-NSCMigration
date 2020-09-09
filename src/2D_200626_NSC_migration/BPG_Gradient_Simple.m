% BIDIRECTOINAL PATH GENERATION

%%_________________
% Using 5 voxel kernerl blur orientation and coherency map 

% Info about 'NSC_atlas_Vikram/WM_M43.tif' 
% conversion 1000 pixel [1444 ?m] radius
% Total size 12.840 mm 7.226 mm [8892  5004 pixel]

% LM-NSC008 cells
% Median distances at 3, 6, and 9 months post-injection were 437 ?m, 1030 ?m, and 902 ?m, respectively. 
% 66, 72, and 59 percent of NSCs were found within the WM (Fig 3D).
% median distances between NSCs and the nearest WM/GM interface were 263 ?m, 118 ?m, and 87 ?m. 

%% Load data 

disp('Gathering data...')

if( ~exist('coh_map','var') ) 
    WM = imread(strcat('NSC_atlas_Vikram/WM_M43.tif'));
    WM = logical(WM.*(WM>10));                                                                              % Threshold image to remove background
    ori_map = 180*imread(strcat('NSC_atlas_Vikram/Orientation_5vox_blur.tif'))/pi;                          % Orientation map converted to degrees (-90 to 90 deg)
    temp = zeros(size(ori_map));
    temp(ori_map>0) = ori_map(ori_map>0);                                                                   % positive angles are ascribed to 0-90 degrees
    temp(ori_map<0) = 180+ori_map(ori_map<0);                                                               % Negative angles are ascribed to 90 - 180 deg
    ori_map = temp; clear temp;
    coh_map = imread(strcat('NSC_atlas_Vikram/Coherency_5vox_blur.tif')).*WM;
    
    clear WM 
end


%% Set parameters         
set_parameters_Simple()

% Coherency limit for terminating path
n_seeds =           1000;                                                                                   % Total number of paths generated
Finaltimestep =     5000;
CONVERT2MICROM =    1.444;

%% Initialize a matrix of same size as the white matter image and initialize a seed point (inj center) around which the seeds for each simulation are placed
[inj_center, seed_ind] = set_initial(coh_map); 

% Storing the coorinate, angle of eigen vector and coherency value for both
% positive and negative direction paths that start from a seed point
p = struct('coord',{},'angle',[],'coherency',[]);             % the three fields 


%% Cancer injection

cancer_center = [1500, 3500];
% cancer_center = [3000,1500];
% cancer_center = [2000,2000];
% cancer_center = [4000,1000];
% cancer_center = [5000,4000];

cancer_sd = 1000; 
[concentration, cgradX, cgradY] = set_cancer( cancer_center, cancer_sd, coh_map ); 


%% 
Tstart = tic;
disp('Running simulation...')


%% For each seed 
for seed_loop = 1:n_seeds
    
    %%%% initialize max step length
    dmax = betainv(rand(1), 1, beta4dist); %(1+beta4dist);

    %%%% initialize first two steps 
    seed = seed_ind(round(rand*length(seed_ind)));                                                          % Choose a random point index to start the simulation
    [seed(1,1),seed(1,2)] = ind2sub(size(coh_map),seed);                                                    % Coordinate of this starting point
    p(seed_loop).coord(:,1) = seed';                                                                        % Starting location stored in data matrix
    ind = sub2ind(size(coh_map),seed(1,1),seed(1,2));                                                       % Index of this point in the white matter image
    % If coherency map data is present then use the angle indicated by the
    % data, else use a random angle
    p(seed_loop).angle(1) = -ori_map(ind);                                                                  % Starting orientation

    %%%% second step 
    i = 2; 
    ind = sub2ind(size(coh_map),round(p(seed_loop).coord(1,i-1)),round(p(seed_loop).coord(2,i-1)));         % Index of previous coordinate point
            
    angle = -ori_map(ind);
    p(seed_loop).angle(i) = angle;            
    p(seed_loop).coherency(i) = coh_map(ind);
    
    eigen_vect = [sind(angle),cosd(angle)];                                                                 % Calculate the eigen vector (EV) direction
    dstep = (rand(1)*2-1)*dmax; %* double(coh_map(ind1));     
    p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*eigen_vect';                                % Find next coordinate point at step d along EV direction
    
    
    for i = 3:Finaltimestep                                                                                 % Run the loop a large number of steps            
        
        %%%% While NSC stays in bounds
        
        seedx = round(p(seed_loop).coord(1,i-1));
        seedy = round(p(seed_loop).coord(2,i-1)); 
        ind = sub2ind(size(coh_map),seedx,seedy);                                                           % Index of previous coordinate point
        
        % If NSC has preferred direction
        if coh_map(ind) > coh_limit                                                                         % Get Eigen vector angle at the previous coordinate point
            angle = -ori_map(ind);
            dstep = dmax * d_w; 
        else
            % Choose rand angle and move with less magnitude if no WM present 
            angle_sampled = -90+2*rand*90;
            angle = rem(angle_sampled + p(seed_loop).angle(i-1) ,360);                                      % Real value of angle if sampled angle + original angle > 360 deg
            dstep = dmax * d_g;
        end 
        
        p(seed_loop).angle(i) = angle; 
        p(seed_loop).coherency(i) = coh_map(ind);
        eigen_vect = [sind(angle),cosd(angle)];                                                             % Calculate the eigen vector (EV) direction

        if modelNum == 1
            % Find next coordinate point at step d along EV direction 
            % Let it move away from the direction it came from. 
            % Find next coordinate point at step d along EV direction
            modelType = 'model1';
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) ... 
                - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * dstep*eigen_vect';            

        elseif modelNum == 2
            modelType = 'model2';
            eigen_vect = - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * eigen_vect';
            chmtx_vect = [cgradY(seedx,seedy); cgradX(seedx,seedy)];

            %  Move along chemotaxis once it is near by cancer
            if has_cancer && (concentration(seedx,seedy) >= chmtx_limit)
                chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
            else
                chmtx_bias = 0;
            end
            
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*eigen_vect + chmtx_bias*chmtx_vect;
            
        end

        % The path is out of bounds or hit max steps then don't proceed
        if ~(all(p(seed_loop).coord(:,i)' - size(coh_map)<=-1)  && all(round(p(seed_loop).coord(:,i)')>1))
            break
        end
          
        
    end 
end



coord = zeros(Finaltimestep,2,n_seeds);
for i = 1:n_seeds 
    if( size( p(i).coord, 2 ) < Finaltimestep ) % if stopped before Finaltimestep
        ind = size( p(i).coord, 2 ) - 1; 
        p(i).coord(1,(ind+1):Finaltimestep) = p(i).coord(1,ind); 
        p(i).coord(2,(ind+1):Finaltimestep) = p(i).coord(2,ind);         
    end
    coord(:,:,i) = p(i).coord';
end

toc( Tstart );


%% Plotting the paths
% figure; imagesc(WM_path); axis equal;
%save 2D_5vox_blur_100k.mat
%edit_indices = [109,93,303,416,133,112,378];

disp('Plotting graphs...')


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
%%% to plot final location 
plot(yfinal, xfinal,'ow', 'markersize', 4, 'markerfacecolor', 'w'); 
savethis('trajectoryAll');



%%% to plot distribution each time  
Tstep = floor( Finaltimestep/4 ); 
figure; 
subplot( 1, 5, 1 ); imagesc(coh_map);   colormap gray;  hold on; title( 'time step 0' )
plot(squeeze( coord(1,2,:)), squeeze( coord(1,1,:)),'ow', 'markersize', 4, 'markerfacecolor', 'w');     
for n = 1:4 
    subplot( 1, 5, n+1 ); imagesc(coh_map);   colormap gray;  hold on; title( int2str(Tstep*n) )  
    plot(squeeze( coord(Tstep*n,2,:)), squeeze( coord(Tstep*n,1,:)),'ow', 'markersize', 4, 'markerfacecolor', 'w'); 
end 
savethis('trajectory');



%%%% plot density only when enough n_seeds 
% if( n_seeds > 100 ) 
%     figure; 
%     for n = 1:5 
%         subplot( 1, 5, n ); 
%         [X,Y] = meshgrid(1:(size(coh_map,2)/150):size(coh_map,2),1:(size(coh_map,1)/150):size(coh_map,1));
%         h = ksdensity( [squeeze( coord(Tstep*n,2,:)), squeeze( coord(Tstep*n,1,:))], [X(:),Y(:)] ); 
%         surfo( X, Y, reshape( h, size(X,1), size(X,2) ) ); 
%         set(gca, 'YDir','reverse')
%     end 
% end 



%%% Boxplots of NSC's distance from injection site
distInit = zeros(n_seeds,Finaltimestep);
for n = 1:length(p) 
%%%% distance from the center of injection site 
    distInit(n,:) = sqrt( sum ( (p(n).coord - inj_center'*ones(1,size(p(n).coord,2))).^2, 1 ) );  
%%%% distance from its Initial injection location 
%     distInit(n,:) = sqrt( sum ( (p(n).coord - p(n).coord(:,1)*ones(1,size(p(n).coord,2))).^2, 1 ) );  
end

figure;     boxplot( distInit(:, [1:1000:4001, Finaltimestep])*CONVERT2MICROM ); ylabel( 'Distance from injection site (\mu m)' ); 
hold on;    plot( ((8+84/(Finaltimestep)):(84/(Finaltimestep)):92), median( distInit )*CONVERT2MICROM )
hold on;    plot( ((8+84/(Finaltimestep)):(84/(Finaltimestep)):92), mean( distInit )*CONVERT2MICROM )
% set(gca,'XTick', [(8):(84/(5)):92], 'XTicklabel',{'step 0', '1000',  '2000',  '3000',  '4000',  '5000'});
savethis('distboxplot');



%%% Determine whether NSC made it within a certain radius to cancer center
if has_cancer
    percAtCancerGraph = zeros(Finaltimestep);
    timeInterval = zeros(Finaltimestep);

    for i = 1:Finaltimestep
        numAtCancer = 0;
        for j = 1:n_seeds
            if (sqrt((cancer_center(1) - p(j).coord(1,i)).^2 + (cancer_center(2) - p(j).coord(2,i)).^2) < cancer_radius)
                numAtCancer = numAtCancer + 1;
            end
        end
        percAtCancerGraph(i) = numAtCancer/(n_seeds)*100;
        timeInterval(i) = i;
    end

    figure()
    plot(timeInterval,percAtCancerGraph,'r')
    xlabel('Time Intervals')
    ylabel('Percent of NSC that Reach Cancer Site')
    savethis('percentArrived');
end



%%% Plot percent of NSC on WM past 3000
acc = 100;
percAtWMGraph = zeros(acc);
timeInterval = zeros(acc);
for i = 1:acc
    k= round(i*Finaltimestep/acc);
    numAtWM = 0;
    for j = 1:n_seeds
        seedx2 = round(p(j).coord(1,k)); 
        seedy2 = round(p(j).coord(2,k));     
        if (coh_map(seedx2,seedy2) > coh_limit && seedy2 > 3000)
            numAtWM = numAtWM + 1;
        end
    end
    percAtWMGraph(i) = numAtWM/(n_seeds)*100;
    timeInterval(i) = k;
end

figure()
% Plot regular graph
plot(timeInterval,percAtWMGraph,'r')
hold on
% Plot deg2 interpolation
p = polyfit(timeInterval,percAtWMGraph,2);
p_ = polyval(p,timeInterval);
plot(timeInterval,p_,'--','Color','r','LineWidth',2)
hold off

xlabel('Time Step')
ylabel('Percent of NSC on WM passed 3000')
savethis('percentAtWM');


disp( 'done' )


%% Functions
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
    inj_center = [2150,1000];
    %inj_center = [5100,730];

    [X,Y] = meshgrid(1:size(coh_map,2),1:size(coh_map,1));

    seed_sd = 250; 
    seed_ROI = sqrt((X-inj_center(2)).^2 + (Y-inj_center(1)).^2)<seed_sd;                                   % All points within 500 pixel distance from inj_center
    seed_ind = find(seed_ROI);  

end 

function savethis(title)
    global modelType d_w d_g chemo_sensitivity alpha4chmtx beta4dist cancer_center FolderName1 FolderName2 has_cancer;
    
    if has_cancer
        saveas( gcf, [pwd strcat( FolderName2, '200626_',modelType,'_',title,'_d',num2str(d_w),'_', num2str(d_g),...
            '_a',int2str(alpha4chmtx),'_b',int2str(beta4dist),'_c',num2str(chemo_sensitivity),'_[',int2str(cancer_center(1)),',',int2str(cancer_center(2)),']', '.fig')] );
    else
        saveas( gcf, [pwd strcat( FolderName1, '200626_',modelType,'_',title,'_d',num2str(d_w),'_', num2str(d_g),...
            '_a',int2str(alpha4chmtx),'_b',int2str(beta4dist),'_c',num2str(chemo_sensitivity),'_[NA]', '.fig')] );
    end
end
