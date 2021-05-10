% BIDIRECTOINAL PATH GENERATION

%%_________________
%Using 5 voxel kernerl blur orientation and coherency map 

%% Load data 
% clear all;
close all;

disp('Gathering data...')
[eigen_map, coh_map] = load_3D(); 


%% Set parameters 
minmax_param()

% Storing the coorinate, angle of eigen vector and coherency value for both
% positive and negative direction paths that start from a seed point
p = struct('coord',{});                                                                                     % The two fields 

%% Initialize a matrix of same size as the white matter image and initialize a seed point (inj center) around which the seeds for each simulation are placed
[seed_ind,seed_sd,ubound,lbound] = set_initial(coh_map);                                                    % Nasal injection initial 


%% Cancer injection
[concentration, cgradX, cgradY, cgradZ, cancer_sd] = set_cancer(coh_map);                                   % Near contralateral corpus callosum


%% For each seed
Tstart = tic; 
disp('Running simulation...')


for seed_loop = 1:n_seeds
    
    %%%% Initialize particular reference step length
    dmax = betainv(rand(1), 1, beta4dist);
    
    %%%% Spawn seed within injection site
    i = 1;
    seed = seed_ind(ceil(rand*length(seed_ind)));                                                           % Choose a random point index to start the simulation
    [seed(1,1),seed(1,2),seed(1,3)] = ind2sub(size(coh_map),seed);
    if seed(3) > ubound(seed(1),seed(2))
        seed(3) = ubound(seed(1),seed(2));
    elseif seed(3) < lbound(seed(1),seed(2))
        seed(3) = lbound(seed(1),seed(2));
    end
    p(seed_loop).coord(:,i) = seed';                                                                        % Coordinate of this starting point

    %%%% First step 
    i = 2; 
    eigen_vect = squeeze( eigen_map( round(p(seed_loop).coord(1,i-1)),round(p(seed_loop).coord(2,i-1)),round(p(seed_loop).coord(3,i-1)), : ));
    dstep = (rand(1)*2-1)*dmax; 
    p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*eigen_vect;                                 % Find next coordinate at step d along EV direction
    
    %%%% Algorithm: while NSC stays in bounds
    for i = 3:Finaltimestep                                                                                 % Run the loop a large number of steps            
        
        seedx = round( p(seed_loop).coord(1,i-1) );                                                         % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
        ind = sub2ind(size(coh_map),seedx,seedy,seedz);

        if coh_map(ind) > coh_limit                                                                         % Get Eigen vector at the previous coordinate point
            eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) ); 
            d = dmax * d_w; 
        else  
            eigen_vect_noise = (rand(3,1)*2-1); 
            eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) ) + eigen_vect_noise/norm(eigen_vect_noise)*norm(eigen_vect); 
            d = dmax * d_g;                                                                                 % Move with less magnitude if no WM present 
        end 

        %%%%%% uniform[0, 1] 
%         dstep = rand(1)*dmax ; 
        %%%%%% beta[1, beta], if beta==1, uniform 
%         dstep = betainv( rand(1), 1, beta4dist )*d ; 
        %%%%%% deterministic
        dstep = d; 

        chmtx_vect = [cgradX(seedx,seedy,seedz); cgradY(seedx,seedy,seedz); cgradZ(seedx,seedy,seedz)];
        chmtx_bias = 0;
        switch modelNum
            
            case 1  % Chemotaxis
                modelType = 'model1';                                                                       % Move along chemotaxis (deterministic)
                dstep = dstep/dmax;                                                                         % Keeps dstep deterministic (dmax is stochastic)
                if has_cancer
                    chmtx_bias = chemo_sensitivity;
                end
                
            case 2  % Chemotaxis, Stochasticity 
                modelType = 'model2';                                                                       %  Move along chemotaxis (stochasticity)
                if has_cancer
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                end
                
            case 3  % Chemotaxis, Chemotax threshold, Stochasticity
                modelType = 'model3';                                                                       %  Move along chemotaxis once near cancer
                if has_cancer && (concentration(seedx,seedy,seedz) >= chmtx_limit)
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                end
                
        end
        
        % take step
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*eigen_vect + chmtx_bias*chmtx_vect;

        % if path is out of bounds then stop
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                               %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...             %Seed is outside of the lower and upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                                %Seed is outside of the side bounds of brain (x,y)
                break
        end


    end    
end

%Fill in rest of coordinates with last valid coordinate
for i = 1:n_seeds 
    if( size( p(i).coord, 2 ) < Finaltimestep )                                                             % If stopped before Finaltimestep
        ind = size( p(i).coord, 2 ) - 1; 
        p(i).coord(1,(ind+1):Finaltimestep) = p(i).coord(1,ind);
        p(i).coord(2,(ind+1):Finaltimestep) = p(i).coord(2,ind);
        p(i).coord(3,(ind+1):Finaltimestep) = p(i).coord(3,ind);
    end
end

toc( Tstart ); 


%% Plot Graphs 
disp('Plotting graphs...')


%% Minimum Path

stepDistance = zeros(n_seeds, Finaltimestep - 1);

for seed_loop = 1:n_seeds
    for i = 2:Finaltimestep
        % record distance traveled each step
        stepDistance(seed_loop, i - 1) = sqrt((p(seed_loop).coord(1,i) - p(seed_loop).coord(1,i-1)).^2 + (p(seed_loop).coord(2,i) - p(seed_loop).coord(2,i-1)).^2 + (p(seed_loop).coord(3,i) - p(seed_loop).coord(3,i-1)).^2);
    end
end

if has_cancer
    dist2cancer = zeros(n_seeds, 1); 
    
    cellind = [ ];
    
    for j = 1:n_seeds
        if (all( abs(cancer_center - p(j).coord(:,end)') < cancer_size ))
            cellind = [cellind, j];
        end
    end 

    % compute cell trajectory lengths 
    stepDistanceTotal = sum( stepDistance, 2 ); 

    %Min Calculations
    % find minimum among cells that are close to cancer 
    minPath = min( stepDistanceTotal(cellind) ); 
    minPath_seed = find( stepDistanceTotal == minPath ); 

    % find cells that are close enough to minimum path 
    mincloseindex = find( abs( stepDistanceTotal-minPath ) < 100 ); 

    % percentage of cells that are near minimum path 
    minclosepercent = length(mincloseindex) / length(stepDistanceTotal);
    
    %Max Calculations
    % find maximum among cells that are close to cancer 
    maxPath = max( stepDistanceTotal(cellind) ); 
    maxPath_seed = find( stepDistanceTotal == maxPath ); 

    % find cells that are close enough to maximum path 
    maxcloseindex = find( abs( stepDistanceTotal-maxPath ) < 100 ); 

    % percentage of cells that are near maximum path 
    maxclosepercent = length(maxcloseindex) / length(stepDistanceTotal);
end

figure; hold on;

%Plot Min Path
plot3(p(minPath_seed).coord(1,:),p(minPath_seed).coord(2,:),p(minPath_seed).coord(3,:),'linewidth',3)

%Plot Max Path
plot3(p(maxPath_seed).coord(1,:),p(maxPath_seed).coord(2,:),p(maxPath_seed).coord(3,:),'linewidth',3)

%Plot White Matter
surf( ubound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );                   %Plot upper bound of mouse brain
surf( lbound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );   

[x,y,z] = sphere;
x = x*seed_sd + inj_center(1);
y = y*seed_sd + inj_center(2);
z = z*seed_sd + inj_center(3);
h = surf(x, y, z);                                                                                          %Plot injection site
set(h,'FaceColor',[1 0 1],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
    
if has_cancer
    [x,y,z] = sphere;
    x = x*cancer_size(1) + cancer_center(1);
    y = y*cancer_size(2) + cancer_center(2);
    z = z*cancer_size(3) + cancer_center(3);
    h = surf(x, y, z);                                                                                      %Plot cancer site
    set(h,'FaceColor',[1 .7 0],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
end

ii = find( coh_map > .7 ); 
[Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
hold on; plot3( X(ii),Y(ii),Z(ii), 'bo');

xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; daspect([1 1 1]); camlight;

%% Plot path of seeds
figure; hold on;

finalCoords = zeros(3,n_seeds);
for i=1: n_seeds
    finalCoords(:,i) = p(i).coord(:, end);
end

%Plot path of seed
for i=1: n_seeds
    plot3(p(i).coord(1,:),p(i).coord(2,:),p(i).coord(3,:),'linewidth',3)  
end

%plot3(p(minPath_seed).coord(1,:),p(minPath_seed).coord(2,:),p(minPath_seed).coord(3,:),'linewidth',3)


plot3(finalCoords(1,:),finalCoords(2,:),finalCoords(3,:),'ok','Markersize', 7,'markerfacecolor','w');       %Plot final coords of each seed
surf( ubound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );                   %Plot upper bound of mouse brain
surf( lbound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );                   %Plot lower bound of mouse brain

% [x,y,z] = sphere;
% x = x*seed_sd + inj_center(1);
% y = y*seed_sd + inj_center(2);
% z = z*seed_sd + inj_center(3);
% h = surf(x, y, z);                                                                                          %Plot injection site
% set(h,'FaceColor',[1 0 1],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
%     
% if has_cancer
%     [x,y,z] = sphere;
%     x = x*cancer_size(1) + cancer_center(1);
%     y = y*cancer_size(2) + cancer_center(2);
%     z = z*cancer_size(3) + cancer_center(3);
%     h = surf(x, y, z);                                                                                      %Plot cancer site
%     set(h,'FaceColor',[1 .7 0],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
% end

ii = find( coh_map > .7 ); 
[Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
hold on; plot3( X(ii),Y(ii),Z(ii),'bo');

xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; daspect([1 1 1]); camlight;
savethis('trajectoryAll');


%% Boxplots of NSC's distance from center of injection site
clear distInit;
distInit = zeros(n_seeds,Finaltimestep);
for n = 1:n_seeds
    distInit(n,:) = sqrt( sum ( (p(n).coord - inj_center'*ones(1,size(p(n).coord,2))).^2, 1 ) );   
end

acc = 10;
figure;     boxplot( distInit(:, (1:acc)*floor(Finaltimestep/acc))*CONVERT2MICRON );
hold on;    plot( (1:Finaltimestep)*(acc/Finaltimestep), median( distInit )*CONVERT2MICRON );
%hold on;    plot( (1:Finaltimestep)*(acc/Finaltimestep), mean( distInit )*CONVERT2MICRON );

xlim([0,7])
ylim([0,4000])
ylabel( 'Distance from injection site (\mu)' );
savethis('distboxplot');




%% Percent of NSC on WM
percOnWM = zeros(Finaltimestep,1);
for i = 1:Finaltimestep
    numAtWM = 0;
    for j = 1:n_seeds
        seedx = round(p(j).coord(1,i)); 
        seedy = round(p(j).coord(2,i)); 
        seedz = round(p(j).coord(3,i)); 

        if coh_map(seedx,seedy,seedz) > coh_limit      
            numAtWM = numAtWM + 1;
        end 
    end
    percOnWM(i,1) = numAtWM;
end
percOnWM = 100*percOnWM/n_seeds;

% Plot regular graph
figure();   plot(percOnWM,'r');
% Plot smooth graph
hold on;   plot(smoothdata(percOnWM),'--','Color','black','LineWidth',2);

xlabel('Time Intervals');   ylabel('Percent of NSC on WM');
savethis('percentOnWM');




%% Determine whether NSC made it within a certain radius to cancer center
if has_cancer
    percAtCancerGraph = zeros(Finaltimestep,1);
    for i = 1:Finaltimestep
        numAtCancer = 0;
        for j = 1:n_seeds
            if all( abs(cancer_center - p(j).coord(:,i)') < cancer_size )
                numAtCancer = numAtCancer + 1;
            end
        end
        percAtCancerGraph(i,1) = numAtCancer;
    end
    percAtCancerGraph = 100*percAtCancerGraph/n_seeds;

    figure();   plot(percAtCancerGraph,'r');
    xlabel('Time Intervals');   ylabel('Percent of NSC that Reach Cancer Site');
    savethis('percentAtCancer');
end




%% Entire procedure finished
disp( 'done' );
%% Functions
function [concentration, cgradX, cgradY, cgradZ, cancer_sd] = set_cancer(coh_map) 
    global cancer_center cancer_size
    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    cancer_sd = cancer_size;
    cancer_sd_X = cancer_sd(1); cancer_sd_Y = cancer_sd(2); cancer_sd_Z = cancer_sd(3);
    
    concen_sd = 25; 
    concen_sd = concen_sd*cancer_sd/cancer_sd(1); 
    concen_param = 2; 
    concentration = 1./(1+(sqrt(((X - cancer_center(1))/concen_sd(1)).^2+((Y - cancer_center(2))/concen_sd(2)).^2+((Z - cancer_center(3))/concen_sd(3)).^2)).^concen_param);    
 
    valcap = 1./(1+(sqrt(((cancer_sd_X)/concen_sd(1)).^2+((0)/concen_sd(2)).^2+((0)/concen_sd(3)).^2)).^concen_param);      %valcap = concentration at border of cancer
    
    concentration(concentration > valcap) = valcap;
    concentration( coh_map==0 ) = 0; 
    
    cgradX = zeros(size(X));
    cgradX(2:end-1,:,:) = concentration(3:end,:,:) - concentration(1:end-2,:,:);

    cgradY = zeros(size(X));
    cgradY(:,2:end-1,:) = concentration(:,3:end,:) - concentration(:,1:end-2,:);

    cgradZ = zeros(size(X));
    cgradZ(:,:,2:end-1) = concentration(:,:,3:end) - concentration(:,:,1:end-2);
    
    % find largest vector magnitude
    L = (cgradX.^2 + cgradY.^2 + cgradZ.^2).^(.5);
    maxL = max(max(L));

    % normalize gradient
    cgradX = cgradX./maxL;
    cgradY = cgradY./maxL;
    cgradZ = cgradZ./maxL;
end

function [seed_ind,seed_sd,indup,inddown] = set_initial(coh_map) 
    global inj_center
    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    seed_sd = 10; 
    seed_ROI = sqrt((X-inj_center(1)).^2 + (Y-inj_center(2)).^2 + (Z-inj_center(3)).^2)<seed_sd;   
    seed_ind = find(seed_ROI);
    
    %Create better upper and lower bounds (z-axis) for mouse brain
    for n = 1:size(coh_map,1)
      for m = 1:size(coh_map,2)
        ind = find( squeeze(coh_map(n,m,:)) ~=0 );
        if( isempty(ind) )
            indup(n,m) = NaN;
            inddown(n,m) = NaN;  
        else
            indup(n,m) = max(ind);
            inddown(n,m) = min(ind);
        end
      end
    end
end

function [EV, FA] = load_3D()
    % OPEN AND LOAD EIGENVECTOR AND FRACTIONAL ANISOTROPY FILE
    EV_file = strcat('../NSC_atlas_and_path_simulation/CTRLP60_avg_vec.img');
    FA_file = strcat('../NSC_atlas_and_path_simulation/CTRLP60_avg_fa.img');
    fid = fopen(EV_file,'r');
    EV = fread(fid,'float32','ieee-le');
    EV = permute(reshape(EV,3,200,280,128),[2,3,4,1]);
    fclose(fid);

    % fractional anisotropy 
    fid = fopen(FA_file,'r');
    FA = fread(fid,'float32','ieee-le');
    FA = reshape(FA,200,280,128);
    fclose(fid);

    % INTERPOLATE TO THE APPROPRIATE SIZE
    limits = [21,180,21,160,36,105];    scale = 2;
    %limits = [1,200,1,280,1,128];    scale = 0;
    EV = EV(limits(1):limits(2),limits(3):limits(4),limits(5):limits(6),:);
    EV_x = interp3(EV(:,:,:,1),scale);  EV_y = interp3(EV(:,:,:,2),scale);  EV_z = interp3(EV(:,:,:,3),scale);
    clear EV; 
    EV(:,:,:,1) = EV_x; EV(:,:,:,2) = EV_y; EV(:,:,:,3) = EV_z;    
    clear EV_x EV_y EV_z

    FA = interp3(FA(limits(1):limits(2),limits(3):limits(4),limits(5):limits(6)),scale);
end 

function savethis(title)
    global modelType d_w d_g chemo_sensitivity alpha4chmtx beta4dist cancer_center FolderName1 FolderName2 has_cancer;
    
    type = '.png';
    if strcmp(title,'trajectoryAll')
        type = '.fig';
    end
    
    if has_cancer
        saveas( gcf, [pwd strcat( FolderName2, '3D_191125_',modelType,'_',title,'_d',num2str(d_w),'_',num2str(d_g),...
            '_a',int2str(alpha4chmtx),'_b',int2str(beta4dist),'_c',num2str(chemo_sensitivity),...
            '_[',int2str(cancer_center(1)),',',int2str(cancer_center(2)),',',int2str(cancer_center(3)),']',type)] );
    else
        saveas( gcf, [pwd strcat( FolderName1, '3D_191125_',modelType,'_',title,'_d',num2str(d_w),'_',num2str(d_g),...
            '_a',int2str(alpha4chmtx),'_b',int2str(beta4dist),'_c',num2str(chemo_sensitivity),'_[NA]',type)] );
    end
end
