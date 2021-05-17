% BIDIRECTOINAL PATH GENERATION

%%_________________
%Using 5 voxel kernerl blur orientation and coherency map 

%% Load data
global eigen_map coh_map seed_sd ubound lbound cancer_sd

% clear all;
close all;

disp('Gathering data...')
[eigen_map, coh_map] = load_3D(); 


%% Set parameters 
set_parameters_newfigure()

% Storing coorinate in struct p
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
            eigen_vect = eigen_vect_noise/norm(eigen_vect_noise); 
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
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) ...
            - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * dstep*eigen_vect + chmtx_bias*chmtx_vect;       

        % if path is out of bounds then put back to previous step
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
            
                break
%             p(seed_loop).coord(:,i) =  p(seed_loop).coord(:,i-1) ;                      
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




%% Plot Graphs, Use smaller p
disp('Plotting graphs...')

if exist('p_original','var') && exist('acc','var')
    for n=1:n_seeds
        if p(n).coord(:,2) ~= p_original(n).coord(:,1+acc)
            p_original = p; break
        end
    end
    p = p_original;
end

acc = 1;                                                                                                  %Use a smaller timestep (record every acc number steps)
p_original = p;
p_new = struct('coord',{});
for i=1:n_seeds
    p_new(i).coord = p(i).coord(:,1:acc:Finaltimestep); 
end
p = p_new; clear p_new;
xvals = (1:acc:Finaltimestep);




%% Plot path of seeds
figure;
plotenvironment(true);

for i=1: n_seeds
    plot3(p(i).coord(1,:),p(i).coord(2,:),p(i).coord(3,:),'linewidth',3);                                       %Plot path of seed
    plot3(p(i).coord(1, end),p(i).coord(2, end),p(i).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w'); %Plot final coords of each seed
end

savethis('trajectoryAll','fig');




%% Boxplots of NSC's distance from center of injection site
clear distInit;
distInit = zeros(n_seeds,Finaltimestep/acc);
for n = 1:n_seeds
    distInit(n,:) = sqrt( sum ( (p(n).coord - inj_center'*ones(1,size(p(n).coord,2))).^2, 1 ) );
end

num_plots = 10;
data = distInit(:, (1:num_plots)*floor(Finaltimestep/(num_plots*acc)))*CONVERT2MICRON;
labels = {num_plots};
for i=1:num_plots
    labels{i} = num2str( floor(Finaltimestep/num_plots)*i );
end
figure;     boxplot( data, labels );
hold on;    plot( (1:acc:Finaltimestep)*(num_plots/Finaltimestep), median( distInit )*CONVERT2MICRON );
% hold on;    plot( (1:acc:Finaltimestep)*(num_plots/Finaltimestep), mean( distInit )*CONVERT2MICRON );

ylim([0 max(max(data))+500]);
xlim([0 num_plots+1]);
xlabel( 'Time Intervals' ); ylabel( 'Distance from injection site (\mu)' );
savethis('distboxplot','png');




%% Percent of NSC on WM
numAtWM = zeros(Finaltimestep/acc,1);
for i = 1:n_seeds
    seedx = round(p(i).coord(1,:)); 
    seedy = round(p(i).coord(2,:)); 
    seedz = round(p(i).coord(3,:)); 
    ind = sub2ind(size(coh_map),seedx, seedy, seedz); 
    numAtWM = numAtWM + double( coh_map(ind) > coh_limit )'; 
end
percOnWM = numAtWM / n_seeds * 100;

% Plot regular graph
figure();   plot(xvals,percOnWM,'r');
% Plot smooth graph
% hold on;   plot(xvals,smoothdata(percOnWM),'--','Color','black','LineWidth',2);

xlabel('Time Intervals');   ylabel('Percent of NSC on WM');
savethis('percentOnWM','png');




%% Determine whether NSC made it within a certain radius to cancer center
if has_cancer
    percAtCancerGraph = zeros(Finaltimestep/acc,1);
    for i = 1:Finaltimestep/acc
        numAtCancer = 0;
        for j = 1:n_seeds
            if all( abs(cancer_center - p(j).coord(:,i)') < cancer_size )
                numAtCancer = numAtCancer + 1;
            end
        end
        percAtCancerGraph(i,1) = numAtCancer;
    end
    percAtCancerGraph = 100*percAtCancerGraph/n_seeds;

    figure();   plot(xvals,percAtCancerGraph,'r');
    xlabel('Time Intervals');   ylabel('Percent of NSC that Reach Cancer Site');
    savethis('percentAtCancer','png');
end




%% Percent of WM in brain
% figure();
% percWM = @(x) 100 * sum(coh_map > x,'all')/sum(coh_map > 0,'all');
% for i=1:50
%     percWMArr(i)=percWM(i/50);
% end
% plot((1:50)/50,percWM((1:50)/50));




%% Plot path of seeds on and off WM
figure; 
plotenvironment(false);
a=1;
plotonoffWM( p(1:a:end) );
savethis('trajectoryOnOffWM','fig');




%% Min/Max Path
figure;
plotenvironment(false);

stepDistance = zeros(n_seeds, Finaltimestep/acc - 1);
for i = 1:n_seeds
    for j = 2:Finaltimestep/acc
        % record distance traveled each step
        stepDistance(i, j - 1) = sqrt((p(i).coord(1,j) - p(i).coord(1,j-1)).^2 + (p(i).coord(2,j) - p(i).coord(2,j-1)).^2 + (p(i).coord(3,j) - p(i).coord(3,j-1)).^2);
    end
end

cellind = [ ];
for j = 1:n_seeds
    if (all( abs(cancer_center - p(j).coord(:,end)') < cancer_size ))
        cellind = [cellind, j];
    end
end 

% compute cell trajectory lengths 
stepDistanceTotal = sum( stepDistance, 2 ); 

% find minimum path among all cells
minPath = min( stepDistanceTotal(:) ); 
minPath_seed = find( stepDistanceTotal == minPath );

% find maximum path among all cells
maxPath = max( stepDistanceTotal(:) ); 
maxPath_seed = find( stepDistanceTotal == maxPath );

if has_cancer
    dist2cancer = zeros(n_seeds, 1); 

    % Min Calculations
    % find minimum path among those that are close to cancer 
    minPath_cancer = min( stepDistanceTotal(cellind) ); 
    minPath_cancer_seed = find( stepDistanceTotal == minPath_cancer ); 

    % find cells that are close enough to minimum path to cancer 
    mincloseindex = find( abs( stepDistanceTotal-minPath ) < 100 ); 

    % percentage of cells that are near minimum path to cancer
    minclosepercent = length(mincloseindex) / length(stepDistanceTotal);
    
    % Max Calculations
    % find maximum path among cells that are close to cancer 
    maxPath_cancer = max( stepDistanceTotal(cellind) ); 
    maxPath_cancer_seed = find( stepDistanceTotal == maxPath_cancer ); 

    % find cells that are close enough to maximum path 
    maxcloseindex = find( abs( stepDistanceTotal-maxPath ) < 100 ); 

    % percentage of cells that are near maximum path 
    maxclosepercent = length(maxcloseindex) / length(stepDistanceTotal);
        
    %Plot cancer Min/Max Paths:
    plot3(p(minPath_cancer_seed).coord(1,:),p(minPath_cancer_seed).coord(2,:),p(minPath_cancer_seed).coord(3,:),'Color','m','linewidth',3);
    plot3(p(minPath_cancer_seed).coord(1, end),p(minPath_cancer_seed).coord(2, end),p(minPath_cancer_seed).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w');

    plot3(p(maxPath_cancer_seed).coord(1,:),p(maxPath_cancer_seed).coord(2,:),p(maxPath_cancer_seed).coord(3,:),'y','linewidth',3);
    plot3(p(maxPath_cancer_seed).coord(1, end),p(maxPath_cancer_seed).coord(2, end),p(maxPath_cancer_seed).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w');
end

%Plot overall Min/Max Paths:
plot3(p(minPath_seed).coord(1,:),p(minPath_seed).coord(2,:),p(minPath_seed).coord(3,:),'Color','m','linewidth',3);
plot3(p(minPath_seed).coord(1, end),p(minPath_seed).coord(2, end),p(minPath_seed).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w');

plot3(p(maxPath_seed).coord(1,:),p(maxPath_seed).coord(2,:),p(maxPath_seed).coord(3,:),'Color','y','linewidth',3);
plot3(p(maxPath_seed).coord(1, end),p(maxPath_seed).coord(2, end),p(maxPath_seed).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w');

savethis('minmaxpath','fig');




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
    EV_file = strcat('NSC_atlas_Vikram/CTRLP60_avg_vec.img');
    FA_file = strcat('NSC_atlas_Vikram/CTRLP60_avg_fa.img');
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

function savethis(title,type)
    global modelType d_w d_g chemo_sensitivity alpha4chmtx beta4dist cancer_center FolderName1 FolderName2 has_cancer;    
    if has_cancer
        saveas( gcf, [pwd strcat( FolderName2, '3D_191125_',modelType,'_',title,'_d',num2str(d_w),'_',num2str(d_g),...
            '_a',int2str(alpha4chmtx),'_b',int2str(beta4dist),'_c',num2str(chemo_sensitivity),...
            '_[',int2str(cancer_center(1)),',',int2str(cancer_center(2)),',',int2str(cancer_center(3)),'].',type)] );
    else
        saveas( gcf, [pwd strcat( FolderName1, '3D_191125_',modelType,'_',title,'_d',num2str(d_w),'_',num2str(d_g),...
            '_a',int2str(alpha4chmtx),'_b',int2str(beta4dist),'_c',num2str(chemo_sensitivity),'_[NA].',type)] );
    end
end

function plotenvironment(plotWM)
    global coh_limit coh_map ubound lbound seed_sd inj_center cancer_size cancer_center has_cancer;
    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    hold on; surf( ubound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );      %Plot upper bound of mouse brain
    hold on; surf( lbound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );      %Plot lower bound of mouse brain

    [x,y,z] = sphere;
    x = x*seed_sd + inj_center(1);
    y = y*seed_sd + inj_center(2);
    z = z*seed_sd + inj_center(3);
    h = surf(x, y, z);                                                                                      %Plot injection site
    set(h,'FaceColor',[1 0 1],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');

    if has_cancer
        [x,y,z] = sphere;
        x = x*cancer_size(1) + cancer_center(1);
        y = y*cancer_size(2) + cancer_center(2);
        z = z*cancer_size(3) + cancer_center(3);
        h = surf(x, y, z);                                                                                  %Plot cancer site
        set(h,'FaceColor',[1 .7 0],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
    end
    if plotWM
        ii = find( coh_map > coh_limit ); 
        hold on; plot3( X(ii(1:200:end)),Y(ii(1:200:end)),Z(ii(1:200:end)),'b.');                            %Plot WM track
    end
    
    xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; daspect([1 1 1]); camlight;
end

function plotonoffWM(p)
    global coh_limit coh_map eigen_map;
    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    for i=1: size(p,2)
        plot3(p(i).coord(1, end),p(i).coord(2, end),p(i).coord(3, end),'ok','Markersize', 7,'markerfacecolor','w'); %Plot final coords of each seed
    end
    
    isOnWM = @(x) logical( coh_map(x) > coh_limit );
    for i=1: size(p,2)
        p_temp = p(i).coord;
        ind = sub2ind( size(coh_map),round(p_temp(1,:)),round(p_temp(2,:)),round(p_temp(3,:)) );

        p_on = p_temp .* isOnWM(ind);
        p_on(p_on == 0) = NaN;
        p_off = p_temp .* ~isOnWM(ind);
        p_off(p_off == 0) = NaN;

        onNaN=0;
        for j=1: size(p_temp,2)-1
            if all( isnan(p_off(:,j+1)) )
                if ~onNaN
                    p_on(:,j)=p_off(:,j);
                    onNaN=1;
                end
            else
                onNaN=0;
            end
        end
        onNaN=0;
        for j=1: size(p_temp,2)-1
            if all( isnan(p_on(:,j+1)) )
                if ~onNaN
                    p_off(:,j)=p_on(:,j);
                    onNaN=1;
                end
            else
                onNaN=0;
            end
        end
        
        hold on; plot3(p_on(1,:),p_on(2,:),p_on(3,:),'Color',[0.6350 0.0780 0.1840],'linewidth',3);         %Plot path of seed on WM
        hold on; plot3(p_off(1,:),p_off(2,:),p_off(3,:),'Color',[0 0.4470 0.7410],'linewidth',3);           %Plot path of seed off WM

        emap_reshaped = reshape( eigen_map, numel(coh_map), 3 ); 
        
        ind_on = sub2ind( size(coh_map),round(p_on(1,:)),round(p_on(2,:)),round(p_on(3,:)) );
        ind_on = ind_on(~isnan(ind_on));
        hold on; quiver3( X(ind_on),Y(ind_on),Z(ind_on),emap_reshaped(ind_on,1)',emap_reshaped(ind_on,2)',emap_reshaped(ind_on,3)','Color',[0.6350 0.0780 0.1840]);     %Plot eigen vects on WM
        ind_off = sub2ind( size(coh_map),round(p_off(1,:)),round(p_off(2,:)),round(p_off(3,:)) );
        ind_off = ind_off(~isnan(ind_off));
        hold on; quiver3( X(ind_off),Y(ind_off),Z(ind_off),emap_reshaped(ind_off,1)',emap_reshaped(ind_off,2)',emap_reshaped(ind_off,3)','Color',[0 0.4470 0.7410]);    %Plot eigen vects off WM
            
        plot3( p_temp(1,:), p_temp(2,:), p_temp(3,:), 'k.' );                                               %Plot trajectory dots 
    end
end