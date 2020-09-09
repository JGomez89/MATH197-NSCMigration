% BIDIRECTOINAL PATH GENERATION

%%_________________
%Using 5 voxel kernerl blur orientation and coherency map 

%% Load data 
clear all

disp('Gathering data...')
[eigen_map, coh_map] = load_3D(); 


%% Set parameters 
set_parameters_newfigure()

% Storing the coorinate, angle of eigen vector and coherency value for both
% positive and negative direction paths that start from a seed point
p = struct('coord',{},'coherency',[]);  % the two fields 

%% Initialize a matrix of same size as the white matter image and initialize a seed point (inj center) around which the seeds for each simulation are placed
[seed_ind,seed_sd] = set_initial(coh_map);                                                                          %Nasal injection initial 


%% Cancer injection
[concentration, cgradX, cgradY, cgradZ, cancer_sd] = set_cancer(coh_map);                  %Near contralateral corpus callosum


%% For each seed
Tstart = tic; 
disp('Running simulation...')


for seed_loop = 1:n_seeds
    
    %%%% initialize particular standard step length
    dmax = betainv(rand(1), 1, beta4dist);
    
    %%%% initialize first two steps 
    seed = seed_ind(round(rand*length(seed_ind)));                                                          % Choose a random point index to start the simulation
    [seed(1,1),seed(1,2),seed(1,3)] = ind2sub(size(coh_map),seed);                                          % Coordinate of this starting point
    p(seed_loop).coord(:,1) = seed';                                                                        % Starting location stored in data matrix
%     ind = sub2ind(size(coh_map),seed(1,1),seed(1,2),seed(1,3));                                           % Index of this point in the white matter image
    % If coherency map data is present then use the angle indicated by the
    % data, else use a random angle

    %%%% second step 
    i = 2; 
%     ind = sub2ind(size(coh_map),round(p(seed_loop).coord(1,i-1)),round(p(seed_loop).coord(2,i-1)),round(p(seed_loop).coord(3,i-1)));       
            
%     angle = -ori_map(ind);
%     p(seed_loop).angle1(i) = angle;            
%     p(seed_loop).coherency1(i) = coh_map(ind);
    
    eigen_vect = squeeze( eigen_map( round(p(seed_loop).coord(1,i-1)),round(p(seed_loop).coord(2,i-1)),round(p(seed_loop).coord(3,i-1)), : ));
    dstep = (rand(1)*2-1)*dmax; 
    p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*eigen_vect;                                 % Find next coordinate point at step d along EV direction
    
%     path = zeros(size(WM));
    
    for i = 3:Finaltimestep                                                                                 % Run the loop a large number of steps            
        
        seedx = round(p(seed_loop).coord(1,i-1)); 
        seedy = round(p(seed_loop).coord(2,i-1)); 
        seedz = round(p(seed_loop).coord(3,i-1)); 

        ind = sub2ind(size(coh_map),seedx,seedy,seedz);                                                     % Index of previous coordinate point

        if coh_map(ind) > coh_limit                                                                         % Get Eigen vector angle at the previous coordinate point
            eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) ); 
            d = dmax * d_w; 
        else  
            eigen_vect_noise = (rand(3,1)*2-1); 
            eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) ) + eigen_vect_noise/norm(eigen_vect_noise)*norm(eigen_vect); 
            % move with less magnitude if no WM present 
            d = dmax * d_g; 
        end 

        %%%%%% uniform[0, 1] 
%             dstep = rand(1)*dmax ; 
        %%%%%% beta[1, beta], if beta==1, uniform 
%         dstep = betainv( rand(1), 1, beta4dist )*d ; 
        %%%%%% deterministic
        dstep = d; 


        switch modelNum

            case 1
                modelType = 'model1';
                % Find next coordinate point at step d along EV direction 
                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) ... 
                    - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * dstep*eigen_vect;
                % let it move away from the direction it came from. 

            case 2
                % add chemotaxis 
                modelType = 'model2';
%                 eigen_vect = - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * eigen_vect';
                chmtx_vect = [cgradZ(seedx,seedy,seedz); cgradY(seedx,seedy,seedz); cgradX(seedx,seedy,seedz)];

                %  Move along chemotaxis
                if has_cancer
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                else
                    chmtx_bias = 0;
                end

                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*eigen_vect + chmtx_bias*chmtx_vect;

            case 3 
                %  move along chemotaxis once it is near by cancer 
                modelType = 'model3';
                eigen_vect = - sign(dot( p(seed_loop).coord(:,i-2) - p(seed_loop).coord(:,i-1), eigen_vect ) + eps ) * eigen_vect';
                chmtx_vect = [cgradZ(seedx,seedy,seedz); cgradY(seedx,seedy,seedz); cgradX(seedx,seedy,seedz)];

                %  Move along chemotaxis once it is near by cancer
                if has_cancer && (concentration(seedx,seedy,seedz) >= chmtx_limit)
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                else
                    chmtx_bias = 0;
                end

                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*eigen_vect + chmtx_bias*chmtx_vect;
        end            

        % the path is out of bounds or hit max steps then don't proceed
        seednow = round(p(seed_loop).coord(:,i));
        if ~(all( seednow' - size(coh_map)<=-1) && all( seednow' >1 ) )
            break
        end
       
        
    end    
end

coord = zeros(Finaltimestep,3,n_seeds);
for i = 1:n_seeds 
    if( size( p(i).coord, 2 ) < Finaltimestep )                                                             % If stopped before Finaltimestep
        ind = size( p(i).coord, 2 ) - 1; 
        p(i).coord(1,(ind+1):Finaltimestep) = p(i).coord(1,ind); 
        p(i).coord(2,(ind+1):Finaltimestep) = p(i).coord(2,ind);
        p(i).coord(3,(ind+1):Finaltimestep) = p(i).coord(3,ind);         
    end
    coord(:,:,i) = p(i).coord';
end

toc( Tstart ); 


%% Plotting the paths 

% gif of 2D slice || implay
% somehow 3D?
disp('Plotting graphs...')


%%% Plot path of seeds
figure; hold on;
ind = 0;
acc = 55;
for i=1:    round(size(coh_map,1)/acc):    round(size(coh_map,1))
    for j=1:    round(size(coh_map,2)/acc):    round(size(coh_map,2))
        for k=1:    round(size(coh_map,3)/acc):    round(size(coh_map,3)) 
            
            if coh_map(i,j,k) > coh_limit
                ind=ind+1;
                X(ind,1) = j;
                Y(ind,1) = i;
                Z(ind,1) = k;
            end
            
        end
    end
end

[x,y,z] = sphere;
x = x*seed_sd + inj_center(2);
y = y*seed_sd + inj_center(1);
z = z*seed_sd + inj_center(3);
h = surf(x, y, z);                                                                                          %Plot injection sphere
set(h,'FaceColor',[0 .3 1],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
daspect([1 1 1]);

if has_cancer
    [x,y,z] = sphere;
    x = x*cancer_sd + cancer_center(2);
    y = y*cancer_sd + cancer_center(1);
    z = z*cancer_sd + cancer_center(3);
    h = surf(x, y, z);                                                                                      %If there's cancer, plot cancer sphere
    set(h,'FaceColor',[1 .7 0],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
end

plot3(X,Y,Z,'om','MarkerSize',2)                                                                            %Plot coherent portion of mouse brain
xfinal = zeros(1,n_seeds); yfinal = zeros(1,n_seeds); zfinal = zeros(1,n_seeds);
for i=1: n_seeds
    plot3(p(i).coord(2,:),p(i).coord(1,:),p(i).coord(3,:),'linewidth',3)                                    %Plot path of seed
    xfinal(i) = p(i).coord(1, end); 
    yfinal(i) = p(i).coord(2, end);
    zfinal(i) = p(i).coord(3, end);
end
plot3(yfinal,xfinal,zfinal,'ok','Markersize', 7, 'markerfacecolor', 'w')                                    %Plot final coords of each seed
xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; daspect([1 1 1]); camlight;

savethis('trajectoryAll');



disp( 'done' ); 


%% Functions
function [concentration, cgradX, cgradY, cgradZ, cancer_sd] = set_cancer(coh_map) 
    global cancer_center
    [X,Y,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    cancer_sd = 20;
    concentration = exp((-(X - cancer_center(2)).^2 - (Y - cancer_center(1)).^2 - (Z - cancer_center(3)).^2 )./cancer_sd^2);

    cgradY = zeros(size(X));
    cgradY(2:end-1,:,:) = concentration(3:end,:,:) - concentration(1:end-2,:,:);

    cgradX = zeros(size(X));
    cgradX(:,2:end-1,:) = concentration(:,3:end,:) - concentration(:,1:end-2,:);

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

function [seed_ind,seed_sd] = set_initial(coh_map) 
    global inj_center
    [X,Y,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    seed_sd = 10; 
    seed_ROI = sqrt((X-inj_center(2)).^2 + (Y-inj_center(1)).^2 + (Z-inj_center(3)).^2)<seed_sd;   
    seed_ind = find(seed_ROI);
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

function savethis(title)
    global modelType d_w d_g chemo_sensitivity alpha4chmtx beta4dist cancer_center FolderName1 FolderName2 has_cancer;
    if has_cancer
        saveas( gcf, [pwd strcat( FolderName2, '3D_191125_',modelType,'_',title,'_d',num2str(d_w),'_', num2str(d_g),...
            '_a',int2str(alpha4chmtx),'_b',int2str(beta4dist),'_c',num2str(chemo_sensitivity),...
            '_[',int2str(cancer_center(1)),',',int2str(cancer_center(2)),',',int2str(cancer_center(3)),']', '.fig')] );
    else
        saveas( gcf, [pwd strcat( FolderName1, '3D_191125_',modelType,'_',title,'_d',num2str(d_w),'_', num2str(d_g),...
            '_a',int2str(alpha4chmtx),'_b',int2str(beta4dist),'_c',num2str(chemo_sensitivity),'_[NA]', '.fig')] );
    end
end
