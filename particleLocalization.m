% Robotics: Estimation and Learning 
% WEEK 4
% 
% Complete this function following the instruction. 
function myPose = particleLocalization(ranges, scanAngles, map, param)

% Number of poses to calculate
N = size(ranges, 2);
% Output format is [x1 x2, ...; y1, y2, ...; z1, z2, ...]
myPose = zeros(3, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Map Parameters 
% 
% % the number of grids for 1 meter.
myResolution = param.resol;
% % the origin of the map in pixels
myOrigin = param.origin; 

% The initial pose is given
myPose(:,1) = param.init_pose;
% You should put the given initial pose into myPose for j=1, ignoring the j=1 ranges. 
% The pose(:,1) should be the pose when ranges(:,j) were measured.
[max_rows, max_cols] = size(map);

% Decide the number of particles, M.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 1000;                       % Please decide a reasonable number of M, 
                               % based on your experiment using the practice data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create M number of particles
P = repmat(myPose(:,1), [1, M]);
%weights = 1.0/M * ones(M, 1);
%scores = zeros(M, 1);
%scanAngles = scanAngles;
%x_vel = 0.85;
%y_vel = 0.8;
%theta_vel = 0.017;
%eff_thresh = 80.0;
sigma_m = 0.01 * [ 1; 1; 2 ];
radius = 0.025;
map_threshold_low  = mode(mode(map))-0.2;
map_threshold_high = mode(mode(map))+0.5;
resample_threshold = 0.8;
for j = 2:N % You will start estimating myPose from j=2 using ranges(:,2).
    W = 1.0/M * ones(1,M);
    corr = zeros(1, M);
    P = repmat(myPose(:,j-1), [1, M]);
    R = radius;
    P(1,1:M) = P(1,1:M) + R.*cos(P(3,1:M));
    P(2,1:M) = P(2,1:M) - R.*sin(P(3,1:M));
    P = P + randn(3, M).*sigma_m;
    
    for i = 1:M
    %   2-1) Find grid cells hit by the rays (in the grid map coordinate frame)
        global_coord(:,1) = ceil(( ranges(:,j).*cos(scanAngles + P(3,i)) + P(1,i))*myResolution + myOrigin(1));
        global_coord(:,2) = ceil((-ranges(:,j).*sin(scanAngles + P(3,i)) + P(2,i))*myResolution + myOrigin(2));
        %fprintf('LIDAR max coords (x, y) = %d, %d\n', max(lidar_global(:,1)), max(lidar_global(:,2)));
        
    %   2-2) For each particle, calculate the correlation scores of the particles
        % for out of range positions, assume the original position on the map, resulting in a low correlation score
        global_coord(global_coord(:,1) < 1, 1) = myOrigin(1);
        global_coord(global_coord(:,1) < 1, 2) = myOrigin(2);
        global_coord(global_coord(:,2) < 1, 1) = myOrigin(1);
        global_coord(global_coord(:,2) < 1, 2) = myOrigin(2);
        global_coord(global_coord(:,1) > size(map,2), 1) = myOrigin(1);
        global_coord(global_coord(:,1) > size(map,2), 2) = myOrigin(2);
        global_coord(global_coord(:,2) > size(map,1), 1) = myOrigin(1);
        global_coord(global_coord(:,2) > size(map,1), 2) = myOrigin(2);
        occ_ind = sub2ind(size(map), global_coord(:,2), global_coord(:,1));
        corr_values = map(occ_ind);
        corr(i) = -3*sum(corr_values<=map_threshold_low)+10*sum(corr_values>=map_threshold_high);
    end
    corr = corr - min(corr);    % make range from 0, inf
    
    %   2-3) Update the particle weights
    W = (W(1:M) .* corr) / sum(corr);
    W = W / sum(W);
    
    %   2-4) Choose the best particle to update the pose
    [ value, idx ] = max(W);    % simply the one with the largest weight
    myPose(:,j) = P(:,idx);
        
        % the world coordinates of the robot
    

    
    % 3) Resample if the effective number of particles is smaller than a threshold
    N_eff = sum(W)^2/sum(W.^2);
    if (N_eff < resample_threshold*M)         
        W_cum = cumsum(W);
         for i = 1:M
             index = find(rand <= W_cum,1);
             P_new(:,i) = P(:,index);
             W_new(i) = W(index);
         end
         P = P_new;
         W = W_new;
             
    end
    
    
    
    

   



end       

        

end

