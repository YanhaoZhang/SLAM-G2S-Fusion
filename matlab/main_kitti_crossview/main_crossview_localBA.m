function [vT_us0_usk, vM_uok_usk] = main_crossview_localBA(vT_us0_usk_org, vT_us0_usk, vT_uo0_uok, T_u0a0, P_b0a0, f_a0_ia0, map_ids, ...
                                                           matched_pts_ids, matched_keyframe_ids, matched_feature, ...
                                                           fine_valid, vM_uok_usk, idcv)
% vT_us0_usk_org: original SLAM trajectory to be updated
% vT_us0_usk: SLAM refined and to be refined
% vT_uo0_uok: oxts trajectory (ground truth just for evaluation)
% T_u0a0: cam0 to cam2. The map points are in cam0 frame0
% P_b0a0: projection matrix (front stereo)
% f_a0_ia0: xyz map points
% map_ids: map point index by ORB-SLAM
% matched_pts_ids: matched 3D point SLAM index
% matched_keyframe_ids, matched keyframe SLAM index 
% matched_feature: left camera 2D features and u of right camera 2D feature
% fine_valid: valid [lat, lon, yaw] by cross-view
% vM_uok_usk: store all valid SLAM to oxts prediction

% notice: we need to update pred_vM_uok_usk everytime when we update the
% trajectory. The prediction is based on the trajectory before refinement


%% rewrite visual odometry
num_frame = size(vT_us0_usk_org,1);

% visual odometry measurement
vo_vM_usl_usk = cell(num_frame-1, 3); 

id = 0;
for k=2:num_frame

T_us0_usl = squeeze(vT_us0_usk_org(k-1,:,:));
T_us0_usk = squeeze(vT_us0_usk_org(k,:,:));

T_usl_us0 = T_us0_usl \ eye(4);


% store data
id=id+1;
vo_vM_usl_usk{id, 1} = T_usl_us0 * T_us0_usk;
vo_vM_usl_usk{id, 2} = k-1;
vo_vM_usl_usk{id, 3} = k;

end


%% load cross-view measurement
valid_id = fine_valid(:,1) & fine_valid(:,2) & fine_valid(:,3);    % consider all [lat, lon, yaw] valid as valid cross-view
num_valid = sum(valid_id);
pred_vM_uok_usk = cell(num_valid, 4);      % relative pose from cross-view, oxts pose, id

id = 0;
for k=1:num_frame

if (~valid_id(k))
    continue
end

if k>idx+1
    break
end

id = id + 1;


% T_uo0_uok = squeeze(oxts_vTu0uk(k,:,:));
T_us0_usk = squeeze(vT_us0_usk(k,:,:));

% T_uok_uo0 = T_uo0_uok \ eye(4);
% T_uok_usk = T_uok_uo0 * T_us0_usk;

% [roll_uok_usk, pitch_uok_usk, yaw_uok_usk] = euler_from_rotation_matrix(T_uok_usk(1:3, 1:3));
% eul = rotm2eul(T_uok_usk(1:3, 1:3));    %yaw pitch roll
% R_uok_usk = euler_to_rotation_matrix(roll_uok_usk, pitch_uok_usk, yaw_uok_usk)


% roll_uok_usk = 0;
% pitch_uok_usk = 0;
% T_uok_usk = eye(4);

% load pred result
% pred_lon = pred_lons(k);
% pred_lat = pred_lats(k);
% pred_yaw = pred_oriens(k);

% pred_yaw_uok_usk = -pred_yaw * pi / 180;
% R_uok_usk = euler_to_rotation_matrix(roll_uok_usk, pitch_uok_usk, pred_yaw_uok_usk);

% M_uok_usk = T_uok_usk;
% M_uok_usk(1,4) = -pred_lon;
% M_uok_usk(2,4) = -pred_lat;
% M_uok_usk(1:3,1:3) = R_uok_usk;

M_uok_usk = squeeze(vM_uok_usk(k,:,:)); 

% store predicted pose
pred_vM_uok_usk{id, 1} = M_uok_usk;
pred_vM_uok_usk{id, 3} = k;
pred_vM_uok_usk{id, 4} = T_us0_usk;    % SLAM pose is actually used
end

% save('/home/users/u1111398/YanhaoZhang/ORB_SLAM3_stereoMultiCam/matlab/debug_.mat', 'pred_vM_uok_usk', 'num_valid', 'valid_id')


%% load map point data
num_pts = size(f_a0_ia0,1);

k = pred_vM_uok_usk{idcv, 3};
pred_M_uok_usk = cell(1,3);
pred_M_uok_usk{1} = pred_vM_uok_usk{idcv, 1};
pred_M_uok_usk{2} = pred_vM_uok_usk{idcv, 2};
pred_M_uok_usk{3} = pred_vM_uok_usk{idcv, 3};
pred_M_uok_usk{4} = pred_vM_uok_usk{idcv, 4};


matched_feat = [matched_pts_ids, matched_keyframe_ids, matched_feature];

% retrive valid map points: check previous 5 frames, and select map points
% that are observed in all num_covi frames
% id_obs = matched_keyframe_ids == (k-1);   % python starts from 0
% num_obs = size(matched_feat,1);
% id_obs = false(num_obs,1);    % initialize valid observation id. 
times_obs = zeros(num_pts,1);    % store the number of keyframes that this map points are observed
num_BAframes = 5;
num_covi = 4;    % minimum number of frames that a points should be observed
for i=1:num_BAframes

if(k-i<0)       % only consider valid frame id 
    break;
end

% if i==1
% id_obs = matched_feat(:,2) == (k-i);   
% Index_pts = unique(matched_feat(id_obs,1));    % find potential map index (SLAM id), remove duplicated index
%     % [~,id_pts] = intersect(map_ids, index_pts);  % find the potential map ids
% else
% id_obs = matched_feat(:,2) == (k-i);
% index_pts = unique(matched_feat(id_obs,1));
% 
% 
% [Index_pts, ~, ~] = intersect(Index_pts, index_pts);    % only consider 
% 
% end



% python starts from 0, therefore we need to minus 1 (i starts from 1)
id_obs_i = matched_feat(:,2) == (k-i);   
index_pts_i = unique(matched_feat(id_obs_i,1));    % find potential map index (SLAM id), remove duplicated index
[~, ~, id_pts_i] = intersect(index_pts_i, map_ids, 'stable'); 

times_obs(id_pts_i) = times_obs(id_pts_i) + 1;   

% id_obs = id_obs | id_obs_i;    % here, I consider all map points. 
% 
% 
% times_obs(id_obs_i) = times_obs(id_obs_i) + 1;
end


% select map points that are observed among all num_covi frames
id_valid_pts = times_obs >= num_covi;       % the id in vector of valid map points
index_valid_pts = map_ids(id_valid_pts);    % the SLAM index of valid map points

% retrive valid observation among all keyframes
% Id_obs = false(num_obs,1);    % initialize valid observation id. 
Id_obs_i = cell(num_BAframes,1);
for i=1:num_BAframes

if(k-i<0)       % only consider valid frame id 
    break;
end

id_obs_i = matched_feat(:,2) == (k-i);  

% select the valid map points
for j=1:length(id_obs_i)

    if(id_obs_i(j))
        index_pts_ij = matched_feat(j,1);
        b_index_valid_pts = index_valid_pts== index_pts_ij;
        if(sum(b_index_valid_pts)==0)   % if this map points are valid, the sum should be larger than 1.
            id_obs_i(j) = 0;     % if this map points are not in valid index, we remove it
        end
    end


% index_pts_i = matched_feat(id_obs_i,1); 

end
Id_obs_i{num_BAframes-(i-1)} = id_obs_i;   % the retrived keyframe is in reverse order, we need to re-order it
% Id_obs = Id_obs | id_obs_i;    % here, I consider all observations of valid map points.
end




%% conduct localBA optimization

% initial map points
p_a0_ia0 = f_a0_ia0(id_valid_pts,:);    % initial map points
% p_a0_ia0 = f_a0_ia0;

% get initial pose and intrinsic matrix
% vT_0k = cell(num_frame, 1); 
% for k=1:num_frame
% vT_0k{k} = squeeze(slam_vTu0uk(k,:,:));
% end



vT_uo0uok = cell(num_BAframes, 1);      % use cross-view pose
frame_ids = zeros(num_BAframes,1);
for k=1:num_BAframes
kk = pred_M_uok_usk{3}-(k-1);
% T_us0_usk = squeeze(slam_vTu0uk(kk,:,:));
T_us0_usk = squeeze(vT_us0_usk(kk,:,:));   % what we want to estimate is T_o0ok, not T_s0sk

% the retrived keyframe is in reverse order, we need to re-order it
vT_uo0uok{num_BAframes-(k-1)} = T_us0_usk;    % initialize using SLAM pose
frame_ids(num_BAframes-(k-1)) = pred_M_uok_usk{3} - k;    % the frame id starts from 0 in python
end

% for visualization
% vT_uo0uok_init = vT_uo0uok;
% vT_uo0uok_gt = cell(num_BAframes, 1);
% for k=1:num_BAframes
% kk = pred_M_uok_usk{3}-(k-1);
% T_uo0_uok = squeeze(oxts_vTu0uk(kk,:,:));
% 
% % the retrived keyframe is in reverse order, we need to re-order it
% vT_uo0uok_gt{num_BAframes-(k-1)} = T_uo0_uok;    % initialize using SLAM pose
% end

% for k=1:num_BAframes
% 
% kk = pred_M_uok_usk{3}-(k-1);
% 
% T_us0_usk = squeeze(slam_vTu0uk(kk,:,:));
% T_a0u0 = T_u0a0 \ eye(4);
% T_a0_ak =  T_a0u0*T_us0_usk*T_u0a0;
% T_ak_a0 = T_a0_ak \ eye(4);
% 
% vT_k0{k} = T_ak_a0;
% frame_ids(k) = pred_M_uok_usk{3} - k;    % the frame id starts from 0 in python
% end

% project map points
% vK = cell(2,1);
% vK{2} = P_b0a0;    % front right
% vK{1} = P_b0a0; K0(1,4) = 0;   % front left


% weights
% w = [1, 1 0.5, 0.5];    %[2d feature weights, CrossView rot, CrossView translation]
w = [1, 1, 10, 10];
% GN iteration
iter=0;
flag = true;
while(flag)

    [F] = calculate_Error(vT_uo0uok, T_u0a0, frame_ids, p_a0_ia0, index_valid_pts, sum(id_valid_pts), matched_feat, Id_obs_i, P_b0a0,  pred_M_uok_usk, w, keyframe_name);
    [J] = calculate_Jacobian(vT_uo0uok, T_u0a0, frame_ids, p_a0_ia0, index_valid_pts, sum(id_valid_pts), matched_feat, Id_obs_i, P_b0a0, pred_M_uok_usk, w, keyframe_name);
    
%     A = full(J'*J);
%     S = svd(A);
%     [V,D,W] = eig(A);
    dX = -(J'*J)\(J'*F);

    % update 
    [vT_uo0uok, p_a0_ia0] = update_state(dX, vT_uo0uok, p_a0_ia0);


    % check stopping points
    aa_dx = dX'*dX;
    aa_FX = F'*F;
    iter = iter+1;
    fprintf('iter, dx, Fx:  %d,  %f,  %f\n', iter, aa_dx, aa_FX)
    if (aa_dx < 1e-3 || iter>=20 )
        flag = false;
    end

end



%% update pose
% initial X
vT_0k = cell(num_frame, 1); 
for k=1:num_frame
% vT_0k{k} = squeeze(slam_vTu0uk(k,:,:));   % what we want to estimate is T_o0ok, not T_s0sk
vT_0k{k} = squeeze(vT_us0_usk(k,:,:));   % what we want to estimate is T_o0ok, not T_s0sk
end

% pose-graph measurement
idcv = [num_BAframes];                  % select the last pose of localBA (the current pose to be used)
fine_M_w_uok = cell(length(idcv),2);
for k=1:length(idcv)
fine_M_w_uok{k,1} = vT_uo0uok{idcv(k)};
fine_M_w_uok{k,2} = frame_ids(idcv(k));
end


% weights
w = [1,0.9,0.04,0.002];   %[VO rot, VO translation, CrossView rot, CrossView translation]
% GN iteration
iter=0;
flag = true;
while(flag)

    [F] = calculate_ErrorPG(vT_0k, vo_vM_usl_usk, fine_M_w_uok, w);
    [J] = calculate_JacobianPG(vT_0k, vo_vM_usl_usk, fine_M_w_uok, w);
    
%     A = full(J'*J);
%     S = svd(A);
%     [V,D,W] = eig(A);
    dX = -(J'*J)\(J'*F);

    % update 
    [vT_0k] = update_statePG(dX, vT_0k);


    % check stopping points
    aa_dx = dX'*dX;
    aa_FX = F'*F;
    iter = iter+1;
    fprintf('iter, dx, Fx:  %d,  %f,  %f\n', iter, aa_dx, aa_FX)
    if (aa_dx < 1e-3 || iter>=5 )
        flag = false;
    end

end













%% conduct global pose-graph optimization
% initial X
vT_0k = cell(num_frame, 1); 
for k=1:num_frame
vT_0k{k} = squeeze(vT_us0_usk(k,:,:));   % what we want to estimate is T_o0ok, not T_s0sk
end

w = [1,0.9,0.04,0.002];   %weights  [VO rot, VO translation, CrossView rot, CrossView translation]
% GN iteration
iter=0;
flag = true;
while(flag)

    [F] = calculate_Error(vT_0k, vo_vM_usl_usk, pred_vM_uok_usk, w);
    [J] = calculate_Jacobian(vT_0k, vo_vM_usl_usk, pred_vM_uok_usk, w);
    
%     A = full(J'*J);
%     S = svd(A);
%     [V,D,W] = eig(A);
    dX = -(J'*J)\(J'*F);

    % update 
    [vT_0k] = update_state(dX, vT_0k);


    % check stopping points
    aa_dx = dX'*dX;
    aa_FX = F'*F;
    iter = iter+1;
    fprintf('iter, dx, Fx:  %d,  %f,  %f\n', iter, aa_dx, aa_FX)
    if (aa_dx < 1e-3 || iter>=5 )
        flag = false;
    end

end




%% update vM_uok_usk (vM_uok_usk is based on vT_us0_usk before updates )
for k=1:num_frame

if k>idx+1
    break
end

M_uok_usk = squeeze(vM_uok_usk(k,:,:));   % cross-view prediction before update
T_us0_usk = squeeze(vT_us0_usk(k,:,:));   % trajectory before update
T_us0_uuk = vT_0k{k};                     % trajectory after update

T_usk_us0 = T_us0_usk \ eye(4);

M_uok_uuk = M_uok_usk * T_usk_us0 * T_us0_uuk;

% save result
vM_uok_usk(k,:,:) = M_uok_uuk;

end




%% update vT_us0_usk
for k=1:num_frame
vT_us0_usk(k,:,:) = vT_0k{k};
end



% fprintf('job done \n')

% a=fine_valid;



end