function [vT_us0_usk, vS_0k, vM_usk_uok] = main_crossview_scalePoseGraph_SLAM2oxts_weights(vT_us0_usk_org, vT_us0_usk, vS_0k, fine_valid, vM_usk_uok, vP_usk_uok, idx)
% vT_us0_usk_org: original SLAM trajectory to calculate visual odometry
% vT_us0_usk: SLAM refined and to be refined
% vS_0k: scale of each pose
% vT_uo0_uok: oxts trajectory (ground truth just for evaluation)
% T_u0a0: cam0 to cam2. The map points are in cam0 frame0
% P_b0a0: projection matrix (front stereo)
% f_a0_ia0: xyz map points
% map_ids: map point index by ORB-SLAM
% fine_valid: valid [lat, lon, yaw] by cross-view
% vM_usk_uok: store all valid SLAM to oxts prediction

% notice: we need to update pred_vM_uok_usk everytime when we update the
% trajectory. The prediction is based on the trajectory before refinement


%% rewrite visual odometry
num_frame = size(vT_us0_usk,1);

% visual odometry measurement
vo_vM_usl_usk = cell(num_frame-1, 3); 

id = 0;
for k=2:num_frame

T_us0_usl = squeeze(vT_us0_usk_org(k-1,:,:));
T_us0_usk = squeeze(vT_us0_usk_org(k,:,:));

% T_us0_usl = squeeze(vT_us0_usk(k-1,:,:));
% T_us0_usk = squeeze(vT_us0_usk(k,:,:));

T_usl_us0 = T_us0_usl \ eye(4);


% store data
id=id+1;
vo_vM_usl_usk{id, 1} = T_usl_us0 * T_us0_usk;
vo_vM_usl_usk{id, 2} = k-1;
vo_vM_usl_usk{id, 3} = k;

end


%% load cross-view measurement
valid_id = fine_valid(:,1) & fine_valid(:,2) & fine_valid(:,3);    % consider all [lon, lat, yaw] valid as valid cross-view

% valid_id = fine_valid(:,2) & fine_valid(:,3);    % only use [ lat, yaw]



num_valid = sum(valid_id);
pred_vM_usk_uok = cell(num_valid, 4);      % relative pose from cross-view, oxts pose, id

weights_fine = zeros(num_valid,1);

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
% T_us0_usk = squeeze(vT_us0_usk(k,:,:));

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



% directly use the updated results from python (w.r.t. the new trajectory)
M_usk_uok = squeeze(vM_usk_uok(k,:,:));   
T_us0_usk = squeeze(vT_us0_usk(k,:,:));     % M_uok_usk is w.r.t. the updated trajectory 



% M_uok_usk = squeeze(vM_uok_usk(k,:,:)); 

% store predicted pose
pred_vM_usk_uok{id, 1} = M_usk_uok;    % new cross-view pose: SLAM to OXTS
pred_vM_usk_uok{id, 2} = T_us0_usk;    % SLAM pose is actually used
pred_vM_usk_uok{id, 3} = k;

% pred_vM_uok_usk{id, 4} = [weights_fine_(k), weights_fine_(k)];    % SLAM pose is actually used



% save weight
corr_p = squeeze(vP_usk_uok(k,:,:));
corr_p_ = corr_p(:);
TF_p_ = islocalmax(corr_p_);
max_p_potential = corr_p_(TF_p_);
w_p = max(max_p_potential)/sum(max_p_potential);

weights_fine(id) = w_p;
end


% reweight
w_normalize = mode(weights_fine);
weights_fine_ = weights_fine/w_normalize;
for k=1:num_valid
pred_vM_usk_uok{k, 4} = [weights_fine_(k), weights_fine_(k)];
end



% save('/home/users/u1111398/YanhaoZhang/ORB_SLAM3_stereoMultiCam/matlab/debug_.mat', 'pred_vM_uok_usk', 'num_valid', 'valid_id')

%% conduct global pose-graph optimization
% initial X
vT_0k = cell(num_frame, 1); 
for k=1:num_frame
vT_0k{k} = squeeze(vT_us0_usk(k,:,:));   % what we want to estimate is T_o0ok, not T_s0sk
end

% w = [1, 0.9, 0.04, 0.002, 1];   %weights  [VO rot, VO translation, CrossView rot, CrossView translation, VO scale]
% w = [1, 0.9, 1, 0.002, 1]; 
w = [1, 0.9, 1.0, 0.004, 0.005, 0.000, 10, 1.0];    %[VO rot, VO translation, CrossView rot, CrossView lon, cross-view lat, cross-view hight, VO scale]
% GN iteration
iter=0;
flag = true;
min_FX_old = 0;
while(flag)


    % scale weight huber
    [F] = calculate_Error_scale_huber(vT_0k, vS_0k, vo_vM_usl_usk, pred_vM_usk_uok, w);
    [J] = calculate_Jacobian_scale_huber(vT_0k, vS_0k, vo_vM_usl_usk, pred_vM_usk_uok, w);
    % 
    % % update 
    dX = -(J'*J)\(J'*F);
    [vT_0k, vS_0k] = update_state_scale(dX, vT_0k, vS_0k);

    % scale weight nohuber
    % [F] = calculate_Error_scale(vT_0k, vS_0k, vo_vM_usl_usk, pred_vM_usk_uok, w);
    % [J] = calculate_Jacobian_scale(vT_0k, vS_0k, vo_vM_usl_usk, pred_vM_usk_uok, w);
    % 
    % % update 
    % dX = -(J'*J)\(J'*F);
    % [vT_0k, vS_0k] = update_state_scale(dX, vT_0k, vS_0k);


    % % no scale no weight no huber
    % [F] = calculate_Error_noScale_huber(vT_0k, vo_vM_usl_usk, pred_vM_usk_uok, w);
    % [J] = calculate_Jacobian_noScale_huber(vT_0k, vo_vM_usl_usk, pred_vM_usk_uok, w);

    % update 
    % dX = -(J'*J)\(J'*F);
    % [vT_0k] = update_state_noScale(dX, vT_0k);


    % check stopping points
    aa_dx = dX'*dX;
    min_FX = F'*F;
    iter = iter+1;
    % fprintf('iter, dx, Fx:  %d,  %f,  %f\n', iter, aa_dx, min_FX)
    if (aa_dx < 1e-3 || iter>=5 || abs(min_FX-min_FX_old)<1e-6 )
        flag = false;
    end
    min_FX_old = min_FX;

end




%% update vM_uok_usk (vM_uok_usk is based on vT_us0_usk before updates )
for k=1:num_frame

if k>idx+1
    break
end

M_usk_uok = squeeze(vM_usk_uok(k,:,:));   % cross-view prediction before update   SLAM to OXTS
T_us0_usk = squeeze(vT_us0_usk(k,:,:));   % SLAM trajectory before update
T_us0_uuk = vT_0k{k};                     % SLAM trajectory after update

T_uuk_us0 = T_us0_uuk \ eye(4);

M_uuk_uok = T_uuk_us0 * T_us0_usk * M_usk_uok;

% save result
vM_usk_uok(k,:,:) = M_uuk_uok;

end



% for k=1:num_frame
% 
% if k>idx+1
%     break
% end
% 
% % for the valid cross-view id, we use cross-view measurement
% if(valid_id(k))
% 
% M_uok_usk = squeeze(vM_uok_usk(k,:,:));   % cross-view prediction before update
% T_us0_usk = squeeze(vT_us0_usk(k,:,:));   % trajectory before update
% T_us0_uuk = vT_0k{k};                     % trajectory after update
% 
% T_usk_us0 = T_us0_usk \ eye(4);
% 
% M_uok_uuk = M_uok_usk * T_usk_us0 * T_us0_uuk;
% 
% % save result
% vM_uok_usk(k,:,:) = M_uok_uuk;
% 
% else
% % for the non-valid cross-view measurement, we directly use the relative
% % change between original SLAM pose to the actual SLAM pose
% 
% 
% 
% end
% 
% 
% end



%% update vT_us0_usk
for k=1:num_frame
vT_us0_usk(k,:,:) = vT_0k{k};
end




% fprintf('job done \n')

% a=fine_valid;



end