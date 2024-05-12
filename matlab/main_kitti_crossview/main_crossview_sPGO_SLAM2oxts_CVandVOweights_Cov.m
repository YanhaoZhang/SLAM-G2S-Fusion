function [vT_us0_usk, vS_0k, vM_usk_uok, Cov_circle] = main_crossview_sPGO_SLAM2oxts_CVandVOweights_Cov(vT_us0_usk_org, L, vT_us0_usk, vS_0k, fine_valid, vM_usk_uok, vP_usk_uok, vW_usk_uok, vMask_usk_uok, idx)
% vT_us0_usk_org: original SLAM trajectory to calculate visual odometry
% L: laplacian matrix whose elements are feature number of covisable frames
% vT_us0_usk: SLAM refined and to be refined
% vS_0k: scale of each pose
% vT_uo0_uok: oxts trajectory (ground truth just for evaluation)
% T_u0a0: cam0 to cam2. The map points are in cam0 frame0
% P_b0a0: projection matrix (front stereo)
% f_a0_ia0: xyz map points
% map_ids: map point index by ORB-SLAM
% fine_valid: valid [lat, lon, yaw] by cross-view
% vM_usk_uok: store all valid SLAM to oxts prediction
% vP_usk_uok: store all cross-view prediction weights


% notice: we need to update pred_vM_uok_usk everytime when we update the
% trajectory. The prediction is based on the trajectory before refinement


%% rewrite visual odometry
num_frame = size(vT_us0_usk,1);

% % visual odometry measurement
% vo_vM_usl_usk = cell(num_frame-1, 3); 
% 
% id = 0;
% for k=2:num_frame
% 
% T_us0_usl = squeeze(vT_us0_usk_org(k-1,:,:));
% T_us0_usk = squeeze(vT_us0_usk_org(k,:,:));
% 
% % T_us0_usl = squeeze(vT_us0_usk(k-1,:,:));
% % T_us0_usk = squeeze(vT_us0_usk(k,:,:));
% 
% T_usl_us0 = T_us0_usl \ eye(4);
% 
% 
% % store data
% id=id+1;
% vo_vM_usl_usk{id, 1} = T_usl_us0 * T_us0_usk;
% vo_vM_usl_usk{id, 2} = k-1;
% vo_vM_usl_usk{id, 3} = k;
% 
% end








% visual odometry measurement
% vo_vM_usl_usk = cell(num_frame*10, 5); 
% weight_vo = zeros(num_frame*10, 1); 
% 
% id = 0;
% num_localEdge = 10;   % skip redundant edges within map
% for k=2:num_frame
% 
% [~, id_edge_k, val_edge_k] = find(L(k-1,:));
% l_next = k-1;
% for idl=1:length(id_edge_k)
% l = id_edge_k(idl);
% 
% if(l>l_next)
% T_us0_usl = squeeze(vT_us0_usk_org(k-1,:,:));
% T_us0_usk = squeeze(vT_us0_usk_org(l,:,:));
% T_usl_us0 = T_us0_usl \ eye(4);
% 
% 
% % store data
% id=id+1;
% vo_vM_usl_usk{id, 1} = T_usl_us0 * T_us0_usk;
% vo_vM_usl_usk{id, 2} = k-1;
% vo_vM_usl_usk{id, 3} = l;
% vo_vM_usl_usk{id, 4} = val_edge_k(idl);
% weight_vo(id) = val_edge_k(idl);
% 
% % next l
% l_next = l+num_localEdge;   % 
% 
% end
% end
% end
% 
% % normalize the weights
% vo_vM_usl_usk = vo_vM_usl_usk(1:id,:);
% weight_vo = weight_vo(1:id,:);
% weight_vo_normalize = weight_vo ./ mean(weight_vo);
% 
% % rewrite the normalized weights
% for k=1:size(vo_vM_usl_usk,1)
% vo_vM_usl_usk{k, 5} = sqrt(weight_vo_normalize(k));
% end


% rewrite laplacian matrix by removing some redundant edges
L_ = logical(L.*(L>20));


for k=2:num_frame
L_(k-1,k) = 1;
end



num_localEdge = num_frame/3;   % skip redundant edges within map
for k=2:num_frame

[~, id_edge_k, val_edge_k] = find(L_(k-1,:));
num_edge_k = length(id_edge_k);

% l_next = k-1;
for idl=1:num_edge_k

    % index and if this element is connected
l = id_edge_k(idl);
val = val_edge_k(idl);

if(val)

    if(l==k)
        L_(k-1,l+(1:num_localEdge)) = 0;
        val_edge_k = L_(k-1,id_edge_k);

    end

    if(l>k+num_localEdge)

        L_((k-2)+(1:num_localEdge),(l-1)+(1:num_localEdge)) = 0;
        L_(k-1,l) = 1;
        val_edge_k = L_(k-1,id_edge_k);
    end

end
end
end

L_ = L_(1:num_frame,1:num_frame);



% visual odometry measurement
vo_vM_usl_usk = cell(num_frame*10, 5); 
weight_vo = zeros(num_frame*10, 1); 

id = 0;
for k=2:num_frame
[~, id_edge_k, ~] = find(L_(k-1,:));
val_edge_k = L(k-1,id_edge_k);

num_edge_k = length(id_edge_k);
for idl=1:num_edge_k
l = id_edge_k(idl);

T_us0_usl = squeeze(vT_us0_usk_org(k-1,:,:));
T_us0_usk = squeeze(vT_us0_usk_org(l,:,:));
T_usl_us0 = T_us0_usl \ eye(4);

% store data
id=id+1;
vo_vM_usl_usk{id, 1} = T_usl_us0 * T_us0_usk;
vo_vM_usl_usk{id, 2} = k-1;
vo_vM_usl_usk{id, 3} = l;
vo_vM_usl_usk{id, 4} = val_edge_k(idl);
weight_vo(id) = val_edge_k(idl);

end
end

vo_vM_usl_usk = vo_vM_usl_usk(1:id,:);
weight_vo = weight_vo(1:id,:);


weight_vo = sqrt(weight_vo);
weight_vo(weight_vo==0)=1;   % sometimes, there will be no edge between two frames

weight_vo_normalize = weight_vo ./ mean(weight_vo);

% weight_vo_normalize(weight_vo_normalize==0)=1;    


for k=1:size(vo_vM_usl_usk,1)
% vo_vM_usl_usk{k, 5} = sqrt(weight_vo_normalize(k));
vo_vM_usl_usk{k, 5} = weight_vo_normalize(k);
end








%% load cross-view measurement
valid_id = fine_valid(:,1) & fine_valid(:,2);    % consider all [lon, lat, yaw] valid as valid cross-view

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
pred_vM_usk_uok{id, 2} = T_us0_usk;    % the previous updated trajectory pose is actually used
pred_vM_usk_uok{id, 3} = k;

% pred_vM_uok_usk{id, 4} = [weights_fine_(k), weights_fine_(k)];    % SLAM pose is actually used



% save weight
corr_p = squeeze(vP_usk_uok(k,:,:));
corr_p_ = corr_p(:);
TF_p_ = islocalmax(corr_p_);
max_p_potential = corr_p_(TF_p_);
w_p = max(max_p_potential)/sum(max_p_potential);




% new weight
% corr_p = squeeze(vP_usk_uok(k,:,:));
% w_k = squeeze(vW_usk_uok(k,:,:));
% corr_m = squeeze(vMask_usk_uok(k,:,:));
% corr_p2 = corr_p .* corr_m;      %todo need to check if we use the masked one or original one
% corr_p2_ = corr_p2(:);
% 
% TF_p_ = islocalmax(corr_p2_);
% max_p_potential = corr_p2_(TF_p_);
% % check if the selected corrolation is the maximum one
% if abs(max(max_p_potential)-w_k)>1e-7
%     max_p_potential = [max_p_potential; w_k];
% end
% w_p = w_k/sum(max_p_potential);


w_k = squeeze(vW_usk_uok(k,:,:));
w_p = w_p * w_k;

weights_fine(id) = w_p;
end


% reweight
weights_fine = sqrt(weights_fine);
w_normalize = mean(weights_fine);
weights_fine_ = weights_fine/w_normalize;
for k=1:num_valid
pred_vM_usk_uok{k, 4} = [weights_fine_(k), weights_fine_(k)]; % [1,1]; % sqrt([weights_fine_(k), weights_fine_(k)]); %
end



%% oxts rotation
valid_id_rot = fine_valid(:,3);    % consider all [lyaw] valid as valid cross-view

num_valid_rot = sum(valid_id_rot);

% num_valid_rot = double(idx+1);

pred_vR_usk_uok = cell(num_valid_rot, 3);      % relative rotation from cross-view, oxts pose, id

id = 0;
for k=1:num_frame

if (~valid_id_rot(k))
    continue
end

if k>idx+1
    break
end

id = id + 1;

% directly use the updated results from python (w.r.t. the new trajectory)
R_usk_uok = squeeze(vM_usk_uok(k,1:3,1:3));   
R_us0_usk = squeeze(vT_us0_usk(k,1:3,1:3));     % M_uok_usk is w.r.t. the updated trajectory 


% store predicted pose
pred_vR_usk_uok{id, 1} = R_usk_uok;    % new cross-view pose: SLAM to OXTS
pred_vR_usk_uok{id, 2} = R_us0_usk;    % the previous updated trajectory pose is actually used
pred_vR_usk_uok{id, 3} = k;

end





% save('/home/users/u1111398/YanhaoZhang/ORB_SLAM3_stereoMultiCam/matlab/debug_.mat', 'pred_vM_uok_usk', 'num_valid', 'valid_id')


% divided the number of measurement
% nVO = size(vo_vM_usl_usk,1)/1000;
% nCV = size(pred_vM_usk_uok,1)/1;
% nCV_r = size(pred_vR_usk_uok,1)/1;

nVO = 1;
nCV = 1;
nCV_r = 1;


%% initial X for GN or cross-view
vT_0k = cell(num_frame, 1); 
for k=1:num_frame
vT_0k{k} = squeeze(vT_us0_usk(k,:,:));   % what we want to estimate is T_o0ok, not T_s0sk
end

% w = [1, 0.9, 0.04, 0.002, 1];   %weights  [VO rot, VO translation, CrossView rot, CrossView translation, VO scale]
% w = [1, 0.9, 1, 0.002, 1]; 
w = [0.85/nVO, 0.9/nVO, 1.0/nCV_r, 0.003/nCV, 0.005/nCV, 0.000, 10, 1.0];    %[VO rot, VO translation, CrossView rot, CrossView lon, cross-view lat, cross-view hight, VO scale, huber threshold]

%% if only calculate the covariance for the initial frame (before trajectory updated)

% for visualization
t_w_fkw = zeros(num_frame,3);
for k=1:num_frame
t_w_fkw(k,:) = vT_0k{k}(1:3,4);
end

% % figure
% subplot(1,2,1)
% plot(t_w_fkw(:,1), t_w_fkw(:,2), '.', 'Color',[0.6350 0.0780 0.1840]); hold on
% axis equal


num_cov_pts = 40;
% cov_next = zeros(2,2);
Cov_circle = zeros(num_frame, 2, num_cov_pts+1);
% Cov_circle = zeros(num_frame, 2, num_cov_pts+1);
% meter_per_pixel = 0.19582850865299833;
if idx==0
    [J] = calculate_Jacobian_scale_huber(vT_0k, vS_0k, vo_vM_usl_usk, pred_vR_usk_uok, pred_vM_usk_uok, w);


    % % calculate covariance of next frame
    % k=double(idx)+2;                       % next frame id in matlab
    % [Cov_scale, normalzier] = covEst_scale(J, k);
    % f_0_i0 = GetCov2d_crossview(Cov_scale{1},   vT_0k{k}(1,4),  vT_0k{k}(2,4), normalzier*1, num_cov_pts);
    % 
    % % transform to local frame
    % [~, ~, yaw_0k] = euler_from_rotation_matrix(vT_0k{k}(1:3, 1:3));
    % R_0k = [cos(yaw_0k), -sin(yaw_0k); sin(yaw_0k), cos(yaw_0k)];
    % t_0_k0 = [vT_0k{k}(1,4);  vT_0k{k}(2,4)];
    % 
    % f_k_ik = R_0k'*(f_0_i0-t_0_k0);
    % 
    % cov_next = f_k_ik;


    [Cov_scale, normalzier] = covEst_scale(J);
    for k=1:num_frame

        f_0_i0 = GetCov2d_crossview(Cov_scale{k},   vT_0k{k}(1,4),  vT_0k{k}(2,4), normalzier*1, num_cov_pts);


        % transform to local frame
        [~, ~, yaw_0k] = euler_from_rotation_matrix(vT_0k{k}(1:3, 1:3));
        R_0k = [cos(yaw_0k), -sin(yaw_0k); sin(yaw_0k), cos(yaw_0k)];
        t_0_k0 = [vT_0k{k}(1,4);  vT_0k{k}(2,4)];

        f_k_ik = R_0k'*(f_0_i0-t_0_k0);

        Cov_circle(k,:,:) = f_k_ik;

        

        % % visualize result
        % meter_per_pixel = 0.19582850865299833;
        % f_k_ik = f_k_ik/meter_per_pixel;
        % f_k_ik(1,:) = -f_k_ik(1,:);
        % 
        % subplot(1,2,1)
        % plot(f_0_i0(1,:),f_0_i0(2,:),'Color','k','LineStyle','-','LineWidth',0.05); hold on;
        % axis equal
        % 
        % subplot(1,2,2)
        % bw_size = 512;
        % f_k_ik = f_k_ik+bw_size/2;
        % BW = poly2mask(f_k_ik(1,:),f_k_ik(2,:),bw_size,bw_size);
        % imshow(BW); hold on
        % 
        % plot(0, 0, '.', 'Color',[0.6350 0.0780 0.1840]); hold on
        % plot(f_k_ik(1,:),f_k_ik(2,:),'r.','MarkerSize',2); hold off;
        % axis equal


    end
    return;
end




%% conduct global pose-graph optimization
% GN iteration
iter=0;
flag = true;
min_FX_old = 0;
while(flag)


    % scale weight huber
    [F] = calculate_Error_scale_huber(vT_0k, vS_0k, vo_vM_usl_usk, pred_vR_usk_uok, pred_vM_usk_uok, w);
    [J] = calculate_Jacobian_scale_huber(vT_0k, vS_0k, vo_vM_usl_usk, pred_vR_usk_uok, pred_vM_usk_uok, w);
    % % 
    % % % update 
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



%% calculate covariance
[J] = calculate_Jacobian_scale_huber(vT_0k, vS_0k, vo_vM_usl_usk, pred_vR_usk_uok, pred_vM_usk_uok, w);

% calculate covariance of next frame
% k=double(idx)+2;                       % next frame id in matlab
% 
% f_0_i0 = GetCov2d_crossview(Cov_scale{1},   vT_0k{k}(1,4),  vT_0k{k}(2,4), normalzier*1, num_cov_pts);

[Cov_scale, normalzier] = covEst_scale(J);
for k=1:num_frame

    if k<idx+1
        continue;   % do not consider original frame since they are not used
    end

    % current estimated pose
    [~, ~, yaw_0k] = euler_from_rotation_matrix(vT_0k{k}(1:3, 1:3));
    R_0k = [cos(yaw_0k), -sin(yaw_0k); sin(yaw_0k), cos(yaw_0k)];
    t_0_k0 = [vT_0k{k}(1,4);  vT_0k{k}(2,4)];

    % covariance circle and transform to local frame
    f_0_i0 = GetCov2d_crossview(Cov_scale{k},   t_0_k0(1),  t_0_k0(2), normalzier*1, num_cov_pts);
    f_k_ik = R_0k'*(f_0_i0-t_0_k0);

    % save data
    Cov_circle(k,:,:) = f_k_ik;
end

% if idx<=num_frame-2
%     [J] = calculate_Jacobian_scale_huber(vT_0k, vS_0k, vo_vM_usl_usk, pred_vR_usk_uok, pred_vM_usk_uok, w);
% 
% 
%     % calculate covariance of next frame
%     k=double(idx)+2;                       % next frame id in matlab
%     [Cov_scale, normalzier] = covEst_scale(J, k);
%     f_0_i0 = GetCov2d_crossview(Cov_scale{1},   vT_0k{k}(1,4),  vT_0k{k}(2,4), normalzier*1, num_cov_pts);
% 
%     % transform to local frame
%     [~, ~, yaw_0k] = euler_from_rotation_matrix(vT_0k{k}(1:3, 1:3));
%     R_0k = [cos(yaw_0k), -sin(yaw_0k); sin(yaw_0k), cos(yaw_0k)];
%     t_0_k0 = [vT_0k{k}(1,4);  vT_0k{k}(2,4)];
% 
%     f_k_ik = R_0k'*(f_0_i0-t_0_k0);
% 
%     cov_next = f_k_ik;
% 
% 
% end
end