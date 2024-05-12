function [d_cv, rmse_all, rmse_vld] = show_crossview_ablation_kitti(vT_0k, slam_vTu0uk, oxts_vTu0uk, coarse_vM_ufk_upk, valid_id)


num_frame = size(slam_vTu0uk,1);
% valid_id = fine_valid(:,1) & fine_valid(:,2) & fine_valid(:,3);    % consider all [lon, lat, yaw] valid as valid cross-view
% valid_id = fine_valid(:,2) & fine_valid(:,3);                    % consider only [lat, yaw] valid as valid cross-view
num_valid = sum(valid_id)

% a2t=[];

remove_id = ~valid_id;


%% rewrite refined SLAM result
t_w_skw = zeros(num_frame,3);
t_w_okw = zeros(num_frame,3);
% t_w_fkw = zeros(num_frame,3);
for k=1:num_frame

t_w_skw(k,:) = slam_vTu0uk(k,1:3,4);
t_w_okw(k,:) = oxts_vTu0uk(k,1:3,4);
% t_w_fkw(k,:) = vT_0k(k,1:3,4);
end



%% rewrite valid cross-view prediction
% diff_t_w_okpk_cv = zeros(num_frame, 3);
t_w_pkw = zeros(num_frame,3);
id = 0;
for k=1:num_frame

id = id + 1;

T_us0_usk = squeeze(slam_vTu0uk(k,:,:));
M_usk_upk = squeeze(coarse_vM_ufk_upk(k,:,:));
M_us0_upk = T_us0_usk * M_usk_upk;

% t_w_skw_ = (T_us0_usk(1:3,4))';
t_w_pkw(k,:) = (M_us0_upk(1:3,4))';

% diff_t_w_okpk_cv(id,:) = p_w_pkw(k,:) - t_w_skw_;
end


% rewrite valid gt of cross-view prediction and refined prediction
diff_t_w_oksk = t_w_okw - t_w_skw;
diff_t_w_pksk = t_w_pkw - t_w_skw;

% diff_fine = t_w_fkw - t_w_okw; 



%% show results
if 0
figure

h_slam = plot(t_w_skw(:,1), t_w_skw(:,2), 'b.','MarkerSize', 2); hold on
h_gt = plot(t_w_okw(:,1), t_w_okw(:,2), '.', 'Color', [0.4660 0.6740 0.1880],'MarkerSize', 2); hold on
% h_cv = plot(p_w_pkw(:,1), p_w_pkw(:,2), '.', 'Color',[0.6350 0.0780 0.1840],'MarkerSize', 2); hold on

h_cv_gt  = quiver(t_w_skw(:,1), t_w_skw(:,2), diff_t_w_oksk(:,1),    diff_t_w_oksk(:,2), 'Color', 'k', 'AutoScale', 'off'); hold on
h_cv_raw = quiver(t_w_skw(remove_id,1), t_w_skw(remove_id,2), diff_t_w_pksk(remove_id,1),    diff_t_w_pksk(remove_id,2), 'Color', 'g', 'AutoScale', 'off'); hold on
h_cv_vld = quiver(t_w_skw(valid_id,1), t_w_skw(valid_id,2), diff_t_w_pksk(valid_id,1),    diff_t_w_pksk(valid_id,2), 'Color', 'm', 'AutoScale', 'off'); hold on
% h_cv_fine = quiver(t_w_okw(:,1), t_w_okw(:,2), diff_fine(:,1), diff_fine(:,2), 'Color', 'm', 'AutoScale', 'off'); hold on
% legend([h_slam, h_gt, h_fine, h_cv_fine],{'SLAM', 'GT', 'SLAM+CrossView', 'diff fine'})
% legend([h_slam, h_gt, h_cv_gt, h_cv_pred],{'SLAM traj.', 'GT traj.', 'diff GT', 'used CV pred'})
legend([h_slam, h_gt, h_cv_gt, h_cv_raw, h_cv_vld],{'SLAM Traj', 'GT Traj', 'GT G2S', 'Raw G2S', 'Selected G2S'})
% legend([h_slam, h_gt, h_cv],{'SLAM', 'GT', 'CrossView'})
axis equal

end



%% calculate lon, lat, yaw error:
err_slam = zeros(num_frame, 3);
err_cv = zeros(num_frame, 3);
% err_fine = zeros(num_frame,3);

t_w_pkw = zeros(num_frame,3);     % store cross-view prediction
vT_w_pk = cell(num_frame,1);
for k=1:num_frame

T_w_sk = squeeze(slam_vTu0uk(k,:,:));    % SLAM
T_w_ok = squeeze(oxts_vTu0uk(k,:,:));   % gt
% T_w_fk = squeeze(vT_0k(k,:,:));                      % fine trajectory


T_ok_sk = (T_w_ok \ eye(4)) * T_w_sk;
% T_ok_fk = (T_w_ok \ eye(4)) * T_w_fk;


% cross-view trajectory
% T_us0_usk = squeeze(vT_us0_usk(k,:,:));
T_sk_pk = squeeze(coarse_vM_ufk_upk(k,:,:));
T_w_pk = T_w_sk * T_sk_pk;
T_ok_pk = (T_w_ok \ eye(4)) * T_w_pk;

% store cross-view trajectory
t_w_pkw(k,:) = (T_w_pk(1:3,4))';
vT_w_pk{k} = T_w_pk;


% rewrite to lon lat and yaw error
[~, ~, yaw_ok_sk] = euler_from_rotation_matrix(T_ok_sk(1:3, 1:3));
% [~, ~, yaw_ok_fk] = euler_from_rotation_matrix(T_ok_fk(1:3, 1:3));
[~, ~, yaw_ok_pk] = euler_from_rotation_matrix(T_ok_pk(1:3, 1:3));



yaw_ok_sk_ = rad2deg(wrapToPi(yaw_ok_sk));
% yaw_ok_fk_ = rad2deg(wrapToPi(yaw_ok_fk));
yaw_ok_pk_ = rad2deg(wrapToPi(yaw_ok_pk));


err_slam(k,:) = abs([T_ok_sk(1,4), T_ok_sk(2,4), yaw_ok_sk_]);
% err_fine(k,:) = abs([T_ok_fk(1,4), T_ok_fk(2,4), yaw_ok_fk_]);
err_cv(k,:) = abs([T_ok_pk(1,4), T_ok_pk(2,4), yaw_ok_pk_]);
end


metrics = [1, 3, 5];
angles = [1, 3, 5];









% print error
% fprintf('mean lon, mean lat, mean angle; max lon, max lat, max angle; lon <1, lon<3, lon<5; lat<1, lat<3, lat<5; angle<1, angle<3, angle<5:');
slam = [mean(err_slam(:,1)), mean(err_slam(:,2)), mean(err_slam(:,3)), ...
        max(err_slam(:,1)), max(err_slam(:,2)), max(err_slam(:,3)), ...
        sum(err_slam(:,1)<metrics(1))/num_frame*100, sum(err_slam(:,1)<metrics(2))/num_frame*100, sum(err_slam(:,1)<metrics(3))/num_frame*100, ...
        sum(err_slam(:,2)<metrics(1))/num_frame*100, sum(err_slam(:,2)<metrics(2))/num_frame*100, sum(err_slam(:,2)<metrics(3))/num_frame*100, ...
        sum(err_slam(:,3)<angles(1))/num_frame*100, sum(err_slam(:,3)<angles(2))/num_frame*100, sum(err_slam(:,3)<angles(3))/num_frame*100];

CV = [mean(err_cv(:,1)), mean(err_cv(:,2)), mean(err_cv(:,3)), ...
        max(err_cv(:,1)), max(err_cv(:,2)), max(err_cv(:,3)), ...
        sum(err_cv(:,1)<metrics(1))/num_frame*100, sum(err_cv(:,1)<metrics(2))/num_frame*100, sum(err_cv(:,1)<metrics(3))/num_frame*100, ...
        sum(err_cv(:,2)<metrics(1))/num_frame*100, sum(err_cv(:,2)<metrics(2))/num_frame*100, sum(err_cv(:,2)<metrics(3))/num_frame*100, ...
        sum(err_cv(:,3)<angles(1))/num_frame*100,  sum(err_cv(:,3)<angles(2))/num_frame*100,  sum(err_cv(:,3)<angles(3))/num_frame*100];

% fine = [mean(err_fine(:,1)), mean(err_fine(:,2)), mean(err_fine(:,3)), ...
%         max(err_fine(:,1)), max(err_fine(:,2)), max(err_fine(:,3)), ...
%         sum(err_fine(:,1)<metrics(1))/num_frame*100, sum(err_fine(:,1)<metrics(2))/num_frame*100, sum(err_fine(:,1)<metrics(3))/num_frame*100, ...
%         sum(err_fine(:,2)<metrics(1))/num_frame*100, sum(err_fine(:,2)<metrics(2))/num_frame*100, sum(err_fine(:,2)<metrics(3))/num_frame*100, ...
%         sum(err_fine(:,3)<angles(1))/num_frame*100,  sum(err_fine(:,3)<angles(2))/num_frame*100,  sum(err_fine(:,3)<angles(3))/num_frame*100];


d_cv = vecnorm(err_cv(:,1:2),2,2);

% 
e_cv_all = zeros(num_frame,1);
e_cv_vld = zeros(num_valid,1);
check_ = 1:2;
for k=1:num_frame
    e_cv_all(k) = err_cv(k,check_)*err_cv(k,check_)';
    if(valid_id(k))
        id = id+1;
        e_cv_vld(id) = err_cv(k,check_)*err_cv(k,check_)';
    end
end


rmse_all = sqrt(sum(e_cv_all)/num_frame);
rmse_vld = sqrt(sum(e_cv_vld)/num_valid);

% d_cv = zeros(num_frame,1);
% for k=1:num_frame
% 
% d_cv(k) = sqrt(err_cv(k,1:2)*err_cv(k,1:2)');
% 
% end


% diff_slam = t_w_skw - t_w_okw;
% diff_fine = t_w_fkw - t_w_okw;
% diff_cv = t_w_pkw - t_w_okw;
% diff_cv = 

% e2_slam_ = diff_slam * diff_slam';
% e2_fine_ = diff_fine * diff_fine';

e_slam = zeros(num_valid,1);
% e_fine = zeros(num_valid,1);
e_cv = zeros(num_valid,1);

check_ = 1:2;
id=0;
for k=1:num_frame

if(valid_id(k))
    id = id+1;
    e_slam(id) = err_slam(k,check_)*err_slam(k,check_)';
    % e_fine(k) = diff_fine(k,check_)*diff_fine(k,check_)';
    e_cv(id) = err_cv(k,check_)*err_cv(k,check_)';
end
end

rmse_slam = sqrt(sum(e_slam)/num_valid);
% rmse_fine = sqrt(sum(e_fine)/num_valid);
rmse_cv = sqrt(sum(e_cv)/num_valid);

a=1;


% e2_slam = diag(diff_slam * diff_slam');
% e2_fine = diag(diff_fine * diff_fine');



% err_slam = vecnorm(diff_slam(:,check_),2,2);
% err_fine = vecnorm(diff_fine(:,check_),2,2);
% err_cv = vecnorm(diff_cv(:,check_),2,2);
% slam = [mean(err_slam), max(err_slam), rmse_slam]
% fine = [mean(err_fine), max(err_fine), rmse_fine]


% rotation error
% diff_r_slam = zeros(num_frame,3);
% diff_r_fine = zeros(num_frame,3);
% diff_r_cv = zeros(num_frame,3);
% e_r_slam = zeros(num_frame,1);
% e_r_fine = zeros(num_frame,1);
% e_r_cv = zeros(num_frame,1);
% for k=1:num_frame
% R_w_fk = squeeze(vT_0k(k,1:3,1:3));                % refined trajectory
% R_w_sk = squeeze(slam_vTu0uk(k,1:3,1:3));   % SLAM
% R_w_ok = squeeze(oxts_vTu0uk(k,1:3,1:3));  % Ground truth
% R_w_pk = vT_w_pk{k}(1:3,1:3);                % cross-view
% 
% 
% er_slam = so3_log(R_w_sk'*R_w_ok);
% er_fine = so3_log(R_w_fk'*R_w_ok);
% er_cv = so3_log(R_w_pk'*R_w_ok);
% 
% % store data
% diff_r_slam(k,:) = er_slam';
% diff_r_fine(k,:) = er_fine';
% diff_r_cv(k,:) = er_cv';
% 
% e_r_slam(k,:) = er_slam'*er_slam;
% e_r_fine(k,:) = er_fine'*er_fine;
% e_r_cv(k,:) = er_cv'*er_cv;
% 
% end
% 
% rmse_r_slam = sqrt(sum(e_r_slam)/num_frame);
% rmse_r_fine = sqrt(sum(e_r_fine)/num_frame);
% rmse_r_cv = sqrt(sum(e_r_cv)/num_frame);
% 
% err_r_slam = vecnorm(diff_r_slam,2,2);
% err_r_fine = vecnorm(diff_r_fine,2,2);
% err_r_cv = vecnorm(diff_r_cv,2,2);
% r_slam = [rad2deg(mean(err_r_slam)), rad2deg(max(err_r_slam)), rad2deg(rmse_r_slam)]
% r_fine = [rad2deg(mean(err_r_fine)), rad2deg(max(err_r_fine)), rad2deg(rmse_r_fine)]



% slam_traj = [mean(err_slam), max(err_slam), rmse_slam, rad2deg(mean(err_r_slam)), rad2deg(max(err_r_slam)), rad2deg(rmse_r_slam)];
% cv_traj = [mean(err_cv), max(err_cv), rmse_cv, rad2deg(mean(err_r_cv)), rad2deg(max(err_r_cv)), rad2deg(rmse_r_cv)];
% fine_traj = [mean(err_fine), max(err_fine), rmse_fine, rad2deg(mean(err_r_fine)), rad2deg(max(err_r_fine)), rad2deg(rmse_r_fine)];



% tb_name = ["mean 2D dist", "max 2D dist", "rmse dist",...
%            "mean 3D rot", "max 3D rot", "rmse 3D rot",...
%            "mean lon", "mean lat", "mean angle", ...
%             "max lon", "max lat", "max angle",...
%             "lon <1", "lon<3", "lon<5",... 
%             "lat<1", "lat<3", "lat<5",...
%             "angle<1", "angle<3", "angle<5"];
% tb_value = [slam_traj, slam; cv_traj, CV; fine_traj, fine];
% % tb_value = array2table([slam; CV; fine]);
% 
% a2t = array2table(tb_value,"VariableNames",tb_name);







%% calculate error
%{
diff_slam = t_w_skw - t_w_okw;
diff_fine = p_w_skw - t_w_okw;
% diff_cv = 

% e2_slam_ = diff_slam * diff_slam';
% e2_fine_ = diff_fine * diff_fine';

e_slam = zeros(num_frame,1);
e_fine = zeros(num_frame,1);

check_ = 1:2;

for k=1:num_frame
e_slam(k) = diff_slam(k,check_)*diff_slam(k,check_)';
e_fine(k) = diff_fine(k,check_)*diff_fine(k,check_)';
end

rmse_slam = sqrt(sum(e_slam)/num_frame);
rmse_fine = sqrt(sum(e_fine)/num_frame);


% e2_slam = diag(diff_slam * diff_slam');
% e2_fine = diag(diff_fine * diff_fine');



err_slam = vecnorm(diff_slam(:,check_),2,2);
err_fine = vecnorm(diff_fine(:,check_),2,2);
slam = [mean(err_slam), max(err_slam), rmse_slam]
fine = [mean(err_fine), max(err_fine), rmse_fine]

% rotation error
diff_r_slam = zeros(num_frame,3);
diff_r_fine = zeros(num_frame,3);
e_r_slam = zeros(num_frame,1);
e_r_fine = zeros(num_frame,1);
for k=1:num_frame
pred_R_w_ok = vT_0k{k}(1:3,1:3);
R_w_sk = squeeze(vT_us0_usk(k,1:3,1:3));
R_w_ok = squeeze(oxts_vTu0uk(k,1:3,1:3));

er_slam = so3_log(R_w_sk'*R_w_ok);
er_fine = so3_log(pred_R_w_ok'*R_w_ok);

% store data
diff_r_slam(k,:) = er_slam';
diff_r_fine(k,:) = er_fine';

e_r_slam(k,:) = er_slam'*er_slam;
e_r_fine(k,:) = er_fine'*er_fine;

end

rmse_r_slam = sqrt(sum(e_r_slam)/num_frame);
rmse_r_fine = sqrt(sum(e_r_fine)/num_frame);

err_r_slam = vecnorm(diff_r_slam,2,2);
err_r_fine = vecnorm(diff_r_fine,2,2);
r_slam = [rad2deg(mean(err_r_slam)), rad2deg(max(err_r_slam)), rad2deg(rmse_r_slam)]
r_fine = [rad2deg(mean(err_r_fine)), rad2deg(max(err_r_fine)), rad2deg(rmse_r_fine)]


% a=1;
%}




end