function [F] = calculate_Error_noScale_huber( vT_0k, vo_vM_usl_usk, pred_vR_usk_uok, pred_vM_usk_uok, w )


% nFrames = size(vT_0k,1);
% num_obs = 0;
dim_r = 3;
dim_t = 3;
% dim_s = 1;
nVO = size(vo_vM_usl_usk,1);    % visual odometry
nCV = size(pred_vM_usk_uok,1);  % loop closure by cross-view
nCV_r = size(pred_vR_usk_uok,1);  % loop closure by cross-view
% nPoses = size(vT_0k,1);
% num_obs = nVO + nCV;
dim_vo = dim_r+dim_t;
% F = zeros(dim_vo*num_obs+ dim_s*(nPoses-1), 1);    % initialize F
% F = zeros(dim_vo*num_obs, 1);    % initialize F
F = zeros( dim_vo*nVO + dim_r*nCV_r + dim_t*nCV, 1);    % initialize F

%% loop over visual odometry
id_o = 1;
for a=1:nVO
    % VO measurement
    Tlk = vo_vM_usl_usk{a,1};    % visual odometry
    Rlk = Tlk(1:3,1:3);
    t_l_kl = Tlk(1:3,4);
    l = vo_vM_usl_usk{a,2};            % id of previous pose
    k = vo_vM_usl_usk{a,3};            % id of current pose
    w_vo = vo_vM_usl_usk{a,5};            % visual odometry weights

    % corresponding poses
    R0l = vT_0k{l}(1:3,1:3);
    t_0_l0 = vT_0k{l}(1:3,4);
    R0k = vT_0k{k}(1:3,1:3);
    t_0_k0 = vT_0k{k}(1:3,4);
    % s_0k = 1;

    % error
    eR = R0l'*R0k*Rlk';
    er = so3_log(eR);
    et = R0l'*(t_0_k0-t_0_l0) - t_l_kl;

    
    % store data
    er = er*w(1)*w_vo;
    et = et*w(2)*w_vo;
    F(dim_vo*(id_o-1)+(1:dim_vo),1) = [er;et];
    id_o = id_o+1;
end


%% loop over cross-view loop: rotation
id_o = 1;
for a=1:nCV_r
    % cross-views measurement
    R_sk_k = pred_vR_usk_uok{a,1};    % cross-view measurement: slam to oxts
    R_0_sk = pred_vR_usk_uok{a,2};    % original SLAM pose
    k = pred_vR_usk_uok{a,3};         % id of current pose
    
    % corresponding poses
    R0k = vT_0k{k}(1:3,1:3);

    % error
    eR = R0k'*R_0_sk*R_sk_k;
    er = so3_log(eR);

    wr = [w(3), 0, 0; 
          0, w(3), 0; 
          0, 0, w(3)];

    % store data
    er = wr*er;


    % et_ = w_huber_*et;


    F(dim_vo*nVO+dim_r*(id_o-1)+(1:dim_r),1) = [er];
    id_o = id_o+1;
end


%% loop over cross-view loop: translation
id_o = 1;
for a=1:nCV
    % cross-views measurement
    T_sk_k = pred_vM_usk_uok{a,1};    % cross-view measurement: slam to oxts
    % R_sk_k = T_sk_k(1:3,1:3);
    t_sk_ksk = T_sk_k(1:3,4);

    T_0_sk = pred_vM_usk_uok{a,2};    % original SLAM pose
    R_0_sk = T_0_sk(1:3,1:3);
    t_0_sk0 = T_0_sk(1:3,4);

    k = pred_vM_usk_uok{a,3};            % id of current pose
    

    % corresponding poses
    % R0k = vT_0k{k}(1:3,1:3);
    t_0_k0 = vT_0k{k}(1:3,4);

    % error
    % eR = R0k'*R_0_sk*R_sk_k;
    % er = so3_log(eR);
    et = R_0_sk'*(t_0_k0-t_0_sk0) - t_sk_ksk;

    % wr = [w(3), 0, 0; 
    %       0, w(3), 0; 
    %       0, 0, w(3)];
    wt = [w(4), 0, 0; 
          0, w(5), 0; 
          0, 0, w(6)];

    % cross-view weights
    weight_cv = pred_vM_usk_uok{a,4};    % lon, lat of cross-view
    w_cv = [weight_cv(1), 0, 0; 
            0, weight_cv(2), 0; 
            0, 0, 1];

    wt = w_cv*wt;

    
    % wt = [w(4), 0, 0; 
    %       0, w(4), 0; 
    %       0, 0, w(4)];

    % huber weights
    % et_norm = sqrt(et'*et);
    % w_huber = 1;
    % if(et_norm>w(8))  
    %     w_huber = w(8)*sign(et_norm)/et_norm;    %todo need to check if this is correct
    % end
    % 
    % w_huber_ = sqrt(w_huber);
    % wt = w_huber_*wt;


    w_huber = [1,1,1];
    if(abs(et(1))>w(8))  
        w_huber(1) = w(8)*sign(et(1))/et(1);    %todo need to check if this is correct
    end
    if(abs(et(2))>w(8))  
        w_huber(2) = w(8)*sign(et(2))/et(2);    %todo need to check if this is correct
    end
    if(abs(et(3))>w(8))  
        w_huber(3) = w(8)*sign(et(3))/et(3);    %todo need to check if this is correct
    end
    w_huber_ = diag(sqrt(w_huber));
    wt = w_huber_*wt;


    

    % store data
    % er = wr*er;
    et = wt*et;
    % F(dim_vo*nVO+dim_vo*(id_o-1)+(1:dim_vo),1) = [et];
    F(dim_vo*nVO+dim_r*nCV_r+dim_t*(id_o-1)+(1:dim_t),1) = [et];
    id_o = id_o+1;
end


% %% scale constraint
% id_o = 1;
% for k=2:nPoses
% 
% l=k-1;
% s_0l = vS_0k(l);
% s_0k = vS_0k(k);
% 
% es = s_0k - s_0l;
% es = w(5) * es;
% 
% 
% F(dim_vo*num_obs+dim_s*(id_o-1)+(1:dim_s),1) = [es];
% id_o = id_o+1;
% 
% 
% end

end