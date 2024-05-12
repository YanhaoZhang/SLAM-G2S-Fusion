function [J] = calculate_Jacobian_noScale_huber(vT_0k, vo_vM_usl_usk, pred_vR_usk_uok, pred_vM_usk_uok, w)



dim_r = 3;
dim_t = 3;
dim_s = 0;
nVO = size(vo_vM_usl_usk,1);    % visual odometry
nCV = size(pred_vM_usk_uok,1);  % loop closure by cross-view
nCV_r = size(pred_vR_usk_uok,1);  % loop closure by cross-view
nPoses = size(vT_0k,1);
% num_obs = nVO + nCV;
dim_vo = dim_r+dim_t;
% J = zeros(dim_vo*num_obs+ dim_s*(nPoses-1),   (dim_r+dim_t+dim_s)*nPoses);    % [r t scale] * nPoses
% J = zeros(dim_vo*num_obs,   (dim_r+dim_t)*nPoses);    % [r t] * nPoses
J = zeros(dim_vo*nVO + dim_r*nCV_r + dim_t*nCV,   (dim_r+dim_t+dim_s)*nPoses);    % [r t scale] * nPoses


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
    % s_0k = vS_0k(k);
    s_0k = 1;

    % error
    % eR = R0l'*R0k*Rlk';
    % er = so3_log(eR);
    % et = s_0k*R0l'*(t_0_k0-t_0_l0) - t_l_kl;

    % store data
    % er = er*w(1);
    % et = et*w(2);

    % Jacobian: 
    Jer_Jrl = -R0l';
    Jer_Jrk = R0l';

    Jet_Jrl = s_0k*R0l'*skew(t_0_k0-t_0_l0);
    Jet_Jtl = s_0k*(-R0l');
    Jet_Jtk = s_0k*R0l';
    % Jet_Jsk = R0l'*(t_0_k0-t_0_l0);


    % store
    Jer_Jrl = Jer_Jrl*w(1)*w_vo;
    Jer_Jrk = Jer_Jrk*w(1)*w_vo;
    Jet_Jrl = Jet_Jrl*w(2)*w_vo;
    Jet_Jtl = Jet_Jtl*w(2)*w_vo;
    Jet_Jtk = Jet_Jtk*w(2)*w_vo;
    % Jet_Jsk = Jet_Jsk*w(2);



    J(dim_vo*(id_o-1)+(1:dim_r),         (dim_vo+dim_s)*(l-1)+(1:dim_r)) =  Jer_Jrl;
    J(dim_vo*(id_o-1)+(1:dim_r),         (dim_vo+dim_s)*(k-1)+(1:dim_r)) =  Jer_Jrk;
    J(dim_vo*(id_o-1)+dim_r+(1:dim_t),   (dim_vo+dim_s)*(l-1)+(1:dim_r)) =  Jet_Jrl;
    J(dim_vo*(id_o-1)+dim_r+(1:dim_t),   (dim_vo+dim_s)*(l-1)+dim_r+(1:dim_t)) =  Jet_Jtl;
    J(dim_vo*(id_o-1)+dim_r+(1:dim_t),   (dim_vo+dim_s)*(k-1)+dim_r+(1:dim_t)) =  Jet_Jtk;
    % J(dim_vo*(id_o-1)+dim_r+(1:dim_t),   (dim_vo+dim_s)*(k-1)+dim_r+dim_t+(1:dim_s)) = Jet_Jsk;



    % F(dim*(id_o-1)+(1:dim),1) = [er;et];
    id_o = id_o+1;
end



%% loop over cross-view loop: rotation
id_o = 1;
for a=1:nCV_r
    % cross-views measurement
    % R_sk_k = pred_vR_usk_uok{a,1};    % cross-view measurement: slam to oxts
    % R_0_sk = pred_vR_usk_uok{a,2};    % original SLAM pose
    k = pred_vR_usk_uok{a,3};         % id of current pose

    % corresponding poses
    R0k = vT_0k{k}(1:3,1:3);

    % error
    % eR = R0k'*R_0_sk*R_sk_k;
    % er = so3_log(eR);

    % Jacobian: 
    Jer_Jrk = -R0k';

    wr = [w(3), 0, 0; 
          0, w(3), 0; 
          0, 0, w(3)];


% fprintf('[')
% fprintf('%f, ', Jer_Jrk(:))
% fprintf(']\n')

    % store data
    Jer_Jrk = wr*Jer_Jrk;

% fprintf('%d, ', dim_vo*nVO+dim_r*(id_o-1)+(1:dim_r))
% fprintf(']\n')
% fprintf('%d, ', (dim_vo+dim_s)*(k-1)+(1:dim_r))
% fprintf
% fprintf(']\n')

    J(dim_vo*nVO+dim_r*(id_o-1)+(1:dim_r),       (dim_vo+dim_s)*(k-1)+(1:dim_r)) =  Jer_Jrk;
    id_o = id_o+1;
end





%% loop over cross-view loop
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

    % store data
    % er = er*w(3);
    % et = et*w(4);

    % Jacobian: 
    % Jer_Jrk = -R0k';
    % Jet_Jrk = R0k'*skew(t_0_sk0-t_0_k0);
    Jet_Jtk = R_0_sk';

    % store
    % Jer_Jrk = Jer_Jrk*w(3);
    % Jet_Jrk = Jet_Jrk*w(4);
    % Jet_Jtk = Jet_Jtk*w(4);

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


    % huber weights
    % et_norm = sqrt(et'*et);
    % w_huber = 1;
    % if(et_norm>w(8))  
    %     w_huber = w(8)*sign(et_norm)/et_norm;    %todo need to check if this is correct
    % end
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




    % wt = [0.001*w(4), 0, 0; 
    %       0, w(4), 0; 
    %       0, 0, w(4)];



    % Jer_Jrk = wr*Jer_Jrk;
    % Jet_Jrk = wt*Jet_Jrk;
    Jet_Jtk = wt*Jet_Jtk;

    % J(dim_vo*nVO+dim_vo*(id_o-1)+(1:dim_r),       (dim_vo)*(k-1)+(1:dim_r)) =  Jer_Jrk;
    % J(dim_vo*nVO+dim_vo*(id_o-1)+dim_r+(1:dim_t), (dim_vo+dim_s)*(k-1)+(1:dim_r)) = Jet_Jrk;
    % J(dim_vo*nVO+dim_vo*(id_o-1)+dim_r+(1:dim_t), (dim_vo)*(k-1)+dim_r+(1:dim_t)) =  Jet_Jtk;
    J(dim_vo*nVO+dim_r*nCV_r+dim_t*(id_o-1)+(1:dim_t), (dim_vo+dim_s)*(k-1)+dim_r+(1:dim_t)) =  Jet_Jtk;



    % F(dim*(id_o-1)+(1:dim),1) = [er;et];
    id_o = id_o+1;
end



% %% scale constraint
% id_o = 1;
% for k=2:nPoses
% 
% % s_0l = vS_0k(k-1);
% % s_0k = vS_0k(k);
% % es = s_0k - s_0l;
% % es = w(5) * es;
% 
% 
% % F(dim_vo*num_obs+dim_s*(id_o-1)+(1:dim_s),1) = [es];
% 
% l = k-1;
% 
% 
% Jes_Jsl = -1;
% Jes_Jsk = 1;
% 
% Jes_Jsl = Jes_Jsl*w(5);
% Jes_Jsk = Jes_Jsk*w(5);
% 
% 
% J(dim_vo*num_obs+dim_s*(id_o-1)+(1:dim_s), (dim_vo+dim_s)*(l-1)+dim_r+dim_t+(1:dim_s)) = Jes_Jsl;
% J(dim_vo*num_obs+dim_s*(id_o-1)+(1:dim_s), (dim_vo+dim_s)*(k-1)+dim_r+dim_t+(1:dim_s)) = Jes_Jsk;
% 
% id_o = id_o+1;
% 
% 
% end



% fix pose 0
J(:,1:6) = [];
J = sparse(J);
% id_f = id_f-1;
% J(:,[1:6, dim_p*num_frames + dim_f*(id_f-1)+(1:dim_f)]) = [];
% J(:,[1,end-6+(1:6)]) = [];
% J = sparse(J);

end




% function [J] = calculate_Jacobian(p_3d, p_2d, R, t, K, w)
% 
% num_obs = size(p_2d,1);
% num_f = size(p_3d,1);
% dim_f = 3;
% dim_o = 2;
% J = zeros(dim_o*num_obs,  6+dim_f*num_f);   % pose 6dof + 3*num_f
% 
% 
% fx = K(1,1);
% fy = K(2,2);
% T = [R,t;0 0 0 1];
% 
% for i=1:num_obs
%     u_i = p_2d(i,:)';
%     p_i = p_3d(i,:)';
%     p_tran = ((R')*p_i+(-R'*t));
% 
%     x = p_tran(1);
%     y = p_tran(2);
%     z = p_tran(3);
% 
%     Jp_Je = [fx/x, 0 -fx*x/(z^2); 0 fy/z -fy*y/(z^2)];
%     Je_JT = -inv(T)*[-skew(p_i), eye(3); zeros(1,6)];
%     Je_JT = Je_JT(1:3,:);
%     Je_Jf = R;
% 
% 
%     J(2*(i-1)+(1:2),1:6) =  (1/w)*(- Jp_Je*Je_JT);
%     J(2*(i-1)+(1:2), 6 + dim_f*(i-1)+(1:3)) =  (1/w)*(-Jp_Je*(Je_Jf));
% end
% 
% 
% J = J(:,1:6);
% 
% 
% end