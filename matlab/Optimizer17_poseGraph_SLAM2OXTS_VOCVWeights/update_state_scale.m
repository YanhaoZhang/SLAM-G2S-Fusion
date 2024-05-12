function [vT_0k, vS_0k] = update_state_scale(dX, vT_0k, vS_0k)


dim_r = 3;
dim_t = 3;
dim_s = 1;
dim_vo = dim_r+dim_t;
nPoses = size(vT_0k,1);

for k=2:nPoses

R0k = vT_0k{k}(1:3,1:3);
t_0_k0 = vT_0k{k}(1:3,4);
s_0k = vS_0k(k);

dr = dX((dim_vo+dim_s)*(k-2)+(1:dim_r));
dt = dX((dim_vo+dim_s)*(k-2)+dim_r+(1:dim_t));
ds = dX((dim_vo+dim_s)*(k-2)+dim_r+dim_t+(1:dim_s));

R0k = so3_exp(dr)*R0k;
t_0_k0 = t_0_k0 + dt;
s_0k = s_0k + ds;


vT_0k{k} = [R0k t_0_k0; 0 0 0 1];
vS_0k(k) = s_0k;

end


% update scale
% ds = dX(1:end-6);
% vsa(2:end) = vsa(2:end) + ds;


% update pose. pose0 is fixed
% dr = dX(end-6 + (1:3));
% dt = dX(end-6 + (4:6));
% Rc1c0 = so3_exp(dr)*Rc1c0;
% t_c1_c0c1 = t_c1_c0c1 + dt;

% vsa(2:end) = vsa(2:end) + dX;
% fprintf('-------------------------\n');
end