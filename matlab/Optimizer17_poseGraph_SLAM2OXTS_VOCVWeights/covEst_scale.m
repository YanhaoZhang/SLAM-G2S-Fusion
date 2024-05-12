function [Cov_pose, normal] = covEst_scale(J)


dim_r = 3;
dim_t = 3;
dim_s = 1;
dim = dim_r+dim_t+dim_s;

H = J'*J;
Cov = H \ eye(size(H,1));

nPoses = size(H,1)/dim +1;


Cov_pose = cell(nPoses,1);
Cov_pose{1} = zeros(2,2);

for k=2:nPoses
    cov_k = Cov(dim*(k-2)+(1:dim), dim*(k-2)+(1:dim));
    Cov_pose{k} = cov_k(dim_r+(1:2), dim_r+(1:2));
end


% frame_id = [idk, 2];    % first stores the covariance of next frame, the second is for normalzie the covariance. idx starts from 0
% Cov_pose = cell(2,1);
% for kk=1:2
%     k = frame_id(kk);
% 
% 
%     cov_k = Cov(dim*(k-2)+(1:dim), dim*(k-2)+(1:dim));
%     Cov_pose{kk} = cov_k(dim_r+(1:2), dim_r+(1:2));
% end


%% rescale the pose

P = Cov_pose{2};
if trace(P)<1e-5
    r1=zeros(2,2);
else
    %P
    r1=real(sqrtm(P)); 
end

S = svd(r1);
normal = mean(S) / 0.2;    % assume the second frame is very accurate, the error is within +-0.01*3m




end