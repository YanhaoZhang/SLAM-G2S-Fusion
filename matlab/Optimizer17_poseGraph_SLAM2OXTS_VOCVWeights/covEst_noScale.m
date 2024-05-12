function [Cov_pose] = covEst_noScale(J)


dim_r = 3;
dim_t = 3;


dim = dim_r+dim_t;

H = J'*J;
Cov = H \ eye(size(H,1));

nPoses = size(H,1)/dim +1;


Cov_pose = cell(nPoses,1);
Cov_pose{1} = zeros(6,6);

for k=2:nPoses
    cov_k = Cov(dim*(k-2)+(1:dim), dim*(k-2)+(1:dim));
    Cov_pose{k} = cov_k(1:(dim_r+dim_t), 1:(dim_r+dim_t));
end

end


