function [R,t,Noise]=Pose_3D_to_3D(p_1,p_2)
%% This code is used to get the relatice rotation and translation of the
% robot poses

% Closed-Form Solution of Absolute Orientation using Orthonormal Matrices


% Copyright (C) 2018 by Yongbo Chen
[N,D]=size(p_1);% N means the number of matching features and D means the dimension of the space
p_1_c=mean(p_1);
p_2_c=mean(p_2);
H=zeros(D,D);
for i=1:1:N
    q_1(i,:)=p_1(i,:)-p_1_c;
    q_2(i,:)=p_2(i,:)-p_2_c;
    H=H+q_1(i,:)'*q_2(i,:);
end
[U,S,V]=svd(H);
X=V*U';
if norm(det(X)-1)<1e-5
    R=X;
elseif norm(det(X)+1)<1e-5
    V(3,:)=-V(3,:);
    R=V*U';
end
if size(R,1)~=3
    t=[];
else
    t=p_2_c'-R*p_1_c';
end
Noise=[];%

end





% rigidTransform3D(points1, points2) Fit a rigid transformation between two sets of 3-D points
%  [R, t] = rigidTransform3D(points1, poitns2) computes a rigid
%  trasnsormation between two sets of corresponding 3-D points. points1 and
%  points2 are M-by-3 arrays of [x,y,z] coordinates. R is a 3-by-3 rotation
%  matrix and t is a 1-by-3 translation vector. 
%
%  If the number of points is greater than 3, the rigidTransform3D computes
%  a least-squares fit.
%
%  Note: the function returns t as a column vector, and R in the
%  post-multiply convention (matrix times column vector). R and t must be
%  transposed to be used in most CVST functions.

% Copyright 2016 MathWorks, Inc.

%#codegen

function [R, t] = rigidTransform3D(p, q)

n = cast(size(p, 1), 'like', p);
m = cast(size(q, 1), 'like', q);

% Find data centroid and deviations from centroid
pmean = sum(p,1)/n;
p2 = bsxfun(@minus, p, pmean);

qmean = sum(q,1)/m;
q2 = bsxfun(@minus, q, qmean);

% Covariance matrix
C = p2'*q2;

[U,~,V] = svd(C);

% Handle the reflection case
R = V*diag([1 1 sign(det(U*V'))])*U';

% Compute the translation
t = qmean' - R*pmean';
end