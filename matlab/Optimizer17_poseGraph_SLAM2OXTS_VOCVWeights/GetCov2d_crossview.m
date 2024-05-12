

function [cv_]=GetCov2d_crossview(P,x,y, n, k)
s=1/n; 

% k=40;   % original: k=20, changed by zhan

if trace(P)<1e-5
    r=zeros(2,2);
else
    %P
    r=real(sqrtm(P)); 
end


j = 1:k+1;
q = 2*pi*(j-1)/k;   % q is the step theta, roll it around

rr = s*3*r;

% r_max = max(rr(1,1), rr(2,2));

cv_ = rr*[cos(q); sin(q)]+[x;y];   %here: it is equal to confert a circle radius=1 to a ellipse. r represents the rotation matrix (sim(2)). since the diag is sigma1 and sigma2, it represent 1 sigma.


% cv = zeros(2,k+1);
% 
% for j=1:k+1
%   q=2*pi*(j-1)/k; 
%   cv(:,j)=s*3*r*[cos(q); sin(q)]+[x;y];  %here: it is equal to confert a circle radius=1 to a ellipse. r represents the rotation matrix (sim(2)). since the diag is sigma1 and sigma2, it represent 1 sigma.
% end

end

% P is 2x2 covariance matrix, s is scale factor, x,y are estimated robot or landmark location

% A typical call is 

% CV=GetCov(sys.P(1:2,1:2,n),sys.X(1,n),sys.X(2,n));

% sys.P(1:2,1:2,n) is robot covariance, sys.X(1,n) sys.X(2,n) is robot location

% To plot it

% set(fig.car.CV,'xdata',CV(1,:),'ydata',CV(2,:));

% fig.car.CV is a graph handle
