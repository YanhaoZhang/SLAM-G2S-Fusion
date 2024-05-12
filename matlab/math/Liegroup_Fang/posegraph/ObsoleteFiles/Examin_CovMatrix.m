clear all
close all
clc


load('covMatrix.mat');

[u, d, v] = svd(A);

%
subplot(2, 3, 1)
spy(cov)
title('cov')
%
subplot(2, 3, 2)
spy(Q)
title('Q')
%
subplot(2, 3, 3)
spy(A)
title('A')
%
subplot(2, 3, 4)
spy(u)
title('u')
%
subplot(2, 3, 5)
spy(d)
title('d')
%
subplot(2, 3, 6)
spy(v)
title('v')
%

max_A_times_cov = max(max( A * cov))


