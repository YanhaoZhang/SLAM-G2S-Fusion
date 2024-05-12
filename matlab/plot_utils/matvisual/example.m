clear, clc, close all

%% generate a test matrix
% A = eye(10);
% A = magic(30);
% A = ones(10);
% A = zeros(10);
% A = NaN(10);
% A = Inf(10);
% A = rand(20, 10, 3);
A = round(10*rand(50, 50));
% A = round(10*rand(10, 5, 5));
% A = abs(fft(1:10));
% A = complex(rand(5,5), rand(5,5));
% A = char(1300*ones(10));
% A = zeros(0, 20);

%% visualize the matrix
matvisual(A, 'annotation')