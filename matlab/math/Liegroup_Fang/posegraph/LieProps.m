clear all
close all
clc

a = rand(6,1);


% SE3.adHat(a) * a = 0;
SE3.adHat(a) * a

SE3.adHat(a)' * a    % SE3.adHat(a) * a != 0;

% SE3.Jr(a) * a = a
% SE3.InvJr(a) * a = a
% SE3.Jl(a) * a = a
% SE3.InvJl(a) * a = a


SE3.Jr(a) * a - a
SE3.InvJr(a) * a - a
SE3.Jl(a) * a - a
SE3.InvJl(a) * a - a