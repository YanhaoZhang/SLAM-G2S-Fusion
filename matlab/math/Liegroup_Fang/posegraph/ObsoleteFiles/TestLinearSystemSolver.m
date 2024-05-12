clear all
close all
clc

load('LinearSystemEquilibrium.mat')
spy(H)



t0 = cputime;

X = H\h;

t1 = cputime;

tDirectSolve = t1 - t0

%

t00 = cputime;

%L = ichol(H);
[x, flag]= pcg(H, h, 1e-3, size(H,1));

t11 = cputime;

tPcgSolve = t11 - t00

%


switch flag
    case 0
disp('pcg converged to the desired tolerance tol within maxit iterations.')
    case 1
disp('pcg iterated maxit times but did not converge.')
    case 2
disp('Preconditioner M was ill-conditioned.')
    case 3
disp('pcg stagnated. (Two consecutive iterates were the same.)')
    case 4
disp('One of the scalar quantities calculated during pcg became too small or too large to continue computing.')    
end