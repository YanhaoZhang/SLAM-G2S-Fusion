%% compute the skew symetric matrix
%{
This is a second version of Teng's code skew.m For applying it to 2-D
rotation. given input 3x1 vector x represents 3-D rotation vector. Given
2*1 vector represents 2-D rotation vector.
%}


function y = skew( x )
%  
% 
    y=[0 -x(3) x(2);...
       x(3) 0 -x(1);...
       -x(2) x(1) 0];
end

