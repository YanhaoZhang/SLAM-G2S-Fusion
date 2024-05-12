function [R_jacob,R_transpose] = SO3_angle2M_Jacobian (theta, axis)
%angle2RotMatrix: angle to rotation matrix
% https://www.learnopencv.com/rotation-matrix-to-euler-angles/


switch axis
    case 'x'
        R_jacob = [0 0 0; 0 -sin(theta) -cos(theta); 0 cos(theta) -sin(theta)];
        R_transpose = [];
    case 'y'
        R_jacob = [-sin(theta) 0 cos(theta); 0 0 0; -cos(theta) 0 -sin(theta)];
        R_transpose = [];
    case 'z'
        R_jacob = [-sin(theta) -cos(theta) 0; cos(theta) -sin(theta) 0; 0 0 0];
        R_transpose = [];
    otherwise
        error('wrong axis along which object rotates')
end