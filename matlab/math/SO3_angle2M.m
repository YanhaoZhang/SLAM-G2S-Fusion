function [R,R_transpose,varargout] = SO3_angle2M (theta, axis)
%angle2RotMatrix: angle to rotation matrix
% https://www.learnopencv.com/rotation-matrix-to-euler-angles/


switch axis
    case 'x'
        R = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
        R_transpose = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        R_jacob = [0 0 0; 0 -sin(theta) -cos(theta); 0 cos(theta) -sin(theta)];
    case 'y'
        R = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
        R_transpose = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
        R_jacob = [-sin(theta) 0 cos(theta); 0 0 0; -cos(theta) 0 -sin(theta)];
    case 'z'
        R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
        R_transpose = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
        R_jacob = [-sin(theta) -cos(theta) 0; cos(theta) -sin(theta) 0; 0 0 0];
    otherwise
        error('wrong axis along which object rotates')

end


% lazy to change the previous code that uses this function
if nargout==3
    varargout{1}=R_jacob;
end
end