% Jr : right-hand Jacobian
% obtain right-hand Jacobian from vector 3by1
function J_r = Jr ( vec3by1 ) 
    theta = norm(vec3by1);  
    vecHat = [           0,        -vec3by1(3),   vec3by1(2);
                        vec3by1(3),         0,           -vec3by1(1);
                       -vec3by1(2),   vec3by1(1),          0         ];             
    if (theta < 1e-6)
        J_r = eye(3) - vecHat/2 + vecHat*vecHat/6;
    else                 
        J_r = eye(3) - ...
            vecHat*(1 - cos(theta))/(theta*theta) + ...
            (vecHat*vecHat)*(theta - sin(theta))/(theta*theta*theta);
    end           
end