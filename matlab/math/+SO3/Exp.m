% Exponential mapping : Encapsulated matrix exponential
% obtain rotation matrix from vector 3by1
function RotMatr = Exp ( vec3by1 )
    theta = norm(vec3by1);
    vecHat = [           0,        -vec3by1(3),   vec3by1(2);
                        vec3by1(3),         0,           -vec3by1(1);
                       -vec3by1(2),   vec3by1(1),          0         ];   
    if (theta < 1e-6)
        RotMatr = eye(3) + vecHat + vecHat*vecHat/2;
    else
        RotMatr = eye(3) + ...
            vecHat*sin(theta)/theta + ...
            (vecHat*vecHat)*(1 - cos(theta))/(theta*theta);
    end
end
