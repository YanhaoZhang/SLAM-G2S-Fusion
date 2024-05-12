% InvJl : inverse left-hand Jacobian
% obtain inversion of left-hand Jacobian from vector 3by1        
function inv_J_l = InvJl ( vec3by1 ) 
    theta = norm(vec3by1);  
    vecHat = [           0,        -vec3by1(3),   vec3by1(2);
                        vec3by1(3),         0,           -vec3by1(1);
                       -vec3by1(2),   vec3by1(1),          0         ];             
    if (theta < 1e-6)
        inv_J_l = eye(3) - vecHat/2 + vecHat*vecHat/12;
    else
        inv_J_l = eye(3)  ...
            - vecHat/2 + ...
            ((1 - theta*sin(theta)/(2*(1-cos(theta))))/(theta*theta)) * (vecHat*vecHat);
    end           
end   