% Jr : right-hand Jacobian
% obtain right-hand Jacobian from vector 6by1        
function J_r = Jr (vec6by1)
    u = vec6by1([1, 2, 3]);
    Theta =  vec6by1([4, 5, 6]);            
    % righ-hand Jacobian for rotation part
    rHat = [  0,            -Theta(3),   Theta(2);
                  Theta(3),       0,        -Theta(1);
                 -Theta(2),  Theta(1),        0         ];
    angle = norm(Theta);       
    % coefficients
    if(angle < 1e-6)
        ca = 1/6;
        cb = -1/24;
        cc = -1/120;
    else
        ca = (angle - sin(angle))/(angle*angle*angle);      
        cb = (1-angle*angle/2 - cos(angle))/(angle*angle*angle*angle);
        cc = (cb - 3* ( ca/(angle*angle) - 1/(6*angle*angle) ))/2;      
    end             
    % compute SO3 right-hand-Jacobian 
    if (angle < 1e-6)
        J_r_rot = eye(3) - rHat/2 + rHat*rHat/6;
    else                 
        J_r_rot = eye(3) - ...
            rHat*(1 - cos(angle))/(angle*angle) + ...
            ca * (rHat*rHat);
    end       
    % compute Q matrix on upper-right corner
    uHat = [ 0,  -u(3),  u(2);
                 u(3),  0,  -u(1);
                -u(2),  u(1),  0 ];
    % compute Q_rh
    Q_rh = -uHat/2 + ca * (rHat*uHat + uHat*rHat - rHat*uHat*rHat) ...
        + cb * (rHat*rHat*uHat + uHat*rHat*rHat - 3*rHat*uHat*rHat) ...
        - cc * (rHat*uHat*rHat*rHat + rHat*rHat*uHat*rHat);
    % construct righ-hand Jacobian for SE3
    J_r = [  J_r_rot,       Q_rh;
               zeros(3,3),  J_r_rot ];           
end