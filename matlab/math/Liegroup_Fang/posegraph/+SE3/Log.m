% Logarithm mapping : 
% obtain vector 6by1 from transformation matrix          
function vec6by1 = Log (Tmatr)
    RotMatr= Tmatr([1,2,3], [1,2,3]);
    t = Tmatr([1,2,3], 4);
    % minimal representation for rotation part
    % compute Logarithm mapping for Rotation
    angle = acos( (trace(RotMatr) -1 )/2 );
    if(~isreal(angle))
%         fprintf(2, ['complex number: angle @ SE3.Log()', ...
%             '      (trace(RotMatr)-1)/2 = ', ...
%             num2str((trace(RotMatr) -1 )/2), ...
%             '      angle = ', num2str(angle), '\n' ]);
        angle = real(angle);
    end
    if(angle < 1e-8)
        axis = [1; 0; 0];
    elseif(pi - angle < 1e-8)  % RotMatr = RotMatr' iff. angle = pi
        axis = - [ sign(RotMatr(2,3))*sqrt((1+RotMatr(1,1))/2); 
                       sign(RotMatr(1,3))*sqrt((1+RotMatr(2,2))/2); 
                       sign(RotMatr(1,2))*sqrt((1+RotMatr(3,3))/2) ];
    else
        axis = [ RotMatr(3,2)-RotMatr(2,3);
                     RotMatr(1,3)-RotMatr(3,1);
                     RotMatr(2,1)-RotMatr(1,2); ];
        axis = axis/(2*sin(angle));
    end       
    Theta = angle*axis;        
    % compute inv(rotation-left-hand-jacobian)
    rHat = [  0,            -Theta(3),   Theta(2);
                  Theta(3),       0,        -Theta(1);
                 -Theta(2),  Theta(1),        0         ];                 
    if (angle < 1e-6)
        J_l_inv = eye(3) - rHat/2 + rHat*rHat/12;
    else 
        J_l_inv = eye(3)  ...
            - rHat/2 + ...
            ((1 - angle*sin(angle)/(2*(1-cos(angle))))/(angle*angle)) * (rHat*rHat);    
    end
    % minimal representation for translation part
    u = J_l_inv * t;
    % rocover minimal representation for T
    vec6by1 = [   u; 
                       Theta ];            
end