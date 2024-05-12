% Logarithm mapping : 
% obtain vector 3by1 from rotation matrix     
function Theta = Log ( RotMatr )
    angle = acos( (trace(RotMatr) -1 )/2 ); % return [pi, 0]
    if(~isreal(angle))
%         fprintf(2, ['complex number: angle @ SO3.Log()', ...
%             '      (trace(RotMatr)-1)/2 = ', ...
%             num2str((trace(RotMatr) -1 )/2), ...
%             '      angle = ', num2str(angle), '\n' ]);
        angle = real(angle);
    end
    if(angle < 1e-8)
        axis = [1; 0; 0];
    elseif(pi -angle < 1e-6) % RotMatr = RotMatr' iff. angle = pi
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
end

