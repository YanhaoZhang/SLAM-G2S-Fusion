classdef SO3

    % Lie group SO(3)
    % Implement all methods as static to avoid passing class objects
    % --------------------------- methods ---------------------------
    % function matr3by3 = Hat (vec3by1)
    % function vec3by1 = Vee (matr3by3)
    % function RotMatr = Exp ( vec3by1 )
    % function Theta = Log ( RotMatr )
    % function J_r = Jr ( vec3by1 )
    % function J_l = Jl ( vec3by1 ) 
    % function inv_J_r = InvJr ( vec3by1 ) 
    % function inv_J_l = InvJl ( vec3by1 ) 
    % function adj = Adj ( RotMatr )
    
    methods (Static)
        
        % Inv : obtain inversion of rotation matrix
        function invRotMatr = Inv (RotMatr)
            invRotMatr = RotMatr';
        end
        
        % Hat : obtain skew-symmetric matrix from vector 3by1
        function matr3by3 = Hat (vec3by1)
            matr3by3 = [           0,        -vec3by1(3),   vec3by1(2);
                                   vec3by1(3),         0,           -vec3by1(1);
                                  -vec3by1(2),   vec3by1(1),          0         ];
        end

        % Vee : otain vector 3by1 from skew-symmetric matrix
        function vec3by1 = Vee (matr3by3)
            vec3by1 = [ matr3by3(3, 2);
                               matr3by3(1, 3);
                               matr3by3(2, 1); ];
        end
       
        % Exponential mapping : Encapsulated matrix exponential
        % obtain rotation matrix from vector 3by1
        function RotMatr = Exp ( vec3by1 )
            theta = norm(vec3by1);
            vecHat = [           0,        -vec3by1(3),   vec3by1(2);
                                vec3by1(3),         0,           -vec3by1(1);
                               -vec3by1(2),   vec3by1(1),          0         ];   
            if (theta < 1e-10)
                RotMatr = eye(3) + vecHat;
            else
                RotMatr = eye(3) + ...
                    vecHat*sin(theta)/theta + ...
                    (vecHat*vecHat)*(1 - cos(theta))/(theta*theta);
            end
        end
   
        % Logarithm mapping : 
        % obtain vector 3by1 from rotation matrix     
        function Theta = Log ( RotMatr )
            angle = acos( (trace(RotMatr) -1 )/2 ); % return [pi, 0]
            if(angle < 1e-10)
                axis = [1; 0; 0];
            elseif(pi -angle < 1e-10) % RotMatr = RotMatr' iff. angle = pi
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
        
        % Jr : right-hand Jacobian
        % obtain right-hand Jacobian from vector 3by1
        function J_r = Jr ( vec3by1 ) 
            theta = norm(vec3by1);  
            vecHat = [           0,        -vec3by1(3),   vec3by1(2);
                                vec3by1(3),         0,           -vec3by1(1);
                               -vec3by1(2),   vec3by1(1),          0         ];             
            if (theta < 1e-10)
                J_r = eye(3) - vecHat/2;
            else                 
                J_r = eye(3) - ...
                    vecHat*(1 - cos(theta))/(theta*theta) + ...
                    (vecHat*vecHat)*(theta - sin(theta))/(theta*theta*theta);
            end           
        end
        
        % Jl : left-hand Jacobian
        % obtain left-hand Jacobian from vector 3by1
        function J_l = Jl ( vec3by1 ) 
            theta = norm(vec3by1);  
            vecHat = [           0,        -vec3by1(3),   vec3by1(2);
                                vec3by1(3),         0,           -vec3by1(1);
                               -vec3by1(2),   vec3by1(1),          0         ];              
            if (theta < 1e-10)
                J_l = eye(3) + vecHat/2;
            else                 
                J_l = eye(3) + ...
                    vecHat*(1 - cos(theta))/(theta*theta) + ...
                    (vecHat*vecHat)*(theta - sin(theta))/(theta*theta*theta);
            end           
        end        
        
        % InvJr : inverse right-hand Jacobian
        % obtain inversion of right-hand Jacobian from vector 3by1        
        function inv_J_r = InvJr ( vec3by1 ) 
            theta = norm(vec3by1);  
            vecHat = [           0,        -vec3by1(3),   vec3by1(2);
                                vec3by1(3),         0,           -vec3by1(1);
                               -vec3by1(2),   vec3by1(1),          0         ];             
            if (theta < 1e-10)
                inv_J_r = eye(3) + vecHat/2;
            else
                inv_J_r = eye(3)  ...
                    + vecHat/2 + ...
                    ((1 - theta*sin(theta)/(2*(1-cos(theta))))/(theta*theta)) * (vecHat*vecHat);
            end           
        end        
        
        % InvJl : inverse left-hand Jacobian
        % obtain inversion of left-hand Jacobian from vector 3by1        
        function inv_J_l = InvJl ( vec3by1 ) 
            theta = norm(vec3by1);  
            vecHat = [           0,        -vec3by1(3),   vec3by1(2);
                                vec3by1(3),         0,           -vec3by1(1);
                               -vec3by1(2),   vec3by1(1),          0         ];             
            if (theta < 1e-10)
                inv_J_l = eye(3) + vecHat/2;
            else
                inv_J_l = eye(3)  ...
                    - vecHat/2 + ...
                    ((1 - theta*sin(theta)/(2*(1-cos(theta))))/(theta*theta)) * (vecHat*vecHat);
            end           
        end   
        
        % Adj: Adjoint of SO(3)
        function adj = Adj ( RotMatr )
            adj = RotMatr;
        end
        
        
    end
    
end
