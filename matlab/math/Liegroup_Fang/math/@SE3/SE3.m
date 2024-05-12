classdef SE3
    
    % Lie group SE(3)
    % Implement all methods as static to avoid passing class objects
    % --------------------------- methods ---------------------------
    % function invT = Inv (Tmatr)
    % function matr4by4 = Hat(vec6by1)
    % function vec6by1 = Vee(matr4by4)   
    % function Tmatr = Exp (vec6by1)
    % function vec6by1 = Log (Tmatr)
    % function J_r = Jr (vec6by1)
    % function J_l = Jl (vec6by1)
    % function inv_J_r = InvJr (vec6by1)
    % function inv_J_l = InvJl (vec6by1)
    % function adj = Adj( Tmatr )
    % function inv_adj = InvAdj (Tmatr)
    % function Q_rh = Qr (vec6by1)
    % function Q_lh = Ql (vec6by1)
    % minimal representation vector 6by1 should be
    % [ translation part 3by1; rotation part 3by1 ];
    
    methods(Static)
        
        % Inv : obtain inversion of transformation matrix        
        function invT = Inv (Tmatr)
            invR= Tmatr([1,2,3], [1,2,3])';
            t = Tmatr([1,2,3], 4);
            invT = [ invR,  -invR*t;
                         0, 0, 0,  1 ];
        end
         
        % Hat : obtain augmented skew-symmetric matrix from vector 6by1        
        function matr4by4 = Hat (vec6by1)
            % vec6by1 = [trans31, rota31]    
            matr4by4 = [  0,  -vec6by1(6),  vec6by1(5),  vec6by1(1);
                                  vec6by1(6),  0,  -vec6by1(4),  vec6by1(2);
                                  -vec6by1(5),  vec6by1(4),  0,  vec6by1(3);
                                  0,   0,   0,   0  ];                      
        end

        % Vee : otain vector 6by1 from augmented skew-symmetric matrix        
        function vec6by1 = Vee (matr4by4)
            vec6by1 = [ matr4by4(1, 4);
                               matr4by4(2, 4);
                               matr4by4(3, 4);
                               matr4by4(3, 2);
                               matr4by4(1, 3);
                               matr4by4(2, 1); ];             
        end
               
        % Exponential mapping : Encapsulated matrix exponential
        % obtain transformation matrix from vector 6by1        
        function Tmatr = Exp (vec6by1)
            u = vec6by1([1, 2, 3]);
            Theta =  vec6by1([4, 5, 6]);
            rHat = [  0,            -Theta(3),   Theta(2);
                          Theta(3),       0,        -Theta(1);
                         -Theta(2),  Theta(1),        0         ];
            angle = norm(Theta);
            if(angle < 1e-10)
                RotMatr = eye(3) + rHat;
                J_l = eye(3) + rHat/2;
            else
                A = sin(angle)/angle;
                B = (1 - cos(angle))/(angle*angle);
                C = (1 - A)/(angle*angle);                
                RotMatr = eye(3) + A * rHat + B * rHat * rHat;
                J_l = eye(3) + B * rHat + C * rHat * rHat;
            end
            Tmatr = [ RotMatr,   J_l * u;
                            0, 0, 0, 1 ];
        end
        
        % Logarithm mapping : 
        % obtain vector 6by1 from transformation matrix          
        function vec6by1 = Log (Tmatr)
            RotMatr= Tmatr([1,2,3], [1,2,3]);
            t = Tmatr([1,2,3], 4);
            % minimal representation for rotation part
            % compute Logarithm mapping for Rotation
            angle = acos( (trace(RotMatr) -1 )/2 );
            if(angle < 1e-10)
                axis = [1; 0; 0];
            elseif(pi - angle < 1e-10)  % RotMatr = RotMatr' iff. angle = pi
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
            if (angle < 1e-10)
                J_l_inv = eye(3) - rHat/2;
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
            ca = (angle - sin(angle))/(angle*angle*angle);      
            cb = (1-angle*angle/2 - cos(angle))/(angle*angle*angle*angle);
            cc = (cb - 3* ( ca/(angle*angle) - 1/(6*angle*angle) ))/2;                
            % compute SO3 right-hand-Jacobian 
            if (angle < 1e-10)
                J_r_rot = eye(3) - rHat/2;
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
        
        % Jl : left-hand Jacobian
        % obtain left-hand Jacobian from vector 6by1        
        function J_l = Jl (vec6by1)
            u = vec6by1([1, 2, 3]);
            Theta =  vec6by1([4, 5, 6]);            
            % righ-hand Jacobian for rotation part
            rHat = [  0,            -Theta(3),   Theta(2);
                          Theta(3),       0,        -Theta(1);
                         -Theta(2),  Theta(1),        0         ];
            angle = norm(Theta);       
            % coefficients
            ca = (angle - sin(angle))/(angle*angle*angle);      
            cb = (1-angle*angle/2 - cos(angle))/(angle*angle*angle*angle);
            cc = (cb - 3* ( ca/(angle*angle) - 1/(6*angle*angle) ))/2;                
            % compute SO3 right-hand-Jacobian 
            if (angle < 1e-10)
                J_l_rot = eye(3) - rHat/2;
            else                 
                J_l_rot = eye(3) + ...
                    rHat*(1 - cos(angle))/(angle*angle) + ...
                    ca * (rHat*rHat);
            end       
            % compute Q matrix on upper-right corner
            uHat = [ 0,  -u(3),  u(2);
                         u(3),  0,  -u(1);
                        -u(2),  u(1),  0; ];    
            % compute Q_lh
            Q_lh = uHat/2 + ca * (rHat*uHat + uHat*rHat + rHat*uHat*rHat) ...
                - cb * (rHat*rHat*uHat + uHat*rHat*rHat - 3*rHat*uHat*rHat) ...
                - cc * (rHat*uHat*rHat*rHat + rHat*rHat*uHat*rHat);
            % construct righ-hand Jacobian for SE3
            J_l = [  J_l_rot,        Q_lh;
                       zeros(3,3),  J_l_rot ];        
        end
        
        % InvJr : inverse right-hand Jacobian
        % obtain inversion of right-hand Jacobian from vector 6by1 
        function inv_J_r = InvJr (vec6by1)
            u = vec6by1([1, 2, 3]);
            Theta =  vec6by1([4, 5, 6]);            
            % righ-hand Jacobian for rotation part
            rHat = [  0,            -Theta(3),   Theta(2);
                          Theta(3),       0,        -Theta(1);
                         -Theta(2),  Theta(1),        0         ];
            angle = norm(Theta);       
            % coefficients
            ca = (angle - sin(angle))/(angle*angle*angle);      
            cb = (1-angle*angle/2 - cos(angle))/(angle*angle*angle*angle);
            cc = (cb - 3* ( ca/(angle*angle) - 1/(6*angle*angle) ))/2;                
            % compute SO3 right-hand-Jacobian 
            if (angle < 1e-10)
                inv_J_r_rot = eye(3) + rHat/2;
            else
                inv_J_r_rot = eye(3)  ...
                    + rHat/2 + ...
                    ((1 - angle*sin(angle)/(2*(1-cos(angle))))/(angle*angle)) * (rHat*rHat);
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
            inv_J_r = [  inv_J_r_rot,  -inv_J_r_rot*Q_rh*inv_J_r_rot;
                             zeros(3,3),    inv_J_r_rot ];           
        end       
        
        % InvJl : inverse left-hand Jacobian
        % obtain inversion of left-hand Jacobian from vector 6by1             
        function inv_J_l = InvJl (vec6by1)
            u = vec6by1([1, 2, 3]);
            Theta =  vec6by1([4, 5, 6]);            
            % righ-hand Jacobian for rotation part
            rHat = [  0,            -Theta(3),   Theta(2);
                          Theta(3),       0,        -Theta(1);
                         -Theta(2),  Theta(1),        0         ];
            angle = norm(Theta);       
            % coefficients
            ca = (angle - sin(angle))/(angle*angle*angle);      
            cb = (1-angle*angle/2 - cos(angle))/(angle*angle*angle*angle);
            cc = (cb - 3* ( ca/(angle*angle) - 1/(6*angle*angle) ))/2;                
            % compute SO3 right-hand-Jacobian 
            if (angle < 1e-10)
                inv_J_l_rot = eye(3) + rHat/2;
            else       
                inv_J_l_rot = eye(3)  ...
                    - rHat/2 + ...
                    ((1 - angle*sin(angle)/(2*(1-cos(angle))))/(angle*angle)) * (rHat*rHat);
            end       
            % compute Q matrix on upper-right corner
            uHat = [ 0,  -u(3),  u(2);
                         u(3),  0,  -u(1);
                        -u(2),  u(1),  0; ];    
            % compute Q_lh
            Q_lh = uHat/2 + ca * (rHat*uHat + uHat*rHat + rHat*uHat*rHat) ...
                - cb * (rHat*rHat*uHat + uHat*rHat*rHat - 3*rHat*uHat*rHat) ...
                - cc * (rHat*uHat*rHat*rHat + rHat*rHat*uHat*rHat);
            % construct righ-hand Jacobian for SE3
            inv_J_l = [  inv_J_l_rot,  -inv_J_l_rot*Q_lh*inv_J_l_rot;
                             zeros(3,3),   inv_J_l_rot ];        
        end       
        
        % Adj: Adjoint of SE(3)        
        function adj = Adj (Tmatr)
            R= Tmatr([1,2,3], [1,2,3]);
            t = Tmatr([1,2,3], 4);
            tHat = [ 0,  -t(3),  t(2);
                         t(3),  0,  -t(1);
                        -t(2),  t(1),  0; ];
            adj = [ R,  tHat*R;
                       zeros(3,3),  R; ];
        end
                
        % InvAdj: inverse Adjoint of SE(3)        
        function inv_adj = InvAdj (Tmatr)
            invR= Tmatr([1,2,3], [1,2,3])';
            t = Tmatr([1,2,3], 4);
            tHat = [ 0,  -t(3),  t(2);
                         t(3),  0,  -t(1);
                        -t(2),  t(1),  0; ];
            inv_adj = [ invR,  -invR*tHat;
                             zeros(3,3),  invR; ];
        end        
        
        % Qr : Q matrix to constitute SE3 right-hand Jacobian 
        function Q_rh = Qr (vec6by1)
            u = vec6by1([1, 2, 3]);
            Theta =  vec6by1([4, 5, 6]);            
            % righ-hand Jacobian for rotation part
            rHat = [  0,            -Theta(3),   Theta(2);
                          Theta(3),       0,        -Theta(1);
                         -Theta(2),  Theta(1),        0         ];
            angle = norm(Theta);       
            % coefficients
            ca = (angle - sin(angle))/(angle*angle*angle);      
            cb = (1-angle*angle/2 - cos(angle))/(angle*angle*angle*angle);
            cc = (cb - 3* ( ca/(angle*angle) - 1/(6*angle*angle) ))/2;              
            % compute Q matrix on upper-right corner
            uHat = [ 0,  -u(3),  u(2);
                         u(3),  0,  -u(1);
                        -u(2),  u(1),  0 ];
            % compute Q_rh
            Q_rh = -uHat/2 + ca * (rHat*uHat + uHat*rHat - rHat*uHat*rHat) ...
                + cb * (rHat*rHat*uHat + uHat*rHat*rHat - 3*rHat*uHat*rHat) ...
                - cc * (rHat*uHat*rHat*rHat + rHat*rHat*uHat*rHat);        
        end
        
        % Ql : Q matrix to constitute SE3 left-hand Jacobian         
        function Q_lh = Ql (vec6by1)
            u = vec6by1([1, 2, 3]);
            Theta =  vec6by1([4, 5, 6]);            
            % righ-hand Jacobian for rotation part
            rHat = [  0,            -Theta(3),   Theta(2);
                          Theta(3),       0,        -Theta(1);
                         -Theta(2),  Theta(1),        0         ];
            angle = norm(Theta);       
            % coefficients
            ca = (angle - sin(angle))/(angle*angle*angle);      
            cb = (1-angle*angle/2 - cos(angle))/(angle*angle*angle*angle);
            cc = (cb - 3* ( ca/(angle*angle) - 1/(6*angle*angle) ))/2;   
            % compute Q matrix on upper-right corner
            uHat = [ 0,  -u(3),  u(2);
                         u(3),  0,  -u(1);
                        -u(2),  u(1),  0; ];    
            % compute Q_lh
            Q_lh = uHat/2 + ca * (rHat*uHat + uHat*rHat + rHat*uHat*rHat) ...
                - cb * (rHat*rHat*uHat + uHat*rHat*rHat - 3*rHat*uHat*rHat) ...
                - cc * (rHat*uHat*rHat*rHat + rHat*rHat*uHat*rHat);            
        end
        
        
        
    end
    
end

