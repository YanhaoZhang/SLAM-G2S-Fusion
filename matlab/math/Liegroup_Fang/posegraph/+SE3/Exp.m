% Exponential mapping : Encapsulated matrix exponential
% obtain transformation matrix from vector 6by1        
function Tmatr = Exp (vec6by1)
    u = vec6by1([1, 2, 3]);
    Theta =  vec6by1([4, 5, 6]);
    rHat = [  0,            -Theta(3),   Theta(2);
                  Theta(3),       0,        -Theta(1);
                 -Theta(2),  Theta(1),        0         ];
    angle = norm(Theta);
    if(angle < 1e-6)
        A = 1;
        B = 1/2;
        C = 1/6;
    else
        A = sin(angle)/angle;
        B = (1 - cos(angle))/(angle*angle);
        C = (1 - A)/(angle*angle); 
    end
    RotMatr = eye(3) + A * rHat + B * rHat * rHat;
    J_l = eye(3) + B * rHat + C * rHat * rHat;
    Tmatr = [ RotMatr,   J_l * u;
                    0, 0, 0, 1 ];
end