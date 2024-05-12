% Exponential mapping : Encapsulated matrix exponential
% obtain transformation matrix from vector 6by1
function Tmatr = se3_exp (vec6by1)
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

