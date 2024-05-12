% Hat : obtain skew-symmetric matrix from vector 3by1
function matr3by3 = Hat (vec3by1)
    matr3by3 = [           0,        -vec3by1(3),   vec3by1(2);
                           vec3by1(3),         0,           -vec3by1(1);
                          -vec3by1(2),   vec3by1(1),          0         ];
end