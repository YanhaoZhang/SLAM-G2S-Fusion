% Vee : otain vector 3by1 from skew-symmetric matrix
function vec3by1 = Vee (matr3by3)
    vec3by1 = [ matr3by3(3, 2);
                       matr3by3(1, 3);
                       matr3by3(2, 1); ];
end
