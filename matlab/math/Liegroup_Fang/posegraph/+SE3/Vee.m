% Vee : otain vector 6by1 from augmented skew-symmetric matrix        
function vec6by1 = Vee (matr4by4)
    vec6by1 = [ matr4by4(1, 4);
                       matr4by4(2, 4);
                       matr4by4(3, 4);
                       matr4by4(3, 2);
                       matr4by4(1, 3);
                       matr4by4(2, 1); ];             
end