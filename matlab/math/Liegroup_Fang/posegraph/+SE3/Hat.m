% Hat : obtain augmented skew-symmetric matrix from vector 6by1        
function matr4by4 = Hat (vec6by1)
    % vec6by1 = [trans31, rota31]    
    matr4by4 = [  0,  -vec6by1(6),  vec6by1(5),  vec6by1(1);
                          vec6by1(6),  0,  -vec6by1(4),  vec6by1(2);
                          -vec6by1(5),  vec6by1(4),  0,  vec6by1(3);
                          0,   0,   0,   0  ];                      
end
