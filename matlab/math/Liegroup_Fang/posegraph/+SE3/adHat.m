function matr = adHat (vec6by1)
    t = vec6by1([1, 2, 3]);
    r =  vec6by1([4, 5, 6]);
    rHat = [  0,            -r(3),   r(2);
                  r(3),       0,        -r(1);
                 -r(2),  r(1),        0         ];  
    tHat = [  0,            -t(3),   t(2);
                  t(3),       0,        -t(1);
                 -t(2),  t(1),        0         ];
    matr = [ rHat,              tHat;
                  zeros(3,3),     rHat ];
end