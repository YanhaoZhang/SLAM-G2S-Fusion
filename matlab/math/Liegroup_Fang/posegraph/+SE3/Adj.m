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