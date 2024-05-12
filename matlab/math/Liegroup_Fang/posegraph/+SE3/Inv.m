% Inv : obtain inversion of transformation matrix        
function invT = Inv (Tmatr)
    invR= Tmatr([1,2,3], [1,2,3])';
    t = Tmatr([1,2,3], 4);
    invT = [ invR,  -invR*t;
                 0, 0, 0,  1 ];
end
