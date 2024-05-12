function [R] = roty(t)

c = cos(t);
s = sin(t);

R = [c 0 s;
     0 1 0;
    -s 0 c];


end