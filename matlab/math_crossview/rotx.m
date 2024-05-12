function [R] = rotx(t)

c = cos(t);
s = sin(t);

R = [1 0  0;
     0 c -s;
     0 s  c];


end