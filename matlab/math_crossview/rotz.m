function [R] = rotz(t)

c = cos(t);
s = sin(t);

R = [c -s 0;
     s  c 0;
     0  0 1];


end