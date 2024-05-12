clear all
close all
clc


num_poses = 20;

cir_radius = 10 / (2 * pi);

theta = 2 * pi/ num_poses;

sa = sin(theta);
ca = cos(theta);

dx = cir_radius * sa;
dy = cir_radius * (1 - ca);

dT = [ ca, -sa, dx;
           sa, ca, dy;
           0, 0 , 1];
 DTT = eye(4);
 DTT(1:2, 1:2) = dT(1:2, 1:2);
 DTT(1:2, 4) = dT(1:2, 3);
 
 v_dtt = SE3.Log(DTT);

 v_dt = v_dtt([1,2,6]);
 
Pc = eye(3);
x= [];
y = [];
o = [];
for i = 0 : 1 : num_poses-1
    x = [x; Pc(1,3)];
    y = [y; Pc(2,3)];
    o = [o; atan2( Pc(2,1), Pc(1,1) )];
    Pc = Pc * dT;
end


figure
axis square
hold on
plot(x,y, 'ro-')
for i = 1 : 1 : num_poses
    dox = 0.2 * cos ( o(i) );
    doy = 0.2 * sin ( o(i) );
    quiver(x(i), y(i), dox, doy);
end
hold off


       
       