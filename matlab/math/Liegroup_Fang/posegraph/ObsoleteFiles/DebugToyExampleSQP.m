clear all
close all
clc

axis = [0;0;1];

pose1 = [ SO3.Exp(0*axis), [0;0; 0;];
                zeros(1,3),  1 ];
pose2 = [ SO3.Exp(pi/4*axis), [1;0; 0;];
                zeros(1,3),  1 ];
pose3 = [ SO3.Exp(-pi/4*axis), [2;0; 0;];
                zeros(1,3),  1 ]; 
pose4 = [ SO3.Exp(pi/4*axis), [2;1; 0;];
                zeros(1,3),  1 ];     
pose5 = [ SO3.Exp(3*pi/4*axis), [2;2; 0;];
                zeros(1,3),  1 ];     
pose6 = [ SO3.Exp(-3*pi/4*axis), [1;2; 0;];
                zeros(1,3),  1 ];      
pose7 = [ SO3.Exp(3*pi/4*axis), [0;2; 0;];
                zeros(1,3),  1 ];    
pose8 = [ SO3.Exp(-3*pi/4*axis), [0;1; 0;];
                zeros(1,3),  1 ];

% relative poses
% add noise 
transSigma = 0.1;

rotaSigma = 0.01;
%
% # odometry
rlvpose12 = SE3.Exp(SE3.Log(SE3.Inv(pose1)*pose2) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose23 = SE3.Exp(SE3.Log(SE3.Inv(pose2)*pose3) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose34 = SE3.Exp(SE3.Log(SE3.Inv(pose3)*pose4) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose45 = SE3.Exp(SE3.Log(SE3.Inv(pose4)*pose5) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose56 = SE3.Exp(SE3.Log(SE3.Inv(pose5)*pose6) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose67 = SE3.Exp(SE3.Log(SE3.Inv(pose6)*pose7) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose78 = SE3.Exp(SE3.Log(SE3.Inv(pose7)*pose8) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
% # global loopclosure
rlvpose18 = SE3.Exp(SE3.Log(SE3.Inv(pose1)*pose8) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
% # local loopclosure
rlvpose24 = SE3.Exp(SE3.Log(SE3.Inv(pose2)*pose4) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose46 = SE3.Exp(SE3.Log(SE3.Inv(pose4)*pose6) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose68 = SE3.Exp(SE3.Log(SE3.Inv(pose6)*pose8) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
rlvpose28 = SE3.Exp(SE3.Log(SE3.Inv(pose2)*pose8) + ...
    [transSigma * rand(2,1); zeros(3,1); rotaSigma*randn] );
% compute odometry;
dpose1 = pose1;
dpose2 = dpose1 * rlvpose12;
dpose3 = dpose2 * rlvpose23;
dpose4 = dpose3 * rlvpose34;
dpose5 = dpose4 * rlvpose45;
dpose6 = dpose5 * rlvpose56;
dpose7 = dpose6 * rlvpose67;
dpose8 = dpose7 * rlvpose78;

% Graph G
G = ConstrainedFormulation;
%--------- add nodes ----------
G = G.addNode(1, dpose1);
G = G.addNode(2, dpose2);
G = G.addNode(3, dpose3);
G = G.addNode(4, dpose4);
G = G.addNode(5, dpose5);
G = G.addNode(6, dpose6);
G = G.addNode(7, dpose7);
G = G.addNode(8, dpose8);
%--------- add edges ----------
transInfo = 1/(transSigma*transSigma);
rotaInfo = 1/(rotaSigma*rotaSigma);
covinv = blkdiag( 1/(transSigma*transSigma), ...
                           1/(transSigma*transSigma), ...
                           1, ...
                           1, ...
                           1, ...
                           1/(rotaSigma*rotaSigma) );
% %
G = G.addEdge(1, 1, 2, rlvpose12, covinv);
G = G.addEdge(2, 2, 3, rlvpose23, covinv);
G = G.addEdge(3, 3, 4, rlvpose34, covinv);
G = G.addEdge(4, 4, 5, rlvpose45, covinv);
G = G.addEdge(5, 5, 6, rlvpose56, covinv);
G = G.addEdge(6, 6, 7, rlvpose67, covinv);
G = G.addEdge(7, 7, 8, rlvpose78, covinv);
% add global loop-closure
G = G.addEdge(8, 1, 8, rlvpose18, covinv);
% add local loop-closure
if(1)
G = G.addEdge(9, 2, 4, rlvpose24, covinv);  
G = G.addEdge(10, 4, 6, rlvpose46, covinv);
G = G.addEdge(11, 6, 8, rlvpose68, covinv);
G = G.addEdge(12, 2, 8, rlvpose28, covinv);    
end
% display Graph Information
G.MaxIters = 50

G.plotTrajectory();

G = G.solveBySQP()
G.plotTrajectory(1000);
