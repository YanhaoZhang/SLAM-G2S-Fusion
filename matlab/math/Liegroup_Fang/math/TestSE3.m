clear all
clc

SE3 = SE3()

TRIAL = 1;

fprintf(2, ['Test SE3 correctness with ', ...
    num2str(TRIAL), ',0000 randomly generated data \n']);

tt0 = cputime;
for j = 1 : TRIAL
for i = 1 : 10000
    % generate valid axis-angle
    vec = 100*rand(3,1);
    axis = vec/norm(vec);
    angle = abs(randn);
    while(angle) > 2*pi
        angle = angle - 2 * pi;
    end
    while (angle > pi)
        angle = 2*pi - angle;
        axis = -axis;
    end
    vecR = angle*axis;
    vect = 100*rand(3,1);
    vecRt = [vect; vecR];
    %
    
    DiffExpLog = vecRt - SE3.Log(SE3.Exp(vecRt));
    DiffAdjInvAdjExp = SE3.Adj(SE3.Exp(vecRt)) * SE3.InvAdj(SE3.Exp(vecRt)) - eye(6);
    DiffJlExpJr = SE3.Jl(vecRt) - SE3.Adj(SE3.Exp(vecRt)) * SE3.Jr(vecRt);    
    DiffJl_nInvJr = SE3.Jl(vecRt)*SE3.InvJr(-vecRt) - eye(6);    
    DiffJr_nInvJl = SE3.Jr(vecRt)*SE3.InvJl(-vecRt) - eye(6);
    DiffQRJ = SE3.Ql(vecRt) - SO3.Exp(vecR) * SE3.Qr(vecRt) ...
        - SO3.Hat(SO3.Jl(vecR)*vect) * SO3.Exp(vecR) * SO3.Jr(vecR);
    
    if(max(DiffExpLog) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffExpLog = ', num2str(max(DiffExpLog)) ]);
        disp('Exception occured: au =! Log(Exp(au))');
    end
    if(max(max(DiffAdjInvAdjExp)) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffAdjInvAdjExp = ', num2str(max(DiffAdjInvAdjExp)) ]);
        disp('Exception occured: Adj(Exp(au)) * InvAdj(Exp(au)) != eye(6)');        
    end
    if(max(max(DiffJlExpJr)) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffExpLog = ', num2str(max(max(DiffJlExpJr))) ]);        
        disp('Exceptoin occured: Jl(au) != Adj(Exp(au))*Jr(au)');
    end
    if(max(max(DiffJl_nInvJr)) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffJl_nInvJr = ', num2str(max(max(DiffJl_nInvJr))) ]);        
        disp('Exceptoin occured: Jl(au) * InvJr(-au) != eye(6)');
    end
    if(max(max(DiffJr_nInvJl)) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffJr_nInvJl = ', num2str(max(max(DiffJr_nInvJl))) ]);        
        disp('Exceptoin occured: Jr(au) * InvJl(-au) != eye(6)');
    end    
    if(max(max(DiffQRJ)) > 1e-10 )
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffQRJ = ', num2str(max(max(DiffQRJ))) ]);        
        disp('Exceptoin occured: Ql - R*Qr + skew(Jl*u) * R * Jr != eye(3)');
    end           
end
end
tt1 = cputime;
disp(['SE3 Implementaton is correct in the test! \n']);
disp(['Time consumed for testing:        ', num2str(tt1 - tt0)]);


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
TRIAL = 100;

fprintf(2, ['Test SE3 Timing with ', ...
    num2str(TRIAL), ',0000 trials for each function \n']);

% generate data
vec = 100*rand(3,1);
axis = vec/norm(vec);
angle = abs(randn);
while(angle) > 2*pi
    angle = angle - 2 * pi;
end
while (angle > pi)
    angle = 2*pi - angle;
    axis = -axis;
end
vecR = angle*axis;
vect = 100*rand(3,1);
vecRt = [vect; vecR];

Tmatr = SE3.Exp(vecRt);
SkewT = SE3.Hat(vecRt);
%

t0 = cputime;

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.Exp(vecRt);
end
end

t1 = cputime;
tExp = t1 - t0

for j = 1 : TRIAL
for i = 1 : 10000
	SE3.Log(Tmatr); 
end
end

t2 = cputime;
tLog = t2 - t1

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.Jr (vecRt);
end
end

t3 = cputime;
tJr = t3 - t2

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.Jl (vecRt);
end
end

t4 = cputime;
tJl = t4 - t3

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.InvJr (vecRt);
end
end

t5 = cputime;
tInvJr = t5 - t4

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.InvJl (vecRt);
end
end

t6 = cputime;
tInvJl = t6 - t5

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.Hat (vecRt);
end
end

t7 = cputime;
tHat = t7 - t6

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.Vee (SkewT);
end
end

t8 = cputime;
tVee = t8 - t7

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.Adj (Tmatr);
end
end

t9 = cputime;
tAdj = t9 - t8

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.InvAdj (Tmatr);
end
end

t10 = cputime;
tInvAdj = t10 - t9

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.Qr (vecRt);
end
end

t11 = cputime;
tQr = t11 - t10

for j = 1 : TRIAL
for i = 1 : 10000
    SE3.Ql (vecRt);
end
end

t12 = cputime;
tQl = t12 - t11


