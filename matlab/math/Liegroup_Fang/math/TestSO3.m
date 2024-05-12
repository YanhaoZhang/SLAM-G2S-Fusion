clc
clear all

SO3 = SO3()

TRIAL = 1;

fprintf(2, ['Test SO3 correctness with ', ...
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
    %
    
    DiffExpLog = vecR - SO3.Log(SO3.Exp(vecR));
    DiffJlExpJr = SO3.Jl(vecR) - SO3.Exp(vecR) * SO3.Jr(vecR);    
    DiffJl_nInvJr = SO3.Jl(vecR)*SO3.InvJr(-vecR) - eye(3);    
    DiffJr_nInvJl = SO3.Jr(vecR)*SO3.InvJl(-vecR) - eye(3);
    if(max(DiffExpLog) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffExpLog = ', num2str(max(DiffExpLog)) ]);
        disp('Exception occured: aa =! Log(Exp(aa))');
    end
    if(max(max(DiffJlExpJr)) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffExpLog = ', num2str(max(max(DiffJlExpJr))) ]);        
        disp('Exceptoin occured: Jl(aa) != Exp(aa)*Jr(aa)');
    end
    if(max(max(DiffJl_nInvJr)) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffJl_nInvJr = ', num2str(max(max(DiffJl_nInvJr))) ]);        
        disp('Exceptoin occured: Jl(aa) * InvJr(-aa) != eye(3)');
    end
    if(max(max(DiffJr_nInvJl)) > 1e-10)
        vecR
        disp(['angle = ', num2str(norm(vecR)), ...
            '    pi - angle = ', num2str(pi-angle), ...
            '    DiffJr_nInvJl = ', num2str(max(max(DiffJr_nInvJl))) ]);        
        disp('Exceptoin occured: Jr(aa) * InvJl(-aa) != eye(3)');
    end    
end
end
tt1 = cputime;
disp(['SO3 Implementaton is correct in the test! \n']);
disp(['Time consumed for testing:        ', num2str(tt1 - tt0)]);

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
TRIAL = 100;

fprintf(2, ['Test SO3 Timing with ', ...
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
MatrR = SO3.Exp(vecR);
SkewM = SO3.Hat(vecR);
%

t0 = cputime;

for j = 1 : TRIAL
for i = 1 : 10000
    SO3.Exp(vecR);
end
end

t1 = cputime;
tExp = t1 - t0

for j = 1 : TRIAL
for i = 1 : 10000
	SO3.Log(MatrR); 
end
end

t2 = cputime;
tLog = t2 - t1

for j = 1 : TRIAL
for i = 1 : 10000
    SO3.Jr (vecR);
end
end

t3 = cputime;
tJr = t3 - t2

for j = 1 : TRIAL
for i = 1 : 10000
    SO3.Jl (vecR);
end
end

t4 = cputime;
tJl = t4 - t3

for j = 1 : TRIAL
for i = 1 : 10000
    SO3.InvJr (vecR);
end
end

t5 = cputime;
tInvJr = t5 - t4

for j = 1 : TRIAL
for i = 1 : 10000
    SO3.InvJl (vecR);
end
end

t6 = cputime;
tInvJl = t6 - t5

for j = 1 : TRIAL
for i = 1 : 10000
    SO3.Hat (vecR);
end
end

t7 = cputime;
tHat = t7 - t6

for j = 1 : TRIAL
for i = 1 : 10000
    SO3.Vee (SkewM);
end
end

t8 = cputime;
tVee = t8 - t7




