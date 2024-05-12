clear all
close all
clc

CHOICES = 4; % 6

for chs = CHOICES
        
    file = ['data_cmd/', 'BatchComandsChs', num2str(chs), '.txt'];
    delete(file)
    diary(file)
    
    G = NLSFormulation;
    G.MaxIters = 20;

    G.IsotropicNoise = 1;
    
    
    switch chs
        case 1
            G.add_data_g2o('g2o/dataset/MITb.g2o')
        case 2
            G.add_data_g2o('g2o/dataset/ring.g2o')
        case 3
            G.add_data_g2o('g2o/dataset/intel.g2o')
        case 4
            G.add_data_g2o('g2o/dataset/M3500.g2o')
        case 5
            G.add_data_g2o('g2o/dataset/ringCity.g2o')
        case 6
            G.add_data_g2o('g2o/dataset/MySphere.g2o')
        case 7
            G.add_data_g2o('g2o/dataset/parking-garage.g2o')
        otherwise
            error('Not a predefined dataset');
    end
    
    %G.DumpData3D('MySphere.txt')
    
    %G.plotTrajectory(1);
    
    G.solveByGaussNewton()
    
    diary off
    %G.plotTrajectory(2);
    
end

