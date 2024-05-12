clear all
close all
clc

G = IncrementalNLSFramework;

dataset_chs = 1;

switch dataset_chs
    case 0
    G.add_data_g2o('input_INTEL_g2o_incremental.g2o')
    case 1
    G.add_data_g2o('g2o/dataset/MITb_incremental.g2o')
    case 2
    G.add_data_g2o('g2o/dataset/ring_incremental.g2o')
    case 3
    G.add_data_g2o('g2o/dataset/intel_incremental.g2o')
    case 4
    G.add_data_g2o('g2o/dataset/M3500_incremental.g2o')
    case 5
    G.add_data_g2o('g2o/dataset/ringCity_incremental.g2o')
    case 6
    G.add_data_g2o('g2o/dataset/MySphere_incremental.g2o')
    case 7 
    G.add_data_g2o('g2o/dataset/parking-garage_incremental.g2o')   
    otherwise
        error('Not a predefined dataset');
end


algorithm_chs = 1;

switch algorithm_chs
    case 1
        algorithm_str = 'Batch';
        file = ['data_cmd/', 'IncrementalComandsDatasetchs', num2str(dataset_chs), algorithm_str, '.txt'];
        delete(file)
        diary(file)
        G.solveByBatchIncrementalGaussNewton ()
        diary off
    case 2
        algorithm_str = 'Online';
        file = ['data_cmd/', 'IncrementalComandsDatasetchs', num2str(dataset_chs), algorithm_str, '.txt'];
        delete(file)
        diary(file)
        G.solveByOnlineIncrementalGaussNewton ()
        diary off
end


