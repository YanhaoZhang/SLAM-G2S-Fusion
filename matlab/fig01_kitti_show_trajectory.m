addpath(genpath('math_crossview'))
addpath(genpath('utils_crossview'))
addpath(genpath('math'))
addpath(genpath('Optimizer17_poseGraph_SLAM2OXTS_Weights'))
addpath(genpath('plot_utils'))
addpath(genpath('utils_fig'))


% clc
clear
close all

data_root = '../../ICRA_result/test_kitti/';
% dataset = '0930_0034';


Datasets = {'1003_0027','1003_0042','1003_0034','0930_0016','0930_0018','0930_0020','0930_0027','0930_0028','0930_0033','0930_0034'};

Ids = [1, 2, 3, 4,5,6,7,8,9,10];

Ids = [6];

data_num = length(Ids);

Distribut_All = cell(length(Ids),2);
num_max_p = 0;

RMSE_All = zeros(data_num, 2);

for k=1:data_num
    id = Ids(k);
    dataset = Datasets{id};
    load([data_root, 'figPlot_poses_', dataset, '.mat'])
    load([data_root, 'paperFig_', dataset,'_ablation_updateTrajectory_now.mat'])  % _slam_updateTrajectory_covBound, 

    % align with trajectory origin
    k
    show_trajectory_kitti(vT_us0_usk, slam_vTu0uk, oxts_vTu0uk, coarse_vM_usk_uok, fine_valid )
    

    % show_trajectory_kitti_globalAlign(vT_us0_usk, slam_vTu0uk, oxts_vTu0uk, coarse_vM_usk_uok, fine_valid )
    
end