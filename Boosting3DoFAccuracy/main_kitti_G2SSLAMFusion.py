
import os

# os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
# os.environ['CUDA_VISIBLE_DEVICES'] = '0'

import torch
import torch.optim as optim
from dataLoader.KITTI_dataset3_outlier import load_test1_data #, dataset #, load_test2_data
from models_kitti import Model
import scipy.io as scio

import ssl

ssl._create_default_https_context = ssl._create_unverified_context  # for downloading pretrained VGG weights

import numpy as np
import os
import argparse

import time


from dataLoader.utils_pykitti import euler_from_rotation_matrix, euler_to_rotation_matrix


import matlab.engine




def test1(net_test, args, dataset_i, save_path, epoch):

    dataset_save = dataset_i[16:18]+dataset_i[19:21]+'_'+dataset_i[28:32]


    # matlab function
    eng = matlab.engine.start_matlab()
    matlab_root = '../matlab'
    eng.addpath(matlab_root + '/main_kitti_crossview')
    eng.addpath(matlab_root + '/Optimizer17_poseGraph_SLAM2OXTS_VOCVWeights')    # use crossview and visual odometry weights
    eng.addpath(matlab_root + '/utils_crossview')
    eng.addpath(matlab_root + '/math')
    eng.addpath(matlab_root + '/math_crossview')


    net_test.eval()

    dataloader = load_test1_data(mini_batch, dataset_i, args.shift_range_lat, args.shift_range_lon, args.rotation_range)
    num_frames = len(dataloader)

    # load laplacian matrix for visual odometry weight
    L = dataloader.load_L_saved()

    # initialize trajectory and scale
    vT_us0_usk = dataloader.get_SLAM_pose()
    vS_0k = np.ones(num_frames)


    # save SLAM and OXTS pose
    vT_us0_usk_org = dataloader.get_SLAM_pose_org()
    vT_uo0_uok = dataloader.get_oxts_pose()
    save_path2 = 'save_path/test_kitti' #
    scio.savemat(os.path.join(save_path2, 'poses_'+dataset_save+'.mat'),
                 {'slam_vTu0uk': vT_us0_usk, 'oxts_vTu0uk': vT_uo0_uok,
                  'keyframe_name': dataloader.keyframe_name
                  })






    
    pred_lons = np.zeros(num_frames)
    pred_lats = np.zeros(num_frames)
    pred_oriens = np.zeros(num_frames)
    
    gt_lons = np.zeros(num_frames)
    gt_lats = np.zeros(num_frames)
    gt_oriens = np.zeros(num_frames)

    fine_valid = np.zeros([num_frames,3], dtype=bool)

    pred_vM_usk_uok = np.zeros([num_frames,4,4])   # store cross-view prediction after updated using the updated SLAM trajectory
    coarse_vM_usk_uok = np.zeros([num_frames,4,4])   # store coarse cross-view prediction: simply zero or use cross-view prediction

    #todo the folowing weights are not used for G2S-SLAM fusion
    weight_vM_usk_uok = np.zeros([num_frames,2])
    corr_vP_usk_uok = np.zeros([num_frames, 155, 155])
    corr_vMask_usk_uok = np.zeros([num_frames, 155, 155], dtype=bool)
    corr_vW_usk_uok = np.zeros([num_frames,1])
    cv_vW_usk_uok = np.zeros([num_frames,1])     # store the valid cross-view weight


    # skip the initial pose
    zero = 0

    b_updated = False
    i_updated = 1

    i_valid = 0
    in_i_valid = 0

    corr_H = 155
    corr_W = 155

    b_l_init = True



    start_time = time.time()
    with torch.no_grad():
        for i in range(num_frames):

            i_updated = i_updated+1

            # skip the first frame
            if i==0:
                pred_lons[i] = zero * args.shift_range_lon
                pred_lats[i] = zero * args.shift_range_lat
                pred_oriens[i] = zero * args.rotation_range

                gt_lons[i]= zero
                gt_lats[i]= zero
                gt_oriens[i]= zero

                fine_valid[i]= [False, False, False]

                pred_vM_usk_uok[i] = np.identity(4)
                coarse_vM_usk_uok[i] = pred_vM_usk_uok[i]

                # save weights
                weight_vM_usk_uok[i, 0] = 1
                weight_vM_usk_uok[i, 1] = 1

                # the following is actually used
                corr_vW_usk_uok[i] = 1
                cv_vW_usk_uok[i] = 1
                corr_vP_usk_uok[i] = np.ones([corr_H, corr_W])

                mask = np.zeros([corr_H, corr_W], dtype=bool)
                th = 2
                u_range = [int(corr_H / 2) - th, int(corr_H / 2) + th + 1]
                v_range = [int(corr_W / 2) - th, int(corr_W / 2) + th + 1]
                mask[v_range[0]:v_range[1], u_range[0]:u_range[1]] = True
                corr_vMask_usk_uok[i] = mask

                # scale pose graph using cross-view and visual odometry weights and output covariance
                _, _, _, Cov_circle = eng.main_crossview_sPGO_SLAM2oxts_CVandVOweights_Cov(
                    vT_us0_usk_org, L, vT_us0_usk, vS_0k,
                    fine_valid, coarse_vM_usk_uok, corr_vP_usk_uok, cv_vW_usk_uok, corr_vMask_usk_uok,
                    i, nargout=4)

                Cov_circle = np.array(Cov_circle)

                continue

            data = dataloader.getitem(i, coarse_vM_usk_uok)
            sat_map, left_camera_k, grd_left_imgs, gt_shift_u, gt_shift_v, gt_heading = [item.unsqueeze(0).to(device) for item in data[:6]]

            if args.proj == 'CrossAttn':
                pred_u, pred_v, pred_orien, corr_u, corr_v, corr_p, corr_w, corr_m \
                = net_test.CVattn_rot_corr(sat_map, grd_left_imgs, left_camera_k, gt_heading=gt_heading, mode='test', cov_i=Cov_circle[i])   # corr_p_u_axis, corr_p_v_axis \
            else:
                raise EOFError

            pred_lons[i] = pred_u.data.cpu().numpy()
            pred_lats[i] = pred_v.data.cpu().numpy()
            pred_oriens[i] = pred_orien.data.cpu().numpy()


            # get gt data
            T_us0_usk_org = dataloader.get_SLAM_pose(i)
            T_uo0_uok = dataloader.get_oxts_pose(i)
            T_usk_uok = np.linalg.inv(T_us0_usk_org) @ T_uo0_uok

            _, _, gt_oriens_ = euler_from_rotation_matrix(T_usk_uok[:3, :3])
            gt_lons[i] = T_usk_uok[0, 3]
            gt_lats[i] = T_usk_uok[1, 3]
            gt_oriens[i] = gt_oriens_ / np.pi * 180



            # save weights   #todo this weights are actually not used for G2S-SLAM fusion
            weight_vM_usk_uok[i ,0] = corr_u
            weight_vM_usk_uok[i, 1] = corr_v
            corr_vP_usk_uok[i] = corr_p[0,:,:].cpu().numpy()
            corr_vW_usk_uok[i] = corr_w[0].cpu().numpy()
            corr_vMask_usk_uok[i] = corr_m[0].cpu().numpy()



            # save raw prediction
            M_usk_uok = np.identity(4)
            raw_roll_usk_uok = 0
            raw_pitch_usk_uok = 0
            raw_yaw_usk_uok= float(pred_oriens[i]) * np.pi / 180
            R_usk_uok = euler_to_rotation_matrix(raw_roll_usk_uok, raw_pitch_usk_uok, raw_yaw_usk_uok)
            M_usk_uok[0, 3] = float(pred_lons[i])
            M_usk_uok[1, 3] = float(pred_lats[i])
            M_usk_uok[:3, :3] = R_usk_uok
            pred_vM_usk_uok[i] = M_usk_uok
            # fine_vM_uok_usk[i] = raw_vM_uok_usk[i]   # initialize prediction


            l = i - 1
            T_us0_usk = dataloader.get_SLAM_pose(i)
            T_us0_usl = dataloader.get_SLAM_pose(l)
            T_usl_usk = np.linalg.inv(T_us0_usl) @ T_us0_usk  # SLAM relative trajectory: use original SLAM visual odometry

            # load the cross-view prediction (after latest update)
            M_usl_uol = coarse_vM_usk_uok[l]
            M_uol_uok = np.linalg.inv(M_usl_uol) @ T_usl_usk @ M_usk_uok



            # calculate error: they are not in the same frame!
            A_uol_uok = M_uol_uok[:3, :3]
            t_ol_okol = M_uol_uok[:3, 3]
            R_usl_usk = T_usl_usk[:3, :3]
            t_sl_sksl = T_usl_usk[:3, 3]

            # R_sl_ol = M_usl_uol[:3, :3]
            # b_sl_okol = R_sl_ol.T @ t_ol_okol
            b_sl_okol = t_ol_okol



            err_roll, err_pitch, err_yaw = euler_from_rotation_matrix(A_uol_uok @ R_usl_usk.T)
            err_t = np.abs(t_sl_sksl - b_sl_okol)



            valid = np.array([False, False, False])
            if(float(corr_vW_usk_uok[l]) == 1.0 and float(corr_vW_usk_uok[i]) == 1.0):
                if (err_t[0] < 0.5 ):
                    valid[0] = True
                if (err_t[1] < 0.5 ):
                    valid[1] = True
                if np.abs(err_yaw) < (0.25 * np.pi / 180):
                    valid[2] = True


            coarse_lons = pred_lons[i]  # gt_lons[i]
            coarse_lats = pred_lats[i]  # gt_lons[i]
            coarse_oriens = pred_oriens[i]
            fine_valid[i] = valid



            coarse_M_usk_uok = np.identity(4)
            coarse_roll_usk_uok = 0
            coarse_pitch_usk_uok = 0
            coarse_yaw_usk_uok = float(coarse_oriens) * np.pi / 180
            coarse_R_usk_uok = euler_to_rotation_matrix(coarse_roll_usk_uok, coarse_pitch_usk_uok, coarse_yaw_usk_uok)
            coarse_M_usk_uok[0, 3] = float(coarse_lons)
            coarse_M_usk_uok[1, 3] = float(coarse_lats)
            coarse_M_usk_uok[:3, :3] = coarse_R_usk_uok
            coarse_vM_usk_uok[i] = coarse_M_usk_uok

            valid_k = fine_valid[i, 0] and fine_valid[i, 1] and fine_valid[i, 2]



            if valid_k:
                cv_vW_usk_uok[i] = 0.5*(1+cv_vW_usk_uok[l])    # this stores the weight of this cross-view measurement


            b_traj_update = False
            if valid[0] and valid[1] and valid[2] :
                i_valid = i_valid+1

                if i_valid % 1 == 0 or i==num_frames-1:
                    b_traj_update = True
                    i_valid = 0

            # update SLAM trajectory
            if valid_k and b_traj_update:

                # debug
                # diff = [pred_lons[i]-gt_lons[i], pred_lats[i]-gt_lats[i], pred_oriens[i]-gt_oriens[i],
                #         pred_lons[i], pred_lats[i], pred_oriens[i],
                #         corr_vW_usk_uok[i], valid_k, cv_vW_usk_uok[i],
                #         # weight_vM_usk_uok[i, 0], weight_vM_usk_uok[i ,1]
                #         # abs(pred_lons[i]-pred_lons[i-1]), abs(pred_lats[i]-pred_lats[i-1]), abs(pred_oriens[i]-pred_oriens[i-1]),
                #         # err_t2
                #         ]
                # print(diff)

                # scale pose graph using cross-view and visual odometry weights and output covariance
                vT_us0_usk, vS_0k, coarse_vM_usk_uok, Cov_circle = eng.main_crossview_sPGO_SLAM2oxts_CVandVOweights_Cov_Ab(
                    vT_us0_usk_org, L, vT_us0_usk, vS_0k,
                    fine_valid, coarse_vM_usk_uok, corr_vP_usk_uok, cv_vW_usk_uok, corr_vMask_usk_uok,
                    i, nargout=4)

                vT_us0_usk = np.array(vT_us0_usk)
                vS_0k = np.array(vS_0k)
                coarse_vM_usk_uok = np.array(coarse_vM_usk_uok)
                Cov_circle = np.array(Cov_circle)
                # update SLAM trajectory
                dataloader.set_SLAM_pose(vT_us0_usk)


            if i % 20 == 0:
                print(i)


    scio.savemat(os.path.join(save_path2, 'paperFig_' + dataset_save + '_ablation_updateTrajectory_now_spatialboundCheck.mat'),
                 {'gt_lons': gt_lons, 'gt_lats': gt_lats, 'gt_oriens': gt_oriens,
                  'pred_lats': pred_lats, 'pred_lons': pred_lons, 'pred_oriens': pred_oriens,
                  # 'fine_lats': fine_lats, 'fine_lons': fine_lons, 'fine_oriens': fine_oriens,
                  'fine_valid': fine_valid, 'vT_us0_usk': vT_us0_usk, 'coarse_vM_usk_uok': coarse_vM_usk_uok,
                  'vS_0k': vS_0k,
                  'weight_vM_usk_uok': weight_vM_usk_uok,
                  'corr_vP_usk_uok': corr_vP_usk_uok,
                  'corr_vW_usk_uok': corr_vW_usk_uok,
                  'cv_vW_usk_uok':cv_vW_usk_uok,
                  'corr_vMask_usk_uok':corr_vMask_usk_uok,
                  # 'corr_vPaxis_usk_uok': corr_vPaxis_usk_uok
                  })
    return





def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--resume', type=int, default=0, help='resume the trained model')
    parser.add_argument('--test', type=int, default=1, help='test with trained model')

    parser.add_argument('--epochs', type=int, default=5, help='number of training epochs')

    parser.add_argument('--lr', type=float, default=1e-4, help='learning rate')  # 1e-2

    parser.add_argument('--rotation_range', type=float, default=10., help='degree')
    parser.add_argument('--shift_range_lat', type=float, default=20., help='meters')
    parser.add_argument('--shift_range_lon', type=float, default=20., help='meters')

    parser.add_argument('--batch_size', type=int, default=1, help='batch size')

    parser.add_argument('--level', type=int, default=3, help='2, 3, 4, -1, -2, -3, -4')
    parser.add_argument('--N_iters', type=int, default=2, help='any integer')

    parser.add_argument('--Optimizer', type=str, default='TransV1G2SP', help='')

    parser.add_argument('--proj', type=str, default='CrossAttn', help='geo, CrossAttn')

    parser.add_argument('--use_uncertainty', type=int, default=1, help='0 or 1')

    args = parser.parse_args()

    return args



def getSavePath(args):
    model_root = '/media/u1111398/yujiao/dataset/Boost3DoFAccuracyModel/'
    save_path = model_root + '/ModelsKitti/3DoF/'\
                + 'lat' + str(args.shift_range_lat) + 'm_lon' + str(args.shift_range_lon) + 'm_rot' + str(
        args.rotation_range) \
                + '_Nit' + str(args.N_iters) + '_' + str(args.Optimizer) + '_' + str(args.proj)

    if args.use_uncertainty:
        save_path = save_path + '_Uncertainty'

    if not os.path.exists(save_path):
        print('not exist save_path:', save_path)
        import sys
        sys.exit()
        # os.makedirs(save_path)

    print('load_model_path:', save_path)

    return save_path


if __name__ == '__main__':

    if torch.cuda.is_available():
        device = torch.device("cuda:0")
    else:
        device = torch.device("cpu")

    np.random.seed(2022)

    args = parse_args()

    mini_batch = args.batch_size

    save_path = getSavePath(args)
    net = Model(args)
    net.to(device)

    if args.test:
        net.load_state_dict(torch.load(os.path.join(save_path, 'model_4.pth')), strict=False)

        Datasets = ['2011_10_03/2011_10_03_drive_0042_sync/']

        for dataset_i in Datasets:
            print('load data: ', dataset_i)
            test1(net, args, dataset_i, save_path, epoch=0)

