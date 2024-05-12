import random

import numpy as np
import os
from PIL import Image
from torch.utils.data import Dataset

import torch
import pandas as pd
import utils
import torchvision.transforms.functional as TF
from torchvision import transforms
import torch.nn.functional as F

from torch.utils.data import DataLoader
from torchvision import transforms

from dataLoader.utils_pykitti import load_oxts_packets_and_poses, euler_from_rotation_matrix, euler_to_rotation_matrix

visualize = False

# from PIL import Image

root_dir = '/media/u1111398/yujiao/dataset/Kitti1'
# dataset = '2011_10_03/2011_10_03_drive_0027_sync/'

test_csv_file_name = 'test.csv'
ignore_csv_file_name = 'ignore.csv'
satmap_dir = 'satmap'
grdimage_dir = 'raw_data'
left_color_camera_dir = 'image_02/data'  # 'image_02\\data' #
right_color_camera_dir = 'image_03/data'  # 'image_03\\data' #
oxts_dir = 'oxts/data'  # 'oxts\\data' #
# depth_dir = 'depth/data_depth_annotated/train/'

GrdImg_H = 256  # 256 # original: 375 #224, 256
GrdImg_W = 1024  # 1024 # original:1242 #1248, 1024
GrdOriImg_H = 375
GrdOriImg_W = 1242
num_thread_workers = 2




        
        



class SatGrdDatasetTest():
    def __init__(self, root, file,
                 transform=None, shift_range_lat=20, shift_range_lon=20, rotation_range=10):
        self.root = root

        self.meter_per_pixel = utils.get_meter_per_pixel(scale=1)
        self.shift_range_meters_lat = shift_range_lat  # in terms of meters
        self.shift_range_meters_lon = shift_range_lon  # in terms of meters
        self.shift_range_pixels_lat = shift_range_lat / self.meter_per_pixel  # shift range is in terms of meters
        self.shift_range_pixels_lon = shift_range_lon / self.meter_per_pixel  # shift range is in terms of meters

        # self.shift_range_meters = shift_range  # in terms of meters

        self.rotation_range = rotation_range  # in terms of degree

        self.skip_in_seq = 2  # skip 2 in sequence: 6,3,1~
        if transform != None:
            self.satmap_transform = transform[0]
            self.grdimage_transform = transform[1]

        self.pro_grdimage_dir = 'raw_data'

        self.satmap_dir = satmap_dir


        # load SLAM Keyframe poses
        self.day_dir = file[:10]
        self.date_dir = file[11:37]
        self.file = file
        self.slam_root = '/home/users/u1111398/YanhaoZhang/tests_kitti_' + self.date_dir
        with open(self.slam_root + '/' + 'keyframe_name.txt', 'r') as f:
            keyframe_name = f.readlines()
        keyframe_name = [p[:-1] for p in keyframe_name]
        self.keyframe_name = keyframe_name
        # keyframe_ids = # the id of each keyframe within each
        num_kf = len(keyframe_name)
        self.num_kf = num_kf

        # let a, c denotes cam0 cam2
        slam_T_a0ak = np.zeros((num_kf, 4,4))     # rectified cam0 frame0
        with open(self.slam_root + '/' + 'keyframe_pose.txt', 'r') as f:
            keyframe_pose = f.readlines()
            for k in range(num_kf):
                pose = keyframe_pose[k]
                content = pose.split(' ')
                for kk in range(12):
                    i = kk // 4
                    j = kk % 4
                    slam_T_a0ak[k, i, j] = float(content[kk])    # first load the SLAM poses in cam0 frame0: T_a0ak
                slam_T_a0ak[k, 3, 3] = 1

        # oxts poses are body pose in world frame: T_wbk
        # the transformation between body frame to cam2 frame are:
        # Eq(8) in https://www.cvlibs.net/publications/Geiger2013IJRR.pdf
        # T_c0bk = T_c0a0 @ T_rect_a0 @ T_a0v @ T_vb @ T_b0bk


        # load camera intrinsics and extrinsics
        # let a b c d refers to cam0 cam1 cam2 cam3
        T_aa0 = np.identity(4)   # the translation must be zero (reference is cam0, only rotation is needed)   # from rectified cam0 to cam0
        T_c0a0 = np.identity(4)  # the rotation must be zero (already rectified)   # from rectified cam2 to rectified cam0
        calib_file_name = os.path.join(self.root, grdimage_dir, self.day_dir, 'calib_cam_to_cam.txt')
        with open(calib_file_name, 'r') as f:
            lines = f.readlines()
            for line in lines:
                # left gray camera rectification (not very useful)
                if 'R_rect_00' in line:
                    # get 3*3 matrix from R_rect_**:
                    items = line.split(':')
                    valus = items[1].strip().split(' ')
                    R = np.identity(3)
                    for kk in range(9):
                        i = kk // 3
                        j = kk % 3
                        R[i, j] = float(valus[kk])
                    T_aa0[:3, :3] = R               # this small rotation is also ignored

                # left gray camera translation
                # if 'P_rect_00' in line:
                #     items = line.split(':')
                #     valus = items[1].strip().split(' ')
                #     P_r0r0 = np.zeros((3,4))
                #     for kk in range(12):
                #         i = kk // 4
                #         j = kk % 4
                #         P_r0r0[i, j] = float(valus[kk])
                #     K = P_r0r0[:3,:3]
                #     T_00[:3,3] = np.linalg.solve(K, P_00[:3,3])

                # left color camera k matrix
                if 'P_rect_02' in line:
                    # get 3*3 matrix from P_rect_**:
                    items = line.split(':')
                    valus = items[1].strip().split(' ')
                    fx = float(valus[0]) * GrdImg_W / GrdOriImg_W
                    cx = float(valus[2]) * GrdImg_W / GrdOriImg_W
                    fy = float(valus[5]) * GrdImg_H / GrdOriImg_H
                    cy = float(valus[6]) * GrdImg_H / GrdOriImg_H
                    left_camera_k = [[fx, 0, cx], [0, fy, cy], [0, 0, 1]]
                    self.cam2_K = torch.from_numpy(np.asarray(left_camera_k, dtype=np.float32))

                    P_c0a0 = np.zeros((3, 4))
                    for kk in range(12):
                        i = kk // 4
                        j = kk % 4
                        P_c0a0[i, j] = float(valus[kk])
                    K = P_c0a0[:3,:3]
                    T_c0a0[:3,3] = np.linalg.solve(K, P_c0a0[:3,3])

                    break

        # transform camera pose from cam0 to cam2
        T_a0c0 = np.linalg.inv(T_c0a0)
        slam_T_c0ck = np.zeros((num_kf, 4, 4))
        for k in range(num_kf):
            slam_T_c0ck[k, :, :] =   T_c0a0 @  slam_T_a0ak[k, :, :]  @ T_a0c0    # Tb0b1 = Tbc * Tc0c1 * Tcb
        # but there is no need to do the transformation, since the rotation is identity.
        # another way to look is that the trajectory is w.r.t. the initial frame. The global transformation does not affect.
        # todo check with other dataset

        # load corresponding oxts GPS data,
        # transform to TUM frame, and transform to keyframe0 frame
        oxts_file_names = []
        for k in range(num_kf):
            oxts_file = os.path.join(self.root, grdimage_dir, file, oxts_dir,
                                     keyframe_name[k].lower().replace('.png', '.txt'))
            oxts_file_names.append(oxts_file)
        self.oxts_vTb0bk, self.oxts_vr_wbk, self.oxts_Twb0 = load_oxts_packets_and_poses(oxts_file_names)    # self.oxts_vr_wbk is just for debug



        # relative pose from body to cam2
        # load relative pose from body to cam0
        T_vb = np.identity(4)    # from imu (body) to velody
        calib_file_name = os.path.join(self.root, grdimage_dir, self.day_dir, 'calib_imu_to_velo.txt')
        with open(calib_file_name, 'r') as f:
            lines = f.readlines()
            for line in lines:
                # left gray camera rectification (not very useful)
                if 'R' in line:
                    # get 3*3 matrix from R_rect_**:
                    items = line.split(':')
                    valus = items[1].strip().split(' ')
                    R = np.identity(3)
                    for kk in range(9):
                        i = kk // 3
                        j = kk % 3
                        R[i, j] = float(valus[kk])
                    T_vb[:3, :3] = R  # this small rotation is also ignored

                if 'T' in line:
                    items = line.split(':')
                    valus = items[1].strip().split(' ')
                    T = np.zeros(3)
                    for kk in range(3):
                        T[kk] = float(valus[kk])
                    T_vb[:3, 3] = T  # this small rotation is also ignored


        T_av = np.identity(4)      # from velody to cam0 (before rectification)
        calib_file_name = os.path.join(self.root, grdimage_dir, self.day_dir, 'calib_velo_to_cam.txt')
        with open(calib_file_name, 'r') as f:
            lines = f.readlines()
            for line in lines:
                # left gray camera rectification (not very useful)
                if 'R' in line:
                    # get 3*3 matrix from R_rect_**:
                    items = line.split(':')
                    valus = items[1].strip().split(' ')
                    R = np.identity(3)
                    for kk in range(9):
                        i = kk // 3
                        j = kk % 3
                        R[i, j] = float(valus[kk])
                    T_av[:3, :3] = R  # this small rotation is also ignored

                if 'T' in line:
                    items = line.split(':')
                    valus = items[1].strip().split(' ')
                    T = np.zeros(3)
                    for kk in range(3):
                        T[kk] = float(valus[kk])
                    T_av[:3, 3] = T  # this small rotation is also ignored



        # transform the oxts from body frame0 to cam2 frame0 (not the trajectory)
        # T_c0bk = T_c0a0 @ T_rect_a0 @ T_a0v @ T_vb @ T_b0bk
        # T_c0a0 = np.linalg.inv(T_a0c0)
        T_a0a = np.linalg.inv(T_aa0)

        # transform the camera frame in the same direction of body frame
        # this is for calculating the yaw angle from rotation matrix. Or the result can be different.
        Tcu = np.array([[0,-1,0,0], [0,0,-1,0], [1,0,0,0], [0,0,0,1]])
        # Tcu = np.identity(4)
        Tuc = np.linalg.inv(Tcu)
        self.slam_vTu0uk = np.zeros((num_kf, 4, 4))
        for k in range(num_kf):
            self.slam_vTu0uk[k, :, :] = Tuc @  slam_T_c0ck[k, :, :] @ Tcu    # Tb0b1 = Tbc * Tc0c1 * Tcb




        Tub = Tuc @ T_c0a0 @ T_a0a @ T_av @ T_vb
        Tbu = np.linalg.inv(Tub)


        # rewrite data. It seems that the network fits Yujiao's configuration
        Tbu[:3,:3] = np.identity(3)
        Tbu[0, 3] = utils.CameraGPS_shift_left[0]
        Tbu[1, 3] = - utils.CameraGPS_shift_left[1]
        Tub = np.linalg.inv(Tbu)


        self.oxts_vTu0uk = np.zeros((num_kf, 4, 4))
        for k in range(num_kf):
            Tb0bk = self.oxts_vTb0bk[k, :, :]
            self.oxts_vTu0uk[k, :, :] = Tub @ Tb0bk @ Tbu    # oxts trajectory of cam2

        self.oxts_Twu0 = self.oxts_Twb0 @ Tbu
        self.Tbu = Tbu
        # self.Tub = np.linalg.inv(Tbu)
        # self.oxts_yaw_wb0 = self.oxts_vr_wbk[0,2]

        ''' load other SLAM data to fine tune SLAM result '''
        # load SLAM projection matrix (cam0 and cam1)


        # the estimated trajectory is T_u0_uk, The original SLAM is in T_a0_ak
        self.T_u0a0 = Tuc @ T_c0a0

        # we need to update the SLAM and map points
        self.slam_vTu0uk_org = self.slam_vTu0uk






    def __len__(self):
        return len(self.keyframe_name)

    def get_file_list(self):
        return self.keyframe_name


    def load_L_from_SLAM(self):
        # Laplacian matrix with elements representing the number of co-visable features
        self.L = np.zeros((self.num_kf, self.num_kf))
        for k in range(self.num_kf):

            # load matched 2D features within a local map
            feat_file_name = self.slam_root + '/' + 'keyframe_featMap_match_' + self.keyframe_name[k] + '_all.txt'
            with open(feat_file_name, 'r') as f:
                lines = f.readlines()

                num_obs_k = len(lines)
                feat_all_k = np.zeros((num_obs_k, 2))  # pts id and keyframe id
                for i in range(num_obs_k):
                    content = lines[i].split(' ')
                    feat_all_k[i, 0] = int(content[0])
                    feat_all_k[i, 1] = int(content[1])

                id_obs_k = np.where(feat_all_k[:, 1] == k)  # valid measurement of current keyframe
                index_pts_k = np.unique(
                    feat_all_k[id_obs_k, 0])  # find potential map index (SLAM id), remove duplicated index

                for i in range(k + 1, self.num_kf):
                    id_obs_i = np.where(feat_all_k[:, 1] == i)  # find the id of covisable keyframe
                    index_pts_i = np.unique(feat_all_k[id_obs_i, 0])
                    val = np.unique(np.intersect1d(index_pts_k, index_pts_i))  # potential pts id

                    self.L[k, i] = len(val)  # store the matched number of features

        return self.L



    def save_L(self):
        L_file_name = self.slam_root + '/' + 'keyframe_matchedFeatNum_L.npy'
        with open(L_file_name, 'wb') as f:
            np.save(f, self.L)

    def load_L_saved(self):
        L_file_name = self.slam_root + '/' + 'keyframe_matchedFeatNum_L.npy'

        # if file not exist, we load it from SLAM
        if not os.path.isfile(L_file_name):
            self.load_L_from_SLAM()
            self.save_L()

        else:
            with open(L_file_name, 'rb') as f:
                self.L = np.load(L_file_name)

        return self.L



    def get_L(self):
        return self.L

    def get_SLAM_pose(self, id=None):
        if id == None:
            return self.slam_vTu0uk
        else:
            return self.slam_vTu0uk[id]

    def set_SLAM_pose(self, vTu0uk):
        self.slam_vTu0uk = vTu0uk


    def get_SLAM_pose_org(self, id=None):
        if id == None:
            return self.slam_vTu0uk_org
        else:
            return self.slam_vTu0uk_org[id]


    def get_oxts_pose(self, id=None):
        if id == None:
            return self.oxts_vTu0uk
        else:
            return self.oxts_vTu0uk[id]




    def get_T_u0a0(self):
        return self.T_u0a0


    def getitem(self, idx, fine_vM_uok_usk):
        '''

        Args:
            idx:
            pred_lons: last timestamp prediction
            pred_lats:
            pred_oriens:

        Returns:

        '''

        # file_name = self.file + self.keyframe_name[idx]
        #
        # day_dir = file_name[:10]
        # drive_dir = file_name[:38]
        # image_no = file_name[38:]

        # =================== read satellite map ===================================
        SatMap_name = os.path.join(self.root, self.satmap_dir, self.file, self.keyframe_name[idx])
        with Image.open(SatMap_name, 'r') as SatMap:
            sat_map = SatMap.convert('RGB')


        # the centre of satellite image is calculated from oxts, (when retriving the satellite image),
        oxts_Twbk = self.oxts_Twb0 @ self.oxts_vTb0bk[idx]
        _, _, oxts_yaw_wbk = euler_from_rotation_matrix(oxts_Twbk[:3, :3])
        sat_rot_wbk = sat_map.rotate(-oxts_yaw_wbk / np.pi * 180)  # Image.rotate(angle): rotate image in counter clockwise direction.
        #
        # # the definition is first translation then rotation!
        sat_align_cam = sat_rot_wbk.transform(sat_rot_wbk.size, Image.AFFINE,
                                             (1, 0, self.Tbu[0, 3] / self.meter_per_pixel,
                                              0, 1, -self.Tbu[1, 3] / self.meter_per_pixel),
                                             resample=Image.BILINEAR)
        _, _, oxts_yaw_bu = euler_from_rotation_matrix(self.Tbu[:3, :3])
        sat_rot = sat_align_cam.rotate(-oxts_yaw_bu / np.pi * 180)  # Image.rotate(angle): rotate image in counter clockwise direction.



        # =================== initialize some required variables ============================
        grd_left_imgs = torch.tensor([])
        left_img_name = os.path.join(self.root, self.pro_grdimage_dir, self.file, left_color_camera_dir, self.keyframe_name[idx])
        with Image.open(left_img_name, 'r') as GrdImg:
            grd_img_left = GrdImg.convert('RGB')
            if self.grdimage_transform is not None:
                grd_img_left = self.grdimage_transform(grd_img_left)

        grd_left_imgs = torch.cat([grd_left_imgs, grd_img_left.unsqueeze(0)], dim=0)



        P_us0_usk = self.slam_vTu0uk[idx]         # todo: we use the updated SLAM trajectory
        T_uo0_uok = self.oxts_vTu0uk[idx]

        # T_uok_usk = T_uok_uol @ T_uol_usl @ T_usl_usk
        T_uok_usk = np.linalg.inv(T_uo0_uok) @ P_us0_usk




        gt_shift_x = T_uok_usk[0, 3]  # --> right as positive, parallel to the heading direction
        gt_shift_y = T_uok_usk[1, 3]  # --> up as positive, vertical to the heading direction
        _, _, gt_theta = euler_from_rotation_matrix(T_uok_usk[:3, :3])
        gt_theta = gt_theta / np.pi * 180


        gt_shift_x = gt_shift_x / self.shift_range_meters_lon
        gt_shift_y = gt_shift_y / self.shift_range_meters_lat
        gt_theta = gt_theta / self.rotation_range




        sat_align_init_shift = sat_rot.transform(sat_rot.size, Image.AFFINE,
                                                (1, 0, gt_shift_x * self.shift_range_pixels_lon,
                                                 0, 1, -gt_shift_y * self.shift_range_pixels_lat), # a negative on y direction
                                                resample=Image.BILINEAR)
        sat_align_init = sat_align_init_shift.rotate(-gt_theta * self.rotation_range)


        sat_map = TF.center_crop(sat_align_init, utils.SatMap_process_sidelength)
        # sat_map = np.array(sat_map, dtype=np.float32)

        # transform
        if self.satmap_transform is not None:
            sat_map = self.satmap_transform(sat_map)

        return sat_map, self.cam2_K, grd_left_imgs[0], \
               torch.tensor(-gt_shift_x, dtype=torch.float32).reshape(1), \
               torch.tensor(-gt_shift_y, dtype=torch.float32).reshape(1), \
               torch.tensor(-gt_theta, dtype=torch.float32).reshape(1), \
               self.keyframe_name[idx]

    


def load_test1_data(batch_size, dataset, shift_range_lat=20, shift_range_lon=20, rotation_range=10):
    SatMap_process_sidelength = utils.get_process_satmap_sidelength()

    satmap_transform = transforms.Compose([
        transforms.Resize(size=[SatMap_process_sidelength, SatMap_process_sidelength]),
        transforms.ToTensor(),
    ])

    Grd_h = GrdImg_H
    Grd_w = GrdImg_W

    grdimage_transform = transforms.Compose([
        transforms.Resize(size=[Grd_h, Grd_w]),
        transforms.ToTensor(),
    ])

    # # Plz keep the following two lines!!! These are for fair test comparison.
    # np.random.seed(2022)
    # torch.manual_seed(2022)

    test1_set = SatGrdDatasetTest(root=root_dir, file=dataset,
                            transform=(satmap_transform, grdimage_transform),
                            shift_range_lat=shift_range_lat,
                            shift_range_lon=shift_range_lon,
                            rotation_range=rotation_range)

    return test1_set









