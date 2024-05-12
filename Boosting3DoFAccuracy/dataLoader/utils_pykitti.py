"""Provides helper methods for loading and parsing KITTI data."""

from collections import namedtuple

import numpy as np
from PIL import Image

__author__ = "Lee Clement"
__email__ = "lee.clement@robotics.utias.utoronto.ca"

# Per dataformat.txt
OxtsPacket = namedtuple('OxtsPacket',
                        'lat, lon, alt, ' +
                        'roll, pitch, yaw, ' +
                        'vn, ve, vf, vl, vu, ' +
                        'ax, ay, az, af, al, au, ' +
                        'wx, wy, wz, wf, wl, wu, ' +
                        'pos_accuracy, vel_accuracy, ' +
                        'navstat, numsats, ' +
                        'posmode, velmode, orimode')

# Bundle into an easy-to-access structure
OxtsData = namedtuple('OxtsData', 'packet, T_w_imu')


def subselect_files(files, indices):
    try:
        files = [files[i] for i in indices]
    except:
        pass
    return files


def rotx(t):
    """Rotation about the x-axis."""
    c = np.cos(t)
    s = np.sin(t)
    return np.array([[1,  0,  0],
                     [0,  c, -s],
                     [0,  s,  c]])


def roty(t):
    """Rotation about the y-axis."""
    c = np.cos(t)
    s = np.sin(t)
    return np.array([[c,  0,  s],
                     [0,  1,  0],
                     [-s, 0,  c]])


def rotz(t):
    """Rotation about the z-axis."""
    c = np.cos(t)
    s = np.sin(t)
    return np.array([[c, -s,  0],
                     [s,  c,  0],
                     [0,  0,  1]])


def transform_from_rot_trans(R, t):
    """Transforation matrix from rotation matrix and translation vector."""
    R = R.reshape(3, 3)
    t = t.reshape(3, 1)
    return np.vstack((np.hstack([R, t]), [0, 0, 0, 1]))


def read_calib_file(filepath):
    """Read in a calibration file and parse into a dictionary."""
    data = {}

    with open(filepath, 'r') as f:
        for line in f.readlines():
            key, value = line.split(':', 1)
            # The only non-float values in these files are dates, which
            # we don't care about anyway
            try:
                data[key] = np.array([float(x) for x in value.split()])
            except ValueError:
                pass

    return data



def pose_inv(T):
    """
    Args:
        T: 4x4 matrix

    Returns: inverse of T
    """

    R = T[:3,:3]
    R_inv = R.T

    T_inv = np.identity(4)

    T_inv[:3,:3] = R_inv
    T_inv[:3, 3] = -R_inv @ T[:3,3]

    # a=1

    return T_inv


def euler_from_rotation_matrix(R):
    """
    Calculate roll, pitch, yaw from a rotation matrix.
    """
    # yaw (z-axis rotation)
    yaw = np.arctan2(R[1, 0], R[0, 0])

    # pitch (y-axis rotation)
    pitch = np.arctan2(-R[2, 0], np.sqrt(R[2, 1]**2 + R[2, 2]**2))

    # roll (x-axis rotation)
    roll = np.arctan2(R[2, 1], R[2, 2])

    return roll, pitch, yaw


def euler_to_rotation_matrix(roll, pitch, yaw):
    # Use the Euler angles to get the rotation matrix
    Rx = rotx(roll)
    Ry = roty(pitch)
    Rz = rotz(yaw)
    R = Rz.dot(Ry.dot(Rx))
    return R

def pose_from_oxts_packet(packet, scale):
    """Helper method to compute a SE(3) pose matrix from an OXTS packet.
    """
    er = 6378137.  # earth radius (approx.) in meters

    # Use a Mercator projection to get the translation vector
    tx = scale * packet.lon * np.pi * er / 180.
    ty = scale * er * \
        np.log(np.tan((90. + packet.lat) * np.pi / 360.))
    tz = packet.alt
    t = np.array([tx, ty, tz])

    # Use the Euler angles to get the rotation matrix
    # Rx = rotx(packet.roll)
    # Ry = roty(packet.pitch)
    # Rz = rotz(packet.yaw)
    # R = Rz.dot(Ry.dot(Rx))
    R = euler_to_rotation_matrix(packet.roll, packet.pitch, packet.yaw)

    # Combine the translation and rotation into a homogeneous transform
    return R, t


def load_oxts_packets_and_poses(oxts_files):
    """Generator to read OXTS ground truth data.

       Poses are given in an East-North-Up coordinate system 
       whose origin is the first GPS position.
    """
    # Scale for Mercator projection (from first lat value)
    scale = None
    # Origin of the global coordinate system (first GPS position)
    Tw0 = None

    num = len(oxts_files)
    oxts_T0i = np.zeros((num, 4,4))
    oxts_r_wi = np.zeros((num, 3))   # roll pitch yaw

    k = 0
    for filename in oxts_files:
        with open(filename, 'r') as f:
            for line in f.readlines():
                line = line.split()
                # Last five entries are flags and counts
                line[:-5] = [float(x) for x in line[:-5]]
                line[-5:] = [int(float(x)) for x in line[-5:]]

                packet = OxtsPacket(*line)
                if scale is None:
                    scale = np.cos(packet.lat * np.pi / 180.)
                R, t = pose_from_oxts_packet(packet, scale)
                Twi = transform_from_rot_trans(R, t)

                if Tw0 is None:
                    Tw0 = Twi
                    T0w = pose_inv(Tw0)

                # transform oxts to keyframe0 frame
                T0i = T0w @ Twi


        # check inverse euler transformation
        # roll, pitch, yaw = euler_from_rotation_matrix(R)
        # print(roll-packet.roll, pitch-packet.pitch, yaw-packet.yaw)


        # save data
        oxts_r_wi[k, :] = line[3:6]
        oxts_T0i[k,:,:] = T0i
        k=k+1

        # T_w_imu = transform_from_rot_trans(R, t - origin)
        # oxts.append(OxtsData(packet, T_w_imu))

    return oxts_T0i, oxts_r_wi, Tw0


# def load_image(file, mode):
#     """Load an image from file."""
#     return Image.open(file).convert(mode)
#
#
# def yield_images(imfiles, mode):
#     """Generator to read image files."""
#     for file in imfiles:
#         yield load_image(file, mode)
#
#
# def load_velo_scan(file):
#     """Load and parse a velodyne binary file."""
#     scan = np.fromfile(file, dtype=np.float32)
#     return scan.reshape((-1, 4))
#
#
# def yield_velo_scans(velo_files):
#     """Generator to parse velodyne binary files into arrays."""
#     for file in velo_files:
#         yield load_velo_scan(file)
