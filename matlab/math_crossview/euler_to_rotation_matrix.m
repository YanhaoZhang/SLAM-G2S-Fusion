function [R] = euler_to_rotation_matrix(roll, pitch, yaw)


Rx = rotx(roll);
Ry = roty(pitch);
Rz = rotz(yaw);
R = Rz * Ry * Rx;

end