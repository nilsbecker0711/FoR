import numpy as np
from numpy import pi
import Becker_Nils_HA1 as ha



def jacobian(frames):
    #Rotation axes: Z, Y, Y, X, Z, X
    #Ji: [Ui-1 x (On - Oi-1), Ui-1]
    n = len(frames) - 1
    on = frames[n][:3, 3]

    #z
    u0 = np.array([0, 0, 1])
    o0 = np.array([0, 0, 0])
    j1 = np.array([np.cross(u0, on - o0), u0]).reshape(6,)

    #y
    u1 = frames[0][:3, 1]
    o1 = frames[0][:3, 3]
    j2 = np.array([np.cross(u1, on - o1), u1]).reshape(6,)

    #y
    u2 = frames[1][:3, 1]
    o2 = frames[1][:3, 3]
    j3 = np.array([np.cross(u2, on - o2), u2]).reshape(6,)

    #x
    u3 = frames[2][:3, 0]
    o3 = frames[2][:3, 3]
    j4 = np.array([np.cross(u3, on - o3), u3]).reshape(6,)

    #z
    u4 = frames[3][:3, 2]
    o4 = frames[3][:3, 3]
    j5 = np.array([np.cross(u4, on - o4), u4]).reshape(6,)

    #x
    u5 = frames[4][:3, 0]
    o5 = frames[4][:3, 3]
    j6 = np.array([np.cross(u5, on - o5), u5]).reshape(6,)

    j = np.transpose(np.matrix([j1, j2, j3, j4, j5, j6]))

    return j

def check_singular(j):
    '''
    The matrix is singular if the determinant is 0.
    :returns: True if j is singular, False otherwise
    '''
    if (np.linalg.det(j)):
        return False
    return True

def cartesian_vel(j, q_dot):
    '''
    :return: j*q_dot
    '''
    return np.transpose(np.matmul(j, q_dot))


#test
'''
#get a forward kinematics solution
l1 = [0, 0, 0.5]
l2 = [0,0, 0.5]
l3 = [0, 0, 0.1]
q = [pi/4, pi/2, pi/2, pi/2, pi/4, pi]
frames = ha.FK_solve(q, "full")
for frame in frames:
    #print(frame)
    pass
j = jacobian(frames)
print(j)
print(check_singular(j))
q_dot = np.array([1,2,3,4,5,6]).reshape(6,)
print(cartesian_vel(j, q_dot))
'''