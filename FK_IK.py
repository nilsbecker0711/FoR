import numpy as np
from numpy import sin, cos, pi


def rotz(theta):
    rot = np.zeros((4,4),dtype=np.float64)
    rot[3,3] = 1
    rot[0,:] = [cos(theta), -sin(theta),0,0]
    rot[1,:] = [sin(theta), cos(theta),0,0]
    rot[2,2] = 1

    return rot

def rotx(theta): 
  rot = np.zeros((4,4),dtype=np.float64)
  rot[3,3] = 1
  rot[0,0] = 1
  rot[1,:] = [0, cos(theta), -sin(theta),0]
  rot[2,:] = [0, sin(theta), cos(theta),0]

  return rot

def trans(vector):
  mat = np.zeros((4,4),dtype=np.float64)
  mat[0:3,3] = vector
  np.fill_diagonal(mat,1)  

  return mat


def FK_solve(q, flag):
  '''
  Solves forward kinematics for given angles.
  :param q: List of input angles
  :param flag: Indicates output
  :returns: If flag = "ee": 4x4 homogenous transformation matrix from base frame to end effector frame 
            If flag = "full": All transformation matrices from base to end effector frame
  '''
  l1 = [0, 0, 0.5]
  l2 = l3 = [0.5 , 0 , 0]
  l4 = [0.1, 0, 0]
  t = [trans(l1), rotz(q[0]), trans(l2), rotx(90), rotz(q[1]), trans(l3), rotz(q[2]), trans(l4), rotz(q[3]), rotx(q[4]), rotz(q[5])]
  start = trans(0)
  for matrix in t:
    start = np.matmul(start, matrix)
  return start
    

def transform_base(trans):
  '''
  Places Robot arm somewhere based on input.
  :param trans: Determines 
  :returns: New Base of robot
  '''
  pass

def IK_solve(base_frame, ee_frame):
  '''
    
  '''
  pass

print(FK_solve([0,0,0,0,0,0], 1))