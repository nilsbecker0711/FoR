import numpy as np
from numpy import arctan, arctan2, sin, cos, pi

l1 = [0,0,0.5]
l2 = [0, 0, 0.5]
l3 = [0, 0, 0.1]

def rotz(theta):
    rot = np.zeros((4,4),dtype=np.float64)
    rot[3,3] = 1
    rot[0,:] = [round(cos(theta), 10), round(-sin(theta), 10),0,0]
    rot[1,:] = [round(sin(theta), 10), round(cos(theta), 10),0,0]
    rot[2,2] = 1

    return rot

def rotx(theta): 
  rot = np.zeros((4,4),dtype=np.float64)
  rot[3,3] = 1
  rot[0,0] = 1
  rot[1,:] = [0, round(cos(theta), 10), round(-sin(theta), 10),0]
  rot[2,:] = [0, round(sin(theta), 10), round(cos(theta), 10),0]

  return rot

def roty(theta): 
  rot = np.zeros((4,4),dtype=np.float64)
  rot[3,3] = 1
  rot[1,1] = 1
  rot[0,:] = [0, round(cos(theta), 10), round(sin(theta), 10),0]
  rot[2,:] = [0, round(-sin(theta), 10), round(cos(theta), 10),0]

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
  print("Input angles:", q)
  #Rotation about z, y, y, x, z, x (all global axis)
  t = [rotz(q[0]), trans(l1), roty(q[1]), trans(l2),  roty(q[2]), trans(l3), rotx(q[3]), rotz(q[4]), rotx(q[5])]
  start = trans(0)
  for matrix in t:
    start = np.matmul(start, matrix)
  return start
    

def transform_base(trans,q):
  '''
  Places Robot arm somewhere based on input.
  :param trans: Vector with translation values for x, y, z
  :returns: New Base of robot
  '''
  newBase = trans(trans)

def IK_solve(base_frame, ee_frame):
  '''
    
  '''
  print(ee_frame)
  print(ee_frame[1][3], ee_frame[0][3])
  theta_1 = arctan2(ee_frame[1][3], ee_frame[0][3])
  print("q1 =",theta_1)
  

IK_solve(0,FK_solve([1,pi/2,2,pi,pi,pi], 1))
#IK_solve(0, FK_solve([0,0,0,0,0,0],0))
