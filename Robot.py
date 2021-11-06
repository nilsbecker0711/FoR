import numpy as np
from numpy import arctan, arctan2, sin, cos, pi, sqrt, arccos

class Robot:
    def __init__(self, q, l1, l2, l3):
        self.q = q
        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.base = trans([0, 0, 0])

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
  rot[0,:] = [round(cos(theta), 10), 0, round(sin(theta), 10),0]
  rot[2,:] = [round(-sin(theta), 10), 0, round(cos(theta), 10),0]

  return rot

def trans(vector):
  mat = np.zeros((4,4),dtype=np.float64)
  mat[0:3,3] = vector
  np.fill_diagonal(mat,1)  

  return mat


def FK_solve(self):
    '''
    Solves FK for input angles and link lenghts of robot.
    :returns: 4x4 matrix of end effector
    '''
    q = self.q
    start = self.base
    start = np.matmul(np.matmul(start, rotz(q[0])), trans(self.l1))
    start = np.matmul(np.matmul(start, roty(q[1])), trans(self.l2))
    start = np.matmul(np.matmul(start, roty(q[2])), trans(self.l3))
    start = np.matmul(start, rotx(q[3]))
    start = np.matmul(start, rotz(q[4]))
    start = np.matmul(start, rotx(q[5]))
    return start
         