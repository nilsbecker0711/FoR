from sympy import *
import sympy as sy
import numpy as np
from numpy import sin, cos

def srotz(theta):
    rot = Matrix([[sy.cos(theta), -sy.sin(theta),0,0],
                [sy.sin(theta), sy.cos(theta),0,0],
                [0,0,1,0],[0,0,0,1]])

    return rot  

def srotx(theta): 
    rot = Matrix([[1,0,0,0],
                [0, sy.cos(theta), -sy.sin(theta),0],
                [0, sy.sin(theta), sy.cos(theta),0],
                [0,0,0,1]])
    return rot

def sroty(theta): 
    rot = Matrix([[sy.cos(theta), 0, sy.sin(theta),0],
                [0,1,0,0],
                [-sy.sin(theta), 0, sy.cos(theta),0],
                [0,0,0,1]])
    return rot

def strans(val):
    mat = Matrix([[1,0,0,val[0]],
                [0,1,0,val[1]],
                [0,0,1,val[2]],
                [0,0,0,1]])
    return mat

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


def FK_solve(q,l):
 
  transformations = []
  
  start = trans(0)
  #transform
  if len(q) > 6:
    start = np.matmul(start, q[6])
    transformations.append(start)
  l1,l2,l3,l4,l5,l6 = l
  #Rotation about z, y, y, x, z, x (all global axis)
  start = np.matmul(np.matmul(start, rotz(q[0])), trans(l1))
  transformations.append(start)
  start = np.matmul(np.matmul(start, roty(q[1])), trans(l2))
  transformations.append(start)
  start = np.matmul(np.matmul(start, roty(q[2])), trans(l3))
  transformations.append(start)
  start = np.matmul(np.matmul(start, roty(q[3])), trans(l4))
  transformations.append(start)
  start = np.matmul(np.matmul(start, roty(q[4])), trans(l5))
  transformations.append(start)
  start = np.matmul(np.matmul(start, roty(q[5])), trans(l6))
  transformations.append(start)

  return transformations

def s_calc_positions():
    '''
    Calculates the CoM symbolically.
    :param lc: List contaning the vectors from joints to CoM
    '''
    q1,q2,q3,q4,q5,q6 = var('q1 q2 q3 q4 q5 q6')
    l1,l2,l3,l4,l4,l6 = var("l1 l2 l3 l4 l5 l6")
    lc1,lc2,lc3,lc4,lc5,lc6= var("lc1 lc2 lc3 lc4 lc5 lc6")
    l = [l1,l2,l3,l4,l4,l6]
    lc = [lc1,lc2,lc3,lc4,lc5,lc6]
    for i in range(len(l)):
        l[i] = strans([0,0,l[i]]).col(-1)
        lc[i] = i = strans([0,0,lc[i]]).col(-1)
      
    start = strans([0,0,0])
    rot = [srotz(q1), sroty(q2), sroty(q3), srotz(q4), srotx(q5), srotz(q6)]
    cm1 = rot[0]*lc[0]
    pos2 = rot[0]*l[0]
    cm2 = rot[0]*rot[1]*strans(pos2)
    pos3 = cm2*l[1]
    cm2 = cm2*lc[1]

    for i in range(3):
        start = start*rot[i]
    start = start * strans(pos3)
    pos4 = start*l[2]
    cm3 = start *lc[2]

    start = strans([0,0,0])
    for i in range(4):
        start = start*rot[i]
    start = start * strans(pos4)
    pos5 = start*l[3]
    cm4 = start *lc[3]

    start = strans([0,0,0])
    for i in range(5):
        start = start*rot[i]
    start = start * strans(pos5)
    pos6 = start*l[4]
    cm5 = start *lc[4]

    start = strans([0,0,0])
    for i in range(6):
        start = start*rot[i]
    start = start * strans(pos6)
    posEE = start*l[5] #Not needed but nice to have
    cm6 = start *lc[5]

    cm = [cm1, cm2, cm3, cm4, cm5, cm6]
    return cm

def s_calc_linear_vel(pos):
    '''
    Calculates the Jv Matrices
    '''
    q1,q2,q3,q4,q5,q6 = var('q1 q2 q3 q4 q5 q6')
    q = [q1,q2,q3,q4,q5,q6]
    all_jac =[]
    jac = []
    der = []
    for i in range(len(pos)):
        for j in range(len(q)):
            for k in range(3):
                der.append(pos[i][k].diff(q[j]))
            jac.append(der)
            der = []
        all_jac.append(jac)
        jac = []
    matrix_jac = []
    for i in range(len(all_jac)):
        matrix_jac.append(Matrix([all_jac[i][0],all_jac[i][1],all_jac[i][2],all_jac[i][3],all_jac[i][4],all_jac[i][5]]))
    return matrix_jac

def euler_lagrange(q,l,lc,m):
    calc_m(q,l,lc,m)

def calc_m(q,l, lc, m):
    linear_vel = s_calc_linear_vel(s_calc_positions())


euler_lagrange(1,1,1,1)