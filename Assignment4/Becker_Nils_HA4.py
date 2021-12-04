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
    l1,l2,l3,l4,l5,l6 = var("l1 l2 l3 l4 l5 l6")
    lc1,lc2,lc3,lc4,lc5,lc6= var("lc1 lc2 lc3 lc4 lc5 lc6")
    l = [l1,l2,l3,l4,l5,l6]
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
                der.append(pos[i][k].diff(q[j])) # get differentials for all positions for x,y,z
            jac.append(der)
            der = []
        all_jac.append(jac)
        jac = []
    matrix_jac = []
    for i in range(len(all_jac)):
        matrix_jac.append(Matrix([all_jac[i][0],all_jac[i][1],all_jac[i][2],all_jac[i][3],all_jac[i][4],all_jac[i][5]]).T)
    return matrix_jac

def s_calc_angular_vel():
    ''''
    Gives a list of all Jw matrices
    '''
    #get symbolic FK solution
    q1,q2,q3,q4,q5,q6 = var('q1 q2 q3 q4 q5 q6')
    l1,l2,l3,l4,l5,l6 = var("l1 l2 l3 l4 l5 l6")
    t1 = srotz(q1)*strans([0,0,l1])
    t2 = t1*sroty(q2)*strans([0,0,l2])
    t3 = t2*sroty(q3)*strans([0,0,l3])
    t4 = t3*srotz(q4)*strans([0,0,l4])
    t5 = t4*srotx(q5)*strans([0,0,l5])
    t6 = t5*srotz(q6)*strans([0,0,l6]) #EE
    transformations = [t1,t2,t3,t4,t5,t6]

    #get symbolic angular velocities
    rotation_axis = [[0,0,1]]
    rotation_axis.append([t2[1],t2[5],t2[9]]) #y
    rotation_axis.append([t3[1],t3[5],t2[9]]) #y
    rotation_axis.append([t4[2],t4[6],t2[10]]) #z
    rotation_axis.append([t5[0],t5[4],t5[8]]) #x
    rotation_axis.append([t6[2],t6[6],t6[10]]) #z
    
    #get angular velocities matrices
    jac = []
    for i in range(len(rotation_axis)):
        matrix = Matrix(rotation_axis[0])
        for j in range(1, len(rotation_axis)):
            if j <= i:
                matrix =matrix.col_insert(j, Matrix(rotation_axis[j]))
            else:
                matrix = matrix.col_insert(j, Matrix([0,0,0]))  
        jac.append(matrix)
    return jac

def inertia(sym):
    '''
    Computes I matrix
    '''
    ixx1, ixy1, ixz1, iyx1,iyy1,iyz1,izx1,izy1,izz1 = var(" ixx1 ixy1 ixz1 iyx1 iyy1 iyz1 izx1 izy1 izz1")
    ixx2, ixy2, ixz2, iyx2,iyy2,iyz2,izx2,izy2,izz2 = var(" ixx2 ixy2 ixz2 iyx2 iyy2 iyz2 izx2 izy2 izz2")
    ixx3, ixy3, ixz3, iyx3,iyy3,iyz3,izx3,izy3,izz3 = var(" ixx3 ixy3 ixz3 iyx3 iyy3 iyz3 izx3 izy3 izz3")
    ixx4, ixy4, ixz4, iyx4,iyy4,iyz4,izx4,izy4,izz4 = var(" ixx4 ixy4 ixz4 iyx4 iyy4 iyz4 izx4 izy4 izz4")
    ixx5, ixy5, ixz5, iyx5,iyy5,iyz5,izx5,izy5,izz5 = var(" ixx5 ixy5 ixz5 iyx5 iyy5 iyz5 izx5 izy5 izz5")
    ixx6, ixy6, ixz6, iyx6,iyy6,iyz6,izx6,izy6,izz6 = var(" ixx6 ixy6 ixz6 iyx6 iyy6 iyz6 izx6 izy6 izz6")

    i1 = [ixx1, ixy1, ixz1, iyx1,iyy1,iyz1,izx1,izy1,izz1]
    i2 = [ixx2, ixy2, ixz2, iyx2,iyy2,iyz2,izx2,izy2,izz2]
    i3 = [ixx3, ixy3, ixz3, iyx3,iyy3,iyz3,izx3,izy3,izz3]
    i4 = [ixx4, ixy4, ixz4, iyx4,iyy4,iyz4,izx4,izy4,izz4]
    i5 = [ixx5, ixy5, ixz5, iyx5,iyy5,iyz5,izx5,izy5,izz5]
    i6 = [ixx6, ixy6, ixz6, iyx6,iyy6,iyz6,izx6,izy6,izz6]

    inertias = [i1,i2,i3,i4,i5,i6]
    matrices = []
    for values in inertias:
        matrix = Matrix([[values[0],0, 0], [0, values[4], 0], [0, 0, values[8]]])
        matrices.append(matrix)
    return matrices

def s_calc_m():
    '''
    Checks if all symblic calculations are correct
    :returns True if m_q is symmetric -> indicates correct calculations, False otherwise
             m_q matrix
    '''
    q1, q2 , q3, q4, q5, q6 = var("q1 q2 q3 q4 q5 q6") 
    m1, m2, m3, m4, m5, m6 = var("m1 m2 m3 m4 m5 m6")
    q = [q1, q2 , q3, q4, q5, q6]
    m = [m1, m2, m3, m4, m5, m6]
    linear_vel = s_calc_linear_vel(s_calc_positions())
    ang_vel =s_calc_angular_vel()
    s_inertias = inertia(True)
    m_q = Matrix([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    for i in range(len(m)):
        if i == 0 or i == 3 or i == 5:
            rot = srotz(q[i])
        elif i == 1 or i == 2:
            rot = sroty(q[i])
        elif i == 4:
            rot = srotx(q[i])
        rot.col_del(3)
        rot.row_del(3)
        m_q = m_q + m[i] * linear_vel[i].T * linear_vel[i] + (ang_vel[i].T*rot) * s_inertias[i] * (rot.T*ang_vel[i])

    #verify m -> Symmetry check
    if (m_q.is_symmetric()):
        print("M(q) is symmetric -> All calculations seem correct :)")
        return (True, m_q)
    else:
        print("False calculation, M(q) is NOT symmetric!")
        return (False, m_q)

def s_c_symb(m,i,j,k):
    q1, q2 , q3, q4, q5, q6 = var("q1 q2 q3 q4 q5 q6") 
    q = [q1, q2 , q3, q4, q5, q6]
    c = m[i+j*6].diff(q[k]) + m[i+k*6].diff(q[j]) - m[j+6*k].diff(q[i])
    return c

def s_calc_c(m):
    '''
    Calculates the c matrix symbolically. 
    :param i: symbolical m matrix
    '''
    c = Matrix([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    for i in range(6):
        for j in range (6):
            runner=0
            for k in range(6):
                runner = runner + s_c_symb(m,i,j,k)
        c[i+j*6] = runner
    return c

def s_calc_g(m):
    '''
    param m: m1-m6
    '''
    lin_vel = s_calc_linear_vel(s_calc_positions())
    gr = Matrix([[0],[0],[-9.81]]) # gravity : +/- -9.81 in opposite of z-direction
    g = []
    for i in range(len(m)):
        start = 0
        jvs = []
        for j in range(6):
            jvs.append(Matrix([lin_vel[i][j], lin_vel[i][j+6], lin_vel[i][j+12]]))
        for k in range(len(jvs)):
            start = start - (jvs[k].T * (m[k]*gr))[0] #[0]-> bc, 1x1 matrix
        g.append(start)
    g_matrix = Matrix([g[0],g[1],g[2],g[3],g[4],g[5]]).T
    print(g_matrix.shape)

def euler_lagrange(q,l,lc,m):
    mq = s_calc_m()
    c = s_calc_c(mq[1]) #this may take a while
    g = s_calc_g(m)
    #Compelty symbolical model
    dq1,dq2,dq3,dq4,dq5,dq6 = var("dq1 dq2 dq3 dq4 dq5 dq6")
    dq = Matrix([dq1,dq2,dq3,dq4,dq5,dq6]).T
    ddq1,ddq2,ddq3,ddq4,ddq5,ddq6 = var("ddq1 ddq2 ddq3 ddq4 ddq5 ddq6")
    ddq = Matrix([ddq1,ddq2,ddq3,ddq4,ddq5,ddq6])
    s_tau = mq[1]*ddq + c*dq + g
    print(s_tau.shape)


euler_lagrange(1,1,1,[1,1,1,1,1,1])