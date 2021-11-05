import numpy as np
from numpy import arctan, arctan2, sin, cos, pi, sqrt, arccos

l1 = [0, 0, 0.5]
l2 = [0,0, 0.5]
l3 = [0, 0, 0.1]
l21 = [0, -0.5, 0]
l32 = [0, -0.1, 0]
lee = [0, 0, 0.01] #length between wrist joints

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


def FK_solve(q, flag):
  '''
  Solves forward kinematics for given angles.
  :param q: List of input angles, if called from transform_base, the transformation from global 
            to robot frame will be in q[6]
  :param flag: Indicates output
  :returns: If flag = "ee": 4x4 homogenous transformation matrix from base frame to end effector frame 
            If flag = "full": All transformation matrices from base to end effector frame
  '''
  print("Input angles:", q)
  
  #Rotation about z, y, y, x, z, x (all global axis)
  t = [np.matmul(rotz(q[0]), trans(l1)), np.matmul(roty(q[1]), trans(l2)), np.matmul(roty(q[2]), trans(l3)), rotx(q[3]), rotz(q[4]), rotz(q[5])]
  t2 = [np.matmul(np.matmul(rotz(q[0]), trans(l1)), rotx(-pi/2)),
       np.matmul(rotz(q[1]), trans(l2)),
       np.matmul(np.matmul(np.matmul(rotz(q[2]), trans(l3)), rotx(pi/2)), roty(pi/2)),
       np.matmul(rotz(q[3]), roty(-pi/2)),
       np.matmul(rotz(q[4]), roty(pi/2)),
       rotz(q[5])]
  transformations = []
  
  start = trans(0)
  #transform
  if len(q) > 6:
    start = np.matmul(start, q[6])
    transformations.append(start)
  
  start = np.matmul(np.matmul(start, rotz(q[0])), trans(l1))
  transformations.append(start)
  start = np.matmul(np.matmul(start, roty(q[1])), trans(l2))
  transformations.append(start)
  start = np.matmul(np.matmul(start, roty(q[2])), trans(l3))
  transformations.append(start)
  start = np.matmul(start, rotx(q[3]))
  transformations.append(start)
  start = np.matmul(start, rotz(q[4]))
  transformations.append(start)
  start = np.matmul(start, rotx(q[5]))
  transformations.append(start)

  if flag == "ee":
    return transformations[len(transformations)-1]
  print(transformations[len(transformations)-1])
  return transformations
    

def transform_base(vec,q):
  '''
  Places Robot arm somewhere based on input.
  :param trans: Vector with translation values for x, y, z
  :param q: The rotation angles of joints 1-6
  :returns: End effector position relative to global frame
  '''
  newBase = trans(vec)
  q.append(newBase)
  return FK_solve(q, "ee")

def IK_solve(base_frame, ee_frame):
  '''
    
  '''
  print(ee_frame)
  px = ee_frame[0][3]
  py = ee_frame[1][3]
  pz = ee_frame[2][3]

  #q1, q2, q3
  theta_1 = []
  theta_2 = []
  theta_3 = []

  theta_11 = arctan2(py, px)
  theta_12  = arctan2(py, px) + pi
  if theta_11 > pi:
    theta_11 = theta_11 - pi
  if theta_12 > pi:
    theta_12 = theta_12 - 2 * pi
  theta_1.append(theta_11)
  theta_1.append(theta_12)
  print("q1 =",theta_1)

  r = sqrt((px ** 2) + (py ** 2))
  print("r =", r)
  pz2 = pz - 0.5 #0.5 = length of l1
  print("pz2 =",pz2)
  s = sqrt((r ** 2) + (pz2 ** 2))
  print("s =",s)
  alpha = arctan2((pz2 ** 2) , r)
  print("alpha =", alpha)
  #print(round(((s ** 2) + (0.5 ** 2) - (0.1 ** 2)) / (2 * 0.5 * s), 10))
  beta =  arccos(round(((s ** 2) + (0.5 ** 2) - (0.1 ** 2)) / (2 * 0.5 * s), 10))
  print("beta =", beta)
  theta_2.append(pi/2 - alpha - beta) 
  theta_2.append(pi/2 - alpha + beta)
  print("q2 =", theta_2)  

  gamma = arccos(round(((0.5 ** 2) + (0.1 ** 2) - (s ** 2)) / (2 * 0.5 * 0.1), 9))
  print("gamma =", gamma)
  theta_31 =  pi - gamma
  theta_32 =  pi + gamma
  #print(theta_31)
  if theta_31 >= pi:
    theta_31 = (theta_31 - pi) * (-1)
  if theta_32 >= pi:
    theta_32 = (theta_32 - pi) * (-1)

  theta_3.append(theta_31)  
  theta_3.append(theta_32)
  print("q3 =", theta_3)

  #q4, q5, q6
  #print(np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(rotz(0), trans(l1)), roty(pi)), trans(l2)), roty(0)), trans(0)))
  #print()
  theta_4 =[]
  theta_5 = []
  theta_6 = []
  for q1 in theta_1:
    for q2 in theta_2:
      for q3 in theta_3:
        r03 = np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(rotz(q1), trans(l1)), roty(pi)), trans(l2)), roty(q3)), trans(q3))
        r03 = np.matrix([r03[0][:3], r03[1][:3], r03[2][:3]])
        r06 = np.matrix([ee_frame[0][:3], ee_frame[1][:3], ee_frame[2][:3]])
        r03t = np.matrix.transpose(r03)
        r36 = np.matmul(r03t, r06)
        r36a = np.squeeze(np.asarray(r36))
        #q5
        c5 = r36a[0][0]
        theta_5.append(arccos(c5))
        #q4
        
  print("q5 =", theta_5)      
        
    



IK_solve(0,FK_solve([pi/4, pi, pi/2 ,pi/2, pi/2,0], "ee"))
#IK_solve(0, transform_base([0,0,0], [pi,pi/2,pi/2,pi/2,pi,pi]))
#IK_solve(0, FK_solve([0,0,0,0,0,0],"ee"))
#print(FK_solve([0,0, pi/2,0,0,0], "ee"))

