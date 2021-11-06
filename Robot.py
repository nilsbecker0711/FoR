import numpy as np
from numpy import arctan, arctan2, sin, cos, pi, sqrt, arccos
import matplotlib.pyplot as plt

class Robot:
    def __init__(self, q, l1, l2, l3):
        self.q = q
        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.transformations = []
        self.base = self.trans([0, 0, 0])

    def rotz(self, theta):
        rot = np.zeros((4,4),dtype=np.float64)
        rot[3,3] = 1
        rot[0,:] = [round(cos(theta), 10), round(-sin(theta), 10),0,0]
        rot[1,:] = [round(sin(theta), 10), round(cos(theta), 10),0,0]
        rot[2,2] = 1

        return rot

    def rotx(self, theta): 
        rot = np.zeros((4,4),dtype=np.float64)
        rot[3,3] = 1
        rot[0,0] = 1
        rot[1,:] = [0, round(cos(theta), 10), round(-sin(theta), 10),0]
        rot[2,:] = [0, round(sin(theta), 10), round(cos(theta), 10),0]

        return rot

    def roty(self, theta): 
        rot = np.zeros((4,4),dtype=np.float64)
        rot[3,3] = 1
        rot[1,1] = 1
        rot[0,:] = [round(cos(theta), 10), 0, round(sin(theta), 10),0]
        rot[2,:] = [round(-sin(theta), 10), 0, round(cos(theta), 10),0]

        return rot

    def trans(self, vector):
        mat = np.zeros((4,4),dtype=np.float64)
        mat[0:3,3] = vector
        np.fill_diagonal(mat,1)  

        return mat


    def FK_solve(self):
        '''
        Solves FK for input angles and link lenghts of robot.
        :returns: 4x4 matrix of end effector
        '''
        self.transformations = []
        q = self.q
        start = self.base
        self.transformations.append(start)
        start = np.matmul(np.matmul(start, self.rotz(q[0])), self.trans(self.l1))
        self.transformations.append(start)
        start = np.matmul(np.matmul(start, self.roty(q[1])), self.trans(self.l2))
        self.transformations.append(start)
        start = np.matmul(np.matmul(start, self.roty(q[2])), self.trans(self.l3))
        self.transformations.append(start)
        start = np.matmul(start, self.rotx(q[3]))
        self.transformations.append(start)
        start = np.matmul(start, self.rotz(q[4]))
        self.transformations.append(start)
        start = np.matmul(start, self.rotx(q[5]))
        self.transformations.append(start)

        return start

    def IK_solve(self, ee_frame, fromFK = False):
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
        

        r = sqrt((px ** 2) + (py ** 2))
        pz2 = pz - self.l1 
        s = sqrt((r ** 2) + (pz2 ** 2))
        alpha = arctan2((pz2 ** 2) , r)
        beta =  arccos(round(((s ** 2) + (self.l2 ** 2) - (self.l3 ** 2)) / (2 * self.l2 * s), 10))
        print("beta =", beta)
        
        theta_2.append(pi/2 - alpha - beta) 
        theta_2.append(pi/2 - alpha + beta)
        

        gamma = arccos(round(((self.l2 ** 2) + (self.l3 ** 2) - (s ** 2)) / (2 * self.l3 * self.l2), 9))
        #gamma = arccos(((px ** 2) + (py ** 2) - (0.5 ** 2) - (0.1 ** 2)) / (0.5 * 0.1))
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

        #q4, q5, q6
        theta_4 =[]
        theta_5 = []
        theta_6 = []
        for q1 in theta_1:
            for q2 in theta_2:
                for q3 in theta_3:
                    r03 = np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(self.rotz(q1), self.trans(self.l1)), self.oty(q2)), self.trans(self.l2)), self.roty(q3)), self.trans(self.q3))
                    r03 = np.matrix([r03[0][:3], r03[1][:3], r03[2][:3]])
                    r06 = np.matrix([ee_frame[0][:3], ee_frame[1][:3], ee_frame[2][:3]])
                    r03t = np.matrix.transpose(r03)
                    r36 = np.matmul(r03t, r06)
                    r36a = np.squeeze(np.asarray(r36))
                    #q5
                    c5 = r36a[0][0]
                    theta_5.append(arccos(c5))
                    #q4
                    s4c5 = r36a[1][0]
                    s4s5 = r36a[2][0]
                
                    theta_4.append(arctan2(s4s5, s4c5))
                
                    #q6
                    s5c6 = r36a[0][1]
                    s5s6 = r36a[0][2]
                
                    theta_6.append(arctan2(s5s6, -s5c6))
                
        print("q1 =", theta_1)
        print("q2 =", theta_2) 
        print("q3 =", theta_3)
        print("q4 =", theta_4)      
        print("q5 =", theta_5)
        print("q6 =", theta_6) 

        all_angles = [theta_1, theta_2, theta_3, theta_4, theta_5, theta_6]

        if not fromFK:
            return all_angles

        else:
            angles =[]
            q = self.q
            #find corresponding angles in q
            i = 0
            for i in range(len(all_angles)):
                angles.append([])
                for angle in all_angles[i]:
                    #rounding, because values change a bit over computations
                    if round(angle, 5) == round(q[i], 5):
                        angles[i] = angle
                        print(f"Found correct angle for q{i}")
                        break
                    elif round(angle, 5) == round(q[i] + pi, 5):
                        angles[i] = angle
                        print(f"Found correct angle for q{i}")
                        break
                    elif round(angle, 5) == round(q[i] - pi, 5):
                        angles[i] = angle
                        print(f"Found correct angle for q{i}")
                        break
            return (all_angles, angles)

    def multiple_FK(self, angles):
        '''
        Solves FK for multiple sets of input angles on robot
        :param angles: List of angle lists
        :returns: ee frame after all transformations are done
        '''
        transformations = len(angles)

    def extract_plot_points(self, serial): 
        xs, ys, zs = [],[],[]
        for trans in serial: 
            x,y,z = trans[0:3,3]
            xs.append(x)
            ys.append(y)
            zs.append(z)

        return xs,ys,zs

    def extract_vectors_from_trans(self, trans):
        x,y,z = trans[0:3,3]
        p = [x,y,z]
        v1 = trans[0:3,0]
        v2 = trans[0:3,1]
        v3 = trans[0:3,2]

        return p, [v1,v2,v3]

    def plot_arrow(self, ax, p,v,color):
        x,y,z = p 
        u,v,w = v 
        ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True,color=color)

    def plot_frame(self,ax, trans): 
        p,vs = self.extract_vectors_from_trans(trans)


        colors = ['r', 'g', 'b']

        for i in range(3): 
            self.plot_arrow(ax, p , vs[i], colors[i])

    def plot(self):
            fk = self.transformations
            fig = plt.figure(figsize=(14,14))
            ax = fig.gca(projection='3d')

            for k in fk: 
                self.plot_frame(ax,k)

            xs,ys,zs = self.extract_plot_points(fk)
            ax.plot(xs,ys,zs, linewidth=1)

            ax.set_xlim3d([-1,1])
            ax.set_ylim3d([-1,1])
            ax.set_zlim3d([-1,1])

            plt.show()


r = Robot([pi/4, pi/4, pi/4, pi/4, pi/4, pi/4], [0, 0, 0.5], [0, 0, 0.5], [0, 0, 0.1])
r.FK_solve()
r.plot()
         