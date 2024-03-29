import numpy as np
from numpy import pi, sqrt
from numpy.core.function_base import linspace
from sympy import symbols, solve, linsolve
from matplotlib.pyplot import *

def trapezoidal_trajectory(q_param, t_param):
    '''
    Method for getting all angles, velocitis and accelerations for one joint at all points during the trajectory.
    :param t_param: [t0, tau, T, tf]
    :param q_param:[q0, qf, dq_m, ddq_m]
    '''

    t0, tau, T, tf = t_param
    tb = T - tau
    q0, qf, dq_m, ddq_m = q_param
    t = np.linspace(t0,tf, int(3E3))

    q = []
    v = []
    a = []

    for i in t: 

        if i <=tb:
            qi = q0 + (0.5*ddq_m*(i-t0)**2)
            q02 = qi
            vi = ddq_m*i
            v02 = vi 
            ai = ddq_m
            
        elif i > tb and i <= T: 
            vi = dq_m
            qi = q02 + v02*(i-tb)
            ai = 0 
            
        elif i > T: 
            vi = ddq_m*(tf-i)
            qi = qf - (0.5*ddq_m*(i-tf)**2)
            ai = -ddq_m

        q.append(qi)
        v.append(vi) 
        a.append(ai)

    return t, q, v, a


def plot_trajecory(q_params, t_params):

    _, tau, T, tf = t_params
    tb = T - tau
    _, qf, dq_m, _ = q_params
    t, q, v, a = trapezoidal_trajectory(q_params, t_params)
    figure(figsize=(16,16))
    subplot(312)
    plot(t,v, linewidth=2, label="v")
    xlabel('t (s)', fontsize=10)
    ylabel(r'v(t) ($\degree$/s)', fontsize=10)
    grid(color='black', linestyle='--', linewidth=1.0, alpha = 0.7)
    grid(True)
    xlim([0,tf])
    ylim([0,max(v)+1])
    hlines(dq_m, 0,tf, linestyles='--', color='r', label=r"$v_{max}$")
    vlines([tb,T], 0, max(q), linestyles='--', linewidth=2)
    legend()

    subplot(311)
    plot(t,q, linewidth=2, label="q")
    xlabel('t (s)', fontsize=10)
    ylabel(r'q(t) ($\degree$)', fontsize=10)
    grid(color='black', linestyle='--', linewidth=1.0, alpha = 0.7)
    grid(True)
    xlim([0,tf])
    ylim([0,max(q)+10])
    hlines(qf, 0,tf, linestyles='--', color='r', label=r"$q_{final}$")
    vlines([tb,T], 0, max(q)+10, linestyles='--', linewidth=2)
    legend()

    subplot(313)
    plot(t,a, linewidth=3, label="a")
    xlabel('t (s)', fontsize=10)
    ylabel(r'a(t) ($\degree^2$/s)', fontsize=10)
    grid(color='black', linestyle='--', linewidth=1.0, alpha = 0.7)
    grid(True)
    xlim([0,tf])
    ylim([min(a)-1,max(a)+1])
    vlines([tb,T], min(a)-1, max(a)+1, linestyles='--', linewidth=2)
    legend()

    show()

def trajectory_time(q_params, t0 = 0):
    '''
    Perform tracetory planning for a single joint.
    :returns: t_params = [t0, tau, T, tf]
    '''
    q0, qf, dq_m, ddq_m = q_params
    dq = qf - q0

    #check if trianglular trajectory
    limit = sqrt(dq * ddq_m)

    if limit <= dq_m: #tringle
        tb = sqrt(dq/ddq_m)
        tf = 2*tb
        T = tf - tb
        tau = T - tb
    else: #trapez
        tb = dq_m / ddq_m
        T = dq / dq_m
        tf = T + tb
        tau = T - tb

    return [t0, tau, T, tf]


def time_sync(q_params_list):
    '''
    Synchronizes motion through all 6 joints
    '''
    
    t0 = 0
    taus = []
    tbs = []

    for p in q_params_list:
        current = trajectory_time(p, t0)
        tau = current[1]
        taus.append(tau)
        tbs.append(current[2] - tau)

    tau = max(taus)
    tb = max(tbs)
    print(f'Rise time {tb}')
    T = tb + tau
    tf = T + tb
    print(f'Total trajectory time: {tf}s')
    t_params = [t0, tau, T, tf]

    #update velocities and acceleration
    i = 1
    new_q_params_list =[]
    for p in q_params_list:
        dq_new = (p[1] - p[0]) / T
        ddq_new = dq_new / tb
        #print(f'Velocity of joint {i} changed form {p[2]} to {dq_new}')
        #print(f'Acceleration of joint {i} changed form {p[3]} to {ddq_new}')
        p[2:] = [dq_new, ddq_new]
        #show new plots for each joint
        #trapezoidal_trajectory(p, t_params)
        new_q_params_list.append(p)
        i += 1
    return (t_params, new_q_params_list)


def numerical_sync(q_params_list, frequency = 0.01):
    '''
    This method synchronizes motion given a frequency (Default t=1/100).
    '''
    synchronized = time_sync(q_params_list)
    
    new_q_params_list = []
    t = synchronized[0]
    tau = (int((t[1] / frequency)) + 1) * frequency
    tb = (int(((t[2] - t[1]) / frequency)) + 1) * frequency
    tf = 2 * tb + tau
    T = tf - tb
    print(f"Total time increased from {t[3]} to {tf}s")
    new_t_params_list = [t[0], tau, T, tf]

    for p in synchronized[1]:
        dq_new = (p[1] - p[0]) / T
        ddq_new = dq_new / tb
        new_q_params_list.append([p[0], p[1], dq_new, ddq_new])
    return (new_t_params_list, new_q_params_list)

 

def three_p_trajecory(q_params, qb, t):
    '''
    Method to move the robot from A to C through B. For explanation, please view the report. This method plots the solution automatically.
    :param q_params: Specifies q0, qf, dq_max and dd/q_max
    :param qb, Specifies the position at point B in degrees
    :param: t: Specifies t0, tb, and tf, -> the times of each position
    '''
    q0, qf, _, _= q_params
     
    a, b, c, d = symbols('a b c d')
    #calculating A->B
    t0 = t[0]     
    a0 = solve(a-q0 + b*t0 +c*t0**2+ d*t0**3, a)
    a1 = solve(b + 2*c*t0+ 3*d*t0**2, b)
    t0 = t[1]
    for a_ in a0:
        for b_ in a1:
            solution, = linsolve([a_-qb + b_ + c*t0**2+ d*t0**3, b_ + 2*c*t0+ 3*d*t0**2], (c, d))
            a2 = solution[0]
            a3 = solution[1]

    x1 = np.linspace(t[0], t[1],int(3E3))

    q_1 = a0[0] + a1[0]*x1 +a2*x1**2 + a3*x1**3
    dq_1 = a1[0] + 2*a2*x1 + 3*a3*x1**2
    ddq_1 = 2*a2 + 6*a3*x1

    #calculating B -> C
    t0 = 0 
    a0 = solve(a-qb + b*t0 +c*t0**2+ d*t0**3, a)
    a1 = solve(b + 2*c*t0+ 3*d*t0**2, b)
    t0 = t[2] - t[1]
    for a in a0:
        for b in a1:
            solution, = linsolve([a-qf + b + c*t0**2+ d*t0**3, b + 2*c*t0+ 3*d*t0**2], (c, d))
            a2 = solution[0]
            a3 = solution[1]
         
    x2 = linspace(t[1], t[2], int(3E3))
    q_2 = a0[0] + (a1[0]*(x2-t[1])) +a2*((x2-t[1])**2) + a3*((x2-t[1])**3)
    dq_2 = (a1[0]) + 2*(a2)*(x2-t[1]) + 3*(a3)*(x2-t[1])**2
    ddq_2 = 2*(a2) + 6*(a3)*(x2-t[1])

    #plot polynomials
    figure(figsize=(16,16))
    subplot(311)
    plot(x1, q_1, 'r', label = "q(t); t_0 <= t <= t_b")
    plot(x2, q_2,'g', label = "q(t); t_b < t <= t_f")
    grid(color='black', linestyle='--', linewidth=1.0, alpha = 0.7)
    grid(True)
    xlabel('t (s)', fontsize=10)
    ylabel(r'q(t) ($\degree$/s)', fontsize=10)
    legend()

    subplot(312)
    plot(x1, dq_1, 'r', label = "v(t); t_0 <= t <= t_b")
    plot(x2, dq_2,'g', label = "v(t); t_b < t <= t_f")
    grid(color='black', linestyle='--', linewidth=1.0, alpha = 0.7)
    grid(True)
    xlabel('t (s)', fontsize=10)
    ylabel(r'v(t) ($\degree$/s)', fontsize=10)
    legend()

    subplot(313)
    plot(x1, ddq_1, 'r', label = "a(t); t_0 <= t <= t_b")
    plot(x2, ddq_2,'g', label = "a(t); t_b < t <= t_f")
    grid(color='black', linestyle='--', linewidth=1.0, alpha = 0.7)
    grid(True)
    xlabel('t (s)', fontsize=10)
    ylabel(r'a(t) ($\degree^2$/s)', fontsize=10)
    legend()
    show()



#test input values
q_params1 = [0, 90, 8, 4]
q_params2 = [10, 170, 7, 5]
q_params3 = [50, 100, 7, 4]
q_params4 = [0, 270, 6, 4]
q_params5 = [180, 359, 6, 5]
q_params6 = [0, 180, 7, 8]


q_params_list = [q_params1, q_params2, q_params3, q_params4, q_params5, q_params6]

#Examples

#-> show trajectory for first joint
#plot_trajecory(q_params1, trajectory_time(q_params1, 0)) 


#-> synchronize motions and plot
#sync1 = time_sync(q_params_list)
#for sync in sync1[1]:
    #plot_trajecory(sync, sync1[0])

#-> synchrinize motions numerically and plot
#sync2 = numerical_sync(q_params_list)
#for sync in sync2[1]:
    #plot_trajecory(sync, sync2[0])

#-> Move joint from position A to C through B
#three_p_trajecory([0,-10,0,0], 15, [0,2,4])
