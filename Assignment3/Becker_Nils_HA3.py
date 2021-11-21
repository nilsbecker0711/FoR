import numpy as np
from numpy import pi, sqrt
from matplotlib.pyplot import *

def draw_plot(t0, tb, tf, T, q_param):
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
            qi =  q02 + v02*(i-tb)
            ai = 0 
            
        elif i > T: 
            vi = ddq_m*(tf-i)
            qi = qf - (0.5*ddq_m*(i-tf)**2)
            ai = -ddq_m

        q.append(qi)
        v.append(vi) 
        a.append(ai)

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

def trajectory_time(q_params, t0):
    '''
    Perform tracetory planning for a single joint.
    :returns: t_params = [t0, tau, T, tf], where tau = dq_max / ddq_max and T = dq / dq_max 
    '''
    q0, qf, dq_m, ddq_m = q_params
    dq = qf - q0

    #check if trianglular trajectory
    limit = sqrt(dq * ddq_m)

    if limit <= dq_m: #tringle
        tb = np.sqrt(dq/ddq_m)
        tf = 2*tb
        T = tf - tb
        tau = T - tb
    else: #trapez
        tb = dq_m/ddq_m
        T = dq/dq_m
        tf = T+tb
        tau = T - tf 

    draw_plot(t0, tb, tf, T, q_params)
    return [t0, tau, T, tf]

q_params1 = [0, 90, 8, 4]
q_params2 = [10, 170, 50, 10]
q_params3 = [50, 100, 5, 2]
q_params4 = [0, 270, 6, 6]
q_params5 = [180, 359, 10, 5]
q_params6 = [0, 180, 7, 8]

t0 = 0
q_params_list = [q_params1, q_params2, q_params3, q_params4, q_params5, q_params6]
for p in q_params_list:
    trajectory_time(p, 0)
