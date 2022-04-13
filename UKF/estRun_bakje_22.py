from curses.ascii import SP
from email.mime import base
from sys import intern
import numpy as np
import scipy as sp

from scipy.linalg import block_diag
from scipy.linalg import cholesky
#NO OTHER IMPORTS ALLOWED (However, you're allowed to import e.g. scipy.linalg)

def q_pred(s, omega, gamma, dt):
    s_p = s
    # x,y,theta, r,B,omega,gamma, w1,w2
    x = s[0]; y = s[1]; theta = s[2]
    r = s[3]; B = s[4]
    vel = 5 * r * omega
    s_p[0] = x + dt * vel * np.cos(theta)
    s_p[1] = y + dt * vel * np.sin(theta)
    s_p[2] = theta + dt * vel / B * np.tan(gamma)

    return s_p

def h(s):
    x = s[0]; y = s[1]; theta = s[2]
    r = s[3]; B = s[4]; w1 = s[5]; w2 = s[6]
    s_z = np.zeros((2,))
    s_z[0] = x + 1/2 * B * np.cos(theta) + w1
    s_z[1] = y + 1/2 * B * np.sin(theta) + w2

    return s_z



def estRun(time, dt, internalStateIn, steeringAngle, pedalSpeed, measurement):
    # In this function you implement your estimator. The function arguments
    # are:
    #  time: current time in [s] 
    #  dt: current time step [s]
    #  internalStateIn: the estimator internal state, definition up to you. 
    #  steeringAngle: the steering angle of the bike, gamma, [rad] 
    #  pedalSpeed: the rotational speed of the pedal, omega, [rad/s] 
    #  measurement: the position measurement valid at the current time step
    #
    # Note: the measurement is a 2D vector, of x-y position measurement.
    #  The measurement sensor may fail to return data, in which case the
    #  measurement is given as NaN (not a number).
    #
    # The function has four outputs:
    #  x: your current best estimate for the bicycle's x-position
    #  y: your current best estimate for the bicycle's y-position
    #  theta: your current best estimate for the bicycle's rotation theta
    #  internalState: the estimator's internal state, in a format that can be understood by the next call to this function

    # Example code only, you'll want to heavily modify this.
    # this internal state needs to correspond to your init function:
    
    # receive state from previous step
    x = internalStateIn[0]
    y = internalStateIn[1]
    theta = internalStateIn[2]
    gamma = steeringAngle
    omega = pedalSpeed
    P_m = internalStateIn[3]
    myColor = internalStateIn[4]

    # NOTE(youngsang): define radius, baseline, gamma, omega uncertainty as process noise
    # since gamma, omega measurement and piecewise constant assumption may not be perfect 
    radius_m = 0.425
    baseline_m = 0.8
    # NOTE(youngsang): assume radius and baseline follow uniform distribution
    # also assume that variance of gamma and omega process noise is 0.5
    radius_var = pow(2 * 0.05 * radius_m, 2) / 12
    baseline_var = pow(2 * 0.1 * baseline_m, 2) / 12
    sig_vv = np.array([[radius_var, 0], [0, baseline_var]])
    process_noise = np.array([radius_m, baseline_m])

    # NOTE(youngsang): measurement noise computed from calibration data (run 0)
    # sig_ww = np.array([[1.08933973, 1.53329122], [1.53329122, 2.98795486]], dtype=object)
    # mean_w = np.array([-0.01891406, 1.62806509], dtype=object)
    sig_ww = np.array([[1.08933973, 1.53329122], [1.53329122, 2.98795486]])
    mean_w = np.array([0.187523, 0.39626958])

    # Defining the sigma points - NOTE(JH) # of sigma points = 2 * (n_x + n_v + n_w) + 1 // n_ksi = n_x + n_v + n_w = 9
    ksi = np.array([x,y,theta, radius_m, baseline_m, mean_w[0], mean_w[1]]) # NOTE(JH) the last two term should contain mean_w??
    # ksi = np.array([x,y,theta, radius_m, baseline_m, omega, gamma, 0, 0]) 

    n_ksi = ksi.shape[0]
    n_points = 2*n_ksi + 1
    sigma_ksi = np.zeros((n_ksi, n_points)) ## 9,19
    
    sigma_ksi[:,0] = ksi

    # print('vel:', omega * 5 * radius_m)

    P_ksi = block_diag(P_m, sig_vv, sig_ww)

    # L_ksi = cholesky(P_ksi, lower=True)
    S, T = sp.linalg.eig(P_ksi)
    L_ksi = T.T @ np.diag(np.sqrt(np.real(S)*n_points))

    for i in range(n_ksi):
        sigma_ksi[:,2*i+1] = ksi + L_ksi[:,i]
        sigma_ksi[:,2*i+2] = ksi - L_ksi[:,i]

    ### Prior update
    sigma_ksi_p = np.zeros((n_ksi, n_points)) ## 9,19
    for i in range(n_points):
        sigma_ksi_p[:,i] = q_pred(sigma_ksi[:,i], omega, gamma, dt)

    Xp_hat = np.sum(sigma_ksi[0:3,:], axis = 1) / n_points

    P_p = np.zeros((3,3))
    for i in range(n_points):
        X = (sigma_ksi[0:3,i] - Xp_hat).reshape(-1,1)
        P_p += 1/n_points * X @ X.T
    
    ### Posterior update
    if not (np.isnan(measurement[0]) or np.isnan(measurement[1])):
        sigma_z = np.zeros((2,n_points))
        for i in range (n_points):
            sigma_z[:,i] = h(sigma_ksi_p[:,i])

        z_hat = np.sum(sigma_z, axis = 1) / n_points
        
        P_zz = np.zeros((2,2))
        for i in range(n_points):
            Z = (sigma_z[:,i] - z_hat).reshape(-1,1)
            P_zz += 1/n_points * Z @ Z.T
        
        P_xz = np.zeros((3,2))
        for i in range(n_points):
            X = (sigma_ksi[0:3,i] - Xp_hat).reshape(-1,1)
            Z = (sigma_z[:,i] - z_hat).reshape(-1,1)
            P_xz += 1/n_points * X @ Z.T
        K = P_xz @ np.linalg.inv(P_zz)
        Xm = Xp_hat + K @ (measurement - z_hat)
        P_m = P_p - K @ P_zz @ K.T
    else:
        Xm = Xp_hat
        P_m = P_p

    Xm = Xm.reshape(-1)


    #we're unreliable about our favourite colour: 
    if myColor == 'green':
        myColor = 'red'
    else:
        myColor = 'green'


    #### OUTPUTS ####
    # Update the internal state (will be passed as an argument to the function
    # at next run), must obviously be compatible with the format of
    # internalStateIn:
    internalStateOut = [Xm[0],
                        Xm[1],
                        Xm[2],
                        P_m, 
                        myColor
                        ]
    
    # NOTE(youngsang): x, y estimated values are position of a rear wheel thus, 
    # to compare with measurement values at last, we add the offset
    x = Xm[0] + 0.5 * baseline_m * np.cos(Xm[2])
    y = Xm[1] + 0.5 * baseline_m * np.sin(Xm[2])
    # x = state_m[0].item()
    # y = state_m[1].item()
    theta = Xm[2].item()

    # DO NOT MODIFY THE OUTPUT FORMAT:
    return x, y, theta, internalStateOut 


