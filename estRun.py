from email.mime import base
from sys import intern
import numpy as np
import scipy as sp
#NO OTHER IMPORTS ALLOWED (However, you're allowed to import e.g. scipy.linalg)

def dqdx(x, v, dt):
    theta = x[2]
    vel = 5 * v[0] * v[2]
    A = np.array([[1, 0, -vel*np.sin(theta)*dt], [0, 1, vel*np.cos(theta)*dt], [0, 0, 1]], dtype=object)
    return A

def dqdv(x, v, dt):
    theta = x[2]
    gamma = x[3]
    omega = x[4]
    L = np.array([[5*omega*np.cos(theta)*dt, 0, 5*v[0]*np.cos(theta)*dt, 0], [5*omega*np.sin(theta)*dt, 0, 5*v[0]*np.sin(theta)*dt, 0], [5*omega/v[1]*np.tan(gamma)*dt, -1/pow(v[1],2)*5*omega*v[0]*np.tan(gamma)*dt, 5*v[0]/v[1]*np.tan(gamma)*dt, 5*v[0]*omega*dt/(v[1]*pow(np.cos(gamma),2))]], dtype=object)
    return L

def dhdx(x):
    theta = x[2]
    baseline_m = 0.8
    H = np.array([[1, 0, -0.5*baseline_m*np.sin(theta)], [0, 1, 0.5*baseline_m*np.cos(theta)]], dtype=object)
    return H

def dhdw(x):
    return np.identity(2)


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
    gamma = internalStateIn[3]
    omega = internalStateIn[4]
    P_m = internalStateIn[5]
    myColor = internalStateIn[6]

    # NOTE(youngsang): define radius, baseline, gamma, omega uncertainty as process noise
    # since gamma, omega measurement and piecewise constant assumption may not be perfect 
    radius_m = 0.425
    baseline_m = 0.8
    # NOTE(youngsang): assume radius and baseline follow uniform distribution
    # also assume that variance of gamma and omega process noise is 0.5
    radius_var = pow(2 * 0.05 * radius_m, 2) / 12
    baseline_var = pow(2 * 0.1 * baseline_m, 2) / 12
    sig_vv = np.array([[radius_var, 0, 0, 0], [0, baseline_var, 0, 0], [0, 0, 0.1, 0], [0, 0, 0, 0.1]], dtype=object)
    process_noise = np.array([radius_m, baseline_m, omega, gamma], dtype=object)

    # NOTE(youngsang): measurement noise computed from calibration data (run 0)
    sig_ww = np.array([[1.08933973, 1.53329122], [1.53329122, 2.98795486]], dtype=object)
    mean_w = np.array([0.18752364358711568, 0.3962695841848831], dtype=object)

    # prior update
    vel = 5 * radius_m * omega 
    x_p = x + vel * np.cos(theta) * dt 
    y_p = y + vel * np.sin(theta) * dt
    theta_p = theta + vel / baseline_m * np.tan(gamma) * dt
    state_p = np.array([x_p, y_p, theta_p])
    
    A = dqdx(internalStateIn, process_noise, dt)
    L = dqdv(internalStateIn, process_noise, dt)

    P_p = A @ P_m @ A.T + L @ sig_vv @ L.T

    # measurement update
    if not (np.isnan(measurement[0]) or np.isnan(measurement[1])):
        # have a valid measurement
        x = measurement[0]
        y = measurement[1]

        H = dhdx(internalStateIn)
        M = dhdw(internalStateIn)
        inv_mat = np.linalg.inv(np.matrix(H @ P_p @ H.T + M @ sig_ww @ M.T, dtype=np.float))
        K = P_p @ H.T @ inv_mat
        z1 = x - (x_p + 0.5 * baseline_m * np.cos(theta_p) + mean_w[0])
        z2 = y - (y_p + 0.5 * baseline_m * np.sin(theta_p) + mean_w[1])
        # z1 = x - (x_p + 0.5 * baseline_m * np.cos(theta_p))
        # z2 = y - (y_p + 0.5 * baseline_m * np.sin(theta_p))
        z = np.array([[z1], [z2]])
        state_m = state_p.reshape(3,1) + K @ z
        P_m = P_p - K @ H @ P_p
                
    else:
        state_m = state_p
        P_m = P_p

    #we're unreliable about our favourite colour: 
    if myColor == 'green':
        myColor = 'red'
    else:
        myColor = 'green'


    #### OUTPUTS ####
    # Update the internal state (will be passed as an argument to the function
    # at next run), must obviously be compatible with the format of
    # internalStateIn:
    internalStateOut = [state_m[0].item(),
                     state_m[1].item(),
                     state_m[2].item(),
                     steeringAngle,
                     pedalSpeed,
                     P_m, 
                     myColor
                     ]
    
    # NOTE(youngsang): x, y estimated values are position of a rear wheel thus, 
    # to compare with measurement values at last, we add the offset
    x = state_m[0].item() + 0.5 * baseline_m * np.cos(state_m[2].item())
    y = state_m[1].item() + 0.5 * baseline_m * np.sin(state_m[2].item())
    # x = state_m[0].item()
    # y = state_m[1].item()
    theta = state_m[2].item()

    # DO NOT MODIFY THE OUTPUT FORMAT:
    return x, y, theta, internalStateOut 


