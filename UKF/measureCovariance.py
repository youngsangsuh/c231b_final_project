import numpy as np
import matplotlib.pyplot as plt

# N = 100
# measErr = np.zeros((N,2))

experimentalRun = 0

experimentalData = np.genfromtxt ('data/run_{0:03d}.csv'.format(experimentalRun), delimiter=',')

numDataPoints = experimentalData.shape[0]

idx = np.all([~np.isnan(experimentalData[:,3]), ~np.isnan(experimentalData[:,4])],axis=0)

meas_x = experimentalData[idx, 3]
meas_y = experimentalData[idx, 4]

true_x = experimentalData[-1,5]
true_y = experimentalData[-1,6]

N = np.sum(idx)
measErr = np.zeros((N,2))
for i in range(N):
    measErr[i,0] = true_x - meas_x[i]
    measErr[i,1] = true_y - meas_y[i]

# print(measErr[1:50,:])




covMatrix = np.cov(measErr.T)
# meanMatrix = np.mean(np.array([measx, measy]), axis=1)
meanMatrix = np.mean(measErr, axis = 0)

finalPosX = experimentalData[-1,5]
finalPosY = experimentalData[-1,6]

print('Covariance matrix of measurement noise: \n', covMatrix)
print('Mean of measurement noise: \n', meanMatrix)
print('Final position x: ', finalPosX)
print('Final position y: ', finalPosY)

