import numpy as np
import matplotlib.pyplot as plt

#provide the index of the experimental run you would like to use.
# Note that using "0" means that you will load the measurement calibration data.
experimentalRun = 0

print('Loading the data file #', experimentalRun)
experimentalData = np.genfromtxt ('data/run_{0:03d}.csv'.format(experimentalRun), delimiter=',')

numDataPoints = experimentalData.shape[0]

idx = np.all([~np.isnan(experimentalData[:,3]), ~np.isnan(experimentalData[:,4])],axis=0)

measx = experimentalData[idx, 3]
measy = experimentalData[idx, 4]

covMatrix = np.cov(np.array([measx, measy]))
meanMatrix = np.mean(np.array([measx, measy]), axis=1)

finalPosX = experimentalData[-1,5]
finalPosY = experimentalData[-1,6]

print('Covariance matrix of measurement noise: \n', covMatrix)
print('Mean of measurement noise: \n', meanMatrix)
print('Final position x: ', finalPosX)
print('Final position y: ', finalPosY)

