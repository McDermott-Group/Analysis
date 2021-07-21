from sklearn import mixture
import dataChest
import importlib
importlib.reload(dataChest)
from dataChest import *
import numpy as np
import matplotlib.pyplot as plt
import general.calibration as calkit
importlib.reload(calkit)
from general.calibration import *

# dc = dataChest(['fluxNoise2','DR1 - 2019-12-17','CorrFar','Q1Q2Q3Q4Corr','General','04-15-20','Charge_resetting','HDF5Data'])
# # dc.openDataset('cwh0705lum_Charge_resetting')
# dc.openDataset('cwh0704ccm_Charge_resetting')
# data = dc.getData(variablesList=['Trial', 'Is SB2', 'Qs SB2'])
# Is = data[:,1].reshape((10,10000))
# Qs = data[:,2].reshape((10,10000))

# o2 = []
# o3_01 = []
# o3_012 = []

# # for i in range(10):
# for i in [0,8]:
    # IQs = np.dstack((Is[i],Qs[i]))[0]
    # scale_factor = np.abs(IQs.max())
    # gmm3 = mixture.GaussianMixture(n_components=3,
                                  # covariance_type='spherical',
                                  # max_iter = 300, tol = 5e-7).fit(IQs/scale_factor)
    # gmm2 = mixture.GaussianMixture(n_components=2,
                                  # covariance_type='spherical',
                                  # max_iter = 300, tol = 5e-7).fit(IQs/scale_factor)
    # bic = np.abs(gmm2.bic(IQs/scale_factor) - gmm3.bic(IQs/scale_factor))/np.abs(gmm2.bic(IQs/scale_factor) + gmm3.bic(IQs/scale_factor))
    # # print i, bic
    # if bic > 0.02:
        # gmm = gmm3
    # else:
        # gmm = gmm2
    # gaussian_means = gmm.means_ * scale_factor
    # gaussian_std = np.sqrt(gmm.covariances_) * scale_factor
    
    # labels2 = gmm2.predict(IQs/scale_factor)
    # labels3 = gmm.predict(IQs/scale_factor)
    # _o2 = np.unique(labels2, return_counts=True)[1]
    # _o3 = np.unique(labels3, return_counts=True)[1]
    # o2 += [1.*_o2[1]/(_o2[0]+_o2[1])]
    # o3_01 += [1.*_o3[1]/(_o3[0]+_o3[1])]
    # o3_012 += [1.*_o3[1]/labels3.size]
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.scatter( Is[i], Qs[i], s=0.5, c=labels3)
    # ax.set_title(i)
    # for state in range(len(gaussian_means)):
        # I0,Q0 = np.mean(Is[i]), np.mean(Qs[i])
        # ax.plot(I0, Q0, 'X')
        # I1,Q1 = gaussian_means[state]
        # ax.plot([I0, I1], [Q0, Q1])
        # ax.add_artist(plt.Circle((I1,Q1), gaussian_std[state], fill=False, zorder=10))
    # plt.draw()
    # plt.pause(0.05)
    
# # fig = plt.figure()
# # ax = fig.add_subplot(111)
# # ax.plot(o2)
# # ax.plot(o3_01)
# # ax.plot(o3_012)
# # plt.draw()
# # plt.pause(0.05)


# # cal = CalibratedStates(plot=True)
# # for i in range(10):
# # for i in [0,8,9]:
# for i in [0,1,5]:
    # cal = CalibratedStates()
    # cal.set_id_fn(cal._id_centers_by_angle)
    # states, SSO, fits = cal.get_single_shot_occupation((Is[i],Qs[i]))
    # plt.title(i)
    # o2 += [SSO[1]]
    # # print cal.gmm.weights_
    # # print SSO
    # plt.pause(0.2)
    
# # fig = plt.figure()
# # ax = fig.add_subplot(111)
# # ax.plot(o2)
# # # ax.plot(o3_01)
# # # ax.plot(o3_012)
# # plt.draw()
# # plt.pause(0.05)


dc = dataChest(['fluxNoise2','DR1 - 2019-12-17','CorrFar','Q2','General','04-22-20','Calibration_split','HDF5Data'])
for file in ['cwo1301fnu_Calibration_split','cwo1317tgv_Calibration_split']:
    dc.openDataset(file)
    data = dc.getData(variablesList=['Is', 'Qs'])
    Is = data[:,0]
    Qs = data[:,1]
    cal = CalibratedStates()
    cal.set_id_fn(cal._id_centers_by_angle)
    states, SSO, fits = cal.get_single_shot_occupation((Is,Qs))
