import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pylab
import scipy
from scipy import special
from scipy import stats

from luminosity_distance_v3 import DL_unitless_fcn, DL_units

def mu_simple(DLPass): # input is in Mpc
    
    muPass1 = 5*np.log10(DLPass) + 25
    
    return muPass1


def mu_fcn(zPass,H0Pass,omega_MPass,omega_lambdaPass): # input is in Mpc
    
    muPass2 = mu_simple(DL_units(DL_unitless_fcn(zPass,omega_MPass,omega_lambdaPass),H0Pass))
    
    return muPass2


def chiSq_fcn(zEmpiricalPass,muEmpiricalPass,H0Pass,omega_MPass,omega_lambdaPass,sigma_mu,sigma_v):
    
    numPts = np.size(zEmpiricalPass) # this is the number of true data points we are using
    answer = 0 # initialize
    for elemNum in range(0,numPts): # if I don't use a for-loop, Python thinks the arrays are lists
        bitToAdd = np.sum(((mu_fcn(zEmpiricalPass[elemNum],H0Pass,omega_MPass,omega_lambdaPass)-\
                      muEmpiricalPass[elemNum]))**2)/((sigma_mu[elemNum]**2)+(sigma_v**2))
        answer += np.copy(bitToAdd)
    
    return answer

'''
# make the arrays to show parameter space

omega_m_array = np.linspace(0.1,2.5,num=25)
omega_lambda_array = np.linspace(0.1,3,num=30)

parameter_space_array = [[0 for x in range(25)] for y in range(30)] 

for colNum in range(0,25): # omega_M
    for rowNum in range(0,30): # omega_lambda
        parameter_space_array[rowNum][colNum] = \
        chiSq_fcn(zVals_empirical,muVals_empirical,H0,omega_m_array[colNum],omega_lambda_array[rowNum],sigma_mu,sigma_v)

plt.contour(omega_m_array,omega_lambda_array,parameter_space_array)
plt.title('Chi-squared', fontsize=20)
plt.xlabel('$\Omega_{M}$', fontsize=20)
plt.ylabel('$\Omega_{\Lambda}$', fontsize=20)
'''
