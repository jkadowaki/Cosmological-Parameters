import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pylab
import scipy
from scipy import special
from scipy import stats

# put everything together into a single set of functions


def kappa(omega_MPass, omega_lambdaPass):

    answer1 = 1.0 - omega_MPass - omega_lambdaPass

    return answer1


def integrand1(omega_MPass, omega_lambdaPass, zPass):

    if (kappa(omega_MPass, omega_lambdaPass) != 0): # not a flat universe
        answer2 = (kappa(omega_MPass, omega_lambdaPass)*(1+zPass)**2 + omega_MPass*(1+zPass)**3 + omega_lambdaPass)**(-0.5)
    else: # flat universe
        answer2 = (omega_MPass*(1+zPass)**3 + omega_lambdaPass)**(-0.5)
        
    return answer2


def DL_unitless_fcn(zPass,omega_MPass,omega_lambdaPass):
    
    delta_z = 0.0005 # step size in the sum
    zMin = 0.0 # our observation point
    zMax = np.copy(zPass) # this is how far out we are looking
    numSteps = (zMax-zMin)/delta_z + 1
    z_array = np.linspace(zMin, zMax, num=numSteps, endpoint=True)
    
    DL_array_sumTerms = integrand1(omega_MPass, omega_lambdaPass, z_array)

    if (kappa(omega_MPass, omega_lambdaPass) < 0): # hyperbolic universe
        DL_unitless = ((1.0+z_array[np.where(z_array==zPass)[0][0]]))/(np.sqrt(np.abs(kappa(omega_MPass, omega_lambdaPass))))*\
            np.sin(\
            np.sqrt(np.abs(kappa(omega_MPass, omega_lambdaPass)))*\
            np.sum(DL_array_sumTerms[0:np.where(z_array==zPass)[0][0]+1.0])*delta_z\
            )
    elif (kappa(omega_MPass, omega_lambdaPass) > 0): # positive curvature universe
        DL_unitless = ((1.0+z_array[np.where(z_array==zPass)[0][0]]))/(np.sqrt(np.abs(kappa(omega_MPass, omega_lambdaPass))))*\
            np.sinh(\
            np.sqrt(np.abs(kappa(omega_MPass, omega_lambdaPass)))*\
            np.sum(DL_array_sumTerms[0:np.where(z_array==zPass)[0][0]+1.0])*delta_z\
            )
    else: # flat universe
        DL_unitless = (1.0+z_array[np.where(z_array==zPass)[0][0]])*\
            np.sum(DL_array_sumTerms[0:np.where(z_array==zPass)[0][0]+1.0])*delta_z
    return DL_unitless


def DL_units(DL_unitlessPass,H0Pass): # this is just for adding in the units
    
    #H0 = 68 # (km s^-1 Mpc^-1)
    c = 3.0E8 # (m s^-1)
    withUnits = DL_unitlessPass*(c/H0Pass)*(1./10**3) # the last bit is to convert km <-> m
    
    return withUnits # (Mpc)

# sanity-check plot
#abcissa = np.linspace(0,5,num=500)
#ordinate = np.copy(abcissa)
#for i in range(0,500):
#    ordinate[i] = DL_unitless_fcn(abcissa[i],0.32,0.68)
#plt.plot(abcissa,ordinate)
