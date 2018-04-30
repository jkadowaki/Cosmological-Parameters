import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pylab
import scipy
from scipy import special
from scipy import stats

# function that returns a 3D cube of chi-squared in H0/omega_M/omega_lambda cosmological parameter space, based on
# empirical data from Betoule+ 2014

# inputs:
# zBinPass- redshift bin, as defined in data.py (use integer 0 or 1)
# MBinPass- host galaxy mass bin, as defined in data.py (use integers 0, 1 or 2)
# omega_M_startPass- lower boundary of omega_M in parameter space (use floats)
# omega_M_stopPass- upper boundary of omega_M in parameter space (use floats)
# omega_M_divisions- number of steps in omega_M space between the upper and lower boundaries of omega_M (use integer)
# omega_lambda_startPass- (analogous to the above, for lambda)
# omega_lambda_stopPass- (analogous to the above)
# omega_lambda_divisionsPass- (analogous to the above)
# omega_H0_startPass- (analogous to the above, for the Hubble constant)
# omega_H0_stopPass- (analogous to the above)
# omega_H0_divisionsPass- (analogous to the above)

# output:
# 3D cube of dimensions [H0 steps][omega_lambda steps][omega_M steps]

def generate_chi_sq_space(zBinPass,MBinPass,\
                          omega_M_startPass,omega_M_stopPass,omega_M_divisionsPass,\
                          omega_lambda_startPass,omega_lambda_stopPass,omega_lambda_divisionsPass,\
                          omega_H0_startPass,omega_H0_stopPass,omega_H0_divisionsPass):
          
    ## get the empirical data
    from data import getdata  # for generating array of distance moduli, array of errors, and array of redshifts

    ## get things to calculate distance modulus based on cosmological parameters
    from luminosity_distance_v3 import kappa, integrand1, DL_unitless_fcn, DL_units
    from mu_and_chiSqd_v3 import mu_simple, mu_fcn, chiSq_fcn
    from no_bigbang_region_v3 import test_no_big_bang

    # set which data to use
    mu_empirical = getdata(zBinPass,MBinPass)[0]
    mu_empirical_err = getdata(zBinPass,MBinPass)[1]
    z_empirical = getdata(zBinPass,MBinPass)[2]
    sigma_v = 0.1

    # set the part of parameter space to use, and how fine of a grid to use
    omega_m_array = np.linspace(omega_M_startPass,omega_M_stopPass,num=omega_M_divisionsPass)
    omega_lambda_array = np.linspace(omega_lambda_startPass,omega_lambda_stopPass,num=omega_lambda_divisionsPass) # 2.2,2.4; 0.1,3.0
    H0_array = np.linspace(omega_H0_startPass,omega_H0_stopPass,num=omega_H0_divisionsPass)
    parameter_space_array = [[[0 for x in range(np.size(omega_m_array))] for y in range(np.size(omega_lambda_array))] \
                             for z in range(np.size(H0_array))]
    ################################################################################################
    ################################################################################################

    # define function for making H0/omega_lambda/omega_m chi-squared space
    def make_chi_sqd_space(zVals_empiricalPass, muVals_empiricalPass, H0Pass, omega_m_arrayPass, \
                           omega_lambda_arrayPass, parameter_space_arrayPass, sigma_muPass, sigma_vPass):
                       
        for depthNum in range(0,np.size(H0_array)): # range in H0
            for colNum in range(0,np.size(omega_m_arrayPass)): # range in omega_M
                for rowNum in range(0,np.size(omega_lambda_arrayPass)): # range in omega_lambda
                    chiSq_val = chiSq_fcn(zVals_empiricalPass,muVals_empiricalPass,H0Pass[depthNum],omega_m_arrayPass[colNum],\
                        omega_lambda_arrayPass[rowNum],sigma_muPass,sigma_vPass)
                    parameter_space_arrayPass[depthNum][rowNum][colNum] = chiSq_val.astype(float)

                    if (test_no_big_bang(omega_m_arrayPass[colNum], omega_lambda_arrayPass[rowNum]) != 1): # remove nonphysical solns (approx)
                        parameter_space_arrayPass[depthNum][rowNum][colNum] = 0
            
        return parameter_space_arrayPass
            
    parameter_space_array = make_chi_sqd_space(z_empirical, mu_empirical, H0_array, omega_m_array, \
        omega_lambda_array, parameter_space_array, mu_empirical_err, sigma_v)

    print('Returning cube for ')
    print('omega_m: '+str(omega_lambda_array))
    print('omega_lambda: '+str(omega_m_array))
    print('H0: '+str(H0_array))
    return parameter_space_array

## make plots
                     
#import matplotlib.cm as cm

## make 'flat universe' line
#abcissa_flatUniv = np.linspace(0.0,1.0,num=10)
#ordinate_flatUniv = 1 - abcissa_flatUniv

# color plot
#plt.imshow(np.log10(parameter_space_array), extent=(omega_m_array.min(), omega_m_array.max(), omega_lambda_array.min(), omega_lambda_array.max()),
#           interpolation='nearest', cmap=cm.gist_rainbow, origin="lower")
#plt.title('Log10 Chi-squared of getdata(0,1)', fontsize=20)
#plt.xlabel('$\Omega_{M}$', fontsize=20)
#plt.ylabel('$\Omega_{\Lambda}$', fontsize=20)
#plt.plot(abcissa_flatUniv, ordinate_flatUniv, '-', color='k') # 'flat universe' line
#plt.colorbar()
#plt.show()

# contour plot                
#plt.contour(omega_m_array,omega_lambda_array,parameter_space_array)
#plt.title('Chi-squared', fontsize=20)
#plt.xlabel('$\Omega_{M}$', fontsize=20)
#plt.ylabel('$\Omega_{\Lambda}$', fontsize=20)
#plt.colorbar()
#plt.show()

