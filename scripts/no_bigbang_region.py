import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pylab
import scipy
from scipy import special
from scipy import stats

# this function returns zero if there is no soln, and 1 if there is a soln
def test_no_big_bang(omega_mPass, omega_lambdaPass):
    zRangePass = np.linspace(1,10,num=10) # range of redshifts to test over
    check = 0 # initialize
    for el in range(0,np.size(zRangePass)):
        limit = omega_mPass*(zRangePass[el]*(1+zRangePass[el])**2)/((1+zRangePass[el])**2 - 1) + \
            (1+zRangePass[el])**2/((1+zRangePass[el])**2 - 1)
        if (omega_lambdaPass > limit):
            check += 1

    if (check != 0): # if a nonphysical soln was found
        result = 0
    else:
        result = 1

    return result
