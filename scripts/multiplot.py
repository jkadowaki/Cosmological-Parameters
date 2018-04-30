#! /usr/bin/env python

################################################################################

import os
import re
import subprocess
import glob as g
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

################################################################################

class ImageFollower(object):
    'update image in response to changes in clim or cmap on another image'
    
    def __init__(self, follower):
        self.follower = follower
    
    def __call__(self, leader):
        self.follower.set_cmap(leader.get_cmap())
        self.follower.set_clim(leader.get_clim())

################################################################################

def plot_OM_OL(arr, fname, peak=-1):
    
    plt.clf()

    plt.rcParams['savefig.facecolor'] = "1."
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size='11')
    
    # Probability Density
    arr[np.where(arr==0.)] = np.nan
    im = plt.imshow(arr, interpolation='None', origin='lower', extent=[0,2.5,-1,3])
    plt.colorbar(im)
    if peak > 0: im.set_clim(vmin=0, vmax=peak)
    
    
    # q0 Lines
    rot = 27.5  # Degrees
    plt.plot([0.0, 2.5],[-0.5,0.8], 'w--')
    plt.plot([0.0, 2.5],[0.0, 1.3], 'w--')
    plt.plot([0.0, 2.5],[0.5, 1.8], 'w--')
    plt.text(2.00, 0.80, '$q_0=-0.5$',              rotation=rot, color='w')
    plt.text(2.00, 1.30, '$q_0=0.0$',               rotation=rot, color='w')
    plt.text(2.00, 1.80, '$q_0=0.5$',               rotation=rot, color='w')
    plt.text(1.25, 0.85, '$\mathrm{Decelerating}$', rotation=rot, color='w')
    plt.text(1.25, 1.03, '$\mathrm{Accelerating}$', rotation=rot, color='w')
    
    # Omega_tot Line
    rot = -46.5  # Degrees
    plt.plot([0.0, 2.0],[1.0, -1.0], 'w-')
    plt.text(1.70, -0.70, '$\Omega_\mathrm{tot}=1$',rotation=rot, color='w')
    plt.text(1.20, -0.37, '$\mathrm{Open}$',        rotation=rot, color='w')
    plt.text(1.20, -0.20, '$\mathrm{Closed}$',      rotation=rot, color='w')
    
    # Omega_Lambda=0 Line
    rot = 5    # Degrees
    plt.plot([0.0, 2.5],[0.0, 0.0], 'w:')
    plt.plot([0.0, 1.35, 1.8, 2.5],[0.0, 0.0, 0.02, 0.11], 'w-')
    plt.text(2.10, -0.12, '$\Omega_\mathrm{\Lambda}=0$', color='w')
    plt.text(1.40, -0.08, '$\mathrm{Recollapes}$', rotation=rot, color='w')
    plt.text(1.40,  0.13, '$\mathrm{Expands \, to \, Infinity}$', rotation=rot, color='w')
    
    # No Big Bang
    plt.text(0.10, 2.75, '$\mathrm{No \, Big \, Bang}$', rotation=65)
    
    # Tick Label Size
    plt.xticks(size=15)
    plt.yticks(size=15)
    
    if prob_tag in fname:
        h0 = re.findall('[-+]?\d+[\.]?\d*', fname)[-1]
        print '\t', fname, '\t', np.max(arr[~np.isnan(arr)]), '\t', np.sum(arr[~np.isnan(arr)])
        plt.title('$H_0 = %s$' % h0,  fontsize=20)
    
    plt.xlabel('$\Omega_\mathrm{M}$', fontsize=20)
    plt.ylabel('$\Omega_\Lambda$',    fontsize=20)
    plt.savefig(fname, bbox_inches='tight')


    Nr = 4
    Nc = 5

    fig = figure()
    cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])

    w = 0.4
    h = 0.22
    ax = []
    images = []
    vmax = -1.
    num = 0

    for i in range(Nr):
        for j in range(Nc):
            pos = [0.075 + j*1.1*w, 0.18 + i*1.2*h, w, h]
            a = fig.add_axes(pos)
            if i > 0:
                a.set_xticklabels([])
            
            # Make some fake data with a range that varies
            # somewhat from one plot to the next.
            data = np.genfromtxt(file[num])
            dd = ravel(data)
        
            # Manually find the max of all colors for
            # use in setting the color scale.
            vmax = max(vmax, amax(dd))
            images.append(a.imshow(data, interpolation='None', origin='lower', extent=[0,2.5,-1,3]))
            
            ax.append(a)
            num += 1

    # Set the first image as the master, with all the others
    # observing it for changes in cmap or norm.

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    for i, im in enumerate(images):
        im.set_norm(norm)
        if i > 0:
            images[0].callbacksSM.connect('changed', ImageFollower(im))

    # The colorbar is also based on this master image.
    fig.colorbar(images[0], cax, orientation='horizontal')


################################################################################

def main():

    global prob_tag, plot_ext
    
    dir      = 'fine_h0_resolution/' # All files to be in directory dir

    text_ext = '.txt'  # Any data/text file to be saved with extention text_ext
    plot_ext = '.0.png'  # Any plots to be saved with extension plot_ext
    mov_ext  = '.gif'
    
    prob_tag  = 'prob_'     # All prob. dist. to start with prefix prob_tag

    
    # Assumes file located in dir & named with prefix cs_tag & extension text_ext
    # Assumes information on z-, m- bin #, and h0 value in filename (in order)
    # e.g. 'chisq_z1_m1_h70.txt'
    # Code currently breaks if h0 is a non-integer value.
    prob_list = sorted(g.glob(os.path.join(dir, prob_tag + '*' + text_ext)))
    z, m, h0 = np.transpose(np.array([re.findall('[-+]?\d+[\.]?\d*', os.path.basename(f)) \
               for f in prob_list]))

    # Iterate through every redshift, stellar mass bin
    for z_param in np.unique(z):
        for m_param in np.unique(m):

            print '\nz: ', z_param, ';\t m: ', m_param, '\n'
            index = np.where((z_param==z) * (m_param==m) == True)[0]
            
            # Iterate through every chisq file
            for num, file in enumerate([prob_list[i] for i in index]):
                

                
            # Normalize P(H0 | mu_0)
            h0_path = os.path.join(dir, h0_tag + 'z' + z_param + '_m' + m_param + text_ext)
            h0_plot            = h0_path.replace(text_ext, plot_ext)
            h0_prob, h0_factor = normalize(np.array(h0_prob))
            np.savetxt(h0_path, h0_prob)
            plot_H0(h0[index], h0_prob, h0_plot)

        
            #"""
            # Create 3D normalized cubes
            prob       = np.array(OmOL_prob)
            norm       = np.sum(prob[~np.isnan(prob)])
            peak_value = np.max(prob[~np.isnan(prob)]/norm)
        
            print 'Prob Sum:\t', norm
            print 'Peak Value:\t', peak_value, '\n'

            for num, frame in enumerate(OmOL_prob):
                np.savetxt(prob_list[num], frame/norm)
                plot_OM_OL(frame/norm, plot_list[num], peak_value)
            #"""
        
        
            # Normalize P(Omega_M, Omega_Lambda | mu_0)
            OM_OL_path = os.path.join(dir, om_ol_tag + 'z' + z_param + '_m' +  m_param + text_ext)
            OM_path = os.path.join(dir, om_tag + 'z' + z_param + '_m' +  m_param + text_ext)
            OL_path = os.path.join(dir, ol_tag + 'z' + z_param + '_m' +  m_param + text_ext)
                    
                    
            OM_OL_plot      = OM_OL_path.replace(text_ext, plot_ext)
            OM_OL_prob, fac = normalize(np.array(OM_OL_prob))
            np.savetxt(OM_OL_path, OM_OL_prob)
            plot_OM_OL(OM_OL_prob, OM_OL_plot)
            
            OM_OL_prob[np.isnan(OM_OL_prob)]=0.
            np.savetxt(OM_path, np.sum(OM_OL_prob, axis=0))
            np.savetxt(OL_path, np.sum(OM_OL_prob, axis=1))

            # Movie Iterating through different H0 in each z,m-bin
            imlist = os.path.join(dir, prob_tag + 'z' + z_param + '_m'        \
                     + m_param + '_h*' + plot_ext)
            mov = imlist.replace(plot_ext, mov_ext).replace('_h*','')
            run('convert -delay 05x100 %s %s' % (imlist, mov))


################################################################################

if __name__ == "__main__":
	main()

################################################################################
