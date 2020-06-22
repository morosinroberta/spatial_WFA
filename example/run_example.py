import wfa_spatial as wfa
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

if __name__ == "__main__":

    #
    # Load data
    #
    d = np.ascontiguousarray(fits.open('CRISP_5173_plage_dat.fits', 'readonly')[0].data, dtype='float64')
    w = np.ascontiguousarray(fits.open('CRISP_5173_plage_wav.fits', 'readonly')[0].data, dtype='float64')
    ny, nx, nStokes, nWav = d.shape
    
        
    #
    # Define weights array, let's use only the three inner points of the array
    #
    intensity_level = np.median(d[:,:,0,0])
    sig = np.zeros((4,w.size), dtype='float64', order='c') 
    sig[:,:] = 3.0e-3 * intensity_level


    #
    # Init line object
    #
    lin = wfa.line(5173)

    
    #
    # Compute Blos, Bhor, Bazi
    #
    alpha_Blos = 10.0
    Blos = wfa.getBlos(w, d, sig, lin, alpha_Blos, mask = (5,6,7))

    alpha_Bhor = 0.1
    Bhor, Bazi = wfa.getBhorAzi(w, d, sig, lin, alpha_Bhor, mask = (5,6,7), vdop = 0.035)

    
    #
    # Plots
    #
    plt.close("all"); plt.ion()
    extent = np.float64((0,nx,0,ny))*2*0.059
    
    f, ax = plt.subplots(nrows=1, ncols=3, figsize=(9,4.5))
    
    im0 = ax[0].imshow(Blos, vmax=800,   vmin=-800, cmap='RdGy',      interpolation='nearest', extent=extent)
    im1 = ax[1].imshow(Bhor, vmax=800,   vmin=0,    cmap='gist_gray', interpolation='nearest', extent=extent)
    im2 = ax[2].imshow(Bazi, vmax=np.pi, vmin=0,    cmap='RdGy',      interpolation='nearest', extent=extent)
    
    names = [r'B$_\parallel$', r'B$_\bot$', r'B$_\chi$']
    f.colorbar(im0, ax=ax[0], orientation='horizontal', label=names[0]+' [G]', )
    f.colorbar(im1, ax=ax[1], orientation='horizontal', label=names[1]+' [G]')
    f.colorbar(im2, ax=ax[2], orientation='horizontal', label=names[2]+' [rad]')

    for ii in range(1,3):
        ax[ii].set_yticklabels([])

    for ii in range(3):
        ax[ii].set_xlabel('x [arcsec]')
    ax[0].set_ylabel('y [arcsec]')
    
    f.set_tight_layout(True)
    f.show()
