"""
Spatially constrained Weak Field Approximation module
Uses Tikhonov regularization to set spatial constraints on the WFA

Reference: Morosin, de la Cruz Rodriguez, Vissers & Yadav (2020)
           https://arxiv.org/abs/2006.14487

Dependences:
- SpatWFA is an external cython/C++ module that contains 
  optimized routines for building and solving the sparse linear system
  as well as for computing the Steffen (1990) high order derivatives


"""
import numpy as np
import spatWFA as spa

# *********************************************************************************************** #

def cder(x, y, Steffen = True):
    """
    function cder computes the derivatives of Stokes I (y)
    
    Input: 
            x: 1D wavelength array
            y: 4D data array (ny, nx, nStokes, nw)
      Steffen: (optional, True) if set, use high order derivatives from Steffen (1990)
               which are implemented in the C++ module. Otherwise use the usual
               centered derivatives formula for non-equidistant grids.
    """
    if(not Steffen):
        ny, nx, nstokes, nlam = y.shape[:]
        yp = np.zeros((ny, nx, nlam), dtype='float32')
        
        odx = x[1]-x[0]; ody = (y[:,:,0,1] - y[:,:,0,0]) / odx
        yp[:,:,0] = ody
        
        for ii in range(1,nlam-1):
            dx = x[ii+1] - x[ii]
            dy = (y[:,:,0,ii+1] - y[:,:,0,ii]) / dx
            
            yp[:,:,ii] = (odx * dy + dx * ody) / (dx + odx)
            
            odx = dx; ody = dy
        
        yp[:,:,-1] = ody    
        return yp
    else:
        return spa.calculate_derivatives(x,y)


# *********************************************************************************************** #

   
class line:
    """
    Class line is used to store the atomic data of spectral lines. We use this
    class as input for the WFA routines below.
    Usage: 
        lin = line(8542)
    """
    def __init__(self, cw=8542):

        self.larm = 4.668645048281451e-13
        
        if(cw == 8542):
            self.j1 = 2.5; self.j2 = 1.5; self.g1 = 1.2; self.g2 = 1.33; self.cw = 8542.091
        elif(cw == 6301):
            self.j1 = 2.0; self.j2 = 2.0; self.g1 = 1.84; self.g2 = 1.50; self.cw = 6301.4995
        elif(cw == 6302):
            self.j1 = 1.0; self.j2 = 0.0; self.g1 = 2.49; self.g2 = 0.0; self.cw = 6302.4931
        elif(cw == 8468):
            self.j1 = 1.0; self.j2 = 1.0; self.g1 = 2.50; self.g2 = 2.49; self.cw = 8468.4059
        elif(cw == 6173):
            self.j1 = 1.0; self.j2 = 0.0; self.g1 = 2.50; self.g2 = 0.0; self.cw = 6173.3340
        elif(cw == 5173):
            self.j1 = 1.0; self.j2 = 1.0; self.g1 = 1.50; self.g2 = 2.0; self.cw = 5172.6843
        elif(cw == 5896):
            self.j1 = 0.5; self.j2 = 0.5; self.g1 = 2.00; self.g2 = 2.0/3.0; self.cw = 5895.9242
        else:
            print("line::init: ERROR, line not implemented")
            self.j1 = 0.0; self.j2 = 0.0; self.g1 = 0.0; self.g2 = 0.0; self.cw = 0.0
            return

        j1 = self.j1; j2 = self.j2; g1 = self.g1; g2 = self.g2
        
        d = j1 * (j1 + 1.0) - j2 * (j2 + 1.0);
        self.geff = 0.5 * (g1 + g2) + 0.25 * (g1 - g2) * d;
        ss = j1 * (j1 + 1.0) + j2 * (j2 + 1.0);
        dd = j1 * (j1 + 1.0) - j2 * (j2 + 1.0);
        gd = g1 - g2;
        self.Gg = (self.geff * self.geff) - (0.0125  * gd * gd * (16.0 * ss - 7.0 * dd * dd - 4.0));

        print("line::init: cw={0}, geff={1}, Gg={2}".format(self.cw, self.geff, self.Gg))

      
# *********************************************************************************************** #
     
def getBlos(w, d, sig, line, alpha, beta = 0.0, mask = None, Bnorm=100.0, nthreads = 1, w0=0, w1=-1, Steffen = False):
    """
    Function getBlos computes Blons with the WFA using spatial constraints
    Usage: Blos = getBlosSpat(w, d, sig, line, alpha, mask = None, Bnorm=100.0, nthreads = 2, w0=0, w1=-1)
    Input:
             w: 1D array containing the wavelength offsets from line center. Dimensions: nw
             d: 4D data array with the dimensions arranged as (ny,nx,nStokes,nw)
           sig: 2D array with the estimate of the noise for each Stokes parameter (nStokes,nw)
         alpha: regularization weight. If the estimate of the noise is correct, it should be around 1.
          beta: low-norm regularization weight. Only use if the noise is extremely high.
          mask: (optional) if provided, only the indexes in mask will be used to compute the WFA
         Bnorm: (optional) Typical norm for the magnetic field gradients.
      nthreads: (not implemented yet)
            w0: (optional) together with w1 allow selecting a wavelength window so compute the WFA
            w1: (optional) together with w0 allow selecting a wavelength window so compute the WFA 
       Steffen: (optional) use Steffen (1990) harmonic centered derivatives, or standard centered ones.

    Reference: Morosin, de la Cruz Rodriguez, Vissers & Yadav (2020)
    """
    ny, nx, ns, nw = d.shape

    # Init constants
    c = -line.larm * line.cw**2 * line.geff; cc = c*c

    # calculate derivatives of the intensity
    der = cder(w, d, Steffen = Steffen)

    # Init tmp storage
    lhs = np.zeros((ny, nx), dtype='float64', order='c')
    rhs = np.zeros((ny, nx), dtype='float64', order='c')

    
    if(mask is None):
        
        if(w1 == -1): w1 = nw
        if(w0 < 0):  w0  = 0
        if(w1 > nw): w1  = nw
        
        inw = w1-w0
        for ii in range(w0,w1):
            isig2 = inw*sig[3,ii]**2
            lhs += cc*der[:,:,ii]**2 / isig2
            rhs += c *der[:,:,ii]*d[:,:,3,ii] / isig2

    else:
        inw = len(mask)
        for ii in mask:
            isig2 = inw*sig[3,ii]**2
            lhs += cc*der[:,:,ii]**2 / isig2
            rhs += c *der[:,:,ii]*d[:,:,3,ii] / isig2
        
    Blos = spa.spatial_constraints_double(ny, nx, alpha/(4*Bnorm**2), beta/(4*Bnorm**2), lhs.flatten(), rhs.flatten(), int(nthreads))

    return Blos

# *********************************************************************************************** #

def getBhorAzi(w, d, sig, lin, alpha, beta=0.0, vdop = 0.05, mask = None, Bnorm=100.0, w0=0, w1=-1, nthreads=2, Steffen = True):
    """
    Function getBhorAzi computes Btrans and Bazi with the WFA using spatial constraints
    Usage: Bhor, Bazi = getBhorAziSpat(w, d, sig, lin, alpha = 0.0, vdop = 0.05, mask = None, Bnorm=100.0, w0=0, w1=-1, nthreads=2)
    Input:
             w: 1D array containing the wavelength offsets from line center. Dimensions: nw
             d: 4D data array with the dimensions arranged as (ny,nx,nStokes,nw)
           sig: 2D array with the estimate of the noise for each Stokes parameter (nStokes,nw)
         alpha: regularization weight. If the estimate of the noise is correct, it should be around 1.
          beta: low-norm regularization weight. Only use if the noise is extremely high.
          vdop: Rough estimate of the Doppler with of the line. Wavelengths below this number will be ignored
          mask: (optional) if provided, only the indexes in mask will be used to compute the WFA
         Bnorm: (optional) Typical norm for the magnetic field gradients.
      nthreads: (not implemented yet)
            w0: (optional) together with w1 allow selecting a wavelength window so compute the WFA
            w1: (optional) together with w0 allow selecting a wavelength window so compute the WFA 
       Steffen: (optional) use Steffen (1990) harmonic centered derivatives, or standard centered ones.

    Reference: Morosin, de la Cruz Rodriguez, Vissers & Yadav (2020)
    """
    ny, nx, ns, nw = d.shape

    # Prepare constants
    c = 0.75 * (lin.larm * lin.cw**2)**2 * lin.Gg; cc=c*c

    # Calculate derivatives
    der = cder(w,d, Steffen = Steffen)

    # Adjust derivatives
    for ii in range(len(w)):
        if(np.abs(w[ii]) >= vdop): scl = 1./w[ii]
        else: scl = 0.0
        der[:,:,ii] *= scl

    # Init tmp storage and results
    lhsQ = np.zeros((ny, nx), dtype='float64', order='c')
    rhsQ = np.zeros((ny, nx), dtype='float64', order='c')
    lhsU = np.zeros((ny, nx), dtype='float64', order='c')
    rhsU = np.zeros((ny, nx), dtype='float64', order='c')

    # prepare LHS and RHS for Stokes Q
    if(mask is None):

        if(w1 == -1): w1 = nw
        if(w0 < 0): w1 = 0
        if(w1 > nw): w1 = nw
        
        inw = w1-w0
        
        for ii in range(w0,w1):
            isig2Q = inw*sig[1,ii]**2
            isig2U = inw*sig[2,ii]**2

            ider2 = der[:,:,ii]**2
            
            lhsQ += cc*ider2 / isig2Q
            lhsU += cc*ider2 / isig2U

            rhsQ += c*d[:,:,1,ii]*der[:,:,ii] / isig2Q
            rhsU += c*d[:,:,2,ii]*der[:,:,ii] / isig2U
    else:
        inw = len(mask)
        for ii in mask:
            isig2Q = inw*sig[1,ii]**2
            isig2U = inw*sig[2,ii]**2

            ider2 = der[:,:,ii]**2
            
            lhsQ += cc*ider2 / isig2Q
            lhsU += cc*ider2 / isig2U

            rhsQ += c*d[:,:,1,ii]*der[:,:,ii] / isig2Q
            rhsU += c*d[:,:,2,ii]*der[:,:,ii] / isig2U

    # Construct and solve linear system for each component of Bhor
    BhorQ_2 = spa.spatial_constraints_double(ny, nx, alpha/(4*Bnorm**4), beta/(4*Bnorm**4), lhsQ.flatten(), rhsQ.flatten(), int(nthreads))
    BhorU_2 = spa.spatial_constraints_double(ny, nx, alpha/(4*Bnorm**4), beta/(4*Bnorm**4), lhsU.flatten(), rhsU.flatten(), int(nthreads))

    # calculate the azimuth
    azimuth = 0.5 * np.arctan2(BhorU_2, BhorQ_2)
    azimuth[np.where(azimuth < 0)] += np.pi
    
    # combine contributions to calculate Btrans
    Bhor = np.sqrt(np.sqrt(BhorQ_2**2 + BhorU_2**2))

    return Bhor, azimuth
