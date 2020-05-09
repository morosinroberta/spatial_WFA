# spatial_WFA
Spatially constrained weak-field approximation routines

The weak-field approximation allows for a quick computation of the
magnetic field vector from spectropolarimetric observations assuming
Zeeman induced polarization.

Tradiationally the WFA has been calculated in an unconstrained manner
assuming that the solution from each pixel in the FOV is independent
from the rest. In these routines we implement the ideas presented in
Morosin et al. (2020) to impose spatial coherency in the reconstruction
of the WFA, which improves the usability and fidelity of the WFA.

## Compilation of the C++ module
These routines require building and solving a sparse linear system of
equations. In order to improve performance we have implemented that
part in a C++ module. To compile it simply use:
```
python3 setup.py build_ext --inplace
```

If everything compiles well, you should see a new file called spatWFA.???.so
that should be copied along with wfa_spatial.py to your PYTHONPATH folder or
to the folder where you want to execute these routines.

## Usage
The input data must be a wavelength array (1D) and a data array (4D) with
dimensions (ny,nx,nStokes,nw). If we assume that we have loaded such data
into two arrays named wav and data and that the observed spectral line is
the Ca II 8542 line:

```
import numpy as np
import wfa_spatial as wa

# get dimensions from the data
ny,nx,nStokes,nw = data.shape

# load the line atomic data
lin = wa.line(8542)

# define the noise estimate arrays. In this case we assume
# that the observations are normalized to the average continuum
# and that the noise level is 10**-3 relative to that level.
sig = np.zeros((nStokes,nw), dtype='float64') + 1.e-3

# Get Blos assuming a regularization weight of 1.0

alpha = 1.0 
Blos  = wa.getBlos(wav, data, sig, lin, alpha)

# Get Btrans and Bazi. Assume a Doppler width of 70 mA to define what is wing
# and what is core (not used)

Btr, Bazi = wa.getBhorAzi(wav, data, sig, lin, alpha, vdop=0.070)

