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
part in a C++ module. This module makes use of the Eigen 3 library,
which should be in your path.

To compile it simply use:
```
python3 setup.py build_ext --inplace
```

If everything compiles well, you should see a new file called spatWFA.???.so
that should be copied along with wfa_spatial.py to your PYTHONPATH folder or
to the folder where you want to execute these routines.

### Using Anaconda python
You can use Anaconda python as a package manager to install all dependencies that are required.
This is particularly convenient in OSX, where the Eigen3 and python are not installed by default.

To do so, we can create a separate environment to install all packages, in this case the environment is called bla but feel free to replace that with a different name:
```
conda create --name bla
conda activate bla
conda install clangxx_osx-64 eigen ipython matplotlib numpy cython scipy astropy llvm-openmp

```
After that, you will be able to compile the binary as explained above. Just remember to load this environment every time you want to use this module.


If you want to use anaconda python in a similar way in Linux, you can follow a very similar approach,
but replacing the compiler packages for gcc:
```
conda create --name bla
conda activate bla
conda install gxx_linux-64 eigen ipython matplotlib numpy cython scipy astropy
```

## Usage
The input data must be a wavelength array (1D) and a data array (4D) with
dimensions (ny,nx,nStokes,nw). If we assume that we have loaded such data
into two arrays named wav and data and that the observed spectral line is
the Ca II 8542 line:

```python
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

```

We have included an example with a real SST/CRISP Mg I 5173 dataset,
which can be found in the example subfolder.

## Citing this method
These routines were developed as part of a publication by
[R. Morosin et al. (2020)](https://arxiv.org/abs/2006.14487).
Please if you use them in your project we would appreciate
it if you cite our publication.

## Acknowledgements
This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (SUNMAG, grant agreement 759548).
