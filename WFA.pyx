"""
CYTHON interface for C++ WFA tools.
Author: R. Morosin, J. de la Cruz Rodriguez, G. Vissers & R. Yadav (ISP-SU, 2020)
https://arxiv.org/abs/2006.14487
"""
cimport numpy as np
from numpy cimport ndarray as ar
from numpy import zeros, abs, sqrt, arctan2, where, pi, float32, float64
from libcpp cimport bool


__author__="R. Morosin, J. de la Cruz Rodriguez, G. Vissers & Yadav (ISP-SU 2020)"
__status__="Developing"


# Expose solver only for float32 and float64 types

cdef extern from "spatial_wfa.hpp":
    cdef void set_spatial_constraints_double "wfa::set_spatial_constraints<double>"(int   ny, int   nx, double  alpha, double beta, const double* lhs, const double* rhs, double* result, int nthreads)
    cdef void set_spatial_constraints_float   "wfa::set_spatial_constraints<float>"(int  ny, int   nx,  float  alpha, float beta, const  float* lhs, const  float* rhs,  float* result, int nthreads)

    cdef void compute_derivatives_float    "wfa::compute_derivatives<float>"(int n,  const float* x,  const float* y,  float* yp)
    cdef void compute_derivatives_double  "wfa::compute_derivatives<double>"(int n, const double* x, const double* y, double* yp)
    cdef void compute_derivatives_many_float "wfa::compute_derivatives_many<float>"(long npix, int n,  const float* x,
				                                                    const float* y, float*  yp, int nthreads)
    cdef void compute_derivatives_many_double "wfa::compute_derivatives_many<double>"(long npix, int n,  const double* x,
				                                                      const double* y, double*  yp, int nthreads)

#
# Constructs the sparse matrix and solves the system (with float64 input)
#
def spatial_constraints_double(int ny, int nx, double alpha, double beta, ar[double,ndim=1] lhs, ar[double,ndim=1] rhs, int nthreads=1):
    cdef ar[double,ndim=2] result = zeros((ny,nx), dtype='float64', order='c')

    set_spatial_constraints_double(<int>ny, <int>nx, <double>alpha, <double>beta, <double*>lhs.data, <double*>rhs.data, <double*>result.data, <int>nthreads)
    
    return result
    

#
# Constructs the sparse matrix and solves the system (with float64 input)
#
def spatial_constraints_float(int ny, int nx, float alpha, float beta, ar[float,ndim=1] lhs, ar[float,ndim=1] rhs, int nthreads=1):
    cdef ar[float,ndim=2] result = zeros((ny,nx), dtype='float32', order='c')

    set_spatial_constraints_float(<int>ny, <int>nx, <float>alpha, <float>beta, <float*>lhs.data, <float*>rhs.data, <float*>result.data, <int>nthreads)
    
    return result



cdef calculate_derivatives_float(ar[float,ndim=1] w, ar[float,ndim=3] d, int nthreads=4):
    cdef int ny = d.shape[0]
    cdef int nx = d.shape[1]
    cdef int nw = d.shape[2]

    cdef ar[float,ndim=3] dp = zeros((ny, nx, nw), dtype='float32')

    cdef long npix = nx*ny
    compute_derivatives_many_float(npix, nw, <float*>w.data, <float*>d.data, <float*>dp.data, nthreads)
    
    return dp



cdef calculate_derivatives_double(ar[double,ndim=1] w, ar[double,ndim=3] d, int nthreads = 4):

    cdef int ny = d.shape[0]
    cdef int nx = d.shape[1]
    cdef int nw = d.shape[2]


    cdef ar[double,ndim=3] dp = zeros((ny, nx, nw), dtype='float64')

    cdef long npix = nx*ny
    compute_derivatives_many_double(npix, nw, <double*>w.data, <double*>d.data, <double*>dp.data, nthreads)
    
    return dp




#
# Fast derivatives
# 
def calculate_derivatives(ar w_in, ar d, int nthreads=1):
    """
    Calculates the derivatives of Stokes I using the high order derivatives proposed by
    Steffen (1990). 
    Input:
         w: 1D array of dimensions (nw) containing the wavelength array
         d: 4D array of dimensions (ny,nx,nStokes,nw) containing the data. Only nStokes=0 will be used
    Returns: 
        dp: 3D array of dimensions (ny,nx,nw) constaining the derivatives of Stokes I
    
    Notes: the data must be continuous in memory. If you have used the "transpose" method without
           also using "ascontiguousarray", it will return garbage because the data won't be 
           contiguous in memory.
    """
    cdef ar[float,ndim=1] wf
    cdef ar[double,ndim=1] wd

    cdef int nw = d.shape[2]

    
    if(d.dtype == 'float32'):
        wf = zeros(nw, dtype='float32')
        wf[:] = w_in[:]

        return calculate_derivatives_float(wf, d, nthreads)

        
    else:
        wd = zeros(nw, dtype='float64')
        wd[:] = w_in[:]
        
        return calculate_derivatives_double(wd, d, nthreads)
