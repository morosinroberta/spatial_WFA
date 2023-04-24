#ifndef WFASPAT
#define WFASPAT

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseLU> 
#include <eigen3/Eigen/SparseCholesky>

#include <chrono>
#include <iostream>
#include <omp.h>
#include <iostream>

//
// C++ routines to solve the spatially coupled WFA problem
// Since the routines are the same for all components of B
// We only need to implement them once
//
// Coded by J. de la Cruz Rodriguez (ISP-SU 2020)
//
// Reference: Morosin, de la Cruz Rodriguez, Vissers & Yadav (2020)
//            https://arxiv.org/abs/2006.14487
//
//
// Modifications:
//
//   2021-05-17, JdlCR: modified penalty function. Only take pixels (x-1) and (y-1),
//                      otherwise we are unnecessarily counting each interval twice.
//                      This should lead to a simpler linear system, faster to solve.
//
//   2021-10-01, JdlCR: reverted to original formula. Counting only for (x-1) and (y-1)
//                      can lead to inbalance of the regularization when alpha is very large
//                      because of pixel (0,0) does not have anything in the diagonal.
//                      Need to think this one. Maybe a low-norm term alone will do.
//  
//

namespace wfa{
  
  // ********************************************************************************* //

  template<class T> inline T signFortran(T const &val)
  {
    return ((val < static_cast<T>(0)) ? static_cast<T>(-1) : static_cast<T>(1));
  }
  
  // ********************************************************************************* //

  template<class T> inline T harmonic_derivative_Steffen_one(T const &xu, T const &xc, T const &xd, T const &yu, T const &y0, T const &yd)
  {
    // ---
    // High order harmonic derivatives
    // Ref: Steffen (1990), A&A..239..443S
    //
    // Arguments:
    //   Assuming three consecutive points of a function yu, y0, yd and the intervals between them odx, and dx:
    //       yu: Upwind point
    //       y0: Central point
    //       yd: downwind point
    //
    // ---

    T const odx = (xc-xu);
    T const dx  = (xd-xc);
    
    T const S0 = (yd - y0) / dx;
    T const Su = (y0 - yu) / odx;
    T const P0 = std::abs((Su*dx + S0*odx) / (odx+dx)) / 2;
    
    return (signFortran(S0) + signFortran(Su)) * std::min<T>(std::abs(Su),std::min<T>(std::abs(S0), P0));
  }
  
  
  // ******************************************************************************************* //

  template<typename T> inline
  void compute_derivatives(int const n, const T* const __restrict__ x, const T* const __restrict__ y, T* const __restrict__ yp)
  {

    yp[0]   = (y[1]-y[0])     / (x[1]-x[0]); 
    yp[n-1] = (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]); 

    int const n1 = n-1;
    
    for(int ii=0; ii<n1; ++ii)
      yp[ii] = harmonic_derivative_Steffen_one(x[ii-1], x[ii], x[ii+1], y[ii-1], y[ii], y[ii+1]);
    
  }
  
  // ******************************************************************************************* //

  template<typename T>
  void compute_derivatives_many(long const npix, int const n,  const T* const __restrict__ x,
				const T* const __restrict__ y, T* __restrict__ yp, int const nthreads)
  {

#pragma omp parallel default(shared) num_threads(nthreads)  
    {
#pragma omp for schedule(static)
      for(long ipix=0; ipix<npix; ++ipix)
	compute_derivatives(n, x, &y[ipix*n], &yp[ipix*n]);
    } // parallel block
  }
  
  // ******************************************************************************************* //

  template<typename T>
  Eigen::SparseMatrix<T,Eigen::RowMajor,int> make_sparse_system(int const ny, int const nx,
								const T* __restrict__ lhs,  T const alpha,
								T const beta, int nthreads)
  {
    
    using namespace Eigen;
    
    int const npix   = ny*nx;
    T   const malpha = - alpha;
    
    // --- Init sparse matrix to store the LHS --- //

    SparseMatrix<T,RowMajor,int> A(npix, npix);


    
    // --- for each pixel count the number of nearest neighbours and store them in an array to set the matrix --- //

    VectorXi nElements_per_row(npix); // 1D vector of integers
    nElements_per_row.setZero();

#pragma omp parallel default(shared) num_threads(nthreads)  
    {
#pragma omp for schedule(static) collapse(2)
    for(int yy=0; yy<ny; ++yy)
      for(int xx=0; xx<nx; ++xx){
	
	int nEl = 1; // diagonal element
	
	// --- count nearest neighbors --- //
	
	nEl += (((xx-1)>=0)? 1 : 0);
	nEl += (((xx+1)<nx)? 1 : 0);
	nEl += (((yy-1)>=0)? 1 : 0);
	nEl += (((yy+1)<ny)? 1 : 0);
	
	nElements_per_row[yy*nx+xx] = nEl;
      }
    
    }
      
    // --- Now set the matrix non-Zero elements of each row --- //
    
    A.reserve(nElements_per_row);


#pragma omp parallel default(shared) num_threads(nthreads)  
    {
    // --- fill the non-zero elements of the matrix --- //
#pragma omp for schedule(static) collapse(2)
    for(int yy=0; yy<ny; ++yy)
      for(int xx=0; xx<nx; ++xx){
	int const ipix = yy*nx + xx;
	
	if((yy-1)>=0) A.insert(ipix,ipix - nx) = malpha;
	if((xx-1)>=0) A.insert(ipix,ipix - 1)  = malpha;
	
	A.insert(ipix,ipix) = lhs[ipix] + alpha*(nElements_per_row[ipix]-1) + beta;

	if((xx+1)<nx) A.insert(ipix,ipix + 1)  = malpha;
	if((yy+1)<ny) A.insert(ipix,ipix + nx) = malpha;
      }
    }

    return A;
  } 

  // ******************************************************************************************* //

  template<typename T>
  void solve_sparse_system(int const npix, Eigen::SparseMatrix<T,Eigen::RowMajor,int> &A,
			   const T* __restrict__ rhs, T* __restrict__ result)
  {
    using namespace Eigen;

    
    // --- Make Eigen maps for the compact arrays --- //
    
    Map<const Matrix<T,Dynamic,1>> B(rhs,    npix);
    Map<      Matrix<T,Dynamic,1>> X(result, npix);


    
    // --- Solve linear system using BIGSTAB/LU decomposition --- //

    BiCGSTAB<SparseMatrix<T,RowMajor,int>> solver(A);
    //SparseLU<SparseMatrix<T, RowMajor, int>> solver; solver.compute(A);

    X = solver.solve(B);
    
	
  }

  // ******************************************************************************************* //

  template<typename T>
  void set_spatial_constraints(int const ny, int const nx, T const alpha, T const beta, const T* __restrict__ lhs,
			       const T* __restrict__ rhs, T* __restrict__ result, int const nthreads = 1)
  {

    if(nthreads > 1){
      Eigen::initParallel();
      Eigen::setNbThreads(nthreads);
    }

    
    // --- construct sparse matrix --- //

    Eigen::SparseMatrix<T,Eigen::RowMajor,int> A =  make_sparse_system<T>(ny, nx, lhs, alpha, beta, nthreads);


    // --- Solve linear system --- //
    
    solve_sparse_system(nx*ny, A, rhs, result);


  }
  
  // ******************************************************************************************* //



  
}


#endif
