    #include "added.h"
#ifndef _Vector_functions_hpp_
#define _Vector_functions_hpp_

//@HEADER
// ************************************************************************
// 
//               HPCCG: Simple Conjugate Gradient Benchmark Code
//                 Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <vector>
#include <sstream>
#include <fstream>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <TypeTraits.hpp>
#include <Vector.hpp>

// Import the Intel MKL Library
#include "mkl.h"

#define MINIFE_MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

namespace miniFE {


template<typename VectorType>
void write_vector(const std::string& filename,
                  const VectorType& vec)
{
  int numprocs = 1, myproc = 0;
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#endif

  std::ostringstream osstr;
  osstr << filename << "." << numprocs << "." << myproc;
  std::string full_name = osstr.str();
  std::ofstream ofs(full_name.c_str());

  typedef typename VectorType::ScalarType ScalarType;

  const std::vector<ScalarType>& coefs = vec.coefs;
  for(int p=0; p<numprocs; ++p) {
    if (p == myproc) {
      if (p == 0) {
        ofs << vec.local_size << std::endl;
      }
  
      typename VectorType::GlobalOrdinalType first = vec.startIndex;
      for(size_t i=0; i<vec.local_size; ++i) {
        ofs << first+i << " " << coefs[i] << std::endl;
      }
    }
#ifdef HAVE_MPI
    t_start = MPI_Wtime(); //added
    MPI_Barrier(MPI_COMM_WORLD);
    t_stop = MPI_Wtime(); //added
    t_latency = (double)((t_stop - t_start) * 1e6);
    sprintf(my_str[numLine], "MPI_Barrier-2: time: %lf \n", t_latency);
    numLine++;
#endif
  }
}

template<typename VectorType>
void sum_into_vector(size_t num_indices,
                     const typename VectorType::GlobalOrdinalType* indices,
                     const typename VectorType::ScalarType* coefs,
                     VectorType& vec)
{
  typedef typename VectorType::GlobalOrdinalType GlobalOrdinal;
  typedef typename VectorType::ScalarType Scalar;

  GlobalOrdinal first = vec.startIndex;
  GlobalOrdinal last = first + vec.local_size - 1;

  std::vector<Scalar>& vec_coefs = vec.coefs;

  for(size_t i=0; i<num_indices; ++i) {
    if (indices[i] < first || indices[i] > last) continue;
    size_t idx = indices[i] - first;
    vec_coefs[idx] += coefs[i];
  }
}

#ifdef MINIFE_HAVE_TBB
template<typename VectorType>
void sum_into_vector(size_t num_indices,
                     const typename VectorType::GlobalOrdinalType* indices,
                     const typename VectorType::ScalarType* coefs,
                     LockingVector<VectorType>& vec)
{
  vec.sum_in(num_indices, indices, coefs);
}
#endif

//------------------------------------------------------------
//Compute the update of a vector with the sum of two scaled vectors where:
//
// w = alpha*x + beta*y
//
// x,y - input vectors
//
// alpha,beta - scalars applied to x and y respectively
//
// w - output vector
//
template<typename VectorType>
void
  waxpby(typename VectorType::ScalarType alpha, const VectorType& x,
         typename VectorType::ScalarType beta, const VectorType& y,
         VectorType& w)
{
  typedef typename VectorType::ScalarType ScalarType;

#ifdef MINIFE_DEBUG
  if (y.local_size < x.local_size || w.local_size < x.local_size) {
    std::cerr << "miniFE::waxpby ERROR, y and w must be at least as long as x." << std::endl;
    return;
  }
#endif

  int n = x.coefs.size();
  const ScalarType* xcoefs = &x.coefs[0];
  const ScalarType* ycoefs = &y.coefs[0];
        ScalarType* wcoefs = &w.coefs[0];
  
  #pragma omp parallel for
  for(int i=0; i<n; ++i) {
    wcoefs[i] = alpha*xcoefs[i] + beta*ycoefs[i];
  }
}

template<typename VectorType>
void
  daxpby(typename VectorType::ScalarType alpha, const VectorType& x,
         typename VectorType::ScalarType beta, VectorType& y)
{
  typedef typename VectorType::ScalarType ScalarType;

  const int n = MINIFE_MIN(x.coefs.size(), y.coefs.size());
  const MINIFE_SCALAR* xcoefs = &x.coefs[0];
        MINIFE_SCALAR* ycoefs = &y.coefs[0];

  #if defined(MINIFE_MKL_DOUBLE)
	cblas_daxpby(
  #elif defined(MINIFE_MKL_FLOAT)
	cblas_saxpby(
  #else
	#error "Unknown MINIFE MKL type"
  #endif
	(const MKL_INT) n,
	alpha,
	xcoefs,
	(const MKL_INT) 1,
	beta,
	ycoefs,
	(const MKL_INT) 1);
}

//Like waxpby above, except operates on two sets of arguments.
//In other words, performs two waxpby operations in one loop.
template<typename VectorType>
void
  fused_waxpby(typename VectorType::ScalarType alpha, const VectorType& x,
         typename VectorType::ScalarType beta, const VectorType& y,
         VectorType& w,
         typename VectorType::ScalarType alpha2, const VectorType& x2,
         typename VectorType::ScalarType beta2, const VectorType& y2,
         VectorType& w2)
{
  typedef typename VectorType::ScalarType ScalarType;

#ifdef MINIFE_DEBUG
  if (y.local_size < x.local_size || w.local_size < x.local_size) {
    std::cerr << "miniFE::waxpby ERROR, y and w must be at least as long as x." << std::endl;
    return;
  }
#endif

  int n = x.coefs.size();
  const ScalarType* xcoefs = &x.coefs[0];
  const ScalarType* ycoefs = &y.coefs[0];
        ScalarType* wcoefs = &w.coefs[0];

  const ScalarType* x2coefs = &x2.coefs[0];
  const ScalarType* y2coefs = &y2.coefs[0];
        ScalarType* w2coefs = &w2.coefs[0];

  for(int i=0; i<n; ++i) {
    wcoefs[i] = alpha*xcoefs[i] + beta*ycoefs[i];
    w2coefs[i] = alpha2*x2coefs[i] + beta2*y2coefs[i];
  }
}

//-----------------------------------------------------------
//Compute the dot product of two vectors where:
//
// x,y - input vectors
//
// result - return-value
//
template<typename Vector>
typename TypeTraits<typename Vector::ScalarType>::magnitude_type
  dot(const Vector& x,
      const Vector& y)
{
  int n = x.coefs.size();

#ifdef MINIFE_DEBUG
  if (y.local_size < n) {
    std::cerr << "miniFE::dot ERROR, y must be at least as long as x."<<std::endl;
    n = y.local_size;
  }
#endif

  typedef typename Vector::ScalarType Scalar;
  typedef typename TypeTraits<typename Vector::ScalarType>::magnitude_type magnitude;

  const Scalar* xcoefs = &x.coefs[0];
  const Scalar* ycoefs = &y.coefs[0];
#if defined(MINIFE_MKL_DOUBLE)
  magnitude result = cblas_ddot(
#elif defined(MINIFE_MKL_FLOAT)
  magnitude result = cblas_sdot(
#else
  #error Unknown MINIFE_SCALAR type.
#endif
	(MKL_INT) n,
	xcoefs,
	(MKL_INT) 1,
	ycoefs,
	(MKL_INT) 1);

#ifdef HAVE_MPI
  magnitude local_dot = result, global_dot = 0;
  MPI_Datatype mpi_dtype = TypeTraits<magnitude>::mpi_type();  
    t_start = MPI_Wtime(); //added
  MPI_Allreduce(&local_dot, &global_dot, 1, mpi_dtype, MPI_SUM, MPI_COMM_WORLD);
    t_stop = MPI_Wtime(); //added
    t_latency = (double)((t_stop - t_start) * 1e6);
    sprintf(my_str[numLine], "MPI_Allreduce-6: time: %lf \n", t_latency);
    numLine++;
    printf("numLine is: %d\n", numLine);
  return global_dot;
#else
  return result;
#endif
}

}//namespace miniFE

#endif

