//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DualReal.h"
#include "MooseError.h"

#include "libmesh/dense_matrix_base_impl.h"
#include "libmesh/dense_matrix_impl.h"
#include "libmesh/int_range.h"
#include "libmesh/print_trace.h"

#include "metaphysicl/dualnumberarray.h"

namespace libMesh
{
template <typename T>
void
DenseMatrix<T>::_multiply_blas(const DenseMatrixBase<T> &, _BLAS_Multiply_Flag)
{
  mooseError("No blas/lapack support for type ", demangle(typeid(T).name()));
}

template <typename T>
void
DenseMatrix<T>::_lu_decompose_lapack()
{
  mooseError("No blas/lapack support for type ", demangle(typeid(T).name()));
}

template <typename T>
void
DenseMatrix<T>::_svd_lapack(DenseVector<Real> & sigma)
{
  // The calling sequence for dgetrf is:
  // DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )

  // JOBU (input)
  //      Specifies options for computing all or part of the matrix U:
  //      = 'A':  all M columns of U are returned in array U:
  //      = 'S':  the first min(m,n) columns of U (the left singular
  //              vectors) are returned in the array U;
  //      = 'O':  the first min(m,n) columns of U (the left singular
  //              vectors) are overwritten on the array A;
  //      = 'N':  no columns of U (no left singular vectors) are
  //              computed.
  char JOBU = 'N';

  // JOBVT (input)
  //       Specifies options for computing all or part of the matrix
  //       V**T:
  //       = 'A':  all N rows of V**T are returned in the array VT;
  //       = 'S':  the first min(m,n) rows of V**T (the right singular
  //               vectors) are returned in the array VT;
  //       = 'O':  the first min(m,n) rows of V**T (the right singular
  //               vectors) are overwritten on the array A;
  //       = 'N':  no rows of V**T (no right singular vectors) are
  //               computed.
  char JOBVT = 'N';

  std::vector<Real> sigma_val;
  std::vector<Number> U_val;
  std::vector<Number> VT_val;

  _svd_helper(JOBU, JOBVT, sigma_val, U_val, VT_val);

  // Copy the singular values into sigma, ignore U_val and VT_val
  sigma.resize(cast_int<unsigned int>(sigma_val.size()));
  for (auto i : IntRange<int>(0, sigma.size()))
    sigma(i) = sigma_val[i];
}

template <typename T>
void
DenseMatrix<T>::_svd_lapack(DenseVector<Real> & sigma,
                            DenseMatrix<Number> & U,
                            DenseMatrix<Number> & VT)
{
  // The calling sequence for dgetrf is:
  // DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )

  // JOBU (input)
  //      Specifies options for computing all or part of the matrix U:
  //      = 'A':  all M columns of U are returned in array U:
  //      = 'S':  the first min(m,n) columns of U (the left singular
  //              vectors) are returned in the array U;
  //      = 'O':  the first min(m,n) columns of U (the left singular
  //              vectors) are overwritten on the array A;
  //      = 'N':  no columns of U (no left singular vectors) are
  //              computed.
  char JOBU = 'S';

  // JOBVT (input)
  //       Specifies options for computing all or part of the matrix
  //       V**T:
  //       = 'A':  all N rows of V**T are returned in the array VT;
  //       = 'S':  the first min(m,n) rows of V**T (the right singular
  //               vectors) are returned in the array VT;
  //       = 'O':  the first min(m,n) rows of V**T (the right singular
  //               vectors) are overwritten on the array A;
  //       = 'N':  no rows of V**T (no right singular vectors) are
  //               computed.
  char JOBVT = 'S';

  // Note: Lapack is going to compute the singular values of A^T.  If
  // A=U * S * V^T, then A^T = V * S * U^T, which means that the
  // values returned in the "U_val" array actually correspond to the
  // entries of the V matrix, and the values returned in the VT_val
  // array actually correspond to the entries of U^T.  Therefore, we
  // pass VT in the place of U and U in the place of VT below!
  std::vector<Real> sigma_val;
  int M = this->n();
  int N = this->m();
  int min_MN = (M < N) ? M : N;

  // Size user-provided storage appropriately. Inside svd_helper:
  // U_val is sized to (M x min_MN)
  // VT_val is sized to (min_MN x N)
  // So, we set up U to have the shape of "VT_val^T", and VT to
  // have the shape of "U_val^T".
  //
  // Finally, since the results are stored in column-major order by
  // Lapack, but we actually want the transpose of what Lapack
  // returns, this means (conveniently) that we don't even have to
  // copy anything after the call to _svd_helper, it should already be
  // in the correct order!
  U.resize(N, min_MN);
  VT.resize(min_MN, M);

  _svd_helper(JOBU, JOBVT, sigma_val, VT.get_values(), U.get_values());

  // Copy the singular values into sigma.
  sigma.resize(cast_int<unsigned int>(sigma_val.size()));
  for (auto i : IntRange<int>(0, sigma.size()))
    sigma(i) = sigma_val[i];
}

template <typename T>
void
DenseMatrix<T>::_svd_helper(
    char, char, std::vector<Real> &, std::vector<Number> &, std::vector<Number> &)
{
  mooseError("No blas/lapack support for type ", demangle(typeid(T).name()));
}

template <typename T>
void
DenseMatrix<T>::_svd_solve_lapack(const DenseVector<T> & /*rhs*/,
                                  DenseVector<T> & /*x*/,
                                  Real /*rcond*/) const
{
  mooseError("No blas/lapack support for type ", demangle(typeid(T).name()));
}

template <typename T>
void
DenseMatrix<T>::_evd_lapack(DenseVector<T> &, DenseVector<T> &, DenseMatrix<T> *, DenseMatrix<T> *)
{
  mooseError("No blas/lapack support for type ", demangle(typeid(T).name()));
}

template <typename T>
void
DenseMatrix<T>::_lu_back_substitute_lapack(const DenseVector<T> &, DenseVector<T> &)
{
  mooseError("No blas/lapack support for type ", demangle(typeid(T).name()));
}

template <typename T>
void
DenseMatrix<T>::_matvec_blas(T, T, DenseVector<T> &, const DenseVector<T> &, bool) const
{
  mooseError("No blas/lapack support for type ", demangle(typeid(T).name()));
}

template class DenseMatrixBase<DualReal>;
template class DenseMatrix<DualReal>;
}
