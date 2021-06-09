#pragma once
#include <universal/blas/matrix.hpp>
#include <norms.hpp>

namespace sw
{
    namespace hprblas
    {
        template <typename Scalar>
        auto qr(const universal::blas::matrix<Scalar> &A,
                universal::blas::matrix<Scalar> &Q,
                universal::blas::matrix<Scalar> &R)
        {
            using size_type = typename universal::blas::matrix<Scalar>::size_type;
            size_type n = num_rows(A), m = num_cols(A);
            for (size_type i = 0; i < m; ++i)
            {
                universal::blas::matrix<Scalar> v; // Store i'th column of A (A[(rows 0->n), i'th column])
                /*
    Need a function which return's a matrix consisting
    of i'th column of A (similar to slice in matlab)
    */
                if (i > 0)
                {
                    /*
      Need a function similar to slice of a matrix in matlab

      R[(rows 0->i), i'th column] = Q[(rows 0->n, columns 0->i)]' * A[rows 0->n,
      i'th column] v= A[(rows 0->n), ith column] - Q[(rows 0->n), (column 0->k)]
      * R[(rows 0->i), i't column]
      */
                }
                R(i, i) = frobenius_norm(v);
                // i'th column of Q = v/ R(i,i)
            }
            return std::make_pair(Q, R);
        }
        template <typename Scalar>
        auto qr(const universal::blas::matrix<Scalar> &A)
        {
            size_type m = num_cols(A), n = num_rows(A);
            universal::blas::matrix<Scalar> Q(n, m), R(n, n);
            qr(A, Q, R);
            return std::make_pair(Q, R);
        }

    } // namespace hprblas
} // namespace sw
