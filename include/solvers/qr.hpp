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
                universal::blas::matrix<Scalar> v = help(A, 0, n, q);

                if (i > 0)
                {
                    help(R, 0, i, i, i) = transpose(help(R, 0, n, 0, i)) * help(A, 0, n, i, i);
                    //Problem => function should change value of the given matrix
                    v = help(A, 0, n, i, i) - help(Q, 0, n, 0, i) * help(R, 0, i, i, i);
                }
                R(i, i) = frobenius_norm(v);
                help(Q, 0, n, i, i) = v / R(i, i);
            }
            return std::make_pair(Q, R);
        }
        template <typename Scalar>
        auto qr(const universal::blas::matrix<Scalar> &A)
        {
            using size_type = typename universal::blas::matrix<Scalar>::size_type;
            size_type m = num_cols(A), n = num_rows(A);
            universal::blas::matrix<Scalar> Q(n, m), R(n, n);
            qr(A, Q, R);
            return std::make_pair(Q, R);
        }
        template <typename Scalar>
        universal::blas::matrix<Scalar> help(const universal::blas::matrix<Scalar> &A, Scalar rb, Scalar re, Scalar cb, Scalar ce)
        {
            using size_type = typename universal::blas::matrix<Scalar>::size_type;
            universal::blas::matrix<Scalar> B(re - rb + 1, ce - cb + 1);
            size_type p = 1, q = 1;
            for (size_type i = rb; i <= re; ++i)
            {
                for (size_type j = cb; j <= ce; ++j)
                {
                    B(p, q) = A(i, j);
                    ++q;
                }
                ++p;
            }
            return B;
        }

    } // namespace hprblas
} // namespace sw
