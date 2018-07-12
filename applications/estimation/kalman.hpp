// kalman.hpp
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <boost/numeric/mtl/mtl.hpp>

template<typename Ty>
class KalmanFilter {
public:

	/**
	* Create a Kalman filter with the specified matrices.
	*   A - System dynamics matrix
	*   C - Output matrix
	*   Q - Process noise covariance
	*   R - Measurement noise covariance
	*   P - Estimate error covariance
	*/
	KalmanFilter(
		double dt,
		const mtl::dense2D<Ty>& A,
		const mtl::dense2D<Ty>& C,
		const mtl::dense2D<Ty>& Q,
		const mtl::dense2D<Ty>& R,
		const mtl::dense2D<Ty>& P
	) : A(A), C(C), Q(Q), R(R), P0(P),
		m(C.rows()), n(A.rows()), dt(dt), initialized(false),
		I(n, n), x_hat(n)
	{
		I.setIdentity();
	}
	KalmanFilter() {}

	/**
	* Initialize the filter with initial states as zero.
	*/
	void init() {
		x_hat.setToZero();
		P = P0;
		t0 = 0;
		t = t0;
	}

	/**
	* Initialize the filter with a guess for initial states.
	*/
	void init(double _t0, const mtl::dense_vector<Ty>& _x0) {
		x_hat = _x0;
		P = P0;
		t0 = _t0;
		t = t0;
	}

	/**
	* Update the estimated state based on measured values. The
	* time step is assumed to remain constant.
	*/
	void update(const mtl::dense_vector<Ty>& y) {
		mtl::dense_vector<Ty> x_hat_new(n);
		x_hat_new = A * x_hat;
		P = A * P * trans(A) + Q;
		K = P * trans(C) * (C * P * trans(C) + R).inverse();
		x_hat_new += K * (y - C * x_hat_new);
		P = (I - K * C) * P;
		x_hat = x_hat_new;
	}

	/**
	* Update the estimated state based on measured values,
	* using the given time step and dynamics matrix.
	*/
	void update(const mtl::dense_vector<Ty>& _y, double _dt, const mtl::dense2D<Ty> _A) {
		A = _A;
		dt = _dt;
		update(_y);
	}

	/**
	* Return the current state and time.
	*/
	mtl::dense_vector<Ty> state() { return x_hat; };
	double time() { return t; };

private:

	// Matrices for computation
	mtl::dense2D<Ty> A, C, Q, R, P, K, P0;

	// System dimensions
	int m, n;

	// Initial and current time
	double t0, t;

	// Discrete time step
	double dt;

	// n-size identity
	mtl::dense2D<Ty> I;

	// Estimated states
	mtl::dense_vector<Ty> x_hat;
};
