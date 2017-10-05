#include "eig.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"
#include "real.h"
#include "settings.h"
#include "algorithms.h"

//void qr_simple_step(const CMatrix_t* A, CMatrix_t* Q, CMatrix_t* R)
//{
//	//function[Q, R] = QRfactor(A)
//	//% код взят из интернета
//	//% https://www.mathworks.com/matlabcentral/answers/169648-qr-factorization-using-householder-transformations
//	//[m, n] = size(A);
//	//R = A; %Start with R = A
//	//Q = eye(m); %Set Q as the identity matrix
//	//for k = 1:m - 1
//	//	x = zeros(m, 1);
//	//	x(k:m, 1) = R(k:m, k);
//	//	g = norm(x);
//	//	v = x; v(k) = x(k) + g;
//	//	%Orthogonal transformation matrix that eliminates one element
//	//		%below the diagonal of the matrix it is post - multiplying:
//	//	s = norm(v);
//	//	if s~= 0, w = v / s;
//	//		u = 2 * R'*w;
//	//		R = R - w*u'; %Product HR
//	//		Q = Q - 2 * Q*(w*w'); %Product QR
//	//	end
//	//end
//
//	int N = A->numCols;
//	int i, j, idx, k;
//
//	CMatrix_t x;
//	REAL_TYPE g, s;
//	CMatrix_t v;
//	CMatrix_t u;
//	CMatrix_t Rt;
//	CMatrix_t ut;
//	CMatrix_t wut;
//	CMatrix_t wt;
//	CMatrix_t wwt;
//	CMatrix_t Qwwt;
//
//	x.numRows = N;
//	x.numCols = 1;
//	x.pDataReal = (REAL_TYPE*)malloc(N * sizeof(REAL_TYPE));
//	x.pDataImag = (REAL_TYPE*)malloc(N * sizeof(REAL_TYPE));
//
//	v.numRows = N;
//	v.numCols = 1;
//	v.pDataReal = (REAL_TYPE*)malloc(N * sizeof(REAL_TYPE));
//	v.pDataImag = (REAL_TYPE*)malloc(N * sizeof(REAL_TYPE));
//
//	u.numRows = N;
//	u.numCols = 1;
//	u.pDataReal = (REAL_TYPE*)malloc(N * sizeof(REAL_TYPE));
//	u.pDataImag = (REAL_TYPE*)malloc(N * sizeof(REAL_TYPE));
//
//	Rt.numCols = N;
//	Rt.numRows = N;
//	Rt.pDataReal = (REAL_TYPE*)malloc(N *  N * sizeof(REAL_TYPE));
//	Rt.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//
//	ut.numRows = 1;
//	ut.numCols = N;
//	ut.pDataReal = (REAL_TYPE*)malloc(N  * sizeof(REAL_TYPE));
//	ut.pDataImag = (REAL_TYPE*)malloc(N  * sizeof(REAL_TYPE));
//
//	wut.numCols = N;
//	wut.numRows = N;
//	wut.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	wut.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//
//	wt.numRows = 1;
//	wt.numCols = N;
//	wt.pDataReal = (REAL_TYPE*)malloc(N * sizeof(REAL_TYPE));
//	wt.pDataImag = (REAL_TYPE*)malloc(N * sizeof(REAL_TYPE));
//
//
//	wwt.numCols = N;
//	wwt.numRows = N;
//	wwt.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	wwt.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//
//	Qwwt.numCols = N;
//	Qwwt.numRows = N;
//	Qwwt.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	Qwwt.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//
//
//	// R = A; %Start with R = A
//	//	Q = eye(m); %Set Q as the identity matrix
//	for (i = 0; i < N; ++i)
//	{
//		for (j = 0; j < N; ++j)
//		{
//			idx = i * A->numCols + j;
//			Q->pDataReal[idx] = (i == j) ? 1.0f : 0.0f;
//			Q->pDataImag[idx] = 0.0;
//			R->pDataReal[idx] = A->pDataReal[idx];
//			R->pDataImag[idx] = A->pDataImag[idx];
//		}
//	}
//
//	//	for k = 1:m - 1
//	for (k = 0; k < N - 1; ++k)
//	{
//		// x = zeros(m, 1);
//		for (idx = 0; idx < N; ++idx)
//		{
//			x.pDataReal[idx] = 0.0;
//			x.pDataImag[idx] = 0.0;
//		}
//
//		// x(k:m, 1) = R(k:m, k);
//		for (idx = k; idx < N; ++idx)
//		{
//			x.pDataReal[idx] = R->pDataReal[idx * N + k];
//			x.pDataImag[idx] = R->pDataImag[idx * N + k];
//		}
//
//		// g = norm(x);
//		cnorm2(x.pDataReal, x.pDataImag, &g, 1, N);
//		g = sq_rt(g);
//
//		// v = x; v(k) = x(k) + g;
//		for (idx = 0; idx < N; ++idx)
//		{
//			v.pDataReal[idx] = x.pDataReal[idx] + g * (idx == k);
//			v.pDataImag[idx] = x.pDataImag[idx];
//		}
//
//		// s = norm(v);
//		cnorm2(v.pDataReal, v.pDataImag, &s, 1, N);
//		s = sq_rt(s);
//
//		// if s~= 0, w = v / s;
//		// здесь w это v
//		if (s < 1e-6)
//			continue;
//		for (idx = 0; idx < N; ++idx)
//		{
//			v.pDataReal[idx] /= s;
//			v.pDataImag[idx] /= s;
//		}
//
//		// u = 2 * R'*w;
//		complex_transp(R, &Rt, 0);
//		complex_matrix_mult(&Rt, &v, 2.0f, 0.0f, &u);
//
//		// R = R - w*u'; %Product HR
//		complex_transp(&u, &ut, 0);
//		complex_matrix_mult(&v, &ut, 1.0f, 0.0f, &wut);
//		for (i = 0; i < N; ++i)
//		{
//			for (j = 0; j < N; ++j)
//			{
//				idx = i * N + j;
//				R->pDataReal[idx] -= wut.pDataReal[idx];
//				R->pDataImag[idx] -= wut.pDataImag[idx];
//			}
//		}
//
//		// Q = Q - 2 * Q*(w*w'); %Product QR
//		complex_transp(&v, &wt, 0);
//		complex_matrix_mult(&v, &wt, 1.0f, 0.0f, &wwt);
//		complex_matrix_mult(Q, &wwt, 2.0f, 0.0f, &Qwwt);
//
//		for (i = 0; i < N; ++i)
//		{
//			for (j = 0; j < N; ++j)
//			{
//				idx = i * N + j;
//				Q->pDataReal[idx] -= Qwwt.pDataReal[idx];
//				Q->pDataImag[idx] -= Qwwt.pDataImag[idx];
//			}
//		}
//	}
//
//	free(x.pDataReal);
//	free(x.pDataImag);
//
//	free(v.pDataReal);
//	free(v.pDataImag);
//
//	free(u.pDataReal);
//	free(u.pDataImag);
//
//	free(Rt.pDataReal);
//	free(Rt.pDataImag);
//
//	free(ut.pDataReal);
//	free(ut.pDataImag);
//
//	free(wut.pDataReal);
//	free(wut.pDataImag);
//
//	free(wt.pDataReal);
//	free(wt.pDataImag);
//
//	free(wwt.pDataReal);
//	free(wwt.pDataImag);
//
//	free(Qwwt.pDataReal);
//	free(Qwwt.pDataImag);
//}
//
//void eig(const CMatrix_t* A, CMatrix_t* e, CMatrix_t* ev)
//{
//	int N = A->numCols, i, j;
//	CMatrix_t acopy;
//	CMatrix_t Q;
//	CMatrix_t R;
//	CMatrix_t Qt;
//	CMatrix_t AQ;
//	CMatrix_t evcopy;
//
//	acopy.numRows = N;
//	acopy.numCols = N;
//	acopy.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	acopy.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	memcpy(acopy.pDataReal, A->pDataReal, N * N * sizeof(REAL_TYPE));
//	memcpy(acopy.pDataImag, A->pDataImag, N * N * sizeof(REAL_TYPE));
//
//	Q.numRows = N;
//	Q.numCols = N;
//	Q.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	Q.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	memset(Q.pDataReal, 0, N * N * sizeof(REAL_TYPE));
//	memset(Q.pDataImag, 0, N * N * sizeof(REAL_TYPE));
//
//	R.numRows = N;
//	R.numCols = N;
//	R.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	R.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	memset(R.pDataReal, 0, N * N * sizeof(REAL_TYPE));
//	memset(R.pDataImag, 0, N * N * sizeof(REAL_TYPE));
//
//	Qt.numRows = N;
//	Qt.numCols = N;
//	Qt.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	Qt.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	memset(Qt.pDataReal, 0, N * N * sizeof(REAL_TYPE));
//	memset(Qt.pDataImag, 0, N * N * sizeof(REAL_TYPE));
//
//	AQ.numRows = N;
//	AQ.numCols = N;
//	AQ.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	AQ.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	memset(AQ.pDataReal, 0, N * N * sizeof(REAL_TYPE));
//	memset(AQ.pDataImag, 0, N * N * sizeof(REAL_TYPE));
//
//	evcopy.numRows = N;
//	evcopy.numCols = N;
//	evcopy.pDataReal = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//	evcopy.pDataImag = (REAL_TYPE*)malloc(N * N * sizeof(REAL_TYPE));
//
//	for (i = 0; i < N; ++i)
//	{
//		for (j = 0; j < N; ++j)
//		{
//			evcopy.pDataReal[i * N + j] = (i == j) ? 1.0f : 0.0f;
//			evcopy.pDataImag[i * N + j] = 0.0f;
//		}
//	}
//
//	/*for i = 1 : Nsteps
//		[Q, ~] = QRfactor_s(A);
//		A = Q' * A * Q;
//	end*/
//	for (i = 0; i < N; ++i)
//	{
//		qr_simple_step(&acopy, &Q, &R);
//		complex_matrix_mult(&acopy, &Q, 1.0f, 0.0f, &AQ);
//		complex_transp(&Q, &Qt, 0);
//		complex_matrix_mult(&Qt, &AQ, 1.0f, 0.0f, &acopy);
//		complex_matrix_mult(&evcopy, &Q, 1.0f, 0.0f, ev);
//		memcpy(evcopy.pDataReal, ev->pDataReal, N * N * sizeof(REAL_TYPE));
//		memcpy(evcopy.pDataImag, ev->pDataImag, N * N * sizeof(REAL_TYPE));
//	}
//
//	for (i = 0; i < N; ++i)
//	{
//		e->pDataReal[i] = A->pDataReal[i * N + i];
//		e->pDataImag[i] = A->pDataImag[i * N + i];
//	}
//
//	free(acopy.pDataReal);
//	free(acopy.pDataImag);
//
//	free(Q.pDataReal);
//	free(Q.pDataImag);
//
//	free(R.pDataReal);
//	free(R.pDataImag);
//
//	free(Qt.pDataReal);
//	free(Qt.pDataImag);
//
//	free(AQ.pDataReal);
//	free(AQ.pDataImag);
//
//	free(evcopy.pDataReal);
//	free(evcopy.pDataImag);
//}

void complex_triag(const CMatrix_t* A, Matrix_t* T, CMatrix_t* P)
// Трехдиагонализация Хаусхолдера
//////////////////////////////
//// Из статьи: Joachim Kopp. Efficient numerical diagonalization of
//// hermitian matrices [http://arXiv.org/abs/physics/0610206v3]
//////////////////////////////////////
//	function [T, P] = triag(A)
//		m = size(A,1);
//		P = eye(m, m, 'like', A);
//		T = A;
//		for i = 2:m
//			x = zeros(m,1,'like',A);
//			x(i:m) = T(i:m,i-1);
//			e = zeros(m,1,'like',A);
//			e(i) = 1;
//
//			u = x - csgn(x(i)) * norm(x) * e;
//
//			w = 1/(u'*u)*(1 + x'*u/(u'*x));
//
//			p = w' * T * u;
//			K = w/2 * u'*p;
//			q = p - K*u;
//
//			T = T - q*u' - u*q';
//
//			Pk = eye(m,m,'like',A);
//			Pk = Pk - w*u*u';
//
//			P = Pk * P;
//		end
//		T = real(T);
// end
{
	int i, create = 1;
	//            m     v  v  v  v    v  v  v  v    m     m    m     m      m    m    m      m
	CMatrix_td TComplex, x, e, u, xt, ut, p, q, qt, Tcopy, qut, uqt, qut_uqt, Pk, uut, Pcopy, Ad, Pd;
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
	double normx, wre, wim, Kre, Kim;

	//		m = size(A,1);
	int m = A->numRows;

#ifdef PREALLOCATION
#ifdef DOUBLE
	double buffer[2 *(MATRIX_MAX_SIZE*8 + VECTOR_MAX_SIZE*8)];
#else // DOUBLE
	double buffer[2 * (MATRIX_MAX_SIZE * 10 + VECTOR_MAX_SIZE * 8)];
#endif // DOUBLE
	TComplex.pDataReal = buffer;
	TComplex.pDataImag = buffer + MATRIX_MAX_SIZE;
	Tcopy.pDataReal = buffer + 2 * MATRIX_MAX_SIZE;
	Tcopy.pDataImag = buffer + 3 * MATRIX_MAX_SIZE;
	qut.pDataReal = buffer + 4 * MATRIX_MAX_SIZE;
	qut.pDataImag = buffer + 5 * MATRIX_MAX_SIZE;
	uqt.pDataReal = buffer + 6 * MATRIX_MAX_SIZE;
	uqt.pDataImag = buffer + 7 * MATRIX_MAX_SIZE;
	qut_uqt.pDataReal = buffer + 8 * MATRIX_MAX_SIZE;
	qut_uqt.pDataImag = buffer + 9 * MATRIX_MAX_SIZE;
	Pk.pDataReal = buffer + 10 * MATRIX_MAX_SIZE;
	Pk.pDataImag = buffer + 11 * MATRIX_MAX_SIZE;
	uut.pDataReal = buffer + 12 * MATRIX_MAX_SIZE;
	uut.pDataImag = buffer + 13 * MATRIX_MAX_SIZE;
	Pcopy.pDataReal = buffer + 14 * MATRIX_MAX_SIZE;
	Pcopy.pDataImag = buffer + 15 * MATRIX_MAX_SIZE;
	x.pDataReal = buffer + 16 * MATRIX_MAX_SIZE;
	x.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + VECTOR_MAX_SIZE;
	e.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 2 * VECTOR_MAX_SIZE;
	e.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 3 * VECTOR_MAX_SIZE;
	u.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 4 * VECTOR_MAX_SIZE;
	u.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 5 * VECTOR_MAX_SIZE;
	xt.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 6 * VECTOR_MAX_SIZE;
	xt.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 7 * VECTOR_MAX_SIZE;
	ut.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 8 * VECTOR_MAX_SIZE;
	ut.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 9 * VECTOR_MAX_SIZE;
	p.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 10 * VECTOR_MAX_SIZE;
	p.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 11 * VECTOR_MAX_SIZE;
	q.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 12 * VECTOR_MAX_SIZE;
	q.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 13 * VECTOR_MAX_SIZE;
	qt.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 14 * VECTOR_MAX_SIZE;
	qt.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 15 * VECTOR_MAX_SIZE;
#ifndef DOUBLE
	Ad.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 16 * VECTOR_MAX_SIZE;
	Ad.pDataImag = buffer + 17 * MATRIX_MAX_SIZE + 16 * VECTOR_MAX_SIZE;
	Pd.pDataReal = buffer + 18 * MATRIX_MAX_SIZE + 16 * VECTOR_MAX_SIZE;
	Pd.pDataImag = buffer + 19 * MATRIX_MAX_SIZE + 16 * VECTOR_MAX_SIZE;
#endif // DOUBLE
	Pd.pDataReal = P->pDataReal;
	Pd.pDataImag = P->pDataImag;
	Pd.numCols = P->numCols;
	Pd.numRows = P->numRows;
#endif // PREALLOCATION

#ifndef DOUBLE
	single_to_double(A, &Ad, 1);
#else // DOUBLE
	Ad.numRows = A->numRows;
	Ad.numCols = A->numCols;
	Ad.pDataReal = A->pDataReal;
	Ad.pDataImag = A->pDataImag;
#endif // DOUBLE

	//		P = eye(m, m, 'like', A);
	complex_eye_d(m, &Pd, ALLOCATE_MATRIX);
	//		T = A;
	complex_clone_d(&Ad, &TComplex, ALLOCATE_MATRIX);

	//		for i = 2:m
	for (i = 1; i <= m-1; ++i)
	{
		if (create > 0)
		{
			complex_new_d(m, 1, &u);
			complex_new_d(m, 1, &p);
			complex_new_d(m, 1, &q);
			complex_new_d(m, m, &qut);
			complex_new_d(m, m, &uqt);
			complex_new_d(m, m, &qut_uqt);
			complex_new_d(m, m, &uut);
		}

		//			x = zeros(m,1,'like',A);
		complex_zeros_d(m, 1, &x, create);
		//			x(i:m) = T(i:m,i-1);
		complex_partial_copy_d(&TComplex, i, i-1, m-i, 1, &x, i, 0);
		//			e = zeros(m,1,'like',A);
		complex_zeros_d(m, 1, &e, create);
		//			e(i) = 1;
		e.pDataReal[i] = 1.0f;

		//			u = x - csgn(x(i)) * norm(x) * e;
		f1 = complex_sign_d(x.pDataReal[i], x.pDataImag[i]);
		cnorm2_d(x.pDataReal, x.pDataImag, &f2, 1, m);
		normx = sqrt(f2);
		complex_matrix_sum_d(&x, &e, -f1*normx, 0.0f, &u);

		//			w = 1/(u'*u)*(1 + x'*u/(u'*x));
		cnorm2_d(u.pDataReal, u.pDataImag, &f3, 1, m);
		f4 = 1.0f/f3;
		complex_transp_d(&u, &ut, create);
		complex_transp_d(&x, &xt, create);
		complex_scal_prod_d(&xt, &u, &f5, &f6);
		complex_scal_prod_d(&x, &ut, &f7, &f8);
		complex_div_d(f5, f6, f7, f8, &f9, &f10);
		wre = f4 * (1.0 + f9);
		wim = f4 * f10;

		//			p = w' * T * u;
		complex_matrix_mult_d(&TComplex, &u, wre, -wim, &p);
		//			K = w/2 * u'*p;
		complex_scal_prod_d(&ut, &p, &f11, &f12);
		complex_mult_d(wre/2.0f, wim/2.0f, f11, f12, &Kre, &Kim);
		//			q = p - K*u;
		complex_matrix_sum_d(&p, &u, -Kre, -Kim, &q);

		//			T = T - q*u' - u*q';
		complex_transp_d(&q, &qt, create);
		complex_matrix_mult_d(&q, &ut, 1.0f, 0.0f, &qut);
		complex_matrix_mult_d(&u, &qt, 1.0f, 0.0f, &uqt);
		complex_matrix_sum_d(&qut, &uqt, 1.0f, 0.0f, &qut_uqt);
		complex_clone_d(&TComplex, &Tcopy, create);
		complex_matrix_sum_d(&TComplex, &qut_uqt, -1.0f, 0.0f, &TComplex);


		//			Pk = eye(m,m,'like',A);
		complex_eye_d(m, &Pk, create);
		//			Pk = Pk - w*u*u';
		complex_matrix_mult_d(&u, &ut, 1.0f, 0.0f, &uut);
		complex_matrix_sum_d(&Pk, &uut, -wre, -wim, &Pk);

		//			P = Pk * P;
		complex_clone_d(&Pd, &Pcopy, create);
		complex_matrix_mult_d(&Pk, &Pcopy, 1.0f, 0.0f, &Pd);

		create = 0;
	}
	//		T = real(T);
#ifdef DOUBLE
	memcpy(T->pData, TComplex.pDataReal, m*m*sizeof(REAL_TYPE));
#else // DOUBLE
	double_complex_to_single_real(&TComplex, T, 0);
	double_complex_to_single_comlpex(&Pd, P, 0);
#endif // DOUBLE

	complex_free_d(&TComplex);
	complex_free_d(&x);
	complex_free_d(&e);
	complex_free_d(&u);
	complex_free_d(&xt);
	complex_free_d(&ut);
	complex_free_d(&p);
	complex_free_d(&q);
	complex_free_d(&qt);
	complex_free_d(&Tcopy);
	complex_free_d(&qut);
	complex_free_d(&uqt);
	complex_free_d(&qut_uqt);
	complex_free_d(&Pk);
	complex_free_d(&uut);
	complex_free_d(&Pcopy);
}

void givens_rotation_real(REAL_TYPE a, REAL_TYPE b, REAL_TYPE* c, REAL_TYPE* s)
//	Вращения Гивенса
//	function [c,s] = giv(a,b)
//		if abs(b) == 0
//			c = 1;
//			s = 0;
//		else
//			if(abs(b) > abs(a))
//				tau = -a/b;
//				s = 1/sqrt(1+tau*tau');
//				c = s*tau;
//			else
//				tau = -b/a;
//				c = 1/sqrt(1+tau*tau');
//				s = c*tau';
//			end
//		end
//	end
{
	REAL_TYPE tau;

	if (is_zero(b) != 0)
	{
		*c = 1.0f;
		*s = 0.0f;
	}
	else
	{
		if (f_abs(b) > f_abs(a))
		{
			tau = -a / b;
			*s = 1.0f / sq_rt(1.0f + tau*tau);
			*c = *s * tau;
		}
		else
		{
			tau = -b / a;
			*c = 1.0f / sq_rt(1.0f + tau*tau);
			*s = *c * tau;
		}
	}
}

void qr_step_symm_real(const Matrix_t* T, Matrix_t* Q, Matrix_t* R)
// Неявный симметричный QR - шаг неприводимой
// трехдиагональной матрицы со сдвигом Уилкинсона
//function [Q, T] = qr_step(T)
//	n = size(T,2);
//	if n == 1
//		Q = 1;
//		return;
//	end
//	d = (T(n-1,n-1) - T(n,n))/2;
//
//	mu = T(n,n) - T(n,n-1)*T(n,n-1)'/(d + d/norm(d)*sqrt(d*d' + T(n,n-1)*T(n,n-1)'));
//
//	x = T(1,1) - mu;
//	z = T(2,1);
//	Q = eye(n,n,'single');
//	for k = 1:n-1
//		[c, s] = giv(x,z);
//		G = eye(n,n,'single');
//		G(k:k+1,k:k+1) = [c s; -s' c];
//		T = G'*T*G;
//		Q = Q*G;
//
//		if k < n-1
//			x = T(k+1,k);
//			z = T(k+2,k);
//		end
//	end
//end
{
	Matrix_t Tcopy, G, Gt, TG, Qcopy;
	REAL_TYPE d, mu, x, z, c, s;
	REAL_TYPE f1;
	int k, create = 1;

	//	n = size(T,2);
	int n = T->numCols;

#ifdef PREALLOCATION
	static REAL_TYPE buffer[MATRIX_MAX_SIZE * 5];
	Tcopy.pData = buffer;
	G.pData = buffer + MATRIX_MAX_SIZE;
	Gt.pData = buffer + 2 * MATRIX_MAX_SIZE;
	TG.pData = buffer + 3 * MATRIX_MAX_SIZE;
	Qcopy.pData = buffer + 4 * MATRIX_MAX_SIZE;
#endif // PREALLOCATION

	//	if n == 1
	//		Q = 1;
	//		return;
	//	end
	if(n == 1)
	{
		*Q->pData = 1.0f;
		return;
	}

	real_clone(T, &Tcopy, ALLOCATE_MATRIX);
	real_new(n, n, &TG);

	//	d = (T(n-1,n-1) - T(n,n))/2;
	d = (Tcopy.pData[(n-2)*n + n-2] - Tcopy.pData[(n-1)*n + n-1]) / 2.0f;

	//	mu = T(n,n) - T(n,n-1)*T(n,n-1)'/(d + d/norm(d)*sqrt(d*d' + T(n,n-1)*T(n,n-1)'));
	f1 = Tcopy.pData[n*n - 2];
	f1 = f1 * f1;
	mu = Tcopy.pData[n*n - 1] - f1/(d + d/f_abs(d)*sq_rt(d*d + f1));

	//	x = T(1,1) - mu;
	x = *Tcopy.pData - mu;
	//	z = T(2,1);
	z = Tcopy.pData[n];
	//	Q = eye(n,n,'single');
	real_eye(n, Q, 0);

	//	for k = 1:n-1
	for (k = 0; k < n-1; ++k)
	{
		//		[c, s] = giv(x,z);
		givens_rotation_real(x, z, &c, &s);
		//		G = eye(n,n,'single');
		real_eye(n, &G, create);
		//		G(k:k+1,k:k+1) = [c s; -s' c];
		G.pData[k*n + k] = c;
		G.pData[k*n + k+1] = s;
		G.pData[(k+1)*n + k] = -s;
		G.pData[(k+1)*n + k+1] = c;

		//		T = G'*T*G;
		real_transp(&G, &Gt, create);
		real_matrix_mult(&Tcopy, &G, 1.0f, &TG);
		real_matrix_mult(&Gt, &TG, 1.0f, &Tcopy);

		//		Q = Q*G;
		real_clone(Q, &Qcopy, create);
		real_matrix_mult(&Qcopy, &G, 1.0f, Q);

		//		if k < n-1
		//			x = T(k+1,k);
		//			z = T(k+2,k);
		//		end
		if (k < n-2)
		{
			x = Tcopy.pData[(k+1)*n + k];
			z = Tcopy.pData[(k+2)*n + k];
		}

		create = 0;
	}

	real_clone(&Tcopy, R, 0);

	real_free(&Tcopy);
	real_free(&G);
	real_free(&Gt);
	real_free(&TG);
	real_free(&Qcopy);
}

void qr_symm_real(const Matrix_t* T, REAL_TYPE tol, Matrix_t* U, Matrix_t* D)
//////////////////////////////////////////////////////////////////////
// Симметричный QR-алгоритм
// Источник: Голуб Дж., Ван Лоун Ч. Матричные вычисления: Пер. с англ. -
//              М.: Мир, 1999. - 548 с. (стр. 380 )
//////////////////////////////////////////////////////////////////////
//function [U,D] = QRsymm(T, tol)
//	n = size(T,1);
//	U = eye(n,n,'single');
//
//	p = 1;
//	q = n;
//	while q ~= 1
//		fl = 0;
//		for i = q-1:-1:1
//			if (abs(T(i+1,i)) <= tol*(abs(T(i,i)) + abs(T(i+1,i+1))))
//				T(i,i+1) = 0;
//				T(i+1,i) = 0;
//
//				if fl == 0
//					q = i;
//				else
//					break;
//				end
//			else
//				fl = 1;
//				p = i;
//			end
//		end
//
//		if q > 1
//			[Qm,Rm] = qr_step(T(p:q,p:q));
//			T(p:q,p:q) = Rm;
//
//			Q = eye(n, n);
//			Q(p:q,p:q) = Qm;
//			U = U*Q;
//		end
//	end
//	D = T;
//end
{
	int n = T->numRows, i, size, create = 1, free_q = 0;
	//	p = 1;
	//	q = n;
	int p = 1, q = n, fl;
	Matrix_t Tcopy, Tpart, Qm, Rm, Q, Ucopy;

#ifdef PREALLOCATION
	static REAL_TYPE buffer[MATRIX_MAX_SIZE * 6];
	Tcopy.pData = buffer;
	Tpart.pData = buffer + MATRIX_MAX_SIZE;
	Qm.pData = buffer + 2 * MATRIX_MAX_SIZE;
	Rm.pData = buffer + 3 * MATRIX_MAX_SIZE;
	Q.pData = buffer + 4 * MATRIX_MAX_SIZE;
	Ucopy.pData = buffer + 5 * MATRIX_MAX_SIZE;
#endif // PREALLOCATION

	real_clone(T, &Tcopy, ALLOCATE_MATRIX);
	real_eye(n, U, 0);

	//	while q ~= 1
	while (q != 1)
	{
		fl = 0;
		//		for i = q-1:-1:1
		for (i = q-1; i > 0; --i) // здесь индексы Matlab
		{
			//			if (abs(T(i+1,i)) <= tol*(abs(T(i,i)) + abs(T(i+1,i+1))))
			if (f_abs(Tcopy.pData[i*n + i-1]) <=
				tol * (f_abs(Tcopy.pData[(i-1)*n + i-1]) + f_abs(Tcopy.pData[i*n + i])))
			{
				//				T(i,i+1) = 0;
				//				T(i+1,i) = 0;
				Tcopy.pData[(i-1)*n + i] = 0.0f;
				Tcopy.pData[i*n + i-1] = 0.0f;

				if (fl == 0)
					q = i;
				else
					break;
			}
			else
			{
				fl = 1;
				p = i;
			}
		}
		//		if q > 1
		if (q > 1)
		{
			//			[Qm,Rm] = qr_step(T(p:q,p:q));
			size = q-p+1;
			real_new(size, size, &Tpart);
			real_new(size, size, &Qm);
			real_new(size, size, &Rm);
			real_partial_copy(&Tcopy, p-1, p-1, size, size, &Tpart, 0, 0);
			qr_step_symm_real(&Tpart, &Qm, &Rm);

			//			T(p:q,p:q) = Rm;
			real_partial_copy(&Rm, 0, 0, size, size, &Tcopy, p-1, p-1);

			//			Q = eye(n, n);
			real_eye(n, &Q, create);
			//			Q(p:q,p:q) = Qm;
			real_partial_copy(&Qm, 0, 0, size, size, &Q, p-1, p-1);
			//			U = U*Q;
			real_clone(U, &Ucopy, create);
			real_matrix_mult(&Ucopy, &Q, 1.0f, U);

			real_free(&Tpart);
			real_free(&Qm);
			real_free(&Rm);
			create = 0;
			free_q = 1;
		}
	}

	//	D = T;
	real_clone(&Tcopy, D, 0);

	real_free(&Tcopy);
	if (free_q)
	{
		real_free(&Q);
		real_free(&Ucopy);
	}
}

void eig_symm_triag(const CMatrix_t* A, REAL_TYPE tol, Matrix_t* S, CMatrix_t* U)
// нахождение собственных значений комплексной матрицы А с точночтью tol
// S - вектор-строка получившихся собственных значений
// U - матрица собственных векторов
//function [U, S] = app_eigen(A)
//	[T, P] = triag(A);
//	[Q, S] = QRsymm(T, tol);
//	U = P'*Q; % U = (Q'*P)';
//end
{
	int n = A->numRows, i;
	Matrix_t T, Q, Ss; // Q не удалять! у неё общий буфер с Qcmpl!
	CMatrix_t P, Qcmpl, Pt;
	int indices[VECTOR_MAX_SIZE];
#ifdef TEST_RESULT
	CMatrix_t US, USUt, test, Scmpl;
#endif // TEST_RESULT

#ifdef PREALLOCATION
	static REAL_TYPE buffer[MATRIX_MAX_SIZE * (3 + 3*2 - 1)]; // матрицы Q и Qcmpl имеют общий буфер
	T.pData = buffer;
	Q.pData = buffer + MATRIX_MAX_SIZE;
	Ss.pData = buffer + 2 * MATRIX_MAX_SIZE;
	P.pDataReal = buffer + 3 * MATRIX_MAX_SIZE;
	P.pDataImag = buffer + 4 * MATRIX_MAX_SIZE;
	Qcmpl.pDataReal = Q.pData;
	Qcmpl.pDataImag = buffer + 5 * MATRIX_MAX_SIZE;
	Pt.pDataReal = buffer + 6 * MATRIX_MAX_SIZE;
	Pt.pDataImag = buffer + 7 * MATRIX_MAX_SIZE;

#ifdef TEST_RESULT
	static REAL_TYPE bufferTest[MATRIX_MAX_SIZE * 2 * 4];
	US.pDataReal = buffer;
	US.pDataImag = buffer + MATRIX_MAX_SIZE;
	USUt.pDataReal = buffer + 2 * MATRIX_MAX_SIZE;
	USUt.pDataImag = buffer + 3 * MATRIX_MAX_SIZE;
	test.pDataReal = buffer + 4 * MATRIX_MAX_SIZE;
	test.pDataImag = buffer + 5 * MATRIX_MAX_SIZE;
	Scmpl.pDataReal = buffer + 6 * MATRIX_MAX_SIZE;
	Scmpl.pDataImag = buffer + 7 * MATRIX_MAX_SIZE;
#endif // TEST_RESULT
#endif // PREALLOCATION

	real_new(n, n, &T);
	real_new(n, n, &Q);
	real_new(n, n, &Ss);
	complex_new(n, n, &P);

	complex_triag(A, &T, &P);

	qr_symm_real(&T, tol, &Q, &Ss);

	// S = diag(Ss);
	for (i = 0; i < n; ++i)
		S->pData[i] = Ss.pData[i*n + i];

	Qcmpl.numCols = n;
	Qcmpl.numRows = n;
#ifndef PREALLOCATION
	Qcmpl.pDataReal = Q.pData;
	Qcmpl.pDataImag = (REAL_TYPE*)malloc(n * n * sizeof(REAL_TYPE));
	Q.pData = NULL; // на всякий случай, если кто-то попытается освободить Q
#endif
	memset(Qcmpl.pDataImag, 0, n * n * sizeof(REAL_TYPE));
	complex_transp(&P, &Pt, 1);

	complex_matrix_mult(&Pt, &Qcmpl, 1.0f, 0.0f, U);

	special_sort_descending(S->pData, S->numCols, indices);
	complex_swap_columns(U, indices);

#ifdef TEST_RESULT
	complex_new(n, n, &US);
	complex_new(n, n, &USUt);
	complex_new(n, n, &test);
	complex_zeros(n, n, &Scmpl, 1);

	for (i = 0; i < S->numCols; ++i)
		Scmpl.pDataReal[i*S->numCols + i] = S->pData[i];

	complex_matrix_mult(U, &Scmpl, 1.0f, 0.0f, &US);
	complex_matrix_mult_right_transp(&US, U, 1.0, 0.0, &USUt);
	complex_matrix_sum(&USUt, A, -1.0f, 0.0f, &test);

	complex_free(&US);
	complex_free(&USUt);
	complex_free(&test);
	complex_free(&Scmpl);
#endif // TEST_RESULT

	real_free(&T);
	complex_free(&Qcmpl);
	complex_free(&Pt);
	real_free(&Ss);
	complex_free(&P);
}

 /*
REAL_TYPE eig_symm_triag_only_one(const CMatrix_t* A, REAL_TYPE tol, CMatrix_t* U)
// нахождение собственных значений комплексной матрицы А с точночтью tol
// S - вектор-строка получившихся собственных значений
// U - матрица собственных векторов
//function [U, S] = app_eigen(A)
//	[T, P] = triag(A);
//	[Q, S] = QRsymm(T, tol);
//	U = (Q'*P)';
//end
{
	int n = A->numRows, i;
	Matrix_t T, Q, Ss, S; // Q не удалять! у неё общий буфер с Qcmpl!
	CMatrix_t P, Qcmpl, Qt, Ut;
	REAL_TYPE maxEig;
	size_t maxPos;

#ifdef PREALLOCATION
	static REAL_TYPE buffer[MATRIX_MAX_SIZE * (3 + 4 * 2 - 1) + VECTOR_MAX_SIZE]; // матрицы Q и Qcmpl имеют общий буфер
	T.pData = buffer;
	Q.pData = buffer + MATRIX_MAX_SIZE;
	Ss.pData = buffer + 2 * MATRIX_MAX_SIZE;
	P.pDataReal = buffer + 3 * MATRIX_MAX_SIZE;
	P.pDataImag = buffer + 4 * MATRIX_MAX_SIZE;
	Qcmpl.pDataReal = Q.pData;
	Qcmpl.pDataImag = buffer + 5 * MATRIX_MAX_SIZE;
	Qt.pDataReal = buffer + 6 * MATRIX_MAX_SIZE;
	Qt.pDataImag = buffer + 7 * MATRIX_MAX_SIZE;
	Ut.pDataReal = buffer + 8 * MATRIX_MAX_SIZE;
	Ut.pDataImag = buffer + 9 * MATRIX_MAX_SIZE;
	S.pData = buffer + 10 * MATRIX_MAX_SIZE;
#endif // PREALLOCATION

	real_new(n, n, &T);
	real_new(n, n, &Q);
	real_new(n, n, &Ss);
	real_new(n, 1, &S);
	complex_new(n, n, &P);

	complex_triag(A, &T, &P);

	qr_symm_real(&T, tol, &Q, &Ss);
	// S = diag(Ss);
	for (i = 0; i < n; ++i)
		S.pData[i] = Ss.pData[i*n + i];

	Qcmpl.numCols = n;
	Qcmpl.numRows = n;
#ifndef PREALLOCATION
	Qcmpl.pDataReal = Q.pData;
	Qcmpl.pDataImag = (REAL_TYPE*)malloc(n * n * sizeof(REAL_TYPE));
	Q.pData = NULL; // на всякий случай, если кто-то попытается освободить Q
#endif
	memset(Qcmpl.pDataImag, 0, n * n * sizeof(REAL_TYPE));
	complex_transp(&Qcmpl, &Qt, 1);

	complex_new(n, n, &Ut);
	complex_matrix_mult(&Qt, &P, 1.0f, 0.0f, &Ut);
	complex_transp(&Ut, U, 0);

	maxPos = find_max(S.pData, S.numRows, &maxEig);
	complex_partial_copy(U, 0, maxPos, U->numCols, 1, U, 0, 0);

	real_free(&T);
	complex_free(&Qcmpl);
	complex_free(&Qt);
	complex_free(&Ut);
	real_free(&Ss);
	complex_free(&P);
	return maxEig;
}
*/
