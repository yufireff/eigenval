#include "eig.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"
#include "real.h"
#include "settings.h"
#include "algorithms.h"
#include "MatrixSpecial.h"
#include "MatrixDsp.h"
#include "dsp.h"
#include "EigDsp.h"

#ifdef TRIAG_DOUBLE
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
	//         m         m      m    m     m       m   m    m      v  v  v  v  v  m   m
	CMatrix_td TComplex, Tcopy, qut, uqt, qut_uqt, Pk, uut, Pcopy, x, e, u, p, q, Ad, Pd;
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
	double normx, wre, wim, Kre, Kim;

	//		m = size(A,1);
	int m = A->numRows;

#ifdef PREALLOCATION
#ifdef DOUBLE
	static double buffer[2 *(MATRIX_MAX_SIZE*8 + VECTOR_MAX_SIZE*5)];
#else // DOUBLE
	static double buffer[2 * (MATRIX_MAX_SIZE * 10 + VECTOR_MAX_SIZE * 5)];
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
	p.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 6 * VECTOR_MAX_SIZE;
	p.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 7 * VECTOR_MAX_SIZE;
	q.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 8 * VECTOR_MAX_SIZE;
	q.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 9 * VECTOR_MAX_SIZE;
#ifndef DOUBLE
	Ad.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 10 * VECTOR_MAX_SIZE;
	Ad.pDataImag = buffer + 17 * MATRIX_MAX_SIZE + 10 * VECTOR_MAX_SIZE;
	Pd.pDataReal = buffer + 18 * MATRIX_MAX_SIZE + 10 * VECTOR_MAX_SIZE;
	Pd.pDataImag = buffer + 19 * MATRIX_MAX_SIZE + 10 * VECTOR_MAX_SIZE;
#else // DOUBLE
	Pd.pDataReal = P->pDataReal;
	Pd.pDataImag = P->pDataImag;
	Pd.numCols = P->numCols;
	Pd.numRows = P->numRows;
#endif // DOUBLE
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
		complex_scal_prod_left_transp_d(&x, &u, &f5, &f6);
		complex_scal_prod_left_transp_d(&u, &x, &f7, &f8);
		complex_div_d(f5, f6, f7, f8, &f9, &f10);
		wre = f4 * (1.0 + f9);
		wim = f4 * f10;

		//			p = w' * T * u;
		complex_matrix_mult_d(&TComplex, &u, wre, -wim, &p);
		//			K = w/2 * u'*p;
		complex_scal_prod_left_transp_d(&u, &p, &f11, &f12);
		complex_mult_d(wre/2.0f, wim/2.0f, f11, f12, &Kre, &Kim);
		//			q = p - K*u;
		complex_matrix_sum_d(&p, &u, -Kre, -Kim, &q);

		//			T = T - q*u' - u*q';
		complex_matrix_mult_right_transp_d(&q, &u, &qut);
		complex_matrix_mult_right_transp_d(&u, &q, &uqt);
		complex_matrix_sum_d(&qut, &uqt, 1.0f, 0.0f, &qut_uqt);
		complex_clone_d(&TComplex, &Tcopy, create);
		complex_matrix_sum_d(&TComplex, &qut_uqt, -1.0f, 0.0f, &TComplex);

		//			Pk = eye(m,m,'like',A);
		complex_eye_d(m, &Pk, create);
		//			Pk = Pk - w*u*u';
		complex_matrix_mult_right_transp_d(&u, &u, &uut);
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

	for (i = 0; i < T->numCols * T->numRows; ++i)
		if (f_abs(T->pData[i]) < 1e-6f)
			T->pData[i] = 0;

	complex_free_d(&TComplex);
	complex_free_d(&x);
	complex_free_d(&e);
	complex_free_d(&u);
	complex_free_d(&p);
	complex_free_d(&q);
	complex_free_d(&Tcopy);
	complex_free_d(&qut);
	complex_free_d(&uqt);
	complex_free_d(&qut_uqt);
	complex_free_d(&Pk);
	complex_free_d(&uut);
	complex_free_d(&Pcopy);
}
#else // TRIAG_DOUBLE
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
	//         m         m      m    m     m       m   m    m      v  v  v  v  v
	CMatrix_t TComplex, Tcopy, qut, uqt, qut_uqt, Pk, uut, Pcopy, x, e, u, p, q;
	float f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
	float normx, wre, wim, Kre, Kim;

	//		m = size(A,1);
	int m = A->numRows;

#ifdef PREALLOCATION
	static float buffer[2 * (MATRIX_MAX_SIZE * 8 + VECTOR_MAX_SIZE * 5)];
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
	p.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 6 * VECTOR_MAX_SIZE;
	p.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 7 * VECTOR_MAX_SIZE;
	q.pDataReal = buffer + 16 * MATRIX_MAX_SIZE + 8 * VECTOR_MAX_SIZE;
	q.pDataImag = buffer + 16 * MATRIX_MAX_SIZE + 9 * VECTOR_MAX_SIZE;
#endif // PREALLOCATION

	//		P = eye(m, m, 'like', A);
	complex_eye(m, P, ALLOCATE_MATRIX);
	//		T = A;
	complex_clone(A, &TComplex, ALLOCATE_MATRIX);

	//		for i = 2:m
	for (i = 1; i <= m - 1; ++i)
	{
		if (create > 0)
		{
			complex_new(m, 1, &u);
			complex_new(m, 1, &p);
			complex_new(m, 1, &q);
			complex_new(m, m, &qut);
			complex_new(m, m, &uqt);
			complex_new(m, m, &qut_uqt);
			complex_new(m, m, &uut);
		}

		//			x = zeros(m,1,'like',A);
		complex_zeros(m, 1, &x, create);
		//			x(i:m) = T(i:m,i-1);
		complex_partial_copy(&TComplex, i, i - 1, m - i, 1, &x, i, 0);
		//			e = zeros(m,1,'like',A);
		complex_zeros(m, 1, &e, create);
		//			e(i) = 1;
		e.pDataReal[i] = 1.0f;

		//			u = x - csgn(x(i)) * norm(x) * e;
		f1 = complex_sign(x.pDataReal[i], x.pDataImag[i]);
		cnorm2(x.pDataReal, x.pDataImag, &f2, 1, m);
		normx = sqrtf(f2);
		complex_matrix_sum(&x, &e, -f1*normx, 0.0f, &u);

		//			w = 1/(u'*u)*(1 + x'*u/(u'*x));
		cnorm2(u.pDataReal, u.pDataImag, &f3, 1, m);
		f4 = 1.0f / f3;
		complex_scal_prod(&x, &u, &f5, &f6);
		complex_scal_prod(&u, &x, &f7, &f8);
		complex_div(f5, f6, f7, f8, &f9, &f10);
		wre = f4 * (1.0f + f9);
		wim = f4 * f10;

		//			p = w' * T * u;
		complex_matrix_mult_dsp(&TComplex, &u, wre, -wim, &p, 0);
		//			K = w/2 * u'*p;
		complex_scal_prod(&u, &p, &f11, &f12);
		complex_mult(wre / 2.0f, wim / 2.0f, f11, f12, &Kre, &Kim);
		//			q = p - K*u;
		complex_matrix_sum(&p, &u, -Kre, -Kim, &q);

		//			T = T - q*u' - u*q';
		complex_matrix_mult_right_transp_dsp(&q, &u, &qut, 0);
		complex_matrix_mult_right_transp_dsp(&u, &q, &uqt, 0);
		complex_matrix_sum(&qut, &uqt, 1.0f, 0.0f, &qut_uqt);
		complex_clone(&TComplex, &Tcopy, create);
		complex_matrix_sum(&TComplex, &qut_uqt, -1.0f, 0.0f, &TComplex);

		//			Pk = eye(m,m,'like',A);
		complex_eye(m, &Pk, create);
		//			Pk = Pk - w*u*u';
		complex_matrix_mult_right_transp_dsp(&u, &u, &uut, 0);
		complex_matrix_sum(&Pk, &uut, -wre, -wim, &Pk);

		//			P = Pk * P;
		//complex_clone(P, &Pcopy, create);
		complex_matrix_mult_dsp(&Pk, P, 1.0f, 0.0f, P, 0);

		create = 0;
	}
	//		T = real(T);
	memcpy(T->pData, TComplex.pDataReal, m*m*sizeof(REAL_TYPE));

	complex_free(&TComplex);
	complex_free(&x);
	complex_free(&e);
	complex_free(&u);
	complex_free(&p);
	complex_free(&q);
	complex_free(&Tcopy);
	complex_free(&qut);
	complex_free(&uqt);
	complex_free(&qut_uqt);
	complex_free(&Pk);
	complex_free(&uut);
	complex_free(&Pcopy);
}
#endif // TRIAG_DOUBLE

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
			float t0 = 1 / b;
			tau = -a / b;
			float t1 = tau*tau;
			float t2 = 1.0f + tau*tau;
			*s = 1.0f / sq_rt(1.0f + tau*tau);
			*c = *s * tau;
		}
		else
		{
			float t0 = 1 / a;
			tau = -b / a;
			float t1 = tau*tau;
			float t2 = 1.0f + tau*tau;
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
	Matrix_t Tcopy, G, Qcopy;
	REAL_TYPE d, mu, x, z, c, s;
	REAL_TYPE f1;
	int k, create = 1;

	//	n = size(T,2);
	int n = T->numCols;

#ifdef PREALLOCATION
	static REAL_TYPE buffer[MATRIX_MAX_SIZE * 3];
	Tcopy.pData = buffer;
	G.pData = buffer + MATRIX_MAX_SIZE;
	Qcopy.pData = buffer + 2 * MATRIX_MAX_SIZE;
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

	//	d = (T(n-1,n-1) - T(n,n))/2;
	d = (Tcopy.pData[(n-2)*n + n-2] - Tcopy.pData[(n-1)*n + n-1]) / 2.0f;

	//	mu = T(n,n) - T(n,n-1)*T(n,n-1)'/(d + d/norm(d)*sqrt(d*d' + T(n,n-1)*T(n,n-1)'));
	f1 = Tcopy.pData[n*n - 2];
	f1 = f1 * f1;
	float t0 = sq_rt(d*d + f1);
	float t1 = d + d / f_abs(d)*sq_rt(d*d + f1);
	float t2 = 1.0 / (d + d / f_abs(d)*sq_rt(d*d + f1));
	float t3 = f1 / (d + d / f_abs(d)*sq_rt(d*d + f1));
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
		real_matrix_mult_at_b_a_dsp(&G, &Tcopy, &Tcopy, 0);

		//		Q = Q*G;
		real_clone(Q, &Qcopy, ALLOCATE_MATRIX);
		real_matrix_mult_dsp(&Qcopy, &G, 1.0f, Q, 0);

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
	real_free(&Qcopy);
	real_free(&G);
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
			real_clone(U, &Ucopy, ALLOCATE_MATRIX);
			real_matrix_mult_dsp(&Ucopy, &Q, 1.0f, U, 0);

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
	real_free(&Ucopy);
	if (free_q)
		real_free(&Q);
}

int eig_symm_triag(const CMatrix_t* A, REAL_TYPE tol, Matrix_t* S, CMatrix_t* U)
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
    unsigned int start, step, prepare, triag, mult, sort, qr, check;
	Matrix_t T, Q, Ss; // Q не удалять! у неё общий буфер с Qcmpl!
	CMatrix_t P, Qcmpl;
	int indices[VECTOR_MAX_SIZE];
	int result = MATRIX_SUCCESS;
#ifdef TEST_RESULT
	CMatrix_t test;
	float norm = 0.0f;
#endif // TEST_RESULT

#ifdef PREALLOCATION
	static REAL_TYPE buffer[MATRIX_MAX_SIZE * (3 + 2*2 - 1)]; // матрицы Q и Qcmpl имеют общий буфер
	T.pData = buffer;
	Q.pData = buffer + MATRIX_MAX_SIZE;
	Ss.pData = buffer + 2 * MATRIX_MAX_SIZE;
	P.pDataReal = buffer + 3 * MATRIX_MAX_SIZE;
	P.pDataImag = buffer + 4 * MATRIX_MAX_SIZE;
	Qcmpl.pDataReal = Q.pData;
	Qcmpl.pDataImag = buffer + 5 * MATRIX_MAX_SIZE;

#ifdef TEST_RESULT
	static REAL_TYPE bufferTest[MATRIX_MAX_SIZE * 2];
	test.pDataReal = bufferTest;
	test.pDataImag = bufferTest + MATRIX_MAX_SIZE;
#endif // TEST_RESULT
#endif // PREALLOCATION

	start =  GetCP0_Count();

	real_new(n, n, &T);
	real_new(n, n, &Q);
	real_new(n, n, &Ss);
	complex_new(n, n, &P);

    step =  GetCP0_Count();
    prepare = step - start;

	complex_triag(A, &T, &P);

    step =  GetCP0_Count();
    triag = step - prepare;

	qr_symm_real_dsp(&T, tol, &Q, &Ss, 0);

    step =  GetCP0_Count();
    qr = step - triag;

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

	complex_matrix_mult_left_transp_dsp(&P, &Qcmpl, U, 0);

    step =  GetCP0_Count();
    mult = step - qr;

	special_sort_descending(S->pData, S->numCols, indices);
	complex_swap_columns(U, indices);

    step =  GetCP0_Count();
    sort = step - mult;

#ifdef TEST_RESULT
	complex_new(n, n, &test);
	complex_a_minus_u_s_ut_dsp(A, U, S, &test, &norm, 0);
	if (norm > tol)
		result = MATRIX_EXCEEDS;
	complex_free(&test);
#endif // TEST_RESULT

    step =  GetCP0_Count();
    check = step - sort;

	real_free(&T);
	complex_free(&Qcmpl);
	real_free(&Ss);
	complex_free(&P);

	return result;
}
