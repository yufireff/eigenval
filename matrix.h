/*
 * Matrix.h
 *
 *  Created on: 09.01.2017
 *      Author: Ivan
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "math.h"
#include "stdint.h"
#include "settings.h"

//#define MATRIX_MAX_DIM	4
#define MATRIX_SUCCESS	(0)
#define MATRIX_EXCEEDS	(-1)
//#define MATRIX_DIM  4

#define SUM_ARIF(a1, an, n) (((a1) + (an))/2 * n) ///< Сумма арифметической прогрессии


#define fsign(x)	(((x) > 0.0f) ? 1.0f : (((x) < 0.0f) ? -1.0f : 0.0f))

typedef enum
{
	LeftMult,
	RightMult,
	LeftRightMult
}MatrixFlag_t;

typedef enum
{
	Common,	// обычная матрица
	Ones,	// матрица, состоящая из единиц
	Eye,	// единичная матрица
	Tri,	//  трехдиагональная матрица
	Vector, //вектор-столбец
	VectorT //вектор-строка
}MatrixType_t;

typedef struct
{
	int32_t numRows;
	int32_t numCols;
	REAL_TYPE * pData;
} Matrix_t;

typedef struct
{
	int32_t numRows;
	int32_t numCols;
	REAL_TYPE * pDataReal;
	REAL_TYPE * pDataImag;
} CMatrix_t;

typedef struct
{
	int32_t numRows;
	int32_t numCols;
	double * pDataReal;
	double * pDataImag;
} CMatrix_td;

#ifndef DSP_OPTIMIZATION_FULL
void dot_prod(REAL_TYPE * a, REAL_TYPE * b,	REAL_TYPE * c, uint32_t dira, uint32_t dirb, uint32_t len);

void cdot_prod(REAL_TYPE * a_real, REAL_TYPE * a_imag, REAL_TYPE * b_real, REAL_TYPE * b_imag,
		REAL_TYPE * c_real, REAL_TYPE * c_imag, uint32_t dira, uint32_t dirb, uint32_t len);

void cnorm2(REAL_TYPE * a_real, REAL_TYPE * a_imag, REAL_TYPE * norm, uint32_t dir, uint32_t len);

int complex_matrix_mult(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res);

int complex_transp(const CMatrix_t* a, CMatrix_t* b, int create);

int inv(const CMatrix_t* a, CMatrix_t* ainv);

void complex_zeros(int rows, int columns, CMatrix_t* a, int create);

void complex_eye(int n, CMatrix_t* a, int create);
#endif // !DSP_OPTIMIZATION_FULL

void complex_free(CMatrix_t* a);

#ifndef DSP_OPTIMIZATION_FULL
void complex_clone(const CMatrix_t* source, CMatrix_t* destination, int create);
#endif // !DSP_OPTIMIZATION_FULL

int complex_new(int rows, int columns, CMatrix_t* a);

#ifndef DSP_OPTIMIZATION_FULL
int complex_partial_copy(const CMatrix_t* source, int source_start_row, int source_start_column,
	int rows_to_copy, int columns_to_copy,
	CMatrix_t* destination, int destination_start_row, int destination_start_column);

int complex_matrix_sum(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res); // res = a + factor * b

int complex_scal_prod(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE* re, REAL_TYPE* im);
#endif // !DSP_OPTIMIZATION_FULL

// -----------
int real_new(int rows, int columns, Matrix_t* a);

void real_free(Matrix_t* a);

#ifndef DSP_OPTIMIZATION_FULL
void real_zeros(int rows, int columns, Matrix_t* a, int create);

void real_eye(int n, Matrix_t* a, int create);

void real_clone(const Matrix_t* source, Matrix_t* destination, int create);

int real_partial_copy(const Matrix_t* source, int source_start_row, int source_start_column,
	int rows_to_copy, int columns_to_copy,
	Matrix_t* destination, int destination_start_row, int destination_start_column);


int real_transp(const Matrix_t* a, Matrix_t* b, int create);

int real_matrix_mult(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res);

int real_to_complex(const Matrix_t* real, CMatrix_t* compl, int create);
#endif // ! DSP_OPTIMIZATION_FULL

void complex_swap_columns(CMatrix_t* m, int* new_indices);

// специальные функции для работы с матрицами double

void single_to_double(const CMatrix_t * S, CMatrix_td *D, int create);

int complex_new_d(int rows, int columns, CMatrix_td* a);

void complex_eye_d(int n, CMatrix_td* a, int create);

void complex_clone_d(const CMatrix_td* source, CMatrix_td* destination, int create);

void complex_zeros_d(int rows, int columns, CMatrix_td* a, int create);

int complex_partial_copy_d(const CMatrix_td* source, int source_start_row, int source_start_column,
	int rows_to_copy, int columns_to_copy,
	CMatrix_td* destination, int destination_start_row, int destination_start_column);

void cnorm2_d(double * a_real, double * a_imag, double * norm2, uint32_t dir, uint32_t len);

int complex_matrix_sum_d(const CMatrix_td* a, const CMatrix_td* b, double factor_re, double factor_im, CMatrix_td* res); // res = a + factor * b

#ifndef DSP_OPTIMIZATION_FULL
int complex_transp_d(const CMatrix_td* a, CMatrix_td* b, int create);

int complex_scal_prod_d(const CMatrix_td* a, const CMatrix_td* b, double* re, double* im);
#endif // !DSP_OPTIMIZATION_FULL

int complex_matrix_mult_d(const CMatrix_td* a, const CMatrix_td* b, double factor_re, double factor_im, CMatrix_td* res);

void double_complex_to_single_real(const CMatrix_td* a, Matrix_t* b, int create);

void double_complex_to_single_complex(const CMatrix_td* a, CMatrix_t* b, int create);

void complex_free_d(CMatrix_td* a);

#endif /* MATRIX_H_ */
