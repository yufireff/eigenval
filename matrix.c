/**
 * \file Matrix.c
 *
 * Формат хранения матрицы:
 * 	последовательно друг за другом в памяти располагаются
 * 	numRows элементов строки матрицы, далее начинается новая строка.
 * 	Количество таких строк numCols.
 * 	Пример:
 * 	A(i,j) = A->pData[i*numCols + j]
 *
 * Функции для работы с матрицами:
 * create_matrix - динамическое выделение памяти для матрицы;
 * dot_prod2 - скалярное произведение векторов;
 * givens - расчет элементов матрицы вращений Гивенса;
 * givensApply - применение матрицы вращений Гивенса;
 * qr_step - Неявный симметричный QR - шаг неприводимой
 * 			трехдиагональной матрицы со сдвигом Уилкинсона;
 * house - расчет вектора Хаусхолдера;
 * houseApply - применение вектора Хаусхолдера;
 * triag - трехдиагонализация Хаусхолдера;
 * qr_symm - симметричный QR-алгоритм.
 * */

#include "settings.h"
#include "matrix.h"
#ifndef PREALLOCATION
#include <stdlib.h>
#endif
#include <string.h>
#include "complex.h"
#include "real.h"
#include "algorithms.h"
#ifdef DSP_OPTIMIZATION
#include "dsp.h"
#endif // DSP_OPTIMIZATION

/**
 * \brief Скалярное произведение
 *
 */
void dot_prod(REAL_TYPE * a, REAL_TYPE * b,	REAL_TYPE * c, uint32_t dira, uint32_t dirb, uint32_t len)
{
	uint32_t i;
	REAL_TYPE * pa, * pb;

	pa = a; pb = b;

	for(i = 0; i < len; i++)
	{
		(*c) = (*c) + (*pa) * (*pb);
		pa += dira;
		pb += dirb;
	}
}

/**
 *  \brief Скалярное произведение векторов на DSP
 *
 * a, b - векторы;
 * dira, dirb - тип вектора (столбец, строка). Определяет на скольно необходимо
 * сместить указатель, чтобы получить следующий элемент;
 * len - длина векторов;
 * c - результат.
 */

void cdot_prod(REAL_TYPE * a_real, REAL_TYPE * a_imag, REAL_TYPE * b_real, REAL_TYPE * b_imag,
		REAL_TYPE * c_real, REAL_TYPE * c_imag, uint32_t dira, uint32_t dirb, uint32_t len)
{
	uint32_t i;

	for(i = 0; i < len; i++)
	{
		(*c_real) += ((*a_real) * (*b_real) + (*a_imag) * (*b_imag));
		(*c_imag) += ((*b_real) * (*a_imag) - (*a_real) * (*b_imag));

		a_real += dira; a_imag += dira;
		b_real += dirb; b_imag += dirb;
	}
}

/**
 * \brief Квадрат нормы комплексного вектора
 *
 * Функция вычисляет скалярное произведение комплексного
 * вектора на самого себя.
 *
 * \param a_real - указатель на действительную часть первлого элемента вектора
 * \param a_imag - указатель на мнимую часть первлого элемента вектора
 * \param norm2 - квадрат нормы
 * \param dir - инкремент адреса
 * \param len - длина вектора
 */
void cnorm2(REAL_TYPE * a_real, REAL_TYPE * a_imag, REAL_TYPE * norm2, uint32_t dir, uint32_t len)
{
	uint32_t i;
	REAL_TYPE *pa_real, *pa_imag;
	pa_real = a_real;
	pa_imag = a_imag;

	*norm2 = 0.0;
	for(i = 0; i < len; i++)
	{
		(*norm2) = (*norm2) + (*pa_real) * (*pa_real) + (*pa_imag) * (*pa_imag);
		pa_real += dir;
		pa_imag += dir;
	}
}


int complex_matrix_mult(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res)
{
	int i, j, k, i1, j1;

	REAL_TYPE real, img;
	REAL_TYPE x, y;

	if (a->numCols != b->numRows || res->numRows != a->numRows || res->numCols != b->numCols)
		return MATRIX_EXCEEDS;

	for (i  = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < b->numCols; ++j)
		{
			real = 0.0;
			img = 0.0;
			for (k = 0; k < a->numCols; ++k)
			{
				i1 = i*a->numCols + k;
				j1 = k*b->numCols + j;
				complex_mult(a->pDataReal[i1], a->pDataImag[i1],
					b->pDataReal[j1], b->pDataImag[j1], &x, &y);
				real += x;
				img += y;
			}
			k = i*res->numCols + j;
			complex_mult(real, img, factor_re, factor_im, res->pDataReal + k, res->pDataImag + k);
		}
	}

	return MATRIX_SUCCESS;
}

int complex_transp(const CMatrix_t* a, CMatrix_t* b, int create)
{
	int i, j, idx, idxt;

	if (create > 0)
		complex_new(a->numCols, a->numRows, b);
	else if (b->numCols != a->numRows || b->numRows != a->numCols)
	{
		b->numCols = a->numRows;
		b->numRows = a->numCols;
#ifndef PREALLOCATION
		b->pDataReal = realloc(b->pDataReal, a->numCols*a->numRows*sizeof(REAL_TYPE));
		b->pDataImag = realloc(b->pDataImag, a->numCols*a->numRows*sizeof(REAL_TYPE));
#endif
	}

	if (a->numRows != b->numCols || a->numCols != b->numRows)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < a->numCols; ++j)
		{
			idx = i * a->numCols + j;
			idxt = j * a->numCols + i;
			b->pDataReal[idx] = a->pDataReal[idxt];
			b->pDataImag[idx] = -1.0f * a->pDataImag[idxt];
		}
	}
	return MATRIX_SUCCESS;
}

int inv(const CMatrix_t* a, CMatrix_t* ainv)
// нахождение обратной матрицы методом элементарных преобразований
{
	int n = a->numCols;
	int i, j, k;
	REAL_TYPE factorReal, factorImag, real, imag;
	CMatrix_t t;

#ifdef PREALLOCATION
	REAL_TYPE buffer[MATRIX_MAX_SIZE * 2 * 2]; // одна комплексная матрица удвоенного размера
	t.pDataReal = buffer;
	t.pDataImag = buffer + MATRIX_MAX_SIZE*2;
#endif // PREALLOCATION

	if (a->numCols != a->numRows || ainv->numCols != a->numCols || ainv->numRows != a->numRows)
		return MATRIX_EXCEEDS;

	// создаём и заполняем удвоенную матрицу, над которой будем проводить элементарные преобразования
	complex_new(n, 2 * n, &t);
	for (i = 0; i < n; ++i)
	{
		memcpy(t.pDataReal + i * 2 * n, a->pDataReal + i * n, n * sizeof(REAL_TYPE));
		memcpy(t.pDataImag + i * 2 * n, a->pDataImag + i * n, n * sizeof(REAL_TYPE));
		memset(t.pDataReal + i * 2 * n + n , 0, n * sizeof(REAL_TYPE));
		memset(t.pDataImag + i * 2 * n + n, 0, n * sizeof(REAL_TYPE));
		t.pDataReal[i * 2 * n + n + i] = 1.0;
	}

	// приведение к верхнетреугольному виду
	for (i = 0; i < n; ++i)
	{
		factorReal = t.pDataReal[i * 2 * n + i];
		factorImag = t.pDataImag[i * 2 * n + i];

		// делаем диагональный элемент равным 1
		for (k = 0; k < 2 * n; ++k)
		{
			complex_div(t.pDataReal[i * 2 * n + k], t.pDataImag[i * 2 * n + k], factorReal, factorImag, &real, &imag);
			t.pDataReal[i * 2 * n + k] = real;
			t.pDataImag[i * 2 * n + k] = imag;
		}

		if (i == n - 1)
			continue;

		for (j = i + 1; j < n; ++j)
		{
			factorReal = t.pDataReal[j * 2 * n + i];
			factorImag = t.pDataImag[j * 2 * n + i];
			for (k = 0; k < 2 * n; ++k)
			// t[j * 2*n + k] -= t[i * 2*n + k] * factor
			{
				complex_mult(t.pDataReal[i * 2 * n + k], t.pDataImag[i * 2 * n + k], factorReal, factorImag, &real, &imag);
				t.pDataReal[j * 2 * n + k] -= real;
				t.pDataImag[j * 2 * n + k] -= imag;
			}
		}
	}

	// приведение левой подматрицы к единичному виду
	for (i = n - 1; i > 0; --i)
	{
		for (j = 0; j < i; ++j)
		{
			factorReal = t.pDataReal[j * 2 * n + i];
			factorImag = t.pDataImag[j * 2 * n + i];

			for (k = 0; k < 2 * n; ++k)
			// t[j * 2*n + k] -= t[i * 2*n + k] * factor
			{
				complex_mult(t.pDataReal[i * 2 * n + k], t.pDataImag[i * 2 * n + k], factorReal, factorImag, &real, &imag);
				t.pDataReal[j * 2 * n + k] -= real;
				t.pDataImag[j * 2 * n + k] -= imag;
			}
		}
	}

	for (i = 0; i < n; ++i)
	{
		memcpy(ainv->pDataReal + i * n, t.pDataReal + i * 2 * n + n, n * sizeof(REAL_TYPE));
		memcpy(ainv->pDataImag + i * n, t.pDataImag + i * 2 * n + n, n * sizeof(REAL_TYPE));
	}
	complex_free(&t);

	return MATRIX_SUCCESS;
}

void complex_zeros(int rows, int columns, CMatrix_t* a, int create)
{
	int i, j;

	if (create > 0)
		complex_new(rows, columns, a);
	else if (a->numCols != columns || a->numRows != rows)
	{
		a->numCols = columns;
		a->numRows = rows;
#ifndef PREALLOCATION
		a->pDataReal = realloc(a->pDataReal, rows*columns*sizeof(REAL_TYPE));
		a->pDataImag = realloc(a->pDataImag, rows*columns*sizeof(REAL_TYPE));
#endif
	}

	for (i = 0; i < rows; ++i)
	{
		for (j = 0; j < columns; ++j)
		{
			a->pDataReal[i*columns + j] = 0.0f;
			a->pDataImag[i*columns + j] = 0.0f;
		}
	}
}

void complex_eye(int n, CMatrix_t* a, int create)
{
	int i, j;

	if (create > 0)
		complex_new(n, n, a);
	else if (a->numCols != n || a->numRows != n)
	{
		a->numCols = n;
		a->numRows = n;
#ifndef PREALLOCATION
		a->pDataReal = realloc(a->pDataReal, n*n*sizeof(REAL_TYPE));
		a->pDataImag = realloc(a->pDataImag, n*n*sizeof(REAL_TYPE));
#endif
	}

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			a->pDataReal[i*n + j] = (i == j) ? 1.0f : 0.0f;
			a->pDataImag[i*n + j] = 0.0f;
		}
	}
}

void complex_free(CMatrix_t* a)
{
	a->numCols = 0;
	a->numRows = 0;
#ifndef PREALLOCATION
	free(a->pDataReal);
	free(a->pDataImag);
	a->pDataReal = NULL;
	a->pDataImag = NULL;

#endif // PREALLOCATION
}

void complex_clone(const CMatrix_t* source, CMatrix_t* destination, int create)
{
	int size;
	if (create > 0)
		size = complex_new(source->numRows, source->numCols, destination);
	else
		size = source->numRows * source->numCols * sizeof(REAL_TYPE);

	if (destination->numRows != source->numRows || destination->numCols != source->numCols)
	{
		destination->numRows = source->numRows;
		destination->numCols = source->numCols;
#ifndef PREALLOCATION
		if (create <= 0)
		{
			destination->pDataReal = realloc(destination->pDataReal, size);
			destination->pDataImag = realloc(destination->pDataImag, size);
		}
#endif
	}

	memcpy(destination->pDataReal, source->pDataReal, size);
	memcpy(destination->pDataImag, source->pDataImag, size);
}

int complex_new(int rows, int columns, CMatrix_t* a)
{
	int size = rows * columns * sizeof(REAL_TYPE);

	a->numRows = rows;
	a->numCols = columns;

#ifndef PREALLOCATION
	a->pDataReal = (REAL_TYPE*)malloc(size);
	a->pDataImag = (REAL_TYPE*)malloc(size);
#endif // PREALLOCATION
	return size;
}

int complex_partial_copy(const CMatrix_t* source, int source_start_row, int source_start_column,
	int rows_to_copy, int columns_to_copy,
	CMatrix_t* destination, int destination_start_row, int destination_start_column)
{
	int i, j, i_source, i_destination;

	if (source_start_row + rows_to_copy > source->numRows
		|| source_start_column + columns_to_copy > source->numCols
		|| destination_start_row + rows_to_copy > destination->numRows
		|| destination_start_column + columns_to_copy > destination->numCols)
		return MATRIX_EXCEEDS;

	for (i = 0; i < rows_to_copy; ++i)
	{
		for (j = 0; j < columns_to_copy; ++j)
		{
			i_source = (i + source_start_row) * source->numCols + j + source_start_column;
			i_destination = (i + destination_start_row) * destination->numCols + j + destination_start_column;
			destination->pDataReal[i_destination] = source->pDataReal[i_source];
			destination->pDataImag[i_destination] = source->pDataImag[i_source];
		}
	}

	return MATRIX_SUCCESS;
}

int complex_matrix_sum(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res)
{
	int i, j, idx;
	REAL_TYPE re, im;

	if (a->numRows != b->numRows || a->numRows != res->numRows
		|| a->numCols != b->numCols || a->numCols != res->numCols)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < a->numCols; ++j)
		{
			idx = i*a->numCols + j;
			complex_mult(factor_re, factor_im, b->pDataReal[idx], b->pDataImag[idx], &re, &im);
			res->pDataReal[idx] = a->pDataReal[idx] + re;
			res->pDataImag[idx] = a->pDataImag[idx] + im;
		}
	}

	return MATRIX_SUCCESS;
}

int complex_scal_prod(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE* re, REAL_TYPE* im)
{
	int i;
	REAL_TYPE tmp_re, tmp_im;

	/*if (a->numRows != 1 || b->numCols != 1 || a->numCols != b->numRows)
		return MATRIX_EXCEEDS;*/

	*re = 0.0f;
	*im = 0.0f;

	for (i = 0; i < int_max(a->numRows, a->numCols); ++i)
	{
		complex_mult(a->pDataReal[i], a->pDataImag[i], b->pDataReal[i], b->pDataImag[i], &tmp_re, &tmp_im);
		*re += tmp_re;
		*im += tmp_im;
	}

	return MATRIX_SUCCESS;
}

int real_new(int rows, int columns, Matrix_t* a)
{
	int size = rows * columns * sizeof(REAL_TYPE);

	a->numRows = rows;
	a->numCols = columns;

#ifndef PREALLOCATION
	a->pData = (REAL_TYPE*)malloc(size);
#endif
	return size;
}

void real_free(Matrix_t* a)
{
	a->numCols = 0;
	a->numRows = 0;
#ifndef PREALLOCATION
	free(a->pData);
	a->pData = NULL;
#endif // PREALLOCATION
}

void real_zeros(int rows, int columns, Matrix_t* a, int create)
{
	int i, j;

	if (create > 0)
		real_new(rows, columns, a);
	else if (a->numCols != columns || a->numRows != rows)
	{
		a->numCols = columns;
		a->numRows = rows;
#ifndef PREALLOCATION
		a->pData = realloc(a->pData, rows*columns*sizeof(REAL_TYPE));
#endif
	}

	for (i = 0; i < rows; ++i)
	{
		for (j = 0; j < columns; ++j)
		{
			a->pData[i*columns + j] = 0.0f;
		}
	}
}

void real_eye(int n, Matrix_t* a, int create)
{
	int i, j;

	if (create > 0)
		real_new(n, n, a);
	else if (a->numCols != n || a->numRows != n)
	{
		a->numCols = n;
		a->numRows = n;
#ifndef PREALLOCATION
		a->pData = realloc(a->pData, n*n*sizeof(REAL_TYPE));
#endif
	}

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			a->pData[i*n + j] = (i == j) ? 1.0f : 0.0f;
		}
	}
}

void real_clone(const Matrix_t* source, Matrix_t* destination, int create)
{
	int size;
	if (create > 0)
		size = real_new(source->numRows, source->numCols, destination);
	else
		size = source->numRows * source->numCols * sizeof(REAL_TYPE);

	if (destination->numRows != source->numRows || destination->numCols != source->numCols)
	{
		destination->numRows = source->numRows;
		destination->numCols = source->numCols;
#ifndef PREALLOCATION
		if (create <= 0)
			destination->pData = realloc(destination->pData, size);
#endif
	}

	memcpy(destination->pData, source->pData, size);
}

int real_partial_copy(const Matrix_t* source, int source_start_row, int source_start_column,
	int rows_to_copy, int columns_to_copy,
	Matrix_t* destination, int destination_start_row, int destination_start_column)
{
	int i, j, i_source, i_destination;

	if (source_start_row + rows_to_copy > source->numRows
		|| source_start_column + columns_to_copy > source->numCols
		|| destination_start_row + rows_to_copy > destination->numRows
		|| destination_start_column + columns_to_copy > destination->numCols)
		return MATRIX_EXCEEDS;

	for (i = 0; i < rows_to_copy; ++i)
	{
		for (j = 0; j < columns_to_copy; ++j)
		{
			i_source = (i + source_start_row) * source->numCols + j + source_start_column;
			i_destination = (i + destination_start_row) * destination->numCols + j + destination_start_column;
			destination->pData[i_destination] = source->pData[i_source];
		}
	}

	return MATRIX_SUCCESS;
}

int real_transp(const Matrix_t* a, Matrix_t* b, int create)
{
	int i, j, idx, idxt;

	if (create > 0)
		real_new(a->numCols, a->numRows, b);
	else if (b->numCols != a->numRows || b->numRows != a->numCols)
	{
		b->numCols = a->numRows;
		b->numRows = a->numCols;
#ifndef PREALLOCATION
		b->pData = realloc(b->pData, a->numCols*a->numRows*sizeof(REAL_TYPE));
#endif
	}


	if (a->numRows != b->numCols || a->numCols != b->numRows)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < a->numCols; ++j)
		{
			idx = i * a->numCols + j;
			idxt = j * a->numCols + i;
			b->pData[idx] = a->pData[idxt];
		}
	}
	return MATRIX_SUCCESS;
}

int real_matrix_mult(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res)
{
	int i, j, k, i1, j1;

	REAL_TYPE sum;

	if (a->numCols != b->numRows || res->numRows != a->numRows || res->numCols != b->numCols)
		return MATRIX_EXCEEDS;

	for (i  = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < b->numCols; ++j)
		{
			sum = 0.0f;
			for (k = 0; k < a->numCols; ++k)
			{
				i1 = i*a->numCols + k;
				j1 = k*b->numCols + j;
				sum += a->pData[i1] * b->pData[j1];
			}
			res->pData[i*res->numCols + j] = factor * sum;
		}
	}

	return MATRIX_SUCCESS;
}

int real_to_complex(const Matrix_t* real, CMatrix_t* compl, int create)
{
	int size;

	if (create > 0)
		size = complex_new(real->numRows, real->numCols, compl);
	else
	{
		if (compl->numCols != real->numCols || compl->numRows != real->numRows)
			return MATRIX_EXCEEDS;
		size = real->numCols * real->numRows * sizeof(REAL_TYPE);
	}

	memcpy(compl->pDataReal, real->pData, size);
	memset(compl->pDataImag, 0, size);
	return MATRIX_SUCCESS;
}

void complex_swap_columns(CMatrix_t* m, int* new_indices)
{
	int i, j;
	CMatrix_t t;
#ifdef PREALLOCATION
	static REAL_TYPE buffer[2 * MATRIX_MAX_SIZE];
	t.pDataReal = buffer;
	t.pDataImag = buffer + MATRIX_MAX_SIZE;
#endif // PREALLOCATION

	//complex_clone(m, &t, ALLOCATE_MATRIX);

	for (i = 0; i < m->numCols; ++i)
	{
		for (j = 0; j < m->numRows; ++j)
		{
			t.pDataReal[j*m->numCols + i] = m->pDataReal[j*m->numCols + new_indices[i]];
			t.pDataImag[j*m->numCols + i] = m->pDataImag[j*m->numCols + new_indices[i]];
		}
	}

	swap_ptr((void**)&t.pDataReal, (void**)&m->pDataReal); // swap корректен и в случае преаллокации
	swap_ptr((void**)&t.pDataImag, (void**)&m->pDataImag); // буферы являются static, их жизнь продолжается до конца программы

	complex_free(&t);
}

void single_to_double(const CMatrix_t * S, CMatrix_td *D, int create)
{
	int i;
	if (create != 0)
		complex_new_d(S->numRows, S->numCols, D);

	for (i = 0; i < D->numRows * D->numCols; ++i)
	{
		D->pDataReal[i] = S->pDataReal[i];
		D->pDataImag[i] = S->pDataImag[i];
	}
}

int complex_new_d(int rows, int columns, CMatrix_td* a)
{
	int size = rows * columns * sizeof(double);

	a->numRows = rows;
	a->numCols = columns;

#ifndef PREALLOCATION
	a->pDataReal = (double*)malloc(size);
	a->pDataImag = (double*)malloc(size);
#endif // PREALLOCATION
	return size;
}

void complex_eye_d(int n, CMatrix_td* a, int create)
{
	int i, j;

	if (create > 0)
		complex_new_d(n, n, a);
	else if (a->numCols != n || a->numRows != n)
	{
		a->numCols = n;
		a->numRows = n;
#ifndef PREALLOCATION
		a->pDataReal = realloc(a->pDataReal, n*n*sizeof(double));
		a->pDataImag = realloc(a->pDataImag, n*n*sizeof(double));
#endif
	}

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			a->pDataReal[i*n + j] = (i == j) ? 1.0f : 0.0f;
			a->pDataImag[i*n + j] = 0.0f;
		}
	}
}

void complex_clone_d(const CMatrix_td* source, CMatrix_td* destination, int create)
{
	int size;
	if (create > 0)
		size = complex_new_d(source->numRows, source->numCols, destination);
	else
		size = source->numRows * source->numCols * sizeof(double);

	if (destination->numRows != source->numRows || destination->numCols != source->numCols)
	{
		destination->numRows = source->numRows;
		destination->numCols = source->numCols;
#ifndef PREALLOCATION
		if (create <= 0)
		{
			destination->pDataReal = realloc(destination->pDataReal, size);
			destination->pDataImag = realloc(destination->pDataImag, size);
		}
#endif
	}

	memcpy(destination->pDataReal, source->pDataReal, size);
	memcpy(destination->pDataImag, source->pDataImag, size);
}

void complex_zeros_d(int rows, int columns, CMatrix_td* a, int create)
{
	int i, j;

	if (create > 0)
		complex_new_d(rows, columns, a);
	else if (a->numCols != columns || a->numRows != rows)
	{
		a->numCols = columns;
		a->numRows = rows;
#ifndef PREALLOCATION
		a->pDataReal = realloc(a->pDataReal, rows*columns*sizeof(double));
		a->pDataImag = realloc(a->pDataImag, rows*columns*sizeof(double));
#endif
	}

	for (i = 0; i < rows; ++i)
	{
		for (j = 0; j < columns; ++j)
		{
			a->pDataReal[i*columns + j] = 0.0f;
			a->pDataImag[i*columns + j] = 0.0f;
		}
	}
}

int complex_partial_copy_d(const CMatrix_td* source, int source_start_row, int source_start_column,
	int rows_to_copy, int columns_to_copy,
	CMatrix_td* destination, int destination_start_row, int destination_start_column)
{
	int i, j, i_source, i_destination;

	if (source_start_row + rows_to_copy > source->numRows
		|| source_start_column + columns_to_copy > source->numCols
		|| destination_start_row + rows_to_copy > destination->numRows
		|| destination_start_column + columns_to_copy > destination->numCols)
		return MATRIX_EXCEEDS;

	for (i = 0; i < rows_to_copy; ++i)
	{
		for (j = 0; j < columns_to_copy; ++j)
		{
			i_source = (i + source_start_row) * source->numCols + j + source_start_column;
			i_destination = (i + destination_start_row) * destination->numCols + j + destination_start_column;
			destination->pDataReal[i_destination] = source->pDataReal[i_source];
			destination->pDataImag[i_destination] = source->pDataImag[i_source];
		}
	}

	return MATRIX_SUCCESS;
}

void cnorm2_d(double * a_real, double * a_imag, double * norm2, uint32_t dir, uint32_t len)
{
	uint32_t i;
	double *pa_real, *pa_imag;
	pa_real = a_real;
	pa_imag = a_imag;

	*norm2 = 0.0;
	for (i = 0; i < len; i++)
	{
		(*norm2) = (*norm2) + (*pa_real) * (*pa_real) + (*pa_imag) * (*pa_imag);
		pa_real += dir;
		pa_imag += dir;
	}
}

int complex_matrix_sum_d(const CMatrix_td* a, const CMatrix_td* b, double factor_re, double factor_im, CMatrix_td* res)
{
	int i, j, idx;
	double re, im;

	if (a->numRows != b->numRows || a->numRows != res->numRows
		|| a->numCols != b->numCols || a->numCols != res->numCols)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < a->numCols; ++j)
		{
			idx = i*a->numCols + j;
			complex_mult_d(factor_re, factor_im, b->pDataReal[idx], b->pDataImag[idx], &re, &im);
			res->pDataReal[idx] = a->pDataReal[idx] + re;
			res->pDataImag[idx] = a->pDataImag[idx] + im;
		}
	}

	return MATRIX_SUCCESS;
}

int complex_transp_d(const CMatrix_td* a, CMatrix_td* b, int create)
{
	int i, j, idx, idxt;

	if (create > 0)
		complex_new_d(a->numCols, a->numRows, b);
	else if (b->numCols != a->numRows || b->numRows != a->numCols)
	{
		b->numCols = a->numRows;
		b->numRows = a->numCols;
#ifndef PREALLOCATION
		b->pDataReal = realloc(b->pDataReal, a->numCols*a->numRows*sizeof(double));
		b->pDataImag = realloc(b->pDataImag, a->numCols*a->numRows*sizeof(double));
#endif
	}

	if (a->numRows != b->numCols || a->numCols != b->numRows)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < a->numCols; ++j)
		{
			idx = i * a->numCols + j;
			idxt = j * a->numCols + i;
			b->pDataReal[idx] = a->pDataReal[idxt];
			b->pDataImag[idx] = -1.0 * a->pDataImag[idxt];
		}
	}
	return MATRIX_SUCCESS;
}

int complex_scal_prod_d(const CMatrix_td* a, const CMatrix_td* b, double* re, double* im)
{
	int i;
	double tmp_re, tmp_im;

	/*if (a->numRows != 1 || b->numCols != 1 || a->numCols != b->numRows)
	return MATRIX_EXCEEDS;*/

	*re = 0.0f;
	*im = 0.0f;

	for (i = 0; i < int_max(a->numRows, a->numCols); ++i)
	{
		complex_mult_d(a->pDataReal[i], a->pDataImag[i], b->pDataReal[i], b->pDataImag[i], &tmp_re, &tmp_im);
		*re += tmp_re;
		*im += tmp_im;
	}

	return MATRIX_SUCCESS;
}

int complex_matrix_mult_d(const CMatrix_td* a, const CMatrix_td* b, double factor_re, double factor_im, CMatrix_td* res)
{
	int i, j, k, i1, j1;

	double real, img;
	double x, y;

	if (a->numCols != b->numRows || res->numRows != a->numRows || res->numCols != b->numCols)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < b->numCols; ++j)
		{
			real = 0.0;
			img = 0.0;
			for (k = 0; k < a->numCols; ++k)
			{
				i1 = i*a->numCols + k;
				j1 = k*b->numCols + j;
				complex_mult_d(a->pDataReal[i1], a->pDataImag[i1],
					b->pDataReal[j1], b->pDataImag[j1], &x, &y);
				real += x;
				img += y;
			}
			k = i*res->numCols + j;
			complex_mult_d(real, img, factor_re, factor_im, res->pDataReal + k, res->pDataImag + k);
		}
	}

	return MATRIX_SUCCESS;
}

void double_complex_to_single_real(const CMatrix_td* D, Matrix_t* S, int create)
{
	int i;
	if (create != 0)
		real_new(D->numRows, D->numCols, S);

	for (i = 0; i < D->numRows * D->numCols; ++i)
		S->pData[i] = (REAL_TYPE)D->pDataReal[i];
}

void complex_free_d(CMatrix_td* a)
{
	a->numCols = 0;
	a->numRows = 0;
#ifndef PREALLOCATION
	free(a->pDataReal);
	free(a->pDataImag);
	a->pDataReal = NULL;
	a->pDataImag = NULL;

#endif // PREALLOCATION
}

void double_complex_to_single_comlpex(const CMatrix_td* D, CMatrix_t* S, int create)
{
	int i;
	if (create != 0)
		complex_new(D->numRows, D->numCols, S);

	for (i = 0; i < D->numRows * D->numCols; ++i)
	{
		S->pDataReal[i] = (REAL_TYPE)D->pDataReal[i];
		S->pDataImag[i] = (REAL_TYPE)D->pDataImag[i];
	}
}

int complex_matrix_mult_right_transp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res)
{
	int i, j, k, i1, j1;

	REAL_TYPE real, img;
	REAL_TYPE x, y;

	if (a->numCols != b->numCols || res->numRows != a->numRows || res->numCols != b->numRows)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < b->numRows; ++j)
		{
			real = 0.0;
			img = 0.0;
			for (k = 0; k < a->numCols; ++k)
			{
				i1 = i*a->numCols + k;
				j1 = j*b->numCols + k;
				complex_mult(a->pDataReal[i1], a->pDataImag[i1],
					b->pDataReal[j1], -b->pDataImag[j1], &x, &y);
				real += x;
				img += y;
			}
			k = i*res->numCols + j;
			complex_mult(real, img, factor_re, factor_im, res->pDataReal + k, res->pDataImag + k);
		}
	}

	return MATRIX_SUCCESS;
}

#ifdef DSP_OPTIMIZATION
int complex_matrix_mult_dsp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pDataReal, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), b->pDataReal, b->numRows * b->numCols * sizeof(float));
	memcpy(&(g_pImag1[0]), a->pDataImag, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pImag2[0]), b->pDataImag, b->numRows * b->numCols * sizeof(float));

	g_Factor = factor_re;
	g_FactorIm = factor_im;
	g_nRows1 = a->numRows;
	g_nColumns1 = a->numCols;
	g_nColumns2 = b->numCols;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(ComplexMatrixMult));
	memcpy(res->pDataReal, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));
	memcpy(res->pDataImag, &(g_pImag3[0]), res->numRows * res->numCols * sizeof(float));

	return MATRIX_SUCCESS;
}

int real_matrix_mult_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pData, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), b->pData, b->numRows * b->numCols * sizeof(float));

	g_Factor = factor;
	g_nRows1 = a->numRows;
	g_nColumns1 = a->numCols;
	g_nColumns2 = b->numCols;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(RealMatrixMult));
	memcpy(res->pData, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));

	return MATRIX_SUCCESS;
}

int inv_dsp(const CMatrix_t* a, CMatrix_t* ainv, int numDsp)
{
    init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pDataReal, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pImag1[0]), a->pDataImag, a->numRows * a->numCols * sizeof(float));

	g_nRows1 = a->numRows;
	run_dsp(numDsp, DSP_ROUTINE_ADDR(ComplexMatrixInv));

	memcpy(ainv->pDataReal, &(g_pReal1[0]), ainv->numRows * ainv->numCols * sizeof(float));
	memcpy(ainv->pDataImag, &(g_pImag1[0]), ainv->numRows * ainv->numCols * sizeof(float));
	return MATRIX_SUCCESS;

}

#endif // DSP_OPTIMIZATION
