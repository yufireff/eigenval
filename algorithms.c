#include "algorithms.h"

void swap_real(REAL_TYPE *a, REAL_TYPE* b)
{
	REAL_TYPE t = *a;
	*a = *b;
	*b = t;
}

void swap_int(int *a, int* b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

void swap_ptr(void** a, void** b)
{
	void* t = *a;
	*a = *b;
	*b = t;
}

void special_sort_ascending(REAL_TYPE *array, int n, int* indices)
// сортировка массива элементов REAL_TYPE пузырьковым методом
// в качестве дополнительного выходного аргумента также получаем массив индексов
{
	int i, j;

	for (i = 0; i < n; ++i)
		indices[i] = i;

	for (i = 0; i < n - 1; ++i)
	{
		for (j = 0; j < n - 1; ++j)
		{
			if (array[j] > array[j + 1])
			{
				swap_real(array + j, array + j + 1);
				swap_int(indices + j, indices + j + 1);
			}
		}
	}
}

void special_sort_descending(REAL_TYPE *array, int n, int* indices)
// сортировка массива элементов REAL_TYPE пузырьковым методом
// в качестве дополнительного выходного аргумента также получаем массив индексов
{
	int i, j;

	for (i = 0; i < n; ++i)
		indices[i] = i;

	for (i = 0; i < n - 1; ++i)
	{
		for (j = 0; j < n - 1; ++j)
		{
			if (array[j] < array[j + 1])
			{
				swap_real(array + j, array + j + 1);
				swap_int(indices + j, indices + j + 1);
			}
		}
	}
}

int find_max(REAL_TYPE *array, int n, REAL_TYPE *max)
{
	int res = 0, i;
	*max = array[0];
	for (i = 1; i < n; ++i)
	{
		if (array[i] > *max)
		{
			*max = array[i];
			res = i;
		}
	}
	return res;
}
