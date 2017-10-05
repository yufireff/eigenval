#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "matrix.h"

void swap_real(REAL_TYPE*a, REAL_TYPE* b);

void swap_int(int*a, int* b);

void swap_ptr(void**a, void** b);

void special_sort_ascending(REAL_TYPE *array, int n, int* indices);

void special_sort_descending(REAL_TYPE *array, int n, int* indices);

int find_max(REAL_TYPE *array, int n, REAL_TYPE *max);

#endif // ALGORITHMS_H

