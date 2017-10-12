#ifndef SETTINGS_H
#define SETTINGS_H

#define VECTOR_MAX_SIZE 20
#define MATRIX_MAX_SIZE VECTOR_MAX_SIZE*VECTOR_MAX_SIZE

#define PREALLOCATION

#ifdef PREALLOCATION
#define ALLOCATE_MATRIX 0
#else
#define ALLOCATE_MATRIX 1
#endif // PREALLOCATION

//#define DOUBLE

#ifdef DOUBLE
#define REAL_TYPE double
#else
#define REAL_TYPE float
#endif // DOUBLE

#ifndef WIN32
//#define DSP_OPTIMIZATION
#endif // WIN32

//#define TEST_RESULT
#endif // SETTINGS_H
