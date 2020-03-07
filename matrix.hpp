#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>



double*     allocateMatrix      (uint64_t n,uint64_t m) ;
double*     allocateVector      (uint64_t n) ;
void        freeMatrix          (double *A);
void        freeVector          (double *v);
void        setMatrixZero       (double *A, uint64_t n, uint64_t m);
void        setMatrixIdentity   (double *A, uint64_t n);
void        copyMatrix          (double *B, double *A, uint64_t n, uint64_t m) ;
void        writeMatrix         (FILE *stream, double *A, uint64_t n, uint64_t m);
void        absMatrix           (double *Aabs,double *A, uint64_t n, uint64_t m);
double      getMaxInMatrix      (double max, double *A, uint64_t n, uint64_t m);
void        matrixSub           (double *S, double *A, double *B, uint64_t n, uint64_t m);
void        matrixAdd           (double *S, double *A, double *B, uint64_t n, uint64_t m);

/* Rajouter les prototypes de vos méthodes ici. Par exemple */


/*
Solves a system of linear equations Ax=b, given a matrix A (size n x n) and vector b(size n)
Using LU factorisation algorithm.
*/
bool  SolveSystemLU_Decomp   (double *x, double *A, double *b, uint64_t n);

/*
 Calculation of the determinant of a square matrix A dimensioned n*n
*/
double det(double *A,uint64_t n);

/*
Performs LU Decompostion to square matrix A dimensioned n*n.
*/
bool decompLU(double *A, uint64_t n);
/* Performs naive multiplication of matrix A (size p x k) by a matrix B (size k x r).
The result matrix S = A*B  is of size (k x r).
We assume that S has already been allocated outside the function.
*/
void        matrixMultiplyNaive (double *S, double *A, double *B, uint64_t p, uint64_t k, uint64_t r);

/* Performs a multiplication of two sqaure matrices A and B (size n x n) by Strassen algorithm.
    We assume that S has already been allocated outside the function.
*/
void        matrixMultiplyStrassen (double *S, double *A, double *B, uint64_t n);

/* 
    Solves a system of linear equations Ax=b for a double-precision matrix A (size n x n).
    Uses iterative ascension algorithm. 
    After the procedure, x contains the solution of Ax=b.
    We assume that x has been allocated outside the function.
*/
void        SolveTriangularSystemUP   (double *x, double *A, double *b, uint64_t n);

/* 
    Performs Gauss elimination for given a matrix A (size n x n) and a vector b (size n).
    Modifies directly matrix A and vector b.
    In the end of the procedure, A is upper truangular and b is modified accordingly.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is impossible to triangularize. 
*/
bool        Triangularize           (double *A, double *b, uint64_t n);

/*
    Solves a system of linear equations Ax=b, given a matrix A (size n x n) and vector b(size n).
    Uses Gauss elimination algorithm based on truangularization and the ascension solving.
    After the procedure, vector x contains the solution to Ax=b.
    We assume that x has been allocated outside the function.
        Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is of rank <n .
*/
bool        SolveSystemGauss        (double *x, double *A, double *b, uint64_t n);
