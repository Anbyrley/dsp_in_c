/** @file linear_algebra.h
*   @brief Contains functions for performing various standard linear algebra manipulations.
*
*
*  @author Alex N. Byrley (anbyrley)
*  @date October 2015
*  @bug No known bugs
*/

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

//================================================================================================//
//========================================INCLUDES================================================//
//================================================================================================//

#include "macros.h"
#include "helper.h"
#include "cblas.h"
#include "lapacke.h"

//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//

//================================================================================================//
/**
* @brief This function forward backward averages a double complex matrix.
*
* NOTE: Can be used to force hermitian symmetry.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void fb_average_matrix_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function forward backward averages a double precision matrix.
*
* NOTE: Can be used to force symmetry.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void fb_average_matrix_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function computes the lower cholesky factor for a symmetric positive def dbl matrix.
*
* @param[in] double* matrix
* @param[in] int order
* @param[out] double* cholesky_factor
*
* @return NONE
*/
//================================================================================================//
void compute_lower_cholesky_factor_dbl(double*,int,double*);


//================================================================================================//
/**
* @brief This function computes the upper cholesky factor for a symmetric positive def dbl matrix.
*
* @param[in] double* matrix
* @param[in] int order
* @param[out] double* cholesky_factor
*
* @return NONE
*/
//================================================================================================//
void compute_upper_cholesky_factor_dbl(double*,int,double*);


//================================================================================================//
/**
* @brief This function performs a singular value decomposition on a double complex matrix.
*
* @param[in] double complex* matrix
* @param[in] int num_rows
* @param[in] int num_cols
* @param[out] double complex* U
* @param[out] double complex* S
* @param[out] double complex* VH
*
* @return NONE
*/
//================================================================================================//
void svd_cmplx(double complex*,int,int,double complex*,double complex*,double complex*);


//================================================================================================//
/**
* @brief This function performs a singular value decomposition on a double precision matrix.
*
* @param[in] double* matrix
* @param[in] int num_rows
* @param[in] int num_cols
* @param[out] double* U
* @param[out] double* S
* @param[out] double* VT
*
* @return NONE
*/
//================================================================================================//
void svd_dbl(double*,int,int,double*,double*,double*);


//================================================================================================//
/**
* @brief This function computes the trace of a double precision matrix using the SVD.
*
* @param[in] double* matrix
* @param[in] int num_rows
* @param[in] int num_cols
*
* @return double trace
*/
//================================================================================================//
double compute_trace_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function computes the determinant of a double precision matrix using the SVD.
*
* @param[in] double* matrix
* @param[in] int num_rows
* @param[in] int num_cols
*
* @return double determinant
*/
//================================================================================================//
double compute_determinant_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function performs a low rank approximation of a double complex matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int rank
*
* @return NONE
*/
//================================================================================================//
void low_rank_approximation_cmplx(double complex*,unsigned int,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function performs a low rank approximation of a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int rank
*
* @return NONE
*/
//================================================================================================//
void low_rank_approximation_dbl(double*,unsigned int,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function performs a low rank partition of a double complex matrix into num_rowsxrank.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int rank
*
* @return NONE
*/
//================================================================================================//
void low_rank_partition_cmplx(double complex*,unsigned int,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function raises a double complex matrix to a power.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double power
*
* @return NONE
*/
//================================================================================================//
void matrix_power_cmplx(double complex*,unsigned int,unsigned int,double);


//================================================================================================//
/**
* @brief This function performs QR decomposition on a double precision matrix.
*
* @param[in] double complex* matrix
* @param[in] int num_rows
* @param[in] int num_cols
* @param[out] double* Q
* @param[out] double* R
*
* @return NONE
*/
//================================================================================================//
void qr_decomposition_dbl(double*,int,int,double*,double*);


//================================================================================================//
/**
* @brief This function compute the inverse of a double complex upper triangular matrix.
*
* @param[in] double complex* matrix
* @param[in] int order
* @param[out] double complex* inverse
*
* @return NONE
*/
//================================================================================================//
void compute_upper_triangular_matrix_inverse_cmplx(double complex*,int,double complex*);


//================================================================================================//
/**
* @brief This function compute the inverse of a double precision upper triangular matrix.
*
* @param[in] double* matrix
* @param[in] int order
* @param[out] double* inverse
*
* @return NONE
*/
//================================================================================================//
void compute_upper_triangular_matrix_inverse_dbl(double*,int,double*);


//================================================================================================//
/**
* @brief This function compute the left inverse of a double complex matrix via QR Decomposition.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double complex* inverse
*
* @return NONE
*/
//================================================================================================//
void compute_left_inverse_qr_cmplx(double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function compute the left inverse of a double precision matrix via QR Decomposition.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double* inverse
*
* @return NONE
*/
//================================================================================================//
void compute_left_inverse_qr_dbl(double*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function performs a rank1 update on a double complex matrix.
*
* NOTE: Computes -- alpha * A + beta * xx^H
*
* @param[in,out] double complex* matrix
* @param[in] double complex alpha
* @param[in] double complex* vector
* @param[in] unsigned int vector_length
* @param[in] double complex beta
*
* @return NONE
*/
//================================================================================================//
void rank1_matrix_update_cmplx(double complex*,double complex,double complex*,unsigned int,double complex);


//================================================================================================//
/**
* @brief This function performs a rank1 update on a double precision matrix.
*
* NOTE: Computes -- alpha * A + beta * xx^T
*
* @param[in,out] double* matrix
* @param[in] double alpha
* @param[in] double* vector
* @param[in] unsigned int vector_length
* @param[in] double beta
*
* @return NONE
*/
//================================================================================================//
void rank1_matrix_update_dbl(double*,double,double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function performs a rank1 update on the inverse of a double complex matrix.
*
* @param[in,out] double complex* inverse
* @param[in] double complex alpha
* @param[in] double complex* vector
* @param[in] double complex beta
* @param[in] unsigned int vector_length
*
* @return NONE
*/
//================================================================================================//
void rank1_matrix_inverse_update_cmplx(double complex*,double complex,double complex*,
									   double complex,unsigned int);


//================================================================================================//
/**
* @brief This function performs a rank1 update on the inverse of a double precision matrix.
*
* @param[in,out] double* inverse
* @param[in] double alpha
* @param[in] double* vector
* @param[in] double beta
* @param[in] int vector_length
*
* @return NONE
*/
//================================================================================================//
void rank1_matrix_inverse_update_dbl(double*,double, double*,double,unsigned int);


//================================================================================================//
/**
* @brief This function forms the projection matrix for a double complex vector or matrix.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double complex* projection_matrix
*
* @return NONE
*/
//================================================================================================//
void form_projection_matrix_cmplx(double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function forms the orthogonal projection matrix for a double complex vector or matrix.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double complex* orthogonal_projection_matrix
*
* @return NONE
*/
//================================================================================================//
void form_orthogonal_projection_matrix_cmplx(double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function computes the outer product between two double complex vectors.
*
* @param[in] double complex* vector1
* @param[in] int vector1_length
* @param[in] double complex* vector2
* @param[in] int vector2_length
* @param[out] double complex* outer_product
*
* @return NONE
*/
//================================================================================================//
void vector_outer_product_cmplx(double complex*,int,double complex*,int,double complex*);


//================================================================================================//
/**
* @brief This function computes the outer product between two double precision vectors.
*
* @param[in] double* vector1
* @param[in] int vector1_length
* @param[in] double* vector2
* @param[in] int vector2_length
* @param[out] double* outer_product
*
* @return NONE
*/
//================================================================================================//
void vector_outer_product_dbl(double*,int,double*,int,double*);


//================================================================================================//
/**
* @brief This function multiplies a double precision vector by a double precision matrix.
*
* @param[in] double* vector
* @param[in] int vector_size
* @param[in] double* matrix
* @param[in] int matrix_rows
* @param[in] int matrix_columns
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void vector_matrix_multiply_dbl(double*,int,double*,int,int,double*);


//================================================================================================//
/**
* @brief This function multiplies a double complex vector by a double complex matrix.
*
* @param[in] double complex* vector
* @param[in] int vector_size
* @param[in] double complex* matrix
* @param[in] int matrix_rows
* @param[in] int matrix_columns
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void vector_matrix_multiply_cmplx(double complex*,int,double complex*,int,int,double complex*);


//================================================================================================//
/**
* @brief This function multiplies a double precision matrix by a double precision vector.
*
* @param[in] double* vector
* @param[in] int vector_size
* @param[in] double* matrix
* @param[in] int matrix_rows
* @param[in] int matrix_columns
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void matrix_vector_multiply_dbl(double*,int,double*,int,int,double*);


//================================================================================================//
/**
* @brief This function multiplies a double precision matrix by a double complex vector.
*
* @param[in] double complex* vector
* @param[in] int vector_size
* @param[in] double complex* matrix
* @param[in] int matrix_rows
* @param[in] int matrix_columns
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void matrix_vector_multiply_cmplx(double complex*,int,double complex*,int,int,double complex*);


//================================================================================================//
/**
* @brief This function multiplies two double complex matrices together.
*
* @param[in] double complex* matrix1
* @param[in] int matrix1_rows
* @param[in] int matrix1_columns
* @param[in] double complex* matrix2
* @param[in] int matrix2_rows
* @param[in] int matrix2_columns
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void matrix_matrix_multiply_cmplx(double complex*,int,int,double complex*,int,int,double complex*);


//================================================================================================//
/**
* @brief This function multiplies two double precision matrices together.
*
* @param[in] double* matrix1
* @param[in] int matrix1_rows
* @param[in] int matrix1_columns
* @param[in] double* matrix2
* @param[in] int matrix2_rows
* @param[in] int matrix2_columns
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void matrix_matrix_multiply_dbl(double*,int,int,double*,int,int,double*);


//================================================================================================//
/**
* @brief This function computes the inverse of a double precision matrix
*
* @param[in] double* matrix
* @param[in] int matrix_order
* @param[out] double* inverse
*
* @return NONE
*/
//================================================================================================//
void compute_matrix_inverse_dbl(double*,int,double*);


//================================================================================================//
/**
* @brief This function computes the inverse of a double complex matrix
*
* @param[in] double complex* matrix
* @param[in] int matrix_order
* @param[out] double complex* inverse
*
* @return NONE
*/
//================================================================================================//
void compute_matrix_inverse_cmplx(double complex*,int,double complex*);


//================================================================================================//
/**
* @brief This function computes the solution of the equation Ax = b.
*
* NOTE: A must be square! Computes the solution using LU factorization.
*
* @param[in] double* A
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double* b
* @param[out] double* x
*
* @return NONE
*/
//================================================================================================//
void solve_linear_system_dbl(double*,unsigned int,unsigned int,double*,double*);


//================================================================================================//
/**
* @brief This function computes the solution of the equation Ax = b.
*
* NOTE: A may be over, under, or exactly determined.
*
* @param[in] double* A
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double* b
* @param[out] double* x
*
* @return NONE
*/
//================================================================================================//
void solve_linear_system_qr_dbl(double*,int,int,double*,double*);


//================================================================================================//
/**
* @brief This function creates a double complex matrix psuedo inverse.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double complex* pinv
*
* @return NONE
*/
//================================================================================================//
void create_psuedo_inverse_cmplx(double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function creates a double precision matrix psuedo inverse via the SVD.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double* pinv
*
* @return NONE
*/
//================================================================================================//
void create_psuedo_inverse_svd_dbl(double*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function creates a double precision matrix psuedo inverse.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double* pinv
*
* @return NONE
*/
//================================================================================================//
void create_psuedo_inverse_dbl(double*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes the solution of the equation A^HAx = A^Hb
*
* @param[in] double complex* A
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double complex* b
* @param[out] double complex* x
*
* @return NONE
*/
//================================================================================================//
void solve_least_squares_system_cmplx(double complex*,unsigned int,unsigned int,
									  double complex*,double complex*);


//================================================================================================//
/**
* @brief This function computes the solution of the equation A^TAx = A^Tb.
*
* @param[in] double* A
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double* b
* @param[out] double* x
*
* @return NONE
*/
//================================================================================================//
void solve_least_squares_system_dbl(double*,unsigned int,unsigned int,double*,double*);


//================================================================================================//
/**
* @brief This function computes (y^H)Ax for complex entries.
*
* @param[in] double complex* A
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double complex* x
* @param[in] double complex* y
*
* @return double complex result
*/
//================================================================================================//
double complex compute_quadratic_form_cmplx(double complex*,unsigned int,unsigned int,
											double complex*,double complex*);


//================================================================================================//
/**
* @brief This function computes (y^T)Ax for complex entries.
*
* @param[in] double* A
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double* x
* @param[in] double* y
*
* @return double result
*/
//================================================================================================//
double compute_quadratic_form_dbl(double*,unsigned int,unsigned int,double*,double*);


//================================================================================================//
/**
* @brief This function computes the eigensystem (values and vectors) for a given matrix.
*
* @param[in] double complex* matrix1
* @param[in] int matrix_order
* @param[in] unsigned int sort
* @param[out] double complex* eigenvalues
* @param[out] double complex* right_eigenvectors
* @param[out] double complex* left_eigenvectors
*
* @return NONE
*/
//================================================================================================//
void compute_eigensystem_cmplx(double complex* matrix,int,unsigned int,double complex*,double complex*,double complex*);


//================================================================================================//
/**
* @brief This function computes the eigensystem (values and vectors) for a given double matrix.
*
* @param[in] double* matrix
* @param[in] int matrix_order
* @param[in] unsigned int sort
* @param[out] double complex* eigenvalues
* @param[out] double complex* right_eigenvectors
* @param[out] double complex* left_eigenvectors
*
* @return NONE
*/
//================================================================================================//
void compute_eigensystem_dbl(double*,int,unsigned int,double complex*,double complex*,double complex*);


//================================================================================================//
/**
* @brief This function computes the eigensystem for a given symmetric double matrix.
*
* @param[in] double* matrix
* @param[in] int matrix_order
* @param[in] unsigned int sort
* @param[out] double* eigenvalues
* @param[out] double* eigenvectors
*
* @return NONE
*/
//================================================================================================//
void compute_symmetric_eigensystem_dbl(double*,int,unsigned int,double*,double*);


//================================================================================================//
/**
* @brief This function fills a matrix's columns with eigenvectors.
*
* @param[in] double complex* eigenvectors
* @param[in] unsigned int order
* @param[in] unsigned int num_eigenvectors
* @param[out] double complex* eigenvector_matrix
*
* @return NONE
*/
//================================================================================================//
void fill_eigenvector_matrix(double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function prints some eigenvectors to stdout.
*
* @param[in] double complex* eigenvectors
* @param[in] unsigned int order
*
* @return NONE
*/
//================================================================================================//
void print_eigenvectors(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the multiplicities of eigenvalues.
*
* @param[in] double complex* eigenvalues
* @param[in] unsigned int num_eigenvalues
* @param[out] unsigned int multiplicities
*
* @return NONE
*/
//================================================================================================//
void find_eigenvalue_multiplicities(double complex*,unsigned int,unsigned int*);



#endif //LINEAR_ALGEBRA_H//
