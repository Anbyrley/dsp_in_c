/** @file helper.h
*   @brief Contains functions for helping other modules.
*
*
*  @author Alex N. Byrley (anbyrley)
*  @date September 2015
*  @bug Too many calls to malloc  (June 16) -- change so that there is a max size of input arguments
*/

#ifndef HELPER_H
#define HELPER_H

//================================================================================================//
//===================================STANDARD INCLUDES============================================//
//================================================================================================//

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#include <time.h>
#include <complex.h>
#include <limits.h>
#include "macros.h"

//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//

//================================================================================================//
/**
* @brief This function computes the least common multiple of two numbers.
*
* NOTE: This is numerically unstable because a*b can become bigger than double can hold.
*
* @param[in] double a
* @param[in] double b
*
* @return double lcm
*/
//================================================================================================//
double lcm_dbl(double,double);


//================================================================================================//
/**
* @brief This function computes the least common multiple of a double precision array.
*
* NOTE: This is numerically unstable because a*b can become bigger than double can hold.
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return double lcm
*/
//================================================================================================//
double lcm_array_dbl(double*,unsigned int);



//================================================================================================//
/**
* @brief This function computes the greatest common divisor of two numbers.
*
* @param[in] double a
* @param[in] double b
*
* @return double gcd
*/
//================================================================================================//
double gcd_dbl(double, double);


//================================================================================================//
/**
* @brief This function computes the greatest common divisor of values in a double precision array.
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return double gcd
*/
//================================================================================================//
double gcd_array_dbl(double*, unsigned int);


//================================================================================================//
/**
* @brief This function computes the radian value from a degree value
*
* @param[in] double degree
*
* @return radian
*/
//================================================================================================//
double degree_to_radian_dbl(double);


//================================================================================================//
/**
* @brief This function computes the diriac delta function.
*
* @param[in] double x
*
* @return diriac(x)
*/
//================================================================================================//
double diriac_dbl(double);


//================================================================================================//
/**
* @brief This function computes the modulus of a with respect to q according to numpy.
*
* If results in NaN, then returns 1
*
* @param[in] double a
* @param[in] double q
*
* @return (a/q - floor(a/q))*q
*/
//================================================================================================//
double modulo_dbl(double, double);


//================================================================================================//
/**
* @brief This function computes sinc(M_PI * x)/(M_PI * x) for double complex x.
*
* If results in NaN, then returns 1
*
* @param[in] double complex x
*
* @return sinc(x)
*/
//================================================================================================//
double complex sinc_cmplx(double complex);


//================================================================================================//
/**
* @brief This function computes sinc(M_PI * x)/(M_PI * x) for double precision x.
*
* If results in NaN, then returns 1
*
* @param[in] double x
*
* @return sinc(x)
*/
//================================================================================================//
double sinc_dbl(double);


//================================================================================================//
/**
* @brief This function computes the value of the natural log of the gamma function G(x).
*
* This function evaluates a polynomial approximation of G(x) as described in Numerical Recipes 
* in C 2nd Edition. Input must be greater than zero!
*
* @param[in] double xx
*
* @return double G(x)
*/
//================================================================================================//
double gamma_ln(double);


//================================================================================================//
/**
* @brief This function computes the value of the gamma function G(x).
*
* This function calls gamma_ln. If x < 0, uses the relation -x G(-x) G(x) = pi/sin(pi *x).
*
* @param[in] double xx
*
* @return double G(x)
*/
//================================================================================================//
double gamma_dbl(double);


//================================================================================================//
/**
* @brief This function returns the coefficients of a Legendre polynomial.
*
* This function returns them in order from greatest power to least power. Amenable to root finding.
*
* @param[in] unsigned int order
* @param[out] double* coeffs
*
* @return NONE
*/
//================================================================================================//
void get_legendre_coeffs(unsigned int, double*);


//================================================================================================//
/**
* @brief This function computes the value of the modified bessel function I0(x).
*
* This function evaluates a polynomial approximation of I0(x) via Horner's Method found in
* Numerical Recipes in C 2nd Edition.
*
* @param[in] double x
*
* @return double I0(x)
*/
//================================================================================================//
double compute_bessel_I0(double);


//================================================================================================//
/**
* @brief This function computes the value of the bessel function J0(x).
*
* This function evaluates a polynomial approximation of J0(x) via Horner's Method found in
* Numerical Recipes in C 2nd Edition.
*
* @param[in] double x
*
* @return double J0(x)
*/
//================================================================================================//
double compute_bessel_J0(double);


//================================================================================================//
/**
* @brief This function compute the value of the chebyschev polynomial Tn(x).
*
* @param[in] int n
* @param[in] double x
*
* @return double Tn(x)
*/
//================================================================================================//
double compute_chebyshev_polynomial(int,double);


//================================================================================================//
/**
* @brief This function returns the sign of a double precision number.
*
* @param[in] double x
*
* @return double sign
*/
//================================================================================================//
double sign_dbl(double);


//================================================================================================//
/**
* @brief This function computes the required smoothing weight given a needed data history.
*
* @param[in] unsigned int history
*
* @return double alpha
*/
//================================================================================================//
double compute_exponential_smoothing_weight(unsigned int);


//================================================================================================//
/**
* @brief This function computes the first order exponential integral for double precision data.
*
* @param[in] double x
*
* @return double Ei(x)
*/
//================================================================================================//
double ei_dbl(double);


//================================================================================================//
/**
* @brief This function computes the standard exponential integral for double precision data.
*
* @param[in] unsigned int n
* @param[in] double x
*
* @return double En(x)
*/
//================================================================================================//
double expint_dbl(int,double);


//================================================================================================//
/**
* @brief This function pauses with an input to stdout;
*
* @return NONE
*/
//================================================================================================//
void pause_program();


//================================================================================================//
/**
* @brief This function exits with exit(1);
*
* @return NONE
*/
//================================================================================================//
void quit();


//================================================================================================//
/**
* @brief This function determines if a double array contains a value or not.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] double value
*
* @return NONE
*/
//================================================================================================//
unsigned int array_contains_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function determines if a double complex array has NaNs or not.
*
* @param[in] double complex* array
* @param[in] unsigned int length
*
* @return unsigned int has_nans
*/
//================================================================================================//
unsigned int array_has_nans_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function determines if a double precision array has NaNs or not.
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return unsigned int has_nans
*/
//================================================================================================//
unsigned int array_has_nans_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function determines if two double precision values are equal.
*
* @param[in] double val1
* @param[in] double val2
*
* @return unsigned int are_equal
*/
//================================================================================================//
unsigned int are_equal_dbl(double,double);


//================================================================================================//
/**
* @brief This function determines if two double precision arrays are equal.
*
* @param[in] double* array1
* @param[in] double* array2
* @param[in] unsigned int length
*
* @return unsigned int are_equal
*/
//================================================================================================//
unsigned int are_arrays_equal_dbl(double*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function determines if two unsigned int arrays are equal.
*
* @param[in] unsigned int* array1
* @param[in] unsigned int* array2
* @param[in] unsigned int length
*
* @return unsigned int are_equal
*/
//================================================================================================//
unsigned int are_arrays_equal_uint(unsigned int*,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function determines if two char arrays are equal.
*
* @param[in] char* array1
* @param[in] char* array2
* @param[in] unsigned int length
*
* @return unsigned int are_equal
*/
//================================================================================================//
unsigned int are_arrays_equal_char(char*,char*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the length of a string by searching for the null terminator.
*
* @param[in] char* filename
*
* @return unsigned int length
*/
//================================================================================================//
unsigned int get_string_length(char*);


//================================================================================================//
/**
* @brief Given a dB gain, this function returns the linear gain constant needed.
*
* @param[in] double dB_gain
*
* @return double linear_gain
*/
//================================================================================================//
double compute_linear_gain_constant_from_dB(double);


//================================================================================================//
/**
* @brief This function returns the next power of two above the given length.
*
* @param[in] unsigned int length
*
* @return unsigned int next_pow_2
*/
//================================================================================================//
unsigned int next_pow_2(unsigned int);


//================================================================================================//
/**
* @brief This function returns checks if the given complex number is real or not.
*
* @param[in] double complex point
*
* @return unsigned int is_real
*/
//================================================================================================//
unsigned int is_real_cmplx( double complex point );


//================================================================================================//
/**
* @brief This function checks if the given uint array is greater than a value or another array.
*
* @param[in] unsigned int* array
* @param[in] unsigned int array_length
* @param[in] unsigned int* compare
* @param[in] unsigned int compare_length
*
* @return unsigned int igreater
*/
//================================================================================================//
unsigned int array_is_greater_than_uint(unsigned int*,unsigned int,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function checks if the given double array is greater than a value or another array.
*
* @param[in] double* array
* @param[in] unsigned int array_length
* @param[in] double* compare
* @param[in] unsigned int compare_length
*
* @return unsigned int igreater
*/
//================================================================================================//
unsigned int array_is_greater_than_dbl(double*,unsigned int,double*,unsigned int);


//================================================================================================//
/**
* @brief This function checks if the given array contains infinite points or not.
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return unsigned int is_valid
*/
//================================================================================================//
unsigned int is_array_valid_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function returns checks if the given array contains NaN points or not.
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return unsigned int has_NaNs
*/
//================================================================================================//
unsigned int array_has_NaNs_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function returns the index where the running power exceeds a threshold in [0,100].
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] double threshold
*
* @return unsigned int power_index
*/
//================================================================================================//
unsigned int find_power_cutoff_index_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function returns a double array consisting of only the unique values.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[out] double* unique
* @param[out] unsigned int* new_length
*
* @return NONE
*/
//================================================================================================//
void find_unique_dbl(double*,unsigned int,double*,unsigned int*);


//================================================================================================//
/**
* @brief This function returns an unsigned int array consisting of only the unique values.
*
* @param[in] unsigned int* array
* @param[in] unsigned int length
* @param[out] unsigned int* unique
* @param[out] unsigned int* new_length
*
* @return NONE
*/
//================================================================================================//
void find_unique_uint(unsigned int*,unsigned int,unsigned int*,unsigned int*);


//================================================================================================//
/**
* @brief This function returns the index of the closest element in a unsigned int array.
*
* @param[in] unsigned int* array
* @param[in] unsigned int length
* @param[in] unsigned int value
*
* @return unsigned int closest_index
*/
//================================================================================================//
unsigned int find_closest_index_uint(unsigned int*,unsigned int, unsigned int);


//================================================================================================//
/**
* @brief This function returns the index of the closest element in a double precision array.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] double value
*
* @return unsigned int closest_index
*/
//================================================================================================//
unsigned int find_closest_index_dbl(double*,unsigned int,double value);


//================================================================================================//
/**
* @brief This function returns the index where the value should be inserted for it to be sorted.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] double value
*
* @return int best_index
*/
//================================================================================================//
int find_insertion_index_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function returns the indices of array values equal to the value.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] double value
* @param[out] unsigned int* num_indices
* @param[out] unsigned int* indices
*
* @return NONE
*/
//================================================================================================//
void find_indices_where_equal_dbl(double*,unsigned int,double,unsigned int*,unsigned int*);


//================================================================================================//
/**
* @brief This function returns the indices of array values equal to the value.
*
* @param[in] unsigned int* array
* @param[in] unsigned int length
* @param[in] unsigned int value
* @param[out] unsigned int* num_indices
* @param[out] unsigned int* indices
*
* @return NONE
*/
//================================================================================================//
void find_indices_where_equal_uint(unsigned int*,unsigned int,unsigned int,
								   unsigned int*,unsigned int*);


//================================================================================================//
/**
* @brief This function returns the indices of array values greater or equal to the value.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] double value
* @param[out] unsigned int* num_indices
* @param[out] unsigned int* indices
*
* @return NONE
*/
//================================================================================================//
void find_indices_where_greater_dbl(double*,unsigned int,double,unsigned int*,unsigned int*);


//================================================================================================//
/**
* @brief This function returns the indices of array values which are non-zero.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[out] unsigned int* num_indices
* @param[out] unsigned int* indices
*
* @return NONE
*/
//================================================================================================//
void find_indices_where_nonzero_dbl(double*,unsigned int,unsigned int*,unsigned int*);


//================================================================================================//
/**
* @brief This function returns the indices of array values less than or equal to the value.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] double value
* @param[out] unsigned int* num_indices
* @param[out] unsigned int* indices
*
* @return NONE
*/
//================================================================================================//
void find_indices_where_less_dbl(double*,unsigned int,double,unsigned int*,unsigned int*);


//================================================================================================//
/**
* @brief This function inserts an element into an array at a specified point and increments length.
*
* @param[in,out] double* array
* @param[in,out] unsigned int* length
* @param[in] double value
* @param[in] unsigned int index
*
* @return NONE
*/
//================================================================================================//
void insert_array_dbl(double*,unsigned int*,double,unsigned int);


//================================================================================================//
/**
* @brief This function removes an element from an array at a specified point and decrements length.
*
* @param[in,out] double complex* array
* @param[in,out] unsigned int* length
* @param[in] unsigned int index
*
* @return NONE
*/
//================================================================================================//
void remove_from_array_cmplx(double complex*,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function removes an element from an array at a specified point and decrements length.
*
* @param[in,out] double* array
* @param[in,out] unsigned int* length
* @param[in] unsigned int index
*
* @return NONE
*/
//================================================================================================//
void remove_from_array_dbl(double*,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function removes an element from an array at a specified point and decrements length.
*
* @param[in,out] unsigned int* array
* @param[in,out] unsigned int* length
* @param[in] unsigned int index
*
* @return NONE
*/
//================================================================================================//
void remove_from_array_uint(unsigned int*,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function sorts a double complex array by magnitude in ascending order.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void sort_array_by_magnitude_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function sorts a double precision array in ascending order.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void sort_array_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function sorts an unsigned int array in ascending order.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void sort_array_uint(unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function returns the indices of a sorted double complex array in ascending mag order.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
* @param[out] unsigned int indices
*
* @return NONE
*/
//================================================================================================//
void sort_array_indices_by_magnitude_cmplx(double complex*,unsigned int,unsigned int*);


//================================================================================================//
/**
* @brief This function returns the indices of a sorted double precision array in ascending order.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[out] unsigned int indices
*
* @return NONE
*/
//================================================================================================//
void sort_array_indices_dbl(double*,unsigned int,unsigned int*);


//================================================================================================//
/**
* @brief This function zeros out every (order)th array element for a double precision array.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[out] unsigned int order
*
* @return NONE
*/
//================================================================================================//
void downsample_array_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function zeros out every (order)th array element for an unsigned int array.
*
* @param[in] unsigned int* array
* @param[in] unsigned int length
* @param[out] unsigned int order
*
* @return NONE
*/
//================================================================================================//
void downsample_array_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function places (order) zeros in between every array element.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[out] unsigned int order
*
* @return NONE
*/
//================================================================================================//
void upsample_array_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function initializes an array of type double complex.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void initialize_array_cmplx(double complex*, unsigned int);


//================================================================================================//
/**
* @brief This function initializes an array of type double precision.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void initialize_array_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function initializes an array of type int.
*
* @param[in,out] int* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void initialize_array_int(int*,unsigned int);


//================================================================================================//
/**
* @brief This function initializes an array of type unsigned int.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void initialize_array_uint(unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function initializes an array of type char.
*
* @param[in,out] char* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void initialize_array_char(char*,unsigned int);


//================================================================================================//
/**
* @brief This function initializes an array of type unsigned char.
*
* @param[in,out] unsigned char* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void initialize_array_uchar(unsigned char*,unsigned int);


//================================================================================================//
/**
* @brief This function initializes an array of type double complex to a constant value.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
* @param[in] double complex constant
*
* @return NONE
*/
//================================================================================================//
void initialize_array_constant_cmplx(double complex*,unsigned int,double complex);


//================================================================================================//
/**
* @brief This function initializes an array of type double precision to a constant value.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double constant
*
* @return NONE
*/
//================================================================================================//
void initialize_array_constant_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function initializes an array of type unsigned integer to a constant value.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int length
* @param[in] unsigned int constant
*
* @return NONE
*/
//================================================================================================//
void initialize_array_constant_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function copies a char array.
*
* @param[in] char* array
* @param[in] unsigned int length
* @param[out] char* copy
*
* @return NONE
*/
//================================================================================================//
void copy_array_char(char*,unsigned int,char*);


//================================================================================================//
/**
* @brief This function copies an int array.
*
* @param[in] int* array
* @param[in] unsigned int length
* @param[out] int* copy
*
* @return NONE
*/
//================================================================================================//
void copy_array_int(int*,unsigned int,int*);


//================================================================================================//
/**
* @brief This function copies an unsigned int array.
*
* @param[in] unsigned int* array
* @param[in] unsigned int length
* @param[out] unsigned int* copy
*
* @return NONE
*/
//================================================================================================//
void copy_array_uint(unsigned int*, unsigned int,unsigned int*);


//================================================================================================//
/**
* @brief This function copies a double precision array.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[out] double* copy
*
* @return NONE
*/
//================================================================================================//
void copy_array_dbl(double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function copies a double precision array.
*
* @param[in] double complex* array
* @param[in] unsigned int length
* @param[out] double complex* copy
*
* @return NONE
*/
//================================================================================================//
void copy_array_cmplx(double complex*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function combines two double precision arrays into a complex array.
*
* @param[in] double* real
* @param[in] double* imag
* @param[in] unsigned int length
* @param[in] double complex* out
*
* @return NONE
*/
//================================================================================================//
void combine_arrays_dbl_to_cmplx(double*,double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function splits a double complex array into its component real and imaginary parts.
*
* @param[in] double complex* array
* @param[in] unsigned int length
* @param[out] double* real
* @param[out] double* imag
*
* @return NONE
*/
//================================================================================================//
void split_array_cmplx(double complex*,unsigned int,double*,double*);


//================================================================================================//
/**
* @brief This function swaps the real and imaginary parts of a a double complex array.
*
* @param[in] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void swap_real_imag_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function zeros the real part of a a double complex array.
*
* @param[in] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void zero_array_real_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function zeros the imaginary part of a a double complex array.
*
* @param[in] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void zero_array_imag_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function returns an array where all the imaginary values are replaced by abs values.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void abs_imaginary_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function returns an array where all the values are replaced by their abs values.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void abs_array_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function returns an array where all the values are replaced by their abs values.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void abs_array_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function returns an array where all the values are replaced by their powers.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
* @param[in] double complex power
*
* @return NONE
*/
//================================================================================================//
void pow_array_cmplx(double complex*,unsigned int,double complex);


//================================================================================================//
/**
* @brief This function returns an array where all the values are replaced by their powers.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double power
*
* @return NONE
*/
//================================================================================================//
void pow_array_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function copies a double complex matrix.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double complex* copy
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_cmplx(double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function copies a double precision matrix.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double* copy
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_dbl(double*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function copies an unsigned int matrix.
*
* @param[in] unsigned int* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] unsigned int* copy
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_uint(unsigned int*,unsigned int,unsigned int,unsigned int*);


//================================================================================================//
/**
* @brief This function normalizes the angle of a double complex array.
*
* @param[in,out] double complex* data
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_angle_cmplx(double complex*, unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the 2norm of a double complex vector.
*
* @param[in,out] double complex* vector
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_vector_2norm_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the 2norm of a double precision vector.
*
* @param[in,out] double* vector
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_vector_2norm_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the maximum value.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_max_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the values of the data to be within [0,1].
*
* @param[in,out] double* data
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_range_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the absolute maximum value.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_abs_max_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the maximum magnitude in a complex array.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_magnitude_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the sum.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_sum_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the power.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_power_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the variance.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_variance_dbl(double*, unsigned int);

//================================================================================================//
/**
* @brief This function makes the vector monic.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void make_array_monic_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function makes the vector monic.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void make_array_monic_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function reverses an unsigned int array in place.
*
* @param[in,out] unsigned int * array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void reverse_array_uint(unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function reverses the vector in place.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void reverse_array_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function reverses the vector in place.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void reverse_array_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function reorders an unsigned int vector in place.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int* new_indices
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void reorder_array_uint(unsigned int*,unsigned int*, unsigned int);


//================================================================================================//
/**
* @brief This function reorders a double precision vector in place.
*
* @param[in,out] double* array
* @param[in] unsigned int* new_indices
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void reorder_array_dbl(double*,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function reorders a double complex vector in place.
*
* @param[in,out] double complex* array
* @param[in] unsigned int* new_indices
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void reorder_array_cmplx(double complex*,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function extracts the real and/or imaginary parts of a complex array.
*
* @param[in] double complex* array
* @param[in] unsigned int length
* @param[out] double* real
* @param[out] double* imag
*
* @return NONE
*/
//================================================================================================//
void extract_array_coordinates_cmplx(double complex*, unsigned int, double*, double*);


//================================================================================================//
/**
* @brief This function fills an array with a uniform spacing of numbers between start and end.
*
* @param[in,out] double linspace
* @param[in] unsigned int length
* @param[in] double start
* @param[in] double end
*
* @return NONE
*/
//================================================================================================//
void initialize_linspace_dbl(double*,unsigned int,double,double);


//================================================================================================//
/**
* @brief This function fills an array with a uniform spacing of numbers between start and end.
*
* @param[in,out] unsigned int* linspace
* @param[in] unsigned int length
* @param[in] unsigned int start
* @param[in] unsigned int end
*
* @return NONE
*/
//================================================================================================//
void initialize_linspace_uint(unsigned int*,unsigned int,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function fills an array with a uniform spacing of numbers between start and end.
*
* @param[in,out] int* linspace
* @param[in] int length
* @param[in] int start
* @param[in] int end
*
* @return NONE
*/
//================================================================================================//
void initialize_linspace_int(int*,int,int,int);


//================================================================================================//
/**
* @brief This function fills an array with a logarithmic spacing of numbers between start and end.
*
* @param[in,out] double* logspace
* @param[in] unsigned int length
* @param[in] double start
* @param[in] double end
*
* @return NONE
*/
//================================================================================================//
void initialize_logspace_dbl(double*,unsigned int,double,double);


//================================================================================================//
/**
* @brief This function fills an array with a range of unsigned ints from start to end.
*
* @param[in,out] unsigned int* range
* @param[in] unsigned int start
* @param[in] unsigned int end
*
* @return NONE
*/
//================================================================================================//
void initialize_range_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function fills an array with a range of doubles from start to end in steps of step.
*
* @param[in,out] unsigned int* array
* @param[in] double start
* @param[in] double end
* @param[in] double step
* @param[out] unsigned int* length
*
* @return NONE
*/
//================================================================================================//
void initialize_range_dbl(double*,double,double,double,unsigned int*);


//================================================================================================//
/**
* @brief This function pads a double complex array with zeros at the end.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
* @param[in] unsigned int num_zeros
*
* @return NONE
*/
//================================================================================================//
void pad_zeros_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function pads a double precision array with zeros at the end.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] unsigned int num_zeros
*
* @return NONE
*/
//================================================================================================//
void pad_zeros_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function multiples a double complex array by a double complex constant.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
* @param[in] double complex constant
*
* @return NONE
*/
//================================================================================================//
void gain_array_constant_cmplx(double complex*,unsigned int,double complex);


//================================================================================================//
/**
* @brief This function multiples a double precision array by a double precision constant.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double constant
*
* @return NONE
*/
//================================================================================================//
void gain_array_constant_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function multiples an unsigned int array by an unsigned int constant.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int length
* @param[in] unsigned int constant
*
* @return NONE
*/
//================================================================================================//
void gain_array_constant_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function adds a constant to an unsigned int array.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int length
* @param[in] unsigned int constant
*
* @return NONE
*/
//================================================================================================//
void add_array_constant_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function adds a constant to an int array.
*
* @param[in,out] int* array
* @param[in] int length
* @param[in] int constant
*
* @return NONE
*/
//================================================================================================//
void add_array_constant_int(int*,int,int);


//================================================================================================//
/**
* @brief This function appends a value to an unsigned int array, as long as it not in there already.
*
* @param[in,out] unsigned int* array
* @param[in,out] unsigned int* array_length
* @param[in] unsigned int value
*
* @return NONE
*/
//================================================================================================//
void append_array_unique_uint(unsigned int*,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function appends values to an array.
*
* @param[in,out] unsigned int* array
* @param[in,out] unsigned int* array_length
* @param[in] double* appendix
* @param[in] unsigned int appendix_length
*
* @return NONE
*/
//================================================================================================//
void append_array_dbl(double*,unsigned int*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function adds a constant to a double precision array.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double constant
*
* @return NONE
*/
//================================================================================================//
void add_array_constant_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function scales an array to the range [start, end].
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double start
* @param[in] double end
*
* @return NONE
*/
//================================================================================================//
void scale_array_dbl(double*,unsigned int,double,double);


//================================================================================================//
/**
* @brief This function right shifts a double precision array by shift_length.
*
* Data is shifted out of the right side and placed back into the left side.
*
* @param[in,out] double* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void right_circular_shift_array_dbl(double*, unsigned int, unsigned int);


//================================================================================================//
/**
* @brief This function left shifts a double complex array by shift_length.
*
* Data is shifted out of the left side and placed back into the right side.
*
* @param[in,out] double complex* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void left_circular_shift_array_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function left shifts a double precision array by shift_length.
*
* Data is shifted out of the left side and placed back into the right side.
*
* @param[in,out] double* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void left_circular_shift_array_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function left shifts an unsigned int array by shift_length.
*
* Data is shifted out of the left side and placed back into the right side.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void left_circular_shift_array_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function right shifts an unsigned int array by delay_length.
*
* Data is shifted out of the right side and lost forever.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int array_length
* @param[in] unsigned int delay_length
*
* @return NONE
*/
//================================================================================================//
void right_shift_array_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function right shifts a double precision array by delay_length.
*
* Data is shifted out of the right side and lost forever.
*
* @param[in,out] double* array
* @param[in] unsigned int array_length
* @param[in] unsigned int delay_length
*
* @return NONE
*/
//================================================================================================//
void right_shift_array_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function left shifts an unsigned int array by shift_length.
*
* Data is shifted out of the left side and lost forever.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void left_shift_array_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function left shifts a double precision array by shift_length.
*
* Data is shifted out of the left side and lost forever.
*
* @param[in,out] double* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void left_shift_array_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function left shifts a double complex array by shift_length.
*
* Data is shifted out of the left side and lost forever.
*
* @param[in,out] double complex* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void left_shift_array_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function removes an element from an unsigned int array.
*
* @param[in,out] unsigned int* array
* @param[out] unsigned int* array_length
* @param[in] unsigned int element
*
* @return NONE
*/
//================================================================================================//
void remove_array_element_uint(unsigned int*,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function repeats the array elements back to back multiple times.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] unsigned int num_repeats
*
* @return NONE
*/
//================================================================================================//
void repeat_array_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function shifts the columns of a double cmplex array to the left.
*
* Data is shifted out of the left side and lost forever.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int num_shifts
*
* @return NONE
*/
//================================================================================================//
void left_shift_matrix_columns_cmplx(double complex*,unsigned int,unsigned int, unsigned int);


//================================================================================================//
/**
* @brief This function removes a column from a double complex matrix and remaps the memory.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in,out] unsigned int* num_cols
* @param[in] unsigned int column_number
*
* @return NONE
*/
//================================================================================================//
void remove_matrix_column_cmplx(double complex*,unsigned int,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function removes a column from a double precision matrix and remaps the memory.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in,out] unsigned int* num_cols
* @param[in] unsigned int column_number
*
* @return NONE
*/
//================================================================================================//
void remove_matrix_column_dbl(double*,unsigned int,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function removes a row from a double precision matrix and remaps the memory.
*
* @param[in,out] double* matrix
* @param[in,out] unsigned int* num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int row_number
*
* @return NONE
*/
//================================================================================================//
void remove_matrix_row_dbl(double*,unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function right shifts a double cmplex array by shift_length.
*
* Data is shifted out of the right side and placed back into the left side.
*
* @param[in,out] double complex* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void right_circular_shift_array_cmplx(double complex*, unsigned int, unsigned int);


//================================================================================================//
/**
* @brief This function right shifts an unsigned int array by shift_length.
*
* Data is shifted out of the right side and placed back into the left side.
*
* @param[in,out] unsigned int* array
* @param[in] unsigned int array_length
* @param[in] unsigned int shift_length
*
* @return NONE
*/
//================================================================================================//
void right_circular_shift_array_uint(unsigned int*,unsigned int, unsigned int);


//================================================================================================//
/**
* @brief This function right shifts a double complex array by delay_length.
*
* Data is shifted out of the right side and lost forever.
*
* @param[in,out] double complex* array
* @param[in] unsigned int array_length
* @param[in] unsigned int delay_length
*
* @return NONE
*/
//================================================================================================//
void right_shift_array_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function compresses the zeros out of an double precision array.
*
* @param[in,out] double array
* @param[in] unsigned int array_length
* @param[out] unsigned int new_length
*
* @return NONE
*/
//================================================================================================//
void eliminate_array_zeros_dbl(double*,unsigned int,unsigned int*);


//================================================================================================//
/**
* @brief This function compresses the zeros out of an unsigned int array.
*
* @param[in,out] unsigned int array
* @param[in] unsigned int array_length
* @param[out] unsigned int new_length
*
* @return NONE
*/
//================================================================================================//
void eliminate_array_zeros_uint(unsigned int*,unsigned int,unsigned int*);


//================================================================================================//
/**
* @brief This function compresses the infs out of a double precision array.
*
* @param[in,out] double array
* @param[in] unsigned int array_length
*
* @return NONE
*/
//================================================================================================//
void eliminate_array_infs_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function appends a column to a matrix and updates the number of columns.
*
* @param[in,out] double matrix
* @param[in] unsigned int num_rows
* @param[in,out] unsigned int* num_cols
* @param[in] double* column
*
* @return NONE
*/
//================================================================================================//
void append_matrix_column_dbl(double*,unsigned int,unsigned int*,double*);


//================================================================================================//
/**
* @brief This function inserts a column into a matrix and updates the number of columns.
*
* @param[in,out] double matrix
* @param[in] unsigned int num_rows
* @param[in,out] unsigned int* num_cols
* @param[in] double* column
* @param[in] unsigned int column_number
*
* @return NONE
*/
//================================================================================================//
void insert_matrix_column_dbl(double*,unsigned int,unsigned int*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the max difference between all elements of a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return double maximum_incremement
*/
//================================================================================================//
double find_maximum_pairwise_difference_dbl(double*, unsigned int);


//================================================================================================//
/**
* @brief This function finds the abs min difference between all elements of a double array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return double abs_minimum_incremement
*/
//================================================================================================//
double find_minimum_abs_pairwise_difference_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum magnitude of a double complex array.
*
* @param[in] double complex* data
* @param[in] unsigned int length
*
* @return double maximum_magnitude
*/
//================================================================================================//
double find_maximum_magnitude_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum angle of a double complex array.
*
* @param[in] double complex* data
* @param[in] unsigned int length
*
* @return double maximum_angle
*/
//================================================================================================//
double find_maximum_angle_cmplx(double complex*, unsigned int);


//================================================================================================//
/**
* @brief This function finds the sum of a double complex array.
*
* @param[in] double complex* data
* @param[in] unsigned int length
*
* @return double complex sum
*/
//================================================================================================//
double complex compute_sum_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the sum of a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return double sum
*/
//================================================================================================//
double compute_sum_dbl(double*, unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum of a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return double maximum
*/
//================================================================================================//
double find_maximum_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum absolute value of a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return double abs_maximum
*/
//================================================================================================//
double find_abs_maximum_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum of an unsigned int array.
*
* @param[in] unsigned int* data
* @param[in] unsigned int length
*
* @return unsigned int maximum
*/
//================================================================================================//
unsigned int find_maximum_uint(unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the median of a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return double median
*/
//================================================================================================//
double find_median_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the median of a double precision array using malloc.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return double median
*/
//================================================================================================//
double find_median_malloc_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds all the local maxima of a double precision array.
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] double* threshold
* @param[out] unsigned int* num_maxima
* @param[out] unsigned int* positions
* @param[out] double* maxima
*
* @return NONE
*/
//================================================================================================//
void find_local_maxima_dbl(double*,unsigned int,double*,unsigned int*,unsigned int*,double*);


//================================================================================================//
/**
* @brief This function finds the index of the maximum of a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return unsigned int index
*/
//================================================================================================//
unsigned int find_maximum_index_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the minimum of a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return double minimum
*/
//================================================================================================//
double find_minimum_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the minimum of an unsigned int array.
*
* @param[in] unsigned int* data
* @param[in] unsigned int length
*
* @return unsigned int
*/
//================================================================================================//
unsigned int find_minimum_uint(unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the index of the minimum of a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
*
* @return unsigned int index
*/
//================================================================================================//
unsigned int find_minimum_index_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the index of the maximum magnitude element in a double complex array.
*
* @param[in] double complex* data
* @param[in] unsigned int length
*
* @return unsigned int index
*/
//================================================================================================//
unsigned int find_maximum_magnitude_index_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the index of the minimum magnitude element in a double complex array.
*
* @param[in] double complex* data
* @param[in] unsigned int length
*
* @return unsigned int index
*/
//================================================================================================//
unsigned int find_minimum_magnitude_index_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function initializes a matrix of type unsigned char.
*
* @param[in,out] unsigned char* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void initialize_matrix_uchar(unsigned char*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function initializes a matrix of type double.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void initialize_matrix_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function initializes a matrix of type double complex.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void initialize_matrix_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum value in a matrix of type unsigned int.
*
* @param[in] unsigned int* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return unsigned int maximum
*/
//================================================================================================//
unsigned int find_matrix_maximum_uint(unsigned int*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum value in a matrix of type double.
*
* @param[in] unsigned int* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return double maximum
*/
//================================================================================================//
double find_matrix_maximum_dbl(double*,unsigned int,unsigned int);

//================================================================================================//
/**
* @brief This function converts a Q1.15 value to an double precision value.
*
* @param[in] int16_t x
*
* @return double y
*/
//================================================================================================//
double q115_to_dbl(int16_t);


//================================================================================================//
/**
* @brief This function converts a double precision value to an Q1.15 value.
*
* @param[in] double x
*
* @return int16_t y
*/
//================================================================================================//
int16_t dbl_to_q115(double);


//================================================================================================//
/**
* @brief This function computes the magnitude of a double complex array.
*
* @param[in] double complex* array
* @param[in] unsigned int length
* @param[out] double* magnitude
*
* @return NONE
*/
//================================================================================================//
void get_magnitude_array_cmplx(double complex*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes the angles around the unit circle of a double complex array.
*
* @param[in] double complex* array
* @param[in] unsigned int length
* @param[out] double* angles
*
* @return NONE
*/
//================================================================================================//
void get_angle_array_cmplx(double complex*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function converts an unsigned int vector to double precision.
*
* @param[in] unsigned int* vector
* @param[in] unsigned int length
* @param[out] double* dbl_vector
*
* @return NONE
*/
//================================================================================================//
void convert_vector_uint_to_dbl(unsigned int*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function converts a matrix of type unsigned int to unsigned char.
*
* @param[in] unsigned int* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] unsigned char* char_matrix
*
* @return NONE
*/
//================================================================================================//
void convert_matrix_uint_to_uchar(unsigned int*,unsigned int,unsigned int,unsigned char*);


//================================================================================================//
/**
* @brief This function conjugates a double complex array.
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void conjugate_array_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function returns the binomial coefficient when n and k are unsigned integers.
*
* @param[in] unsigned int n
* @param[in] unsigned int k
*
* @return unsigned int n_choose_k
*/
//================================================================================================//
unsigned int n_choose_k_uint(unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function returns the binomial coefficient when n and k are double precision.
*
* NOTE: Uses the relation that (n k) = gamma(n+1)/(gamma(k+1)*gamma(n-k+1))
*
* @param[in] double n
* @param[in] double k
*
* @return double n_choose_k
*/
//================================================================================================//
double n_choose_k_dbl(double,double);


//================================================================================================//
/**
* @brief This function computes the arithmetic mean of a double complex vector array
*
* @param[in] double complex* array
* @param[in] unsigned int length
*
* @return double mean
*/
//================================================================================================//
double complex compute_mean_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the arithmetic mean of a double precision vector array
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return double mean
*/
//================================================================================================//
double compute_mean_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the arithmetic mean of an unsigned int vector array
*
* @param[in] unsigned int* array
* @param[in] unsigned int length
*
* @return unsigned int mean
*/
//================================================================================================//
unsigned int compute_mean_uint(unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the geometric mean of a double precision vector array
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return double mean
*/
//================================================================================================//
double compute_geometric_mean_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the arithmetic mean of a double precision vector via recursion
*
* This implements: 	xk+1 = xk + 1/(k+1) * (yk+1 - xk)
*
* @param[in] double* array
* @param[in] unsigned int length
* @param[in] unsigned int iteration
* @param[in,out] double* mean
*
* @return NONE
*/
//================================================================================================//
void recurse_mean_dbl(double*, unsigned int, unsigned int, double*);


//================================================================================================//
/**
* @brief This function subtracts out the arithmetic mean of a vector array
*
* @param[in,out] double complex* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void remove_mean_cmplx(double complex*,unsigned int length );


//================================================================================================//
/**
* @brief This function subtracts out the arithmetic mean of a vector array
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void remove_mean_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function adds a constant to every element of a double precision array.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double constant
*
* @return NONE
*/
//================================================================================================//
void add_constant_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function computes the variance of a vector array
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return double variance
*/
//================================================================================================//
double compute_variance_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the sample skewness of a double precision vector array
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return double sample_skewness
*/
//================================================================================================//
double compute_sample_skewness_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the sample kurtosis of a double precision vector array
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return double sample_kurtosis
*/
//================================================================================================//
double compute_sample_kurtosis_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the arithmetic mean of a double complex vector via recursion
*
* This implements: 	xk+1 = xk + 1/(k+1) * (yk+1 - xk)
*
* @param[in] double complex* array
* @param[in] unsigned int length
* @param[in] unsigned int iteration
* @param[in,out] double complex* mean
*
* @return NONE
*/
//================================================================================================//
void recurse_mean_cmplx(double complex*, unsigned int, unsigned int, double complex*);


//================================================================================================//
/**
* @brief This function normalizes a vector to have zero mean and unit variance.
*
* @param[in,out] double* array
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void normalize_standard_dbl(double*, unsigned int);


//================================================================================================//
/**
* @brief This function thresholds a double precision array.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double threshold
*
* @return NONE
*/
//================================================================================================//
void threshold_array_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function thresholds a double precision array.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double threshold
*
* @return NONE
*/
//================================================================================================//
void threshold_array_abs_dbl(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function clips a double precision array to within an interval.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] double minimum
* @param[in] double maximum
*
* @return NONE
*/
//================================================================================================//
void clip_array_dbl(double*,unsigned int,double,double);


//================================================================================================//
/**
* @brief This function computes the mean squared error between two integer vectors.
*
* @param[in] int* vector1
* @param[in] int* vector2
* @param[in] unsigned int length
*
* @return int error
*/
//================================================================================================//
int compute_mean_error_int(int*,int*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the mean squared error between two double complex vectors.
*
* @param[in] double complex* vector1
* @param[in] double complex* vector2
* @param[in] unsigned int length
*
* @return double error
*/
//================================================================================================//
double complex compute_mean_squared_error_cmplx(double complex*,double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the mean squared error between two double precision vectors.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
*
* @return double error
*/
//================================================================================================//
double compute_mean_squared_error_dbl(double*, double*, unsigned int);


//================================================================================================//
/**
* @brief This function computes the sum of the absolute value of the differences between two
* double precision vectors.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
*
* @return double distance
*/
//================================================================================================//
double compute_manhattan_distance_dbl(double*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the sum of the absolute value of the differences between two
* double precision vectors raised to a power, then the sum is rooted.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
* @param[in] double power
*
* @return double distance
*/
//================================================================================================//
double compute_minkowski_distance_dbl(double*,double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function computes the maximum value of the absolute distance between the elements.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
*
* @return double distance
*/
//================================================================================================//
double compute_chebyshev_distance_dbl(double*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the maximum error between two double complex vectors.
*
* @param[in] double complex* vector1
* @param[in] double complex* vector2
* @param[in] unsigned int length
*
* @return double error
*/
//================================================================================================//
double compute_maximum_error_cmplx(double complex*,double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the maximum error between two double precision vectors.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
*
* @return double error
*/
//================================================================================================//
double compute_maximum_error_dbl(double*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the power in a double complex array.
*
* @param[in] double complex* array
* @param[in] unsigned int length
*
* @return double power
*/
//================================================================================================//
double compute_power_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the power in a double precision array.
*
* @param[in] double* array
* @param[in] unsigned int length
*
* @return double power
*/
//================================================================================================//
double compute_power_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the 2norm of a double complex array.
*
* @param[in] double complex* vector
* @param[in] unsigned int length
*
* @return double 2norm
*/
//================================================================================================//
double compute_vector_2norm_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the 2norm of a double precision array.
*
* @param[in] double* vector
* @param[in] unsigned int length
*
* @return double 2norm
*/
//================================================================================================//
double compute_vector_2norm_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes the 1norm of a double precision array.
*
* @param[in] double* vector
* @param[in] unsigned int length
*
* @return double 1norm
*/
//================================================================================================//
double compute_vector_1norm_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum absolute column sum of a double complex matrix.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return double 1norm
*/
//================================================================================================//
double compute_matrix_1norm_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum absolute column sum of a double precision matrix.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return double 1norm
*/
//================================================================================================//
double compute_matrix_1norm_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function finds the maximum absolute row sum of a double precision matrix.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return double infnorm
*/
//================================================================================================//
double compute_matrix_infnorm_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function compute the frobenius norm for a double precision matrix.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return double fnorm
*/
//================================================================================================//
double compute_matrix_fnorm_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function computes result = alpha * array + (1-alpha) * update for alpha in [0,1].
*
* @param[in] double* array
* @param[in] double* update
* @param[in] unsigned int length
* @param[in] double alpha
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void exponentially_smooth_array_dbl(double*,double*,unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function computes the hadamard product for double complex vectors.
*
* @param[in] double complex* vector1
* @param[in] double complex* vector2
* @param[in] unsigned int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void hadamard_product_cmplx(double complex*,double complex*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function computes the hadamard product for double precision vectors.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void hadamard_product_dbl(double*,double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes and returns the inner product for two double complex vectors.
*
* @param[in] double complex* vector1
* @param[in] double complex* vector2
* @param[in] unsigned int length
*
* @return double complex inner_product
*/
//================================================================================================//
double complex compute_inner_product_cmplx(double complex*,double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function computes and returns the inner product for two double precision vectors.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
*
* @return double inner_product
*/
//================================================================================================//
double compute_inner_product_dbl(double*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function computes and returns the p-distance for two double complex vectors.
*
* @param[in] double complex* vector1
* @param[in] double complex* vector2
* @param[in] unsigned int length
* @param[in] double p
*
* @return double distance
*/
//================================================================================================//
double compute_vector_distance_cmplx(double complex*,double complex*,unsigned int,double);


//================================================================================================//
/**
* @brief This function computes and returns the p-distance for two double precision vectors.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
* @param[in] double p
*
* @return double distance
*/
//================================================================================================//
double compute_vector_distance_dbl(double*,double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function subtracts two double complex vectors.
*
* NOTE: Order is vector1 - vector2
*
* @param[in] double complex* vector1
* @param[in] double complex* vector2
* @param[in] unsigned int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void subtract_vectors_cmplx(double complex*,double complex*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function subtracts two double precision vectors.
*
* NOTE: Order is vector1 - vector2
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void subtract_vectors_dbl(double*,double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function adds two double complex vectors.
*
* @param[in] double complex* vector1
* @param[in] double complex* vector2
* @param[in] unsigned int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void add_vectors_cmplx(double complex*,double complex*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function adds two double precision vectors.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void add_vectors_dbl(double*,double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function adds two double complex matrices.
*
* @param[in] double complex* matrix1
* @param[in] double complex* matrix2
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void add_matrices_cmplx(double complex*,double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function subtracts two double complex matrices.
*
* @param[in] double complex* matrix1
* @param[in] double complex* matrix2
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void subtract_matrices_cmplx(double complex*,double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function adds two double precision matrices.
*
* @param[in] double* matrix1
* @param[in] double* matrix2
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void add_matrices_dbl(double*,double*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function negates a double precision vector.
*
* @param[in,out] double* vector
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void negate_vector_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function inverts the entries of a double complex vector.
*
* @param[in,out] double complex* vector
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void invert_vector_cmplx(double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function inverts the entries of a double precision vector.
*
* @param[in,out] double* vector
* @param[in] unsigned int length
*
* @return NONE
*/
//================================================================================================//
void invert_vector_dbl(double*,unsigned int);


//================================================================================================//
/**
* @brief This function divides the entries of two double precision vectors.
*
* @param[in] double* vector1
* @param[in] double* vector2
* @param[in] unsigned int length
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void divide_vectors_dbl(double*,double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function sets the diagonal elements of a double complex matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_cols
* @param[in] double complex* elements
* @param[in] unsigned int num_elements
*
* @return NONE
*/
//================================================================================================//
void set_matrix_diagonal_cmplx(double complex*,unsigned int,double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function sets the diagonal elements of a matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_cols
* @param[in] double* elements
* @param[in] unsigned int num_elements
*
* @return NONE
*/
//================================================================================================//
void set_matrix_diagonal_dbl(double*,unsigned int,double*,unsigned int);


//================================================================================================//
/**
* @brief This function sets the off-diagonal elements of a matrix.
*
* NOTE: The off diagonal to set is indicated by +-index from the main diagonal.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double* elements
* @param[in] unsigned int num_elements
* @param[in] int off_diagonal_index
*
* @return NONE
*/
//================================================================================================//
void set_matrix_offdiagonal_dbl(double*,unsigned int,unsigned int,double*,unsigned int,int);


//================================================================================================//
/**
* @brief This function copies the diagonal elements of a double complex matrix.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_cols
* @param[in] unsigned int num_rows
* @param[out] double complex* diagonal
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_diagonal_cmplx(double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function copies the diagonal elements of a double precision matrix.
*
* @param[in] double* matrix
* @param[in] unsigned int num_cols
* @param[in] unsigned int num_rows
* @param[out] double* diagonal
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_diagonal_dbl(double*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function copies the off-diagonal elements of a double precision matrix.
*
* NOTE: The off diagonal to set is indicated by +-index from the main diagonal.
*
* @param[in] double* matrix
* @param[in] unsigned int num_cols
* @param[in] unsigned int num_rows
* @param[in] int off_diagonal_index
* @param[out] double* diagonal
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_offdiagonal_dbl(double*,unsigned int,unsigned int,int,double*);


//================================================================================================//
/**
* @brief This function copies a row of a double complex matrix.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_cols
* @param[in] unsigned int row_number
* @param[out] double complex* row
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_row_cmplx(double complex*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function copies a row of a double precision matrix.
*
* @param[in] double* matrix
* @param[in] unsigned int num_cols
* @param[in] unsigned int row_number
* @param[out] double* row
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_row_dbl(double*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function copies a row of an unsigned int matrix.
*
* @param[in] unsigned int* matrix
* @param[in] unsigned int num_cols
* @param[in] unsigned int row_number
* @param[out] unsigned int* row
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_row_uint(unsigned int*,unsigned int,unsigned int,unsigned int*);


//================================================================================================//
/**
* @brief This function copies a column of a double complex matrix.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int column_number
* @param[out] double complex* column
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_column_cmplx(double complex*,unsigned int,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function copies a column of a double precision matrix.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int column_number
* @param[out] double* column
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_column_dbl(double*,unsigned int,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function copies a column of an unsigned int matrix.
*
* @param[in] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int column_number
* @param[out] double* column
*
* @return NONE
*/
//================================================================================================//
void copy_matrix_column_uint(unsigned int*,unsigned int,unsigned int,unsigned int,unsigned int*);


//================================================================================================//
/**
* @brief This function replaces a column of a double complex matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double complex* column
* @param[in] unsigned int column_number
*
* @return NONE
*/
//================================================================================================//
void replace_matrix_column_cmplx(double complex*,unsigned int,unsigned int,double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function replaces a column of a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] double* column
* @param[in] unsigned int column_number
*
* @return NONE
*/
//================================================================================================//
void replace_matrix_column_dbl(double*,unsigned int,unsigned int,double*,unsigned int);


//================================================================================================//
/**
* @brief This function replaces a row of a double complex matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_cols
* @param[in] double complex* row
* @param[in] unsigned int row_number
*
* @return NONE
*/
//================================================================================================//
void replace_matrix_row_cmplx(double complex*,unsigned int,double complex*,unsigned int);


//================================================================================================//
/**
* @brief This function replaces a row of a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_cols
* @param[in] double* row
* @param[in] unsigned int row_number
*
* @return NONE
*/
//================================================================================================//
void replace_matrix_row_dbl(double*,unsigned int,double*,unsigned int);


//================================================================================================//
/**
* @brief This function replaces a row of an unsigned int matrix.
*
* @param[in,out] unsigned int* matrix
* @param[in] unsigned int num_cols
* @param[in] unsigned int* row
* @param[in] unsigned int row_number
*
* @return NONE
*/
//================================================================================================//
void replace_matrix_row_uint(unsigned int*,unsigned int,unsigned int*,unsigned int);


//================================================================================================//
/**
* @brief This function partitions a double complex matrix into two matrices at the column break.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int column_break
* @param[out] double complex* block1
* @param[out] double complex* block2
*
* @return NONE
*/
//================================================================================================//
void partition_matrix_columns_cmplx(double complex*,unsigned int,unsigned int,unsigned int,
									double complex*,double complex*);


//================================================================================================//
/**
* @brief This function partitions a double complex matrix into two matrices at the row break.
*
* @param[in] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int row_break
* @param[out] double complex* block1
* @param[out] double complex* block2
*
* @return NONE
*/
//================================================================================================//
void partition_matrix_rows_cmplx(double complex*,unsigned int,unsigned int,unsigned int,
								 double complex*,double complex*);


//================================================================================================//
/**
* @brief This function transposes a double complex matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void transpose_matrix_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function transposes a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void transpose_matrix_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function transposes a double complex matrix and then conjugates it.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void hermitian_transpose_matrix_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function appends rows and columns of zeros to a double complex matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int extra_rows
* @param[in] unsigned int extra_cols
*
* @return NONE
*/
//================================================================================================//
void append_zeros_matrix_cmplx(double complex*,unsigned int,unsigned int,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function appends rows and columns of zeros to a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
* @param[in] unsigned int extra_rows
* @param[in] unsigned int extra_cols
*
* @return NONE
*/
//================================================================================================//
void append_zeros_matrix_dbl(double*,unsigned int,unsigned int,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the columns of a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void normalize_matrix_columns_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the rows of a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void normalize_matrix_rows_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the sum of the rows of a double complex matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void normalize_matrix_rows_sum_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the sum of the rows of a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void normalize_matrix_rows_sum_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the frobenius norm of a double precision matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void normalize_matrix_frobenius_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the maximum absolute column sum of a double complex matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void normalize_matrix_1norm_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function normalizes the 2norm of the rows of a double precision matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void normalize_matrix_rows_2norm_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function creates a double complex identity matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void initialize_identity_matrix_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function creates a double complex exchange matrix.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void initialize_exchange_matrix_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function creates a double precision identity matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void initialize_identity_matrix_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function creates a double precision exchange matrix.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void initialize_exchange_matrix_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function creates a double precision circulant matrix from data.
*
* NOTE: Can be used to create the autocorrelation matrix from an autocorrelation array.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double* data
*
* @return NONE
*/
//================================================================================================//
void form_circulant_matrix_dbl(double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function creates a unit vector.
*
* @param[in,out] double* array
* @param[in] unsigned int length
* @param[in] unsigned int coordinate
*
* @return NONE
*/
//================================================================================================//
void initialize_unit_vector_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function prints a newline to stdout
*
* @return NONE
*/
//================================================================================================//
void newline();


//================================================================================================//
/**
* @brief This function prints a Hello World to stdout
*
* @return NONE
*/
//================================================================================================//
void hello_world();


//================================================================================================//
/**
* @brief This function prints a uint row vector.
*
* @param[in] unsigned int* vector
* @param[in] unsigned int length
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_row_vector_uint(unsigned int*,unsigned int,FILE*);


//================================================================================================//
/**
* @brief This function prints a uint vector.
*
* @param[in] unsigned int* vector
* @param[in] unsigned int length
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_vector_uint( unsigned int*, unsigned int, FILE*);


//================================================================================================//
/**
* @brief This function prints a int row vector.
*
* @param[in] int* vector
* @param[in] int length
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_row_vector_int(int* vector, int vector_size, FILE*);


//================================================================================================//
/**
* @brief This function prints an int vector.
*
* @param[in] int* vector
* @param[in] int length
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_vector_int(int*, unsigned int, FILE*);


//================================================================================================//
/**
* @brief This function prints a row vector.
*
* @param[in] double* vector
* @param[in] unsigned int length
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_row_vector_dbl(double*,unsigned int,FILE*);


//================================================================================================//
/**
* @brief This function prints a row vector double complex.
*
* @param[in] double complex* vector
* @param[in] unsigned int length
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_row_vector_cmplx(double complex*,unsigned int,FILE*);


//================================================================================================//
/**
* @brief This function prints a vector.
*
* @param[in] double* vector
* @param[in] unsigned int length
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_vector_dbl(double*,unsigned int,FILE*);


//================================================================================================//
/**
* @brief This function prints a vector.
*
* @param[in] double complex* vector
* @param[in] unsigned int length
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_vector_cmplx(double complex*, unsigned int, FILE*);


//================================================================================================//
/**
* @brief This function prints an unsigned int matrix.
*
* @param[in] unsigned int* matrix 
* @param[in] int num_rows
* @param[in] int num_cols
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_matrix_uint(unsigned int*,int,int,FILE*);


//================================================================================================//
/**
* @brief This function prints a double precision matrix.
*
* @param[in] double* matrix 
* @param[in] int num_rows
* @param[in] int num_cols
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_matrix_dbl(double*,int,int,FILE*);


//================================================================================================//
/**
* @brief This function prints a matrix.
*
* @param[in] double complex* matrix 
* @param[in] int num_rows
* @param[in] int num_cols
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_matrix_cmplx(double complex*,int,int,FILE*);


//================================================================================================//
/**
* @brief This function reads complex data from a .dat file.
*
* NOTE: Data in the data file must be expressed as two double precision numbers ie. [real imag]
*
* @param[in] char* filename
* @param[out] unsigned int* data_length
* @param[out] double complex* data_out
*
* @return NONE
*/
//================================================================================================//
void read_dat_file_cmplx(char*,unsigned int*,double complex*);


//================================================================================================//
/**
* @brief This function returns the number of double precision elements in a .dat file.
*
* @param[in] char* filename
*
* @return unsigned int length
*/
//================================================================================================//
unsigned int get_dat_file_length_dbl(char*);


//================================================================================================//
/**
* @brief This function reads double precision data from a .dat file.
*
* @param[in] char* filename
* @param[out] unsigned int* data_length
* @param[out] double* data_out
*
* @return NONE
*/
//================================================================================================//
void read_dat_file_dbl(char*,unsigned int*,double*);


//================================================================================================//
/**
* @brief This function reads unsigned int data from a .dat file.
*
* @param[in] char* filename
* @param[out] unsigned int* data_length
* @param[out] unsigned int* data_out
*
* @return NONE
*/
//================================================================================================//
void read_dat_file_uint(char*,unsigned int*,unsigned int*);


//================================================================================================//
/**
* @brief This function prints a double precision matrix stored in column major format.
*
* @param[in] double* matrix 
* @param[in] int num_rows
* @param[in] int num_cols
* @param[in] FILE* file
*
* @return NONE
*/
//================================================================================================//
void print_column_major_matrix_dbl(double*,int,int,FILE*);


//================================================================================================//
/**
* @brief This function converts a dbl cmplx matrix from row major storage to column major storage.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void row_to_column_major_matrix_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function converts a dbl matrix from row major storage to column major storage.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void row_to_column_major_matrix_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function converts a dbl cmplx matrix from column major storage to row major storage.
*
* @param[in,out] double complex* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void column_to_row_major_matrix_cmplx(double complex*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function converts a dbl matrix from column major storage to row major storage.
*
* @param[in,out] double* matrix
* @param[in] unsigned int num_rows
* @param[in] unsigned int num_cols
*
* @return NONE
*/
//================================================================================================//
void column_to_row_major_matrix_dbl(double*,unsigned int,unsigned int);


//================================================================================================//
/**
* @brief This function creates a korobov sequence of lattice sampled doubles from [0,1]
*
* @param[in] unsigned int base
* @param[in] unsigned int num_samples
* @param[out] double* sequence
*
* @return NONE
*/
//================================================================================================//
void make_korobov_sequence(unsigned int, unsigned int, double* sequence);

#endif //HELPER_H//
