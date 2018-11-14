/** @file polynomial.h
*   @brief Contains functions for working with polynomials
*
*	This unit requires the helper.o unit.
*
*  @author Alex N. Byrley (anbyrley)
*  @date September 2015
*  @bug No known bugs
*/

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

//================================================================================================//
//===================================STANDARD INCLUDES============================================//
//================================================================================================//

#include "macros.h"
#include "helper.h"


//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//

//================================================================================================//
/**
* @brief This function creates the companion matrix of a polynomial.
*
* @param[in] double complex* coefficients
* @param[in] unsigned int length
* @param[out] double complex* companion_matrix
*
* @return NONE
*/
//================================================================================================//
void initialize_companion_matrix(double complex*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function finds roots of a polynomial represented by a double coefficient array.
*
* @param[in] double* coefficients
* @param[in] unsigned int length
* @param[out] double complex* roots
*
* @return NONE
*/
//================================================================================================//
void find_roots_dbl(double*, unsigned int, double complex*);


//================================================================================================//
/**
* @brief This function finds roots of a polynomial represented by a double complex coefficient array.
*
* @param[in] double complex* coefficients
* @param[in] unsigned int length
* @param[out] double complex* roots
*
* @return NONE
*/
//================================================================================================//
void find_roots_cmplx(double complex*, unsigned int, double complex*);


//================================================================================================//
/**
* @brief This function converts an array of roots to the coefficients they represent.
*
* @param[in] double complex* roots
* @param[in] unsigned int num_roots
* @param[out] double complex* coefficients
*
* @return NONE
*/
//================================================================================================//
void convert_roots_to_coefficients(double complex*, unsigned int, double complex*);


//================================================================================================//
/**
* @brief This function evaluates a double complex polynomial.
*
* @param[in] double complex* coefficients
* @param[in] unsigned int length
* @param[in] double complex x
*
* @return double complex result
*/
//================================================================================================//
double complex evaluate_polynomial_cmplx(double complex*, unsigned int, double complex);


//================================================================================================//
/**
* @brief This function evaluates a double precision polynomial.
*
* @param[in] double* coefficients
* @param[in] unsigned int length
* @param[in] double x
*
* @return double result
*/
//================================================================================================//
double evaluate_polynomial_dbl(double*, unsigned int, double);



#endif //POLYNOMIAL_H//
