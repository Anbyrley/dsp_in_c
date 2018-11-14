#include "polynomial.h"
#include "lapacke.h"


void initialize_companion_matrix( double complex* coefficients, 
								  unsigned int length,
								  double complex* companion_matrix )
{

	int i,r,position, order;
	double complex coefficients_copy[2*(MAX_POLYNOMIAL_ORDER+1)];

	//===Copy Coefficients===//
	for (i=0; i<(int)length; i++){
		coefficients_copy[i] = coefficients[i];
	}

	//===Prepare For Root Finding===//
	make_array_monic_cmplx(coefficients_copy, length);
	reverse_array_cmplx(coefficients_copy, length);

	//===Fill Companion Matrix===//
	order = length - 1;
	position = -1;
	for (i=0; i<order; i++){
		for (r=0; r<order; r++){
	
			if (i==0){
				if (r < (order-1) ){
					companion_matrix[r +i*order] = 0;
				}
				else{
					companion_matrix[r +i*order] = -coefficients_copy[i];
				}
			} 
			else{
				if ( r < (order-1) ){
					if ((int)r==position){
						companion_matrix[r +i*order] = 1;
					}
					else{
						companion_matrix[r +i*order] = 0;
					}	
					
				} //If its the last column, place the row's coeff
				else{
					companion_matrix[r +i*order] = -coefficients_copy[i];
				}
			}
		}
		
		position++;
	}

	return;


}

void find_roots_dbl( double* coefficients,
					 unsigned int length,
					 double complex* roots )
{

	unsigned int i;
	int order, ok, c1, c2, c3;
	char c4;
	double complex la_dummy[2*(MAX_POLYNOMIAL_ORDER+1)]; //maybe 4?
	double work[4*(MAX_POLYNOMIAL_ORDER+1)];
	double complex la_work[4*(MAX_POLYNOMIAL_ORDER+1)];
	double complex coefficients_cmplx[2*(MAX_POLYNOMIAL_ORDER+1)];
	double complex* companion_matrix;

	//===Create Complex Coefficients===//
	for (i=0; i<length; i++){
		coefficients_cmplx[i] = coefficients[i] + I*0.0;
	}

	//===Create Companion Matrix===//
	companion_matrix = malloc(length*length*sizeof(double complex));
	initialize_companion_matrix(coefficients_cmplx, length, companion_matrix);

	//===Set Up For Function Call===//
	order = length - 1;
	c1=order;			 
	c2=2*order;    			
	c3=1;				
	c4='N';

	//===Find Eigenvalues Of Companion Matrix===//
	if (companion_matrix == NULL){
		fprintf(stderr, "Error: Companion Matrix Is NULL!\n");
	}
	zgeev_(&c4, &c4, &c1, companion_matrix, &c1, roots, la_dummy, &c3,
		   la_dummy, &c3, la_work, &c2, work, &ok);   

	//===Clean Up===//
	free(companion_matrix);

	return;
}

void find_roots_cmplx( double complex* coefficients,
				 	   unsigned int length,
				 	   double complex* roots )
{

	int order, ok, c1, c2, c3;
	char c4;
	double complex la_dummy[2*(MAX_POLYNOMIAL_ORDER+1)]; //maybe 4?
	double work[4*(MAX_POLYNOMIAL_ORDER+1)];
	double complex la_work[4*(MAX_POLYNOMIAL_ORDER+1)];
	double complex* companion_matrix;

	//===Create Companion Matrix===//
	companion_matrix = malloc(length*length*sizeof(double complex));
	initialize_companion_matrix(coefficients, length, companion_matrix);

	//===Set Up For Function Call===//
	order = length - 1;
	c1=order;			 
	c2=2*order;    			
	c3=1;				
	c4='N';

	//===Find Eigenvalues Of Companion Matrix===//
	if (companion_matrix == NULL){
		fprintf(stderr, "Error: Companion Matrix Is NULL!\n");
	}
	zgeev_(&c4, &c4, &c1, companion_matrix, &c1, roots, la_dummy, &c3,
		   la_dummy, &c3, la_work, &c2, work, &ok);   

	//===Clean Up===//
	free(companion_matrix);

	return;

}

void convert_roots_to_coefficients( double complex* roots,
									unsigned int num_roots,
									double complex* coefficients)
{

	unsigned int i, j;
	unsigned int prev_len, cur_len;
	double complex ak[MAX_POLYNOMIAL_ORDER+1]; 	

	//===Initialize Algorithm===//
	prev_len = 2; cur_len = 2;
	ak[0] = -roots[0];
	ak[1] = 1;

	//===Run===//
	if (num_roots > 1){
		for (i=1; i<num_roots; i++){

			//===Compute New Coeffs===//
			cur_len = prev_len + 1;
			for (j=0; j<cur_len; j++){
				if (j == 0){
					coefficients[j] = -roots[i] * ak[j];
				}
				else if (j == cur_len - 1){
					coefficients[j] = 1;
				}
				else{
					coefficients[j] = ak[j-1] + (-roots[i] * ak[j]);
				}	
			}

			//===Save Old Coeffs===//
			for (j=0; j<cur_len; j++){
				ak[j] = coefficients[j];
			}
			prev_len = cur_len;

		}
	}
	else{
		copy_array_cmplx(ak, 2, coefficients);
	}

	reverse_array_cmplx(coefficients, num_roots+1);

	
	return;
}

double complex evaluate_polynomial_cmplx( double complex* coefficients, 
										  unsigned int length,
										  double complex x )
{
		
	//===Locals===//
	unsigned int i;	
	double complex result;

	//===Initialize result (b_n=a_n)===//
	result = coefficients[0];  
	
	//===Evaluate Polynomial===//
	for (i=1; i<length; i++){
		result = result*x + coefficients[i]; 
	}

	return result;
}

double evaluate_polynomial_dbl( double* coefficients,
								unsigned int length, 
								double x )
{
		
	//===Locals===//
	unsigned int i;	
	double result;

	//===Initialize result (b_n=a_n)===//
	result = coefficients[0];  
	
	//===Evaluate Polynomial===//
	for (i=1; i<length; i++){
		result = result*x + coefficients[i]; 
	}

	return result;
}


