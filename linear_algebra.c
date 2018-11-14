#include "linear_algebra.h"

//================================================================================================//
//================================SYMMETRIC MATRIX ROUTINES=======================================//
//================================================================================================//

void fb_average_matrix_cmplx( double complex* matrix,
							  unsigned int num_rows,
							  unsigned int num_cols )
{
	unsigned int maxx;
	double complex exchange[MAX_FRAME_LENGTH];
	double complex temp[MAX_FRAME_LENGTH];
	double complex temp2[MAX_FRAME_LENGTH];

	//===Mallocs===//
	maxx = MAX(num_rows, num_cols);
	if (maxx*maxx > MAX_FRAME_LENGTH){
		fprintf(stderr, "Error:: Matrix Is Too Large! In Function -- fb_average_matrix_cmplx!\n");
		quit();
	}

	//===Make Exchange Matrix===//
	initialize_exchange_matrix_cmplx(exchange, num_rows, num_cols);

	//===Multiply (A*)J===//	
	copy_array_cmplx(matrix, num_rows*num_cols, temp2);
	conjugate_array_cmplx(temp2, num_rows*num_cols);
	matrix_matrix_multiply_cmplx(temp2, num_rows, num_cols, exchange, num_rows, num_cols, temp); 	

	//===Multiply JAJ===//
	initialize_matrix_cmplx(temp2, num_rows, num_cols);
	matrix_matrix_multiply_cmplx(exchange, num_rows, num_cols, temp, num_rows, num_cols, temp2);

	//===Add===//
	add_matrices_cmplx(matrix, temp2, num_rows, num_cols, matrix);

	//===Gain===//
	gain_array_constant_cmplx(matrix, num_rows*num_cols, 0.5); 	

	return;
}

void fb_average_matrix_dbl( double* matrix,
							unsigned int num_rows,
							unsigned int num_cols )
{
	unsigned int maxx;
	double *exchange, *temp1, *temp2;

	//===Mallocs===//
	maxx = MAX(num_rows, num_cols);
	exchange = malloc(maxx*maxx*sizeof(double));
	temp1 = malloc(maxx*maxx*sizeof(double));
	temp2 = malloc(maxx*maxx*sizeof(double));


	//===Make Exchange Matrix J===//
	initialize_exchange_matrix_dbl(exchange, num_rows, num_cols);

	//===Multiply (A*)J===//	
	matrix_matrix_multiply_dbl(matrix, num_rows, num_cols,
							   exchange, num_rows, num_cols, temp1); 	

	//===Multiply JAJ===//
	initialize_matrix_dbl(temp2, num_rows, num_cols);
	matrix_matrix_multiply_dbl(exchange, num_rows, num_cols,
							     temp1, num_rows, num_cols, temp2);

	//===Add===//
	add_matrices_dbl(matrix, temp2, num_rows, num_cols, matrix);

	//===Gain===//
	gain_array_constant_dbl(matrix, num_rows*num_cols, 0.5); 	

	//===Clean Up===//
	free(temp2);
	free(temp1);
	free(exchange);

	return;
}	

void compute_lower_cholesky_factor_dbl( double* matrix,
										int order,
										double* cholesky_factor )
{
	char uplow;
	int success;
	double copy[MAX_FRAME_LENGTH*MAX_FRAME_LENGTH];

	//===Sanity Check===//
	if (order > MAX_FRAME_LENGTH){
		fprintf(stderr, "Error: Matrix Is Too Large! In Function -- compute_lower_cholesky_factor_dbl!\n");
		quit();
	}

	//===Copy Matrix===//
	copy_matrix_dbl(matrix, order, order, copy);

	//===Run Decomposition===//
	uplow = 'L';
	dpotrf_(&uplow, &order, copy, &order, &success);

	//===Copy Over===//
	copy_matrix_dbl(copy, order, order, cholesky_factor);

	return;
}

void compute_upper_cholesky_factor_dbl( double* matrix,
										int order,
										double* cholesky_factor )
{
	char uplow;
	int success;
	double copy[MAX_FRAME_LENGTH*MAX_FRAME_LENGTH];

	//===Sanity Check===//
	if (order > MAX_FRAME_LENGTH){
		fprintf(stderr, "Error: Matrix Is Too Large! In Function -- compute_upper_cholesky_factor_dbl!\n");
		quit();
	}

	//===Copy Matrix===//
	copy_matrix_dbl(matrix, order, order, copy);

	//===Run Decomposition===//
	uplow = 'U';
	dpotrf_(&uplow, &order, copy, &order, &success);

	//===Copy Over===//
	copy_matrix_dbl(copy, order, order, cholesky_factor);

	return;
}


//================================================================================================//
//===========================SINGULAR VALUE DECOMPOSITION ROUTINES================================//
//================================================================================================//

void svd_cmplx( double complex* matrix,
			    int num_rows,
				int num_cols,
				double complex* U,
				double complex* S,
				double complex* VH )
{
	char job;
	int lwork, success;
	double *S_dbl;
	double *rwork;
	double complex *work, *copy;
	double *svals;

	//===Allocate Workspaces===//
	lwork = 10*num_rows*num_cols;
	work = malloc(lwork * sizeof(double complex));
	copy = malloc(num_rows*num_cols * sizeof(double complex));
	rwork = malloc(lwork * sizeof(double));
	S_dbl = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));
	svals = malloc(MAX(num_rows,num_cols) * sizeof(double));

	//===Copy And Transpose For Computation===//
	copy_array_cmplx(matrix, num_rows * num_cols, copy);
	transpose_matrix_cmplx(matrix, num_rows, num_cols);

	//===Run SVD===//	
	job = 'A';
	zgesvd_(&job, &job, &num_rows, &num_cols, matrix, &num_rows, svals, U, &num_rows, VH, &num_cols, work, &lwork, rwork, &success); 	
	
	//===Copy Singular Values===//
	initialize_array_dbl(S_dbl, num_rows * num_cols);
	set_matrix_diagonal_dbl(S_dbl, num_cols, svals, MIN(num_rows,num_cols));
	combine_arrays_dbl_to_cmplx(S_dbl, NULL, num_rows*num_cols, S);

	//===Transpose Results For Return===//
	transpose_matrix_cmplx(U, num_rows, num_rows);
	transpose_matrix_cmplx(VH, num_cols, num_cols);       	

	//===Copy Back===//
	copy_array_cmplx(copy, num_rows * num_cols, matrix);

	//===Clean Up===//
	free(svals);
	free(S_dbl);
	free(rwork);
	free(copy);
	free(work);

	return;
}

void svd_dbl( double* matrix,
			  int num_rows,
			  int num_cols,
			  double* U,
			  double* S,
			  double* VT )
{
	char job;
	int lwork, success;
	double *S_dbl;
	double *work, *copy;
	double* svals;

	//===Allocate Workspaces===//
	lwork = 10*num_rows*num_cols;
	work = malloc(lwork * sizeof(double));
	copy = malloc(num_rows*num_cols * sizeof(double));
	S_dbl = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));
	svals = malloc(MAX(num_rows,num_cols) * sizeof(double));

	//===Copy And Transpose For Computation===//
	copy_array_dbl(matrix, num_rows * num_cols, copy);
	transpose_matrix_dbl(matrix, num_rows, num_cols);

	//===Run SVD===//	
	job = 'A';
	dgesvd_(&job, &job, &num_rows, &num_cols, matrix, &num_rows, svals, U, &num_rows, VT, &num_cols, work, &lwork, &success); 	
	
	//===Copy Singular Values===//
	initialize_array_dbl(S_dbl, num_rows * num_cols);
	set_matrix_diagonal_dbl(S_dbl, num_cols, svals, MIN(num_rows,num_cols));
	copy_array_dbl(S_dbl, num_rows * num_cols, S);

	//===Transpose Results For Return===//
	transpose_matrix_dbl(U, num_rows, num_rows);
	transpose_matrix_dbl(VT, num_cols, num_cols);       	

	//===Copy Back===//
	copy_array_dbl(copy, num_rows * num_cols, matrix);

	//===Clean Up===//
	free(svals);
	free(S_dbl);
	free(copy);
	free(work);

	return;
}

double compute_trace_dbl( double* matrix,
						  unsigned int num_rows,
						  unsigned int num_cols )
{
	double trace;
	double eigenvalues[MAX_FRAME_LENGTH];
	double *U, *S, *VT;

	//===Mallocs===//
	S = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));
	VT = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));
	U = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));

	//===Get SVD===//
	svd_dbl(matrix,num_rows,num_cols,U,S,VT);

	//===Compute Trace===//
	copy_matrix_diagonal_dbl(S, num_rows, num_cols, eigenvalues);
	trace = compute_sum_dbl(eigenvalues, MIN(num_rows,num_cols));

	//===Clean Up===//
	free(U);
	free(VT);
	free(S);

	return trace;
}

double compute_determinant_dbl( double* matrix,
						  		unsigned int num_rows,
						  		unsigned int num_cols )
{
	unsigned int i;
	double determinant;
	double eigenvalues[MAX_FRAME_LENGTH];
	double *U, *S, *VT;

	//===Mallocs===//
	S = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));
	VT = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));
	U = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));

	//===Get SVD===//
	svd_dbl(matrix,num_rows,num_cols,U,S,VT);

	//===Compute Trace===//
	initialize_array_dbl(eigenvalues, MAX_FRAME_LENGTH);
	copy_matrix_diagonal_dbl(S, num_rows, num_cols, eigenvalues);
	determinant = 1;
	for (i=0; i<MAX(num_rows,num_cols); i++){
		determinant *= eigenvalues[i];
	}

	//===Clean Up===//
	free(U);
	free(VT);
	free(S);

	return determinant;
}


void low_rank_approximation_cmplx( double complex* matrix,
								   unsigned int num_rows,
								   unsigned int num_cols,
								   unsigned int rank )
{
	unsigned int i;
	double complex *U, *S, *VH, *temp, *diagonal;

	if (rank > num_cols){
		fprintf(stderr, "Error:: Rank > Input Dimension Of Matrix (Number of Columns)! In Function -- low_rank_approximation_cmplx!\n");
		quit();
	}
	if (rank == 0){
		fprintf(stderr, "Error:: Rank Is Zero! In Function -- low_rank_approximation_cmplx!\n");
		quit();
	}
	
	//===Allocations===//
	diagonal = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double complex));
	temp = malloc(num_rows * num_cols *sizeof(double complex));
	U = malloc(num_rows * num_rows * sizeof(double complex));
	S = malloc(num_rows * num_cols * sizeof(double complex));
	VH = malloc(num_cols * num_cols *sizeof(double complex));

	//===Initializations===//
	initialize_array_cmplx(diagonal, MAX(num_rows,num_cols) * MAX(num_rows,num_cols));
	initialize_array_cmplx(temp, num_rows * num_cols);
	initialize_array_cmplx(U, num_rows * num_rows);
	initialize_array_cmplx(S, num_rows * num_cols);
	initialize_array_cmplx(VH, num_cols * num_cols);

	//===SVD===//
	svd_cmplx(matrix, num_rows, num_cols, U, S, VH);

	//===Take First 'rank' Singular Values===//
	copy_matrix_diagonal_cmplx(S,num_rows,num_cols,diagonal);
	for (i=0; i<MIN(num_rows,num_cols); i++){
		if (i < rank){
			diagonal[i] = diagonal[i];
		}
		else{
			diagonal[i] = 0.0;
		}
	}
	set_matrix_diagonal_cmplx(S, num_cols, diagonal, MIN(num_rows,num_cols));

	//===Reconstruct Matrix===//
	matrix_matrix_multiply_cmplx(S, num_rows, num_cols, VH, num_cols, num_cols, temp);
	matrix_matrix_multiply_cmplx(U, num_rows, num_rows, temp, num_rows, num_cols, matrix);

	//===Clean Up===//
	free(VH);
	free(S);
	free(U);
	free(temp);
	free(diagonal);

	return;
}	

void low_rank_approximation_dbl( double* matrix,
								 unsigned int num_rows,
								 unsigned int num_cols,
								 unsigned int rank )
{
	unsigned int i;
	double *U, *S, *VT, *temp, *diagonal;

	if (rank > num_cols){
		fprintf(stderr, "Error:: Rank > Input Dimension Of Matrix (Number of Columns)! In Function -- low_rank_approximation_dbl!\n");
		quit();
	}
	if (rank == 0){
		fprintf(stderr, "Error:: Rank Is Zero! In Function -- low_rank_approximation_dbl!\n");
		quit();
	}
	
	//===Allocations===//
	diagonal = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double));
	temp = malloc(num_rows * num_cols *sizeof(double));
	U = malloc(num_rows * num_rows * sizeof(double));
	S = malloc(num_rows * num_cols * sizeof(double));
	VT = malloc(num_cols * num_cols *sizeof(double));

	//===Initializations===//
	initialize_array_dbl(diagonal, MAX(num_rows,num_cols) * MAX(num_rows,num_cols));
	initialize_array_dbl(temp, num_rows * num_cols);
	initialize_array_dbl(U, num_rows * num_rows);
	initialize_array_dbl(S, num_rows * num_cols);
	initialize_array_dbl(VT, num_cols * num_cols);

	//===SVD===//
	svd_dbl(matrix, num_rows, num_cols, U, S, VT);

	//===Take First 'rank' Singular Values===//
	copy_matrix_diagonal_dbl(S,num_rows,num_cols,diagonal);
	for (i=0; i<MIN(num_rows,num_cols); i++){
		if (i < rank){
			diagonal[i] = diagonal[i];
		}
		else{
			diagonal[i] = 0.0;
		}
	}
	set_matrix_diagonal_dbl(S, num_cols, diagonal, MIN(num_rows,num_cols));

	//===Reconstruct Matrix===//
	matrix_matrix_multiply_dbl(S, num_rows, num_cols, VT, num_cols, num_cols, temp);
	matrix_matrix_multiply_dbl(U, num_rows, num_rows, temp, num_rows, num_cols, matrix);

	//===Clean Up===//
	free(VT);
	free(S);
	free(U);
	free(temp);
	free(diagonal);

	return;
}	


void low_rank_partition_cmplx( double complex* matrix,
							   unsigned int num_rows,
							   unsigned int num_cols,
							   unsigned int rank )
{
		
	double complex *temp, *temp2;
	double complex *U, *S, *VH;

	//===Allocate Temps===//	
	temp = malloc( MAX(num_rows ,num_cols) * MAX(num_rows , num_cols) * sizeof(double complex));	
	temp2 = malloc( MAX(num_rows , num_cols) * MAX(num_rows , num_cols) * sizeof(double complex));	

	//===Allocate SVDs===//
	U = malloc(num_rows * num_rows * sizeof(double complex));
	S = malloc(num_rows * num_cols * sizeof(double complex));
	VH = malloc(num_cols * num_cols * sizeof(double complex));

	//===Take SVD===//
	svd_cmplx(matrix, num_rows, num_cols, U, S, VH);

	//===Parition U===//
	initialize_array_cmplx(temp, MAX(num_rows ,num_cols) * MAX(num_rows ,num_cols));
	partition_matrix_columns_cmplx(U, num_rows, num_rows, rank, temp, NULL);
	initialize_array_cmplx(U, num_rows * num_rows);
	copy_array_cmplx(temp, num_rows * rank, U); 
	
	//===Parition S===//
	initialize_array_cmplx(temp, MAX(num_rows ,num_cols)*MAX(num_rows ,num_cols));
	partition_matrix_rows_cmplx(S, num_rows, num_cols, rank, temp, NULL);
	initialize_array_cmplx(S, num_rows * num_cols);
	copy_array_cmplx(temp, rank * num_cols, S); 	

	//===Parition VH===//
	initialize_array_cmplx(temp, MAX(num_rows ,num_cols)*MAX(num_rows ,num_cols));
	partition_matrix_columns_cmplx(VH, num_cols, num_cols, rank, temp, NULL);
	initialize_array_cmplx(VH, num_cols * num_cols);
	copy_array_cmplx(temp, rank * num_cols, VH); 	

	//===Recreate Matrix===//
	initialize_array_cmplx(temp, MAX(num_rows ,num_cols)*MAX(num_rows ,num_cols));
	initialize_array_cmplx(temp2, MAX(num_rows ,num_cols)*MAX(num_rows ,num_cols));
	matrix_matrix_multiply_cmplx(S, rank, num_cols, VH, num_cols, rank, temp);	
	matrix_matrix_multiply_cmplx(U, num_rows, rank, temp, rank, rank, temp2);	
	initialize_array_cmplx(matrix, num_rows * num_cols);
	copy_array_cmplx(temp2, num_rows * rank, matrix);

	//===Clean Up===//
	free(VH);
	free(S);
	free(U);
	free(temp2);
	free(temp);

	return;
}

void matrix_power_cmplx( double complex* matrix,
					 	 unsigned int num_rows,
						 unsigned int num_cols,
						 double power )
{
	double swap;
	double complex *U, *S, *VH, *temp;
	double complex diagonal[MAX_POLYNOMIAL_ORDER+1];

	//===Mallocs===//
	temp = malloc(num_rows * num_cols *sizeof(double complex));
	U = malloc(num_rows * num_rows * sizeof(double complex));
	S = malloc(num_rows * num_cols * sizeof(double complex));
	VH = malloc(num_cols * num_cols *sizeof(double complex));

	//===Take Power===//	
	if (power < 0){

		//===Take Psuedo Inverse===//
		initialize_array_cmplx(temp, num_rows*num_cols);
		create_psuedo_inverse_cmplx(matrix, num_rows, num_cols, temp);

		//===Swap Rows and Columns===//
		swap = num_rows;
		num_rows = num_cols;
		num_cols = swap;

		//===Copy Over===//
		initialize_array_cmplx(matrix, num_rows*num_cols);
		copy_array_cmplx(temp, num_rows * num_cols, matrix);

	}

	//===Take SVD===//
	svd_cmplx(matrix, num_rows, num_cols, U, S, VH);

	//===Take Power Of Diagonal Entries===//
	copy_matrix_diagonal_cmplx(S, num_rows, num_cols, diagonal);
	pow_array_cmplx(diagonal, MIN(num_rows, num_cols), (double complex)fabs(power));

	//===Set Diagonal===//
	set_matrix_diagonal_cmplx(S, num_cols, diagonal, MIN(num_rows, num_cols));

	//===Run Multiplication===//
	initialize_array_cmplx(temp, num_rows*num_cols);
	matrix_matrix_multiply_cmplx(S, num_rows, num_cols, VH, num_cols, num_cols, temp);
	matrix_matrix_multiply_cmplx(U, num_rows, num_rows, temp, num_rows, num_cols, matrix);		

	//===Clean Up===//
	free(VH);
	free(S);
	free(U);
	free(temp);

	return;
}

//================================================================================================//
//=================================QR DECOMPOSITON ROUTINES=======================================//
//================================================================================================//

void qr_decomposition_cmplx( double complex* matrix,
						   	 int num_rows,
						   	 int num_cols,
						   	 double complex* Q,
						   	 double complex* R )
{
	double complex* copy;
	int i ,lda, lwork, success;
	double complex r[MAX_POLYNOMIAL_ORDER+1];
	double complex work[MAX_POLYNOMIAL_ORDER+1];
	double complex tau[MAX_POLYNOMIAL_ORDER+1];

	//===Error Pre-Check===//
	if (num_rows < num_cols){
		fprintf(stderr, "Error:: Underdetermined System. Input Transpose Instead. In Function -- qr_decomposition_dbl!\n");
		return;
	}

	//===Copy Matrix For Local Work===//
	copy = malloc(num_rows * num_cols * sizeof(double complex));
	copy_array_cmplx(matrix, num_rows*num_cols, copy);	

	//===Factor A===//
	lda = num_rows; lwork = MAX(1, 10*num_cols); 
	initialize_array_cmplx(work, MAX_POLYNOMIAL_ORDER+1);
	initialize_array_cmplx(tau, MAX_POLYNOMIAL_ORDER+1);
	row_to_column_major_matrix_cmplx(copy, num_rows, num_cols);
	zgeqrf_(&num_rows, &num_cols, copy, &lda, tau, work, &lwork, &success);

	//===Make Q===//
	copy_array_cmplx(copy, num_rows*num_cols, Q);
	zungqr_(&num_rows, &num_rows, &num_cols, Q, &lda, tau, work, &lwork, &success);
	column_to_row_major_matrix_cmplx(Q, num_rows, num_rows);

	//===Make R===//
	column_to_row_major_matrix_cmplx(copy, num_rows, num_cols);
	for (i=0; i<num_cols; i++){
		initialize_array_cmplx(r, MAX_POLYNOMIAL_ORDER+1);
		copy_matrix_column_cmplx(copy, num_rows, num_cols, i, r);
		initialize_array_cmplx(r+(i+1), num_rows-(i+1));
		replace_matrix_column_cmplx(R, num_rows, num_cols, r, i);
	}		

	//===Clean Up===//
	free(copy);

	return;
}

void qr_decomposition_dbl( double* matrix,
						   int num_rows,
						   int num_cols,
						   double* Q,
						   double* R )
{
	double* copy;
	int i ,lda, lwork, success;
	double r[MAX_POLYNOMIAL_ORDER+1];
	double work[MAX_POLYNOMIAL_ORDER+1];
	double tau[MAX_POLYNOMIAL_ORDER+1];

	//===Error Pre-Check===//
	if (num_rows < num_cols){
		fprintf(stderr, "Error:: Underdetermined System. Input Transpose Instead. In Function -- qr_decomposition_dbl!\n");
		return;
	}

	//===Copy Matrix For Local Work===//
	copy = malloc(num_rows * num_cols * sizeof(double));
	copy_array_dbl(matrix, num_rows*num_cols, copy);	

	//===Factor A===//
	lda = num_rows; lwork = MAX(1, 10*num_cols); 
	initialize_array_dbl(work, MAX_POLYNOMIAL_ORDER+1);
	initialize_array_dbl(tau, MAX_POLYNOMIAL_ORDER+1);
	row_to_column_major_matrix_dbl(copy, num_rows, num_cols);
	dgeqrf_(&num_rows, &num_cols, copy, &lda, tau, work, &lwork, &success);

	//===Make Q===//
	copy_array_dbl(copy, num_rows*num_cols, Q);
	dorgqr_(&num_rows, &num_rows, &num_cols, Q, &lda, tau, work, &lwork, &success);
	column_to_row_major_matrix_dbl(Q, num_rows, num_rows);

	//===Make R===//
	column_to_row_major_matrix_dbl(copy, num_rows, num_cols);
	for (i=0; i<num_cols; i++){
		initialize_array_dbl(r, MAX_POLYNOMIAL_ORDER+1);
		copy_matrix_column_dbl(copy, num_rows, num_cols, i, r);
		initialize_array_dbl(r+(i+1), num_rows-(i+1));
		replace_matrix_column_dbl(R, num_rows, num_cols, r, i);
	}		

	//===Clean Up===//
	free(copy);

	return;
}

void compute_upper_triangular_matrix_inverse_cmplx( double complex* matrix, 
										 		    int order,
										 		    double complex* inverse )
{
	char uplo, diag;
	int success;
	double complex* copy;

	//===Make Copy===//
	copy = malloc(order*order*sizeof(double complex));
	copy_array_cmplx(matrix, order*order, copy);

	//===Run Inverse On Copy===//
	uplo = 'U'; diag = 'N';	 
	row_to_column_major_matrix_cmplx(copy, order, order);
	ztrtri_(&uplo, &diag, &order, copy, &order, &success);
	column_to_row_major_matrix_cmplx(copy, order, order);

	//===Copy To Inverse===//
	copy_array_cmplx(copy, order*order, inverse);

	//===Clean Up===//
	free(copy);

	return;
}

void compute_upper_triangular_matrix_inverse_dbl( double* matrix, 
										 		  int order,
										 		  double* inverse )
{
	char uplo, diag;
	int success;
	double* copy;

	//===Make Copy===//
	copy = malloc(order*order*sizeof(double));
	copy_array_dbl(matrix, order*order, copy);

	//===Run Inverse On Copy===//
	uplo = 'U'; diag = 'N';	 
	row_to_column_major_matrix_dbl(copy, order, order);
	dtrtri_(&uplo, &diag, &order, copy, &order, &success);
	column_to_row_major_matrix_dbl(copy, order, order);

	//===Copy To Inverse===//
	copy_array_dbl(copy, order*order, inverse);

	//===Clean Up===//
	free(copy);

	return;
}

void compute_left_inverse_qr_cmplx( double complex* matrix,
								    unsigned int num_rows,
								    unsigned int num_cols,
								    double complex* inverse )
{
	double complex *R, *R_Inverse, *Q;
	double complex* copy;

	//===Mallocs===//
	copy = malloc(num_rows * num_cols * sizeof(double complex));
	R = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double complex));
	R_Inverse = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double complex));
	Q = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double complex));

	//===Run QR Decomposition===//
	copy_array_cmplx(matrix, num_rows*num_cols, copy);
	qr_decomposition_cmplx(copy, num_rows, num_cols, Q, R);

	//===Invert R===//
	compute_upper_triangular_matrix_inverse_cmplx(R,num_cols,R_Inverse);
	append_zeros_matrix_cmplx(R_Inverse, num_cols, num_cols, 0, num_rows-num_cols);

	//===Flip Q===//
	hermitian_transpose_matrix_cmplx(Q, num_rows, num_rows);

	//===Multiply===//
	matrix_matrix_multiply_cmplx(R_Inverse, num_cols, num_rows, Q, num_rows, num_rows, inverse);

	//===Clean Up===//
	free(Q);
	free(R_Inverse);
	free(R);
	free(copy);

	return;
}


void compute_left_inverse_qr_dbl( double* matrix,
								  unsigned int num_rows,
								  unsigned int num_cols,
								  double* inverse )
{
	double R[3*(MAX_POLYNOMIAL_ORDER+1)];
	double R_Inverse[3*(MAX_POLYNOMIAL_ORDER+1)];
	double Q[3*(MAX_POLYNOMIAL_ORDER+1)];

	//===Run QR Decomposition===//
	qr_decomposition_dbl(matrix, num_rows, num_cols, Q, R);

	//===Invert R===//
	compute_upper_triangular_matrix_inverse_dbl(R,num_cols,R_Inverse);
	append_zeros_matrix_dbl(R_Inverse, num_cols, num_cols, 0, num_rows-num_cols);

	//===Flip Q===//
	transpose_matrix_dbl(Q, num_rows, num_rows);

	//===Multiply===//
	matrix_matrix_multiply_dbl(R_Inverse, num_cols, num_rows, Q, num_rows, num_rows, inverse);

	return;
}



//================================================================================================//
//===================================MATRIX UPDATE ROUTINES=======================================//
//================================================================================================//

void rank1_matrix_update_cmplx( double complex* matrix,
							    double complex alpha,
						  	    double complex* vector,
						  	    unsigned int vector_length,
							    double complex beta )
{

	double complex update[MAX_FRAME_LENGTH];
	double complex conjugated[MAX_FRAME_LENGTH];

	//===Sanity Checks===//
	if (vector_length * vector_length > MAX_FRAME_LENGTH){
		fprintf(stderr, "Error:: Vector Length Is Too Long! In Function -- rank1_matrix_update_cmplx!\n");
		quit();
	}

	//===Copy Vector===//
	copy_array_cmplx(vector, vector_length, conjugated);

	//===Compute Rank1 Matrix===//
	vector_outer_product_cmplx(vector, vector_length, conjugated, vector_length, update);
	gain_array_constant_cmplx(update, vector_length * vector_length, beta);

	//===Add For Update===//
	gain_array_constant_cmplx(matrix, vector_length * vector_length, alpha);
	add_matrices_cmplx(matrix, update, vector_length, vector_length, matrix);

	return;
}

void rank1_matrix_update_dbl( double* matrix,
							  double alpha,
						  	  double* vector,
						  	  unsigned int vector_length,
							  double beta )
{

	double* update;
	update = malloc(vector_length * vector_length * sizeof(*vector));

	//===Compute Rank1 Matrix===//
	vector_outer_product_dbl(vector, vector_length, vector, vector_length, update);
	gain_array_constant_dbl(update, vector_length * vector_length, beta);

	//===Add For Update===//
	gain_array_constant_dbl(matrix, vector_length * vector_length, alpha);
	add_matrices_dbl(matrix, update, vector_length, vector_length, matrix);

	//===Clean Up===//
	free(update);		

	return;
}

void rank1_matrix_inverse_update_cmplx( double complex* inverse,
									    double complex alpha,
									    double complex* vector,
									    double complex beta,
									    unsigned int vector_length )
{
	double complex denominator;
	double complex *temp;
	double complex *temp2;
	double complex *temp3;
	
	temp = malloc(vector_length*vector_length*sizeof(*temp));
	temp2 = malloc(vector_length*vector_length*sizeof(*temp2));
	temp3 = malloc(vector_length*vector_length*sizeof(*temp3));

	//===Compute Denominator===//
	initialize_array_cmplx(temp, vector_length*vector_length);
	matrix_vector_multiply_cmplx(vector, vector_length, inverse, vector_length, vector_length, temp);
	denominator = (beta/alpha) * compute_inner_product_cmplx(vector, temp, vector_length) + 1.0;
		
	//===Compute Numerator===//
	initialize_array_cmplx(temp, vector_length*vector_length);
	copy_array_cmplx(vector, vector_length, temp3);
	vector_outer_product_cmplx(vector, vector_length, temp3, vector_length, temp);	
	matrix_matrix_multiply_cmplx(inverse, vector_length, vector_length,
							     temp, vector_length, vector_length, temp2);
	matrix_matrix_multiply_cmplx(temp2,vector_length, vector_length,
							     inverse, vector_length, vector_length, temp);


	//===Gain RHS===//
	gain_array_constant_cmplx(temp, vector_length*vector_length, -(beta/alpha)/denominator);

	//===Update Inverse===//	
	add_matrices_cmplx(inverse, temp, vector_length, vector_length, inverse);

	//===Gain RHS===//
	gain_array_constant_cmplx(inverse, vector_length*vector_length, 1.0/alpha);

	//===Clean Up===//
	free(temp3);
	free(temp2);
	free(temp);

	return;
}


void rank1_matrix_inverse_update_dbl( double* inverse,
									  double alpha,
									  double* vector,
									  double beta,
									  unsigned int vector_length )
{
	double denominator;
	double* temp;
	double* temp2;
	temp = malloc(vector_length*vector_length*sizeof(double));
	temp2 = malloc(vector_length*vector_length*sizeof(double));

	//===Compute Denominator===//
	initialize_array_dbl(temp, vector_length*vector_length);
	matrix_vector_multiply_dbl(vector, vector_length, inverse, vector_length, vector_length, temp);
	denominator = (beta/alpha) * compute_inner_product_dbl(vector, temp, vector_length) + 1.0;
		
	//===Compute Numerator===//
	initialize_array_dbl(temp, vector_length*vector_length);
	vector_outer_product_dbl(vector, vector_length, vector, vector_length, temp);	
	matrix_matrix_multiply_dbl(inverse, vector_length, vector_length,
							   temp, vector_length, vector_length, temp2);
	matrix_matrix_multiply_dbl(temp2,vector_length, vector_length,
							   inverse, vector_length, vector_length, temp);

	//===Gain RHS===//
	gain_array_constant_dbl(temp, vector_length*vector_length, -(beta/alpha)/denominator);

	//===Update Inverse===//	
	add_matrices_dbl(inverse, temp, vector_length, vector_length, inverse);

	//===Gain RHS===//
	gain_array_constant_dbl(inverse, vector_length*vector_length, 1.0/alpha);

	//===Clean Up===//
	free(temp2);
	free(temp);

	return;
}

//================================================================================================//
//==================================PROJECTION ROUTINES===========================================//
//================================================================================================//

void form_projection_matrix_cmplx( double complex* matrix,
								   unsigned int num_rows,
								   unsigned int num_cols,
								   double complex* projection_matrix )
{

	double complex* matrix_H;
	double complex* temp;
	double complex* matrix_H_matrix_1;

	//===Mallocs===//
	matrix_H = malloc(num_rows*num_cols*sizeof(double complex));
	temp = malloc(num_rows*num_cols*sizeof(double complex));
	matrix_H_matrix_1 = malloc(num_rows*num_cols*sizeof(double complex));

	//===Process A Vector===//
	if (num_rows == 1 || num_cols == 1){

		
		//===Form (vH.v)^-1===//
		copy_array_cmplx(matrix, num_rows*num_cols, matrix_H);
		temp[0] = compute_inner_product_cmplx(matrix_H, matrix, num_rows*num_cols);
		matrix_H_matrix_1[0] = 1.0/temp[0];
				
		//===Gain vH===//
		gain_array_constant_cmplx(matrix_H, num_rows*num_cols, matrix_H_matrix_1[0]);

		//===Outer Product===//
		vector_outer_product_cmplx(matrix, num_rows*num_cols, matrix_H, num_rows*num_cols, projection_matrix);
		
	}
	else{ //===Process A Matrix===//

		//===Form C^H===//
		copy_array_cmplx(matrix, num_rows*num_cols, matrix_H);
		hermitian_transpose_matrix_cmplx(matrix_H, num_rows, num_cols);
		
		//===Form (C^H C)^-1===//
		matrix_matrix_multiply_cmplx(matrix_H, num_cols, num_rows, matrix, num_rows, num_cols, temp);
		compute_left_inverse_qr_cmplx(temp, num_cols, num_cols, matrix_H_matrix_1);

		//===Form (C^H C)^-1 C^H===//
		matrix_matrix_multiply_cmplx(matrix_H_matrix_1, num_cols, num_cols, matrix_H, num_cols, num_rows, temp);

		//===Form C [(C^H C)^-1 C^H]===///
		matrix_matrix_multiply_cmplx(matrix, num_rows, num_cols, temp, num_cols, num_rows, projection_matrix);
	
	}


	//===Clean Up===//
	free(matrix_H_matrix_1);
	free(temp);
	free(matrix_H);

	return;
}

//to make I - 11^H/N, make 1 -> 1/(sqrt(N)
void form_orthogonal_projection_matrix_cmplx( double complex* matrix,
								   			  unsigned int num_rows,
								   			  unsigned int num_cols,
								   			  double complex* orthogonal_projection_matrix )
{
	unsigned int dim;
	double complex* projection_matrix;
	double complex* identity_matrix;

	//===Mallocs===//
	dim = MAX(num_rows, num_cols);
	projection_matrix = malloc(dim*dim*sizeof(double complex));
	identity_matrix = malloc(dim*dim*sizeof(double complex));

	//===Process A Vector===//
	if (num_rows == 1 || num_cols == 1){

		//===Get Projection Matrix===//
		form_projection_matrix_cmplx(matrix, num_rows, num_cols, projection_matrix);

		//===Get Identity Matrix===//
		initialize_identity_matrix_cmplx(identity_matrix, dim, dim);

		//===Subtract===//
		subtract_matrices_cmplx(identity_matrix, projection_matrix, dim, dim, orthogonal_projection_matrix);
		
	}
	else{ //===Process A Matrix===//


		//===Get Projection Matrix===//
		form_projection_matrix_cmplx(matrix, num_rows, num_cols, projection_matrix);

		//===Get Identity Matrix===//
		initialize_identity_matrix_cmplx(identity_matrix, num_rows, num_rows);

		//===Subtract===//
		subtract_matrices_cmplx(identity_matrix, projection_matrix, num_rows, num_rows, orthogonal_projection_matrix);
			
	}


	//===Clean Up===//
	free(identity_matrix);
	free(projection_matrix);

	return;
}



//================================================================================================//
//==================================MATRIX OPERATION ROUTINES=====================================//
//================================================================================================//

void vector_outer_product_cmplx( double complex* vector1,
							     int vector1_length,
							     double complex* vector2,
							     int vector2_length,
							     double complex* outer_product )
{
	conjugate_array_cmplx(vector2, vector2_length);
	matrix_matrix_multiply_cmplx( vector1, vector1_length, 1, 
								  vector2, 1, vector2_length, 
								  outer_product );
	return;
}

void vector_outer_product_dbl( double* vector1,
							   int vector1_length,
							   double* vector2,
							   int vector2_length,
							   double* outer_product )
{
	matrix_matrix_multiply_dbl( vector1, vector1_length, 1, 
								vector2, 1, vector2_length, 
								outer_product );
	return;
}


void vector_matrix_multiply_dbl( double* vector,
							 	 int vector_size,
							 	 double* matrix,
							 	 int matrix_rows,
							 	 int matrix_columns,
							 	 double* result )
{

	if (matrix_rows != vector_size){
		fprintf(stderr, "Vector-Matrix Sizes Are Incompatible In Function -- vector_matrix_multiply\n");
		return;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
				1, matrix_columns, vector_size, 
				1.0, vector, vector_size, 
				matrix, matrix_columns, 
				0.0, result, matrix_columns);
	
	return;
}

void vector_matrix_multiply_cmplx( double complex* vector,
							 	   int vector_size,
							 	   double complex* matrix,
							 	   int matrix_rows,
							 	   int matrix_columns,
							 	   double complex* result )
{

	if (matrix_rows != vector_size){
		fprintf(stderr, "Vector-Matrix Sizes Are Incompatible In Function -- vector_matrix_multiply_cmplx\n");
		return;
	}
	double alpha, beta;
	alpha = 1.0; beta = 0.0;

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
				1, matrix_columns, vector_size, 
				&alpha, vector, vector_size, 
				matrix, matrix_columns, 
				&beta, result, matrix_columns);
	
	return;
}



void matrix_vector_multiply_dbl( double* vector,
							 	 int vector_size,
							 	 double* matrix,
							 	 int matrix_rows,
							 	 int matrix_columns,
							 	 double* result )
{


	if (matrix_columns != vector_size){
		fprintf(stderr, "Matrix-Vector Sizes Are Incompatible In Function -- matrix_vector_multiply_dbl\n");
		return;
	}
 
	cblas_dgemv(CblasRowMajor,
		         CblasNoTrans, matrix_rows, matrix_columns,
		         1.0, matrix, matrix_columns,
		         vector, 1, 0.0,
		         result, 1);

	return;
}

void matrix_vector_multiply_cmplx( double complex* vector,
							 	   int vector_size,
							 	   double complex* matrix,
							 	   int matrix_rows,
							 	   int matrix_columns,
							 	   double complex* result )
{

	double alpha, beta;
 	alpha = 1.0; beta = 0.0;

	if (matrix_columns != vector_size){
		fprintf(stderr, "Matrix-Vector Sizes Are Incompatible In Function -- matrix_vector_multiply_cmplx\n");
		return;
	}
	
	cblas_zgemv(CblasRowMajor,
		         CblasNoTrans, matrix_rows, matrix_columns,
		         &alpha, matrix, matrix_columns,
		         vector, 1, &beta,
		         result, 1);

	return;
}


void matrix_matrix_multiply_cmplx( double complex* matrix1,
								   int matrix1_rows,
							 	   int matrix1_columns,
							 	   double complex* matrix2,
							 	   int matrix2_rows,
							 	   int matrix2_columns,
							 	   double complex* result )
{
	double alpha, beta;
	alpha = 1.0; beta = 0.0;

	if (matrix1_columns != matrix2_rows){
		fprintf(stderr, "Matrix-Matrix Sizes Are Incompatible In Function -- matrix_matrix_multiply_cmplx\n");
		return;
	}

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
				matrix1_rows, matrix2_columns, matrix1_columns, 
				&alpha, matrix1, matrix1_columns, 
				matrix2, matrix2_columns, 
				&beta, result, matrix2_columns);	

	return;
}

void matrix_matrix_multiply_dbl( double* matrix1,
							 	 int matrix1_rows,
							 	 int matrix1_columns,
							 	 double* matrix2,
							 	 int matrix2_rows,
							 	 int matrix2_columns,
							 	 double* result )
{

	if (matrix1_columns != matrix2_rows){
		fprintf(stderr, "Matrix-Matrix Sizes Are Incompatible In Function -- matrix_matrix_multiply\n");
		return;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
				matrix1_rows, matrix2_columns, matrix1_columns, 
				1.0, matrix1, matrix1_columns, 
				matrix2, matrix2_columns, 
				0.0, result, matrix2_columns);
	
	return;
}


void compute_matrix_inverse_dbl( double* matrix,
							 	 int order,
							 	 double * inverse )
{

	int N, lwork;
	int success;
	int *pivot;
	double* workspace;

	//===Allocate Space===//
	pivot = malloc(order * order * order * sizeof(*pivot));
	workspace = malloc(order * order * sizeof(*workspace));

	//===Run Setup===//
	N = order;
	copy_array_dbl(matrix, order*order, inverse);
	lwork = order*order;

	//===Factor Matrix===//
	dgetrf_(&N,&N,inverse,&N,pivot,&success);

	//===Compute Inverse===//
	dgetri_(&N, inverse, &N, pivot, workspace, &lwork, &success);

	//===Clean Up===//
	free(workspace);
	free(pivot);

	return;
}

void compute_matrix_inverse_cmplx( double complex* matrix,
							 	   int order,
							 	   double complex* inverse )
{

	int N, lwork;
	int success;
	int *pivot;
	double complex* workspace;

	//===Allocate Space===//
	lwork = order*order;
	pivot = malloc(order * order * order * sizeof(*pivot));
	workspace = malloc(lwork * sizeof(*workspace));

	//===Run Setup===//
	N = order;
	copy_matrix_cmplx(matrix, order, order, inverse);

	//===Factor Matrix===//
	zgetrf_(&N,&N,inverse,&N,pivot,&success);

	//===Compute Inverse===//
	zgetri_(&N, inverse, &N, pivot, workspace, &lwork, &success);

	//===Clean Up===//
	free(workspace);
	free(pivot);

	return;
}

void solve_linear_system_dbl( double* A,
							  unsigned int num_rows,
							  unsigned int num_cols,
							  double* b,
							  double* x )
{
	
	char t;
	int m, n, lda, ldb, nrhs, success;
	int piv[MAX_POLYNOMIAL_ORDER];

	if (num_rows != num_cols){
		fprintf(stderr, "Error:: Matrix Must Be Square! In Function -- solve_linear_system_dbl!\n");
		return;
	}

	m = num_rows; n = num_cols; lda = m; ldb = n; success = 1; nrhs = 1; t = 'T';

	//===Factor Matrix===//
	dgetrf_(&m, &n, A, &lda, piv, &success);

	//===Solve System===//
	dgetrs_(&t, &n, &nrhs, A, &lda, piv, b, &ldb, &success); 

	//===Copy Over===//
	copy_array_dbl(b, num_cols, x);

	return;
} 

void solve_linear_system_qr_dbl( double* A,
								 int num_rows,
								 int num_cols,
								 double* b,
								 double* x )
{
	int lda, ldb, lwork, nrhs, success;
	char trans;
	double* copy;
	double work[MAX_POLYNOMIAL_ORDER+1];
	double b_saved[MAX_POLYNOMIAL_ORDER+1];

	//===Copy Matrix For Local Work===//
	copy = malloc(num_rows * num_cols * sizeof(double));
	copy_array_dbl(A, num_rows * num_cols, copy);
	copy_array_dbl(b, num_cols, b_saved);

	//===Solve Equation===//
	trans = 'N'; nrhs = 1;
	lda = num_rows; ldb = num_rows; lwork = 10*num_rows;
	if (num_cols > num_rows){
		ldb = num_cols;
	}
	row_to_column_major_matrix_dbl(copy, num_rows, num_cols);
	dgels_(&trans, &num_rows, &num_cols, &nrhs, copy, &lda, b, &ldb, work, &lwork, &success);

	//===Copy Back===//
	copy_array_dbl(b, num_cols, x);
	copy_array_dbl(b_saved, num_cols, b);
	
	//===Clean Up===//
	free(copy);

	return;
}

void solve_lower_triangular_system_cmplx( double complex* A,
										  unsigned int num_rows,
										  unsigned int num_cols,
										  double complex* b,
										  double complex* x )
{
	char uplo, trans, diag;
	int N, nrhs, lda, ldb, success;
	double complex b_copy[MAX_POLYNOMIAL_ORDER+1];
	double complex* A_copy;

	//===Mallocs===//
	A_copy = malloc(num_rows*num_cols*sizeof(double complex));
	
	//===Set Up===//
	uplo = 'L'; trans = 'N'; diag = 'N';
	N = (int)num_cols; nrhs = 1;
	lda = (int)num_rows; ldb = (int)num_rows;

	//===Solve System===//
	copy_array_cmplx(b, num_rows, b_copy);
	copy_matrix_cmplx(A, num_rows, num_cols, A_copy);
	row_to_column_major_matrix_cmplx(A_copy, num_rows, num_cols);
	ztrtrs_(&uplo, &trans, &diag, &N, &nrhs, A_copy, &lda, b_copy, &ldb, &success);

	//===Copy Over===//
	copy_array_cmplx(b_copy, num_rows, x);

	//===Clean Up===//
	free(A_copy);

	return;
}

void solve_lower_triangular_system_dbl( double* A,
										unsigned int num_rows,
										unsigned int num_cols,
										double* b,
										double* x )
{
	char uplo, trans, diag;
	int N, nrhs, lda, ldb, success;
	double b_copy[MAX_POLYNOMIAL_ORDER+1];
	double* A_copy;

	//===Mallocs===//
	A_copy = malloc(num_rows*num_cols*sizeof(double));
	
	//===Set Up===//
	uplo = 'L'; trans = 'N'; diag = 'N';
	N = (int)num_cols; nrhs = 1;
	lda = (int)num_rows; ldb = (int)num_rows;

	//===Solve System===//
	copy_array_dbl(b, num_rows, b_copy);
	copy_matrix_dbl(A, num_rows, num_cols, A_copy);
	row_to_column_major_matrix_dbl(A_copy, num_rows, num_cols);
	dtrtrs_(&uplo, &trans, &diag, &N, &nrhs, A_copy, &lda, b_copy, &ldb, &success);

	//===Copy Over===//
	copy_array_dbl(b_copy, num_rows, x);

	//===Clean Up===//
	free(A_copy);

	return;
}



void create_psuedo_inverse_cmplx( double complex* matrix,
								  unsigned int num_rows,
								  unsigned int num_cols,
								  double complex* pinv )
{	
	double complex *temp, *temp2;
	double complex *U, *S, *VH;
	double complex *diagonal;

	//===Mallocs===//
	diagonal = malloc(MAX(num_rows,num_cols) * sizeof(double complex));
	U = malloc(num_rows * num_rows * sizeof(double complex));
	S = malloc(num_rows * num_cols * sizeof(double complex));
	VH = malloc(num_cols * num_cols * sizeof(double complex));
	temp = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double complex));
	temp2 = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double complex));

	//===Get SVD===//
	svd_cmplx(matrix, num_rows, num_cols, U, S, VH);

	//===Invert S===//
	copy_matrix_diagonal_cmplx(S, num_rows, num_cols, diagonal);
	invert_vector_cmplx(diagonal, MIN(num_rows, num_cols));	
	initialize_array_cmplx(S, num_cols * num_rows);
	set_matrix_diagonal_cmplx(S, num_rows, diagonal, MIN(num_rows,num_cols));

	//===Tranpose U and VH===//
	transpose_matrix_cmplx(U,num_rows,num_rows);
	transpose_matrix_cmplx(VH,num_cols,num_cols);

	//===Multiply S*UH===//
	matrix_matrix_multiply_cmplx(S, num_cols, num_rows, U, num_rows, num_rows, temp);
	
	//===Multiply V*(S*UH)===//
	matrix_matrix_multiply_cmplx(VH, num_cols, num_cols, temp, num_cols, num_rows, temp2);
	
	//===Copy Over===//
	initialize_array_cmplx(pinv, num_rows*num_cols);
	copy_array_cmplx(temp2, num_cols*num_rows, pinv);

	//===Clean Up===//
	free(temp2);
	free(temp);
	free(VH);
	free(S);
	free(U);
	free(diagonal);

	return;
}

void create_psuedo_inverse_svd_dbl( double* matrix,
									unsigned int num_rows,
									unsigned int num_cols,
									double* pinv )
{	
	double complex *temp, *temp2;
	double complex *U, *S, *VH, *copy;
	double complex* diagonal;

	//===Mallocs===//
	diagonal = malloc(MAX(num_rows,num_cols) * sizeof(double complex));
	U = malloc(num_rows * num_rows * sizeof(double complex));
	S = malloc(num_rows * num_cols * sizeof(double complex));
	VH = malloc(num_cols * num_cols * sizeof(double complex));
	temp = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double complex));
	temp2 = malloc(MAX(num_rows,num_cols) * MAX(num_rows,num_cols) * sizeof(double complex));
	copy = malloc(num_rows * num_cols * sizeof(double complex));

	//===Get SVD===//
	combine_arrays_dbl_to_cmplx(matrix, NULL, num_rows*num_cols, copy);
	svd_cmplx(copy, num_rows, num_cols, U, S, VH);

	//===Invert S===//
	copy_matrix_diagonal_cmplx(S, num_rows, num_cols, diagonal);
	invert_vector_cmplx(diagonal, MIN(num_rows, num_cols));	
	initialize_array_cmplx(S, num_cols * num_rows);
	set_matrix_diagonal_cmplx(S, num_rows, diagonal, MIN(num_rows,num_cols));

	//===Tranpose U and VH===//
	transpose_matrix_cmplx(U,num_rows,num_rows);
	transpose_matrix_cmplx(VH,num_cols,num_cols);

	//===Multiply S*UH===//
	matrix_matrix_multiply_cmplx(S, num_cols, num_rows, U, num_rows, num_rows, temp);
	
	//===Multiply V*(S*UH)===//
	matrix_matrix_multiply_cmplx(VH, num_cols, num_cols, temp, num_cols, num_rows, temp2);
	
	//===Copy Over===//
	initialize_array_dbl(pinv, num_rows*num_cols);
	split_array_cmplx(temp2, num_cols*num_rows, pinv, NULL);

	//===Clean Up===//
	free(copy);
	free(temp2);
	free(temp);
	free(VH);
	free(S);
	free(U);
	free(diagonal);

	return;
}


void create_psuedo_inverse_dbl( double* matrix,
								unsigned int num_rows,
								unsigned int num_cols,
								double* pinv )
{
	double* matrix_T, *temp;

	matrix_T = malloc(MAX(num_rows,num_cols)*MAX(num_rows,num_cols)*sizeof(double));
	temp = malloc(MAX(num_rows,num_cols)*MAX(num_rows,num_cols)*sizeof(double));

	//===Form A^T===//
	copy_array_dbl(matrix, num_rows * num_cols, matrix_T);
	transpose_matrix_dbl(matrix_T, num_rows, num_cols);
	
	//===Form (A^T*A)^-1===//
	matrix_matrix_multiply_dbl(matrix_T, num_cols, num_rows, matrix, num_rows, num_cols, temp);
	compute_matrix_inverse_dbl(temp, num_cols, temp);

	//===Form (A^T*A)^-1 * A^T===//
	matrix_matrix_multiply_dbl(temp, num_cols, num_cols, matrix_T, num_cols, num_rows, pinv);

	//===Clean Up===//
	free(temp); free(matrix_T);
	
	return;
}

void solve_least_squares_system_cmplx( double complex* A,
								       unsigned int num_rows,
								       unsigned int num_cols,
								       double complex* b,
								       double complex* x )
{	
	double complex* A_plus;
	A_plus = malloc(num_rows * num_cols * sizeof(double complex));

	//===Create A+ and Solve===//
	create_psuedo_inverse_cmplx(A, num_rows, num_cols, A_plus);
	matrix_vector_multiply_cmplx(b, num_rows, A_plus, num_cols, num_rows, x);

	//===Clean Up===//
	free(A_plus);

	return;
} 	

void solve_least_squares_system_dbl( double* A,
								     unsigned int num_rows,
								     unsigned int num_cols,
								     double* b,
								     double* x )
{	
	double* A_plus;
	A_plus = malloc(num_rows * num_cols * sizeof(double));

	//===Create A+ and Solve===//
	create_psuedo_inverse_dbl(A, num_rows, num_cols, A_plus);
	matrix_vector_multiply_dbl(b, num_rows, A_plus, num_cols, num_rows, x);

	//===Clean Up===//
	free(A_plus);

	return;
} 	

double complex compute_quadratic_form_cmplx( double complex* A,
											 unsigned int num_rows,
											 unsigned int num_cols,
											 double complex* x,
											 double complex* y )
{
	double complex result;
	double complex temp[MAX_POLYNOMIAL_ORDER+1];

	//===Compute Ax===//
	matrix_vector_multiply_cmplx(x, num_cols, A, num_rows, num_cols, temp);

	//===Compute y^H * Ax===//	
	result = compute_inner_product_cmplx(y, temp, num_cols);

	return result;
}

double compute_quadratic_form_dbl( double* A,
								   unsigned int num_rows,
								   unsigned int num_cols,
								   double* x,
								   double* y )
{
	double result;
	double temp[MAX_POLYNOMIAL_ORDER+1];

	//===Compute Ax===//
	matrix_vector_multiply_dbl(x, num_cols, A, num_rows, num_cols, temp);

	//===Compute y^T * Ax===//	
	result = compute_inner_product_dbl(y, temp, num_cols);

	return result;
}


//================================================================================================//
//=====================================EIGENVALUE ROUTINES========================================//
//================================================================================================//

void compute_eigensystem_cmplx(  double complex* matrix,
					   	   		 int matrix_order,
								 unsigned int sort,
					   	   		 double complex* eigenvalues,
					   	   		 double complex* right_eigenvectors,
					   	   		 double complex* left_eigenvectors )
{
	char compute_left, compute_right;
	int i, did_not_fail, dbl_order;
	unsigned int* indices;
	double complex* copy;
	double* rwork;
	double complex *work;
	double complex *temp;
	
	//===Set Up===//
	compute_left = 'V';
	compute_right = 'V';
	if (right_eigenvectors == NULL){
		compute_right = 'N';
	}
	if (left_eigenvectors == NULL){
		compute_left = 'N';
	}
	dbl_order = 2*matrix_order*matrix_order;

	work = malloc(2*matrix_order*matrix_order*sizeof(double complex));
	rwork = malloc(2*matrix_order*matrix_order*sizeof(double));
	indices = malloc(matrix_order * sizeof(unsigned int));
	temp = malloc(matrix_order*matrix_order*sizeof(double complex));
	copy = malloc(matrix_order*matrix_order*sizeof(double complex));

	//===Copy Matrix===//
	copy_matrix_cmplx(matrix, matrix_order, matrix_order, copy);

	//===Compute The Eigensystem===//
	zgeev_(&compute_left, &compute_right, &matrix_order, matrix, &matrix_order, eigenvalues, 
		    left_eigenvectors, &matrix_order, right_eigenvectors, &matrix_order, 
			work, &dbl_order, rwork, &did_not_fail);


	if (sort){
		//===Sort By Index===//
		sort_array_indices_by_magnitude_cmplx(eigenvalues, matrix_order, indices);

		//===Rearrange Eigenvalues===//
		for (i=0; i<matrix_order; i++){
			temp[i] = eigenvalues[indices[i]];
		}
		copy_array_cmplx(temp, matrix_order, eigenvalues);
	
		//===Rearrange Eigenvectors===//
		if (right_eigenvectors != NULL){
			initialize_array_cmplx(temp, matrix_order*matrix_order);
			for (i=0; i<matrix_order; i++){
				copy_array_cmplx(right_eigenvectors + indices[i]*matrix_order, matrix_order, temp + i*matrix_order);
			}
			copy_array_cmplx(temp, matrix_order*matrix_order, right_eigenvectors);
		}
		if (left_eigenvectors != NULL){
			initialize_array_cmplx(temp, matrix_order*matrix_order);
			for (i=0; i<matrix_order; i++){
				copy_array_cmplx(left_eigenvectors + indices[i]*matrix_order, matrix_order, temp + i*matrix_order);
			}
			copy_array_cmplx(temp, matrix_order*matrix_order, left_eigenvectors);
		}
	}

	//===Copy Over Matrix===//
	copy_matrix_cmplx(copy, matrix_order, matrix_order, matrix);

	//===Clean Up===//
	free(copy);
	free(temp);
	free(indices);
	free(rwork);
	free(work);

	return;
}

void compute_eigensystem_dbl(  double* matrix,
					   	   	   int matrix_order,
							   unsigned int sort,
					   	   	   double complex* eigenvalues,
					   	   	   double complex* right_eigenvectors,
					   	   	   double complex* left_eigenvectors )
{
	double complex* matrix_cmplx;

	//===Mallocs===//
	matrix_cmplx = malloc(matrix_order*matrix_order*sizeof(double complex));

	//===Compute System===//
	combine_arrays_dbl_to_cmplx(matrix, NULL, matrix_order*matrix_order, matrix_cmplx);
	compute_eigensystem_cmplx(matrix_cmplx, matrix_order, sort, eigenvalues,
					   	   	  right_eigenvectors, left_eigenvectors );

	//===Clean Up===//
	free(matrix_cmplx);

	return;
}

void compute_symmetric_eigensystem_dbl( double* matrix,
					   	   	   			int matrix_order,
										unsigned int sort,
					   	   	   			double* eigenvalues,
					   	   	   			double* eigenvectors )
{
	int i, lda, lwork, did_not_fail;
	char both, uplo;
	unsigned int *indices;
	double *work, *copy, *temp;

	//===Mallocs===//
	work = malloc(2*matrix_order*matrix_order*sizeof(double));
	copy = malloc(matrix_order*matrix_order*sizeof(double));
	indices = malloc(matrix_order * sizeof(unsigned int));
	temp = malloc(matrix_order*matrix_order*sizeof(double));


	//===Copy Over===//
	copy_array_dbl(matrix, matrix_order*matrix_order, copy);

	//===Set Locals===//
	both = 'V';
	uplo = 'U';
	lda = matrix_order;
	lwork = lda*lda*lda; 

	//===Run Eigendecomposition===//
	dsyev_(&both, &uplo, &matrix_order, copy, &lda, eigenvalues, work, &lwork, &did_not_fail);

	//===Copy Over Eigenvectors===//
	copy_array_dbl(copy, matrix_order*matrix_order, eigenvectors);

	//===Sort Eigenvalues===//
	if (sort){

		//===Sort By Index===//
		sort_array_indices_dbl(eigenvalues, matrix_order, indices);
		reverse_array_uint(indices, matrix_order);

		//===Rearrange Eigenvalues===//
		for (i=0; i<matrix_order; i++){
			temp[i] = eigenvalues[indices[i]];
		}
		copy_array_dbl(temp, matrix_order, eigenvalues);
	
		//===Rearrange Eigenvectors===//
		if (eigenvectors != NULL){
			initialize_array_dbl(temp, matrix_order*matrix_order);
			for (i=0; i<matrix_order; i++){
				copy_array_dbl(eigenvectors + indices[i]*matrix_order, matrix_order, temp + i*matrix_order);
			}
			copy_array_dbl(temp, matrix_order*matrix_order, eigenvectors);
		}
	}

	//===Clean Up===//
	free(temp);
	free(indices);
	free(copy);
	free(work);

	return;
}


void fill_eigenvector_matrix( double complex* eigenvectors,
							  unsigned int order,
							  unsigned int num_eigenvalues,
							  double complex* eigenvector_matrix )
{
	unsigned int c,r;
	for (c=0; c<num_eigenvalues; c++){
		for (r=0; r<order; r++){
			eigenvector_matrix[c + r*num_eigenvalues] = eigenvectors[r+c*order];
		}
	}

	return;
}

void print_eigenvectors( double complex* eigenvectors,
						 unsigned int order )
{
	unsigned int i, j;
	//===Print Eigenvectors===//
	for (i=0; i<order; i++){
		fprintf(stdout, "Eigenvector %d: ", i);
		for (j=0; j<order; j++){
			fprintf(stdout, "%+lf%+lfj ", creal(eigenvectors[j+i*order]), cimag(eigenvectors[j+i*order]));
		}
		fprintf(stdout, "\n");
	}	
	fprintf(stdout, "\n");
	return;
}

void print_eigenvectors_dbl( double* eigenvectors,
						     unsigned int order )
{
	unsigned int i, j;
	//===Print Eigenvectors===//
	for (i=0; i<order; i++){
		fprintf(stdout, "Eigenvector %d: ", i);
		for (j=0; j<order; j++){
			fprintf(stdout, "%+lf ", eigenvectors[j+i*order]);
		}
		fprintf(stdout, "\n");
	}	
	fprintf(stdout, "\n");
	return;
}


void find_eigenvalue_multiplicities( double complex* eigenvalues,
									 unsigned int num_eigenvalues,
									 unsigned int* multiplicities )
{

	unsigned int i, j;
	double distance;

	for (i=0; i<num_eigenvalues; i++){
		multiplicities[i] = 0;
		for (j=0; j<num_eigenvalues; j++){
			distance = compute_vector_distance_cmplx(eigenvalues + i, eigenvalues + j, 1, 2.0);
			if (distance < 10e-16){
				multiplicities[i] += 1;
			}
		}
	}

	return;
}

//================================================================================================//
//========================================TEST ROUTINES===========================================//
//================================================================================================//

void test_eigenvalue_routines()
{
	double complex eigenvalues[3];
	double complex *A, *left_eigenvectors, *right_eigenvectors, *eigenvector_matrix;
	unsigned int i, j, order, min_index;

	//===Allocations===//
	A = calloc(4*4, sizeof(double complex));
	left_eigenvectors = calloc(4*4, sizeof(double complex));
	right_eigenvectors = calloc(4*4, sizeof(double complex));
	eigenvector_matrix = calloc(4*4, sizeof(double complex));

	//===Set Data===//
	order = 3;
	A[0] = 1; A[1] = 0; A[2] = 0; 
	A[3] = 0; A[4] = 2; A[5] = 0;
	A[6] = 0; A[7] = 0; A[8] = 3;

	order = 4;
	A[0] = 1; A[1] = 1; A[2] = 1; A[3] = 1;
	A[4] = 1; A[5] = 1; A[6] = 1; A[7] = 1;
	A[8] = 1; A[9] = 1; A[10] = 1; A[11] = 1;
	A[12] = 1; A[13] = 1; A[14] = 1; A[15] = 1;

	//===Print Matrix===//
	fprintf(stdout, "Matrix: \n");
	print_matrix_cmplx(A, order, order, stdout);
	newline();

	//===Compute Eigenvalues And Eigenvectors===//
	compute_eigensystem_cmplx(A, order, 1, eigenvalues, right_eigenvectors, left_eigenvectors);
	
	newline();
	for (i=0; i<order*order; i++){
		fprintf(stdout, "Eigenvector Value %d: %+lf%+lfj\n", i, creal(right_eigenvectors[i]), cimag(right_eigenvectors[i]));
	}
	quit();

	//===Print Eigenvalues===//
	fprintf(stdout, "\n");
	for (i=0; i<order; i++){
		fprintf(stdout, "Eigenvalue %d: %+lf%+lfj\n", i, creal(eigenvalues[i]), cimag(eigenvalues[i]));
	}
	fprintf(stdout, "\n");

	//===Print Eigenvectors===//
	for (i=0; i<order; i++){
		fprintf(stdout, "Right Eigenvector %d: ", i);
		for (j=0; j<order; j++){
			fprintf(stdout, "%+lf%+lfj ", creal(right_eigenvectors[j+i*order]), cimag(right_eigenvectors[j+i*order]));
		}
		fprintf(stdout, "\n");
	}	
	fprintf(stdout, "\n");
	print_eigenvectors(right_eigenvectors, order);
	

	//===Find Minimum===//
	min_index = find_minimum_magnitude_index_cmplx(eigenvalues, order);
	fprintf(stdout, "\nMin Eigenvalue: %+lf%+lfj\n", creal(eigenvalues[min_index]), cimag(eigenvalues[min_index]));
	fprintf(stdout, "Min Eigenvector: ");
	for (j=0; j<order; j++){
		fprintf(stdout, "%+lf%+lfj ", creal(right_eigenvectors[j+min_index*order]), cimag(right_eigenvectors[j+min_index*order]));
	}
	fprintf(stdout, "\n");
	
	fill_eigenvector_matrix(right_eigenvectors,order,order,eigenvector_matrix);
	fprintf(stdout, "\nEigenvector Matrix:\n");
	print_matrix_cmplx(eigenvector_matrix,order,order,stdout);

	//===Clean Up===//
	free(eigenvector_matrix);
	free(right_eigenvectors);
	free(left_eigenvectors);
	free(A);

	return;
}

void test_linear_system_solver()
{
	unsigned int m;
	double A[3*3];
	double b[3];
	double x[3], truth[3];

	A[0] = 1; A[1] = 1; A[2] = 1;
	A[3] = 0; A[4] = 2; A[5] = 5;
	A[6] = 2; A[7] = 5; A[8] = -1;
	
	b[0] = 6;
	b[1] = -4;
	b[2] = 27;	

	truth[0] = 5.0; truth[1] = 3.0; truth[2] = -2.0;

	//===Solve System===//
	solve_linear_system_dbl(A, 3, 3, b, x);
	fprintf(stdout, "\n");
	for (m=0; m<3; m++){
		fprintf(stdout, "Found: %+lf  True: %+lf\n", x[m], truth[m]);
	}
	

	return;
}


void test_svd(void)
{
	double complex U[100];
	double complex S[100];
	double complex VH[100];
	double complex T[100];
	double complex A[12] = {0.0, 1.0, 2.0, 3.0, 0.0, 1.0 + I*0.3, 2.0, 3.0, 0.0, 1.0-3.3*I, 2.0, 3.0};
	int num_rows, num_cols;	

	//===Run SVD===//
	num_rows = 4; num_cols = 3;
	svd_cmplx(A, num_rows, num_cols, U, S, VH);

	//===Print Results===//
	print_matrix_cmplx(U, num_rows, num_rows, stdout);	
	newline();
	print_matrix_cmplx(S, num_rows, num_cols, stdout);
	newline();
	print_matrix_cmplx(VH, num_cols, num_cols, stdout); 
	newline();
	print_matrix_cmplx(A, num_rows, num_cols, stdout);
	newline();

	//===Run Reconstruction===//
	matrix_matrix_multiply_cmplx(U,num_rows, num_rows,
							 	 S, num_rows, num_cols,
							 	 T );
	matrix_matrix_multiply_cmplx(T, num_rows, num_cols,
							 	 VH, num_cols, num_cols,
							 	 A );
	print_matrix_cmplx(A, num_rows, num_cols, stdout);
	newline();

	return;
}

void test_svd_dbl()
{

	double U[100];
	double S[100];
	double VH[100];
	double T[100];
	double A[6] = {3.0, 2.0, 2.0, 2.0, 3.0, -2.0};
	int num_rows, num_cols;	

	//===Run SVD===//
	num_rows = 2; num_cols = 3;
	svd_dbl(A, num_rows, num_cols, U, S, VH);

	//===Print Results===//
	fprintf(stdout, "U:\n");
	print_matrix_dbl(U, num_rows, num_rows, stdout);	
	newline();
	fprintf(stdout, "S:\n");
	print_matrix_dbl(S, num_rows, num_cols, stdout);
	newline();
	fprintf(stdout, "VH:\n");
	print_matrix_dbl(VH, num_cols, num_cols, stdout); 
	newline();
	fprintf(stdout, "Matrix:\n");
	print_matrix_dbl(A, num_rows, num_cols, stdout);
	newline();

	//===Run Reconstruction===//
	matrix_matrix_multiply_dbl(U,num_rows, num_rows,
							 	 S, num_rows, num_cols,
							 	 T );
	matrix_matrix_multiply_dbl(T, num_rows, num_cols,
							 	 VH, num_cols, num_cols,
							 	 A );
	fprintf(stdout, "Reconstructed From SVD:\n");
	print_matrix_dbl(A, num_rows, num_cols, stdout);
	newline();


	//===Low Rank Approximation===//
	fprintf(stdout, "Low Rank Approximation:\n");
	low_rank_approximation_dbl(A, num_rows, num_cols, 1);
	svd_dbl(A, num_rows, num_cols, U, S, VH);

	//===Print Results===//
	fprintf(stdout, "U:\n");
	print_matrix_dbl(U, num_rows, num_rows, stdout);	
	newline();
	fprintf(stdout, "S:\n");
	print_matrix_dbl(S, num_rows, num_cols, stdout);
	newline();
	fprintf(stdout, "VT:\n");
	print_matrix_dbl(VH, num_cols, num_cols, stdout); 
	newline();
	fprintf(stdout, "Matrix:\n");
	print_matrix_dbl(A, num_rows, num_cols, stdout);
	newline();



	return;
}

void test_outer_product(void)
{
	double vector[3];
	double Rxx[9];

	vector[0] = 3.0; vector[1] = 2.0; vector[2] = 2.0;

	vector_outer_product_dbl(vector, 3, vector, 3, Rxx);

	newline();	
	print_matrix_dbl(Rxx, 3, 3, stdout);
	newline();	

	return;
}

void test_QR_decomposition()
{
	int num_rows, num_cols;
	double error;
	double A[100];
	double A_recon[100];
	double R[MAX_POLYNOMIAL_ORDER+1];
	double Q[MAX_POLYNOMIAL_ORDER+1];

	//===Make A Matrix===//
	num_rows = 6; num_cols = 4;
	A[0] = -0.57; A[1] = -1.28; A[2] = -0.39; A[3] = 0.25;
	A[4] = -1.93; A[5] = 1.08; A[6] = -0.31; A[7] = -2.14;
	A[8] = 2.30; A[9] = 0.24; A[10] = 0.40; A[11] = -0.35;
	A[12] = -1.93; A[13] = 0.64; A[14] = -0.66; A[15] = 0.08;
	A[16] = 0.15; A[17] = 0.30; A[18] = 0.15; A[19] = -2.13;
	A[20] = -0.02; A[21] = 1.03; A[22] = -1.43; A[23] = 0.50; 	

	//===Run QR Decomposition===//
	qr_decomposition_dbl(A, num_rows, num_cols, Q, R);

	//===Print Original===//
	newline();
	fprintf(stdout, "Original Matrix\n");
	print_matrix_dbl(A, num_rows, num_cols, stdout);
	newline();

	//===Reconstruct===//
	matrix_matrix_multiply_dbl(Q, num_rows, num_rows , R, num_rows, num_cols, A_recon);	

	//===Print Reconstructed===//
	newline();
	fprintf(stdout, "Reconstructed Matrix\n");
	print_matrix_dbl(A_recon, num_rows, num_cols, stdout);
	newline();

	//===Compute Error===//
	error = compute_mean_squared_error_dbl(A, A_recon, num_rows*num_cols);
	fprintf(stdout, "Error: %+lf\n", error);

	return;
}

void test_qr_solver()
{
	int num_rows, num_cols;
	double error;
	double A[100];
	double b[100];
	double x[100];
	double x_true[100];

	
	fprintf(stdout, "|======Square System=======|");

	num_rows = 3; num_cols = 3;
	A[0] = 3; A[1] = -2; A[2] = 5;
	A[3] = 4; A[4] = -7; A[5] = -1;
	A[6] = 5; A[7] = -6; A[8] = 4;

	b[0] = 2;	x_true[0] = 1.0; 
	b[1] = 19;  x_true[1] = -2.0;
	b[2] = 13;	x_true[2] = -1.0;

	solve_linear_system_qr_dbl(A, num_rows, num_cols, b, x);

	newline();

	fprintf(stdout, "System:\n");
	print_matrix_dbl(A, num_rows, num_cols, stdout);
	newline();

	fprintf(stdout, "Found Answer:\n");
	print_vector_dbl(x, num_cols, stdout);
	newline();

	fprintf(stdout, "True Answer:\n");
	print_vector_dbl(x_true, num_cols, stdout);
	newline();
	
	fprintf(stdout, "Error:\n");
	error = compute_mean_squared_error_dbl(x, x_true, num_cols);
	fprintf(stdout, "%+lf\n", error);
	newline();


	fprintf(stdout, "|======Overdetermined System=======|");

	num_rows = 3; num_cols = 2;
	A[0] = 2; A[1] = 0;
	A[2] = -1; A[3] = 1;
	A[4] = 0; A[5] = 2;

	b[0] = 1;	x_true[0] = 1.0/3.0; 
	b[1] = 0;   x_true[1] = -1.0/3.0;
	b[2] = -1;

	solve_linear_system_qr_dbl(A, num_rows, num_cols, b, x);

	newline();
	fprintf(stdout, "Found Answer:\n");
	print_vector_dbl(x, num_cols, stdout);
	newline();

	fprintf(stdout, "True Answer:\n");
	print_vector_dbl(x_true, num_cols, stdout);
	newline();
	
	fprintf(stdout, "Error:\n");
	error = compute_mean_squared_error_dbl(x, x_true, num_cols);
	fprintf(stdout, "%+lf\n", error);
	newline();

	

	fprintf(stdout, "|======Underdetermined System=======|");

	num_rows = 2; num_cols = 3;
	A[0] = 1; A[1] = 1; A[2] = 1;
	A[3] = -1; A[4] = -1; A[5] = 1;
	
	b[0] = 1;	x_true[0] = 1.0/4.0; 
	b[1] = 0;   x_true[1] = 1.0/4.0;
			    x_true[2] = 1.0/2.0;

	solve_linear_system_qr_dbl(A, num_rows, num_cols, b, x);

	newline();
	fprintf(stdout, "Found Answer:\n");
	print_vector_dbl(x, num_cols, stdout);
	newline();

	fprintf(stdout, "True Answer:\n");
	print_vector_dbl(x_true, num_cols, stdout);
	newline();
	
	fprintf(stdout, "Error:\n");
	error = compute_mean_squared_error_dbl(x, x_true, num_cols);
	fprintf(stdout, "%+lf\n", error);
	newline();

	return;
}

void test_invert_upper_triangular_matrix_dbl()
{
	int num_rows, num_cols;		
	double error;
	double R[MAX_POLYNOMIAL_ORDER+1];
	double inverse[MAX_POLYNOMIAL_ORDER+1];
	double R1[MAX_POLYNOMIAL_ORDER+1];
	
	//===Make Matrix===//
	num_rows = 5; num_cols = 5;
	R[0] = 1.0; R[1] = 3.0; R[2] = 4.0; R[3] = 5.0; R[4] = 6.0;
	R[5] = 0.0; R[6] = 2.0; R[7] = 8.0; R[8] = 9.0; R[9] = 1.0;
	R[10] =0.0; R[11]= 0.0; R[12]= 4.0; R[13]= 8.0; R[14] =4.0;
	R[15] =0.0; R[16]= 0.0; R[17]= 0.0; R[18]=-2.0; R[19] =6.0; 
	R[20] =0.0; R[21]= 0.0; R[22]= 0.0; R[23]= 0.0; R[24]=-1.0; 

	//===Compute Inverse===//
	compute_upper_triangular_matrix_inverse_dbl(R, num_rows, inverse);

	//===Print Found Inverse===//
	newline();
	fprintf(stdout, "Found Inverse\n");
	print_matrix_dbl(inverse, num_rows, num_cols, stdout);
	newline();

	//===Make True Inverse===//
	R1[0] = 1.0; R1[1] =-1.5; R1[2] = 2.0; R1[3] =3.75; R1[4] = 35.0;
	R1[5] = 0.0; R1[6] = 0.5; R1[7] =-1.0; R1[8] =-1.75; R1[9] = -14.0;
	R1[10] =0.0; R1[11]= 0.0; R1[12]=0.25; R1[13]= 1.0; R1[14] =7.0;
	R1[15] =0.0; R1[16]= 0.0; R1[17]= 0.0; R1[18]=-0.5; R1[19] =-3.0; 
	R1[20] =0.0; R1[21]= 0.0; R1[22]= 0.0; R1[23]= 0.0; R1[24]=-1.0; 

	//===Print True Inverse===//
	newline();
	fprintf(stdout, "True Inverse\n");
	print_matrix_dbl(R1, num_rows, num_cols, stdout);
	newline();

	//===Compute And Print Error===//
	fprintf(stdout, "Error:\n");
	error = compute_mean_squared_error_dbl(inverse, R1, num_cols*num_rows);
	fprintf(stdout, "%+lf\n", error);
	newline();


	return;
}

void test_qr_inverse()
{
	int num_rows, num_cols;
	double A[100];
	double A_Inverse[100];
	double I_Recon[100];
	double complex A_cmplx[100];
	double complex A_Inverse_cmplx[100];
	double complex I_Recon_cmplx[100];

	//===Make A Matrix===//
	num_rows = 6; num_cols = 4;
	A[0] = -0.57; A[1] = -1.28; A[2] = -0.39; A[3] = 0.25;
	A[4] = -1.93; A[5] = 1.08; A[6] = -0.31; A[7] = -2.14;
	A[8] = 2.30; A[9] = 0.24; A[10] = 0.40; A[11] = -0.35;
	A[12] = -1.93; A[13] = 0.64; A[14] = -0.66; A[15] = 0.08;
	A[16] = 0.15; A[17] = 0.30; A[18] = 0.15; A[19] = -2.13;
	A[20] = -0.02; A[21] = 1.03; A[22] = -1.43; A[23] = 0.50; 	


	//===Print Matrix===//
	newline();
	fprintf(stdout, "Matrix: \n");
	print_matrix_dbl(A, num_rows, num_cols, stdout);
	newline();

	//===Make Left Inverse===//
	compute_left_inverse_qr_dbl(A, num_rows, num_cols, A_Inverse);
	
	//===Print Left Inverse===//
	newline();
	fprintf(stdout, "Left Inverse: \n");
	print_matrix_dbl(A_Inverse, num_cols, num_rows, stdout);
	newline();
	
	//===Make Identity===//
	matrix_matrix_multiply_dbl(A_Inverse, num_cols, num_rows, A, num_rows, num_cols, I_Recon);

	//===Print Identity===//
	newline();
	fprintf(stdout, "A^-1 * A: \n");
	print_matrix_dbl(I_Recon, num_cols, num_cols, stdout);
	newline();


	//===Make A Matrix===//
	num_rows = 6; num_cols = 4;
	A_cmplx[0] = -0.57; A_cmplx[1] = -1.28; A_cmplx[2] = -0.39; A_cmplx[3] = 0.25;
	A_cmplx[4] = -1.93; A_cmplx[5] = 1.08; A_cmplx[6] = -0.31; A_cmplx[7] = -2.14;
	A_cmplx[8] = 2.30; A_cmplx[9] = 0.24; A_cmplx[10] = 0.40; A_cmplx[11] = -0.35;
	A_cmplx[12] = -1.93; A_cmplx[13] = 0.64; A_cmplx[14] = -0.66; A_cmplx[15] = 0.08;
	A_cmplx[16] = 0.15; A_cmplx[17] = 0.30; A_cmplx[18] = 0.15; A_cmplx[19] = -2.13;
	A_cmplx[20] = -0.02; A_cmplx[21] = 1.03; A_cmplx[22] = -1.43; A_cmplx[23] = 0.50; 	


	//===Print Matrix===//
	newline();
	fprintf(stdout, "Matrix Cmplx: \n");
	print_matrix_cmplx(A_cmplx, num_rows, num_cols, stdout);
	newline();

	//===Make Left Inverse===//
	compute_left_inverse_qr_cmplx(A_cmplx, num_rows, num_cols, A_Inverse_cmplx);
	
	//===Print Left Inverse===//
	newline();
	fprintf(stdout, "Left Inverse Cmplx: \n");
	print_matrix_cmplx(A_Inverse_cmplx, num_cols, num_rows, stdout);
	newline();
	
	//===Make Identity===//
	matrix_matrix_multiply_cmplx(A_Inverse_cmplx, num_cols, num_rows, A_cmplx, num_rows, num_cols, I_Recon_cmplx);

	//===Print Identity===//
	newline();
	fprintf(stdout, "A^-1 * A Cmplx: \n");
	print_matrix_cmplx(I_Recon_cmplx, num_cols, num_cols, stdout);
	newline();


	return;
}

void test_svd_inverse()
{

	int num_rows, num_cols;
	double A[100];
	double A_Inverse[100];
	double I_Recon[100];

	//===Make A Matrix===//
	num_rows = 6; num_cols = 4;
	A[0] = -0.57; A[1] = -1.28; A[2] = -0.39; A[3] = 0.25;
	A[4] = -1.93; A[5] = 1.08; A[6] = -0.31; A[7] = -2.14;
	A[8] = 2.30; A[9] = 0.24; A[10] = 0.40; A[11] = -0.35;
	A[12] = -1.93; A[13] = 0.64; A[14] = -0.66; A[15] = 0.08;
	A[16] = 0.15; A[17] = 0.30; A[18] = 0.15; A[19] = -2.13;
	A[20] = -0.02; A[21] = 1.03; A[22] = -1.43; A[23] = 0.50; 	


	//===Print Matrix===//
	newline();
	fprintf(stdout, "Matrix: \n");
	print_matrix_dbl(A, num_rows, num_cols, stdout);
	newline();

	//===Create Inverse===//
	create_psuedo_inverse_svd_dbl(A, num_rows, num_cols, A_Inverse);

	//===Print Left Inverse===//
	newline();
	fprintf(stdout, "SVD Left Inverse: \n");
	print_matrix_dbl(A_Inverse, num_cols, num_rows, stdout);
	newline();
	
	//===Make Identity===//
	matrix_matrix_multiply_dbl(A_Inverse, num_cols, num_rows, A, num_rows, num_cols, I_Recon);

	//===Print Identity===//
	newline();
	fprintf(stdout, "A^-1 * A: \n");
	print_matrix_dbl(I_Recon, num_cols, num_cols, stdout);
	newline();


	return;
}


void test_lower_triangular_solver()
{
	unsigned int num_rows, num_cols;
	double complex A[4*4];
	double complex b[4];
	double complex x_true[4];
	double complex x[4];

	//===Make A Matrix===//
	num_rows = 4; num_cols = 4;
	A[0] = 3.0; A[1] = 0.0; A[2] = 0.0; A[3] = 0.0;
	A[4] = -1.0; A[5] = 1.0; A[6] = 0.0; A[7] = 0.0;
	A[8] = 3.0; A[9] = -2.0; A[10] = -1.0; A[11] = 0.0;
	A[12] = 1.0; A[13] = -2.0; A[14] = 6.0; A[15] = 2.0;

	//===Set Up RHS and LHS===//
	b[0] = 5.0;	x_true[0] = 5.0/3.0; 
	b[1] = 6.0; x_true[1] = 23.0/3.0;
	b[2] = 4.0; x_true[2] = -43.0/3.0;
	b[3] = 2.0;	x_true[3] = 305.0/6.0;

	//===Print===//
	fprintf(stdout, "True Solution CMPLX: \n");
	print_vector_cmplx(x_true, 4, stdout);

	//===Solve===//
	solve_lower_triangular_system_cmplx(A, num_rows, num_cols, b, x);

	//===Print===//
	fprintf(stdout, "Computed Solution CMPLX: \n");
	print_vector_cmplx(x, 4, stdout);



	double A_dbl[4*4];
	double b_dbl[4];
	double x_true_dbl[4];
	double x_dbl[4];

	//===Make A Matrix===//
	num_rows = 4; num_cols = 4;
	A_dbl[0] = 3.0; A_dbl[1] = 0.0; A_dbl[2] = 0.0; A_dbl[3] = 0.0;
	A_dbl[4] = -1.0; A_dbl[5] = 1.0; A_dbl[6] = 0.0; A_dbl[7] = 0.0;
	A_dbl[8] = 3.0; A_dbl[9] = -2.0; A_dbl[10] = -1.0; A_dbl[11] = 0.0;
	A_dbl[12] = 1.0; A_dbl[13] = -2.0; A_dbl[14] = 6.0; A_dbl[15] = 2.0;

	//===Set Up RHS and LHS===//
	b_dbl[0] = 5.0;	x_true_dbl[0] = 5.0/3.0; 
	b_dbl[1] = 6.0; x_true_dbl[1] = 23.0/3.0;
	b_dbl[2] = 4.0; x_true_dbl[2] = -43.0/3.0;
	b_dbl[3] = 2.0;	x_true_dbl[3] = 305.0/6.0;

	//===Print===//
	fprintf(stdout, "True Solution DBL: \n");
	print_vector_dbl(x_true_dbl, 4, stdout);

	//===Solve===//
	solve_lower_triangular_system_dbl(A_dbl, num_rows, num_cols, b_dbl, x_dbl);

	//===Print===//
	fprintf(stdout, "Computed Solution DBL: \n");
	print_vector_dbl(x_dbl, 4, stdout);


	return;
}
