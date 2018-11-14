#include "dsp.h"

//================================================================================================//
//=============================FAST FOURIER TRANSFORM FUNCTIONS===================================//
//================================================================================================//

void fft_41( double* data, 
		      int nn, 
		      int isign )
{

    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
	double swap_temp;
    
    n = nn << 1;
    j = 1;
    for (i=1; i<n; i+=2){
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
		    j -= m;
		    m >>= 1;
		}
		j += m;
    }
    mmax = 2;
    while (n > mmax){
		istep = 2*mmax;
		theta = TWOPI/(-isign*mmax);  //reversed sign for python compliance
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1; m<mmax; m+=2) {
		    for (i=m; i<=n; i+=istep){
				j = i+mmax;
				tempr = wr*data[j] - wi*data[j+1];
				tempi = wr*data[j+1] + wi*data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
		    }
		    wr = (wtemp = wr)*wpr - wi*wpi + wr;
		    wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
    }
	return;
}

void fft_41_real( double* data, 
				int length,
				int isign )
{
	int i, i1, i2, i3, i4, np3;
	double c1, c2, h1r, h1i, h2r, h2i;
	double wr, wi, wpr, wpi, wtemp, theta;

	//===Init===//
	c1 = 0.5;
	theta = -1.0*M_PI/((double)(length >> 1)); //compensate for reversed sign in fft_41 for python compliance
	if (isign == 1){
		c2 = -0.5;
		fft_41(data, length >> 1, 1);
	}
	else{
		c2 = 0.5;
		theta *= -1.0;
	}
	
	//===Run===//
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0+wpr;
	wi = wpi;
	np3 = length+3;
	for (i=2; i<=(length>>2);i++){
		i4 = 1+(i3=np3-(i2=1+(i1=i+i-1))); //setting i1, then i2 wrt i1, then i3 wrt i2 etc
		h1r = c1*(data[i1]+data[i3]);
		h1i = c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i = c2*(data[i1]-data[i3]);
		data[i1] = h1r+wr*h2r-wi*h2i;
		data[i2] = h1i+wr*h2i+wi*h2r;
		data[i3] = h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr = (wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1){
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	}
	else{
		data[1] = c1*((h1r=data[1])+data[2]);
		data[2] = c1*(h1r-data[2]);
		fft_41(data, length>>1, -1);
	}

	return;
}

//================================================================================================//
//=================================FFT INTERFACE FUNCTIONS========================================//
//================================================================================================//

void fft_dbl( double* data,
		   	  unsigned int length,
		   	  double complex* result )
{
	unsigned int i;
	double XFFT[2*(length+1)];

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- fft_dbl\n");
		quit();
	}									

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = data[i];							
		XFFT[2*i+2] = 0;							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, 1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1] + I*XFFT[2*i+2]; 
	}			

	return;
}

void ifft_dbl( double complex *data, 
		   	   int length,
			   double* result )
{
	int i;
	double XFFT[2*length+1];

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- ifft_dbl\n");
		quit();
	}									

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = creal(data[i]);							
		XFFT[2*i+2] = cimag(data[i]);							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, -1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1];
		result[i] /= (double)length; 
	}			

	return;
}

void fft_cmplx( double complex *data, 
		   		int length, 
		   		double complex* result )
{
	int i;
	double XFFT[2*length+1];

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- fft_cmplx\n");
		quit();
	}									

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = creal(data[i]);							
		XFFT[2*i+2] = cimag(data[i]);							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, 1);								

	//===Convert To Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1] + I*XFFT[2*i+2];
	}			

	return;
}

void ifft_cmplx( double complex *data, 
		   		 int length,
				 double complex *result )
{
	int i;
	double XFFT[2*length+1];

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- ifft_cmplx\n");
		quit();
	}									

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = creal(data[i]);							
		XFFT[2*i+2] = cimag(data[i]);							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, -1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1] + I*XFFT[2*i+2];
		result[i] /= (double complex) length; 
	}			

	return;
}

void fft_normalized_dbl( double* data,
		   	  			 unsigned int length,
		   	  			 double complex* result )
{
	unsigned int i;
	double XFFT[2*(length+1)];

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- fft_normalized_dbl\n");
		quit();
	}									

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = data[i];							
		XFFT[2*i+2] = 0;							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, 1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1] + I*XFFT[2*i+2]; 
		result[i] /= ((double)length);
	}			

	return;
}

void ifft_normalized_dbl( double complex *data, 
		   	   			  int length,
			   			  double* result )
{
	int i;
	double XFFT[2*length+1];

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- ifft_normalized_dbl\n");
		quit();
	}									

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = creal(data[i]);							
		XFFT[2*i+2] = cimag(data[i]);							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, -1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1];
	}			

	return;
}

void fft_normalized_malloc_dbl( double* data,
		   	  			 	    unsigned int length,
		   	  			 	    double complex* result )
{
	unsigned int i;
	double* XFFT;

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- fft_normalized_malloc_dbl\n");
		quit();
	}									

	//===Malloc===//
	XFFT = malloc(2*(length+1)*sizeof(double));

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = data[i];							
		XFFT[2*i+2] = 0;							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, 1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1] + I*XFFT[2*i+2]; 
		result[i] /= ((double)length);
	}			

	//===Clean Up===//
	free(XFFT);

	return;
}

void fft_malloc_dbl( double* data,
		   	  	     unsigned int length,
		   	  	     double complex* result )
{
	unsigned int i;
	double *XFFT;

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- fft_malloc_dbl\n");
		quit();
	}			

	//===Malloc===//
	XFFT = malloc(2*(length+1)*sizeof(double));						

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = data[i];							
		XFFT[2*i+2] = 0;							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, 1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1] + I*XFFT[2*i+2]; 
	}			

	//===Clean Up===//
	free(XFFT);

	return;
}

void ifft_malloc_dbl( double complex *data, 
		   	   		  int length,
			   		  double* result )
{
	int i;
	double* XFFT;

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- ifft_malloc_dbl\n");
		quit();
	}									

	//===Malloc===//
	XFFT = malloc(2*(length+1)*sizeof(double));						

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = creal(data[i]);							
		XFFT[2*i+2] = cimag(data[i]);							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, -1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1];
		result[i] /= (double)length; 
	}			

	//===Clean Up===//
	free(XFFT);

	return;
}

void fft_malloc_cmplx( double complex* data,
		   	  	       unsigned int length,
		   	  	       double complex* result )
{
	unsigned int i;
	double *XFFT;

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- fft_malloc_cmplx\n");
		quit();
	}			

	//===Malloc===//
	XFFT = malloc(2*(length+1)*sizeof(double));						

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = creal(data[i]);							
		XFFT[2*i+2] = cimag(data[i]);							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, 1);								

	//===Convert To Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1] + I*XFFT[2*i+2];
	}			

	//===Clean Up===//
	free(XFFT);

	return;
}

void ifft_malloc_cmplx( double complex *data, 
		   	   		    int length,
			   		    double complex* result )
{
	int i;
	double* XFFT;

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- ifft_malloc_cmplx\n");
		quit();
	}									

	//===Malloc===//
	XFFT = malloc(2*(length+1)*sizeof(double));						

	//===Convert to a Double Array===//
	for(i=0; i<length; i++){
		XFFT[2*i+1] = creal(data[i]);							
		XFFT[2*i+2] = cimag(data[i]);							
	}

	//===Calculate the FFT or IFFT===//
	fft_41(XFFT, length, -1);								

	//===Convert to Complex===//
	for(i=0; i<length; i++){
		result[i] = XFFT[2*i+1] + I*XFFT[2*i+2];
		result[i] /= (double complex) length; 
	}			

	//===Clean Up===//
	free(XFFT);

	return;
}

void real_fft_malloc_dbl( double* data,
				   		  unsigned int length,
				   		  double complex* result )
{
	unsigned int i;
	double* XFFT;
	double complex* temp_result;

	//===Mallocs===//
	XFFT = malloc((length+1)*sizeof(double));
	temp_result = malloc((length+1)*sizeof(double complex));

	//===Init===//
	initialize_array_dbl(XFFT, length+1);
	copy_array_dbl(data, length, XFFT + 1);

	//===Take FFT===//
	fft_41_real(XFFT, length, 1);

	//===Save===//
	temp_result[0] = XFFT[1];
	temp_result[length/2] = XFFT[2];
	for (i=1; i<length/2; i++){
		if (i == length/4){
			temp_result[i] = XFFT[2*i+1] - I*XFFT[2*i+2]; //compensate for isign for python compliance
		}
		else{
			temp_result[i] = XFFT[2*i+1] + I*XFFT[2*i+2]; 
		}
	}

	//===Make Full===//
	convert_real_fft_to_cmplx_malloc(temp_result, length, result);

	//===Clean Up===//
	free(XFFT);

	return;
}


void chirp_fft_dbl( double* signal,
					unsigned int signal_length,
					double reference_frequency,
					double frequency_spacing,
					unsigned int num_frequency_nodes,
					double complex* chirp_fft )
{

	//===Make Impulse Response===//
	unsigned int n, M, N;
	double delta_w, n_dbl, N_dbl, omega_0;
	double complex W;
	double complex* h1;
	double complex* g;
	double complex* output;

	//===Mallocs===//
	h1 = malloc(next_pow_2(2*(2*signal_length+num_frequency_nodes-1))*sizeof(double complex));
	g = malloc(next_pow_2(2*(2*signal_length+num_frequency_nodes-1))*sizeof(double complex));
	output = malloc(next_pow_2(2*(2*signal_length+num_frequency_nodes-1))*sizeof(double complex));

	//===Set Up===//
	delta_w = frequency_spacing; W = cexp(-I*delta_w); omega_0 = reference_frequency;
	M = num_frequency_nodes; N = signal_length; N_dbl = (double)N;

	//===Create Chirp Impulse Response===//
	initialize_array_cmplx(h1, next_pow_2(2*(M+N-1)));
	for (n=0; n<M+N-1; n++){
		n_dbl = (double)n; 
		h1[n] = cpow(W, -0.5*(n_dbl-N_dbl+1.0)*(n_dbl-N_dbl+1.0));
	}

	//===Make Input Vector===//
	initialize_array_cmplx(g, next_pow_2(2*(M+N-1)));
	for (n=0; n<N; n++){
		n_dbl = (double)n;
		g[n] = signal[n] * cexp(-I*omega_0*n) * cpow(W, 0.5*n_dbl*n_dbl);
	}

	//===Convolve===//
	initialize_array_cmplx(output, next_pow_2(2*(M+N-1)));
	convolve_fft_malloc_cmplx(h1, M+N-1, g, N, output);

	//===Weight===//
	for (n=0; n<M+N-1; n++){
		output[n] *= cpow(W, 0.5*(n-N_dbl+1.0)*(n-N_dbl+1.0));
	}

	//===Copy Over===//
	copy_array_cmplx(output + N - 1, M, chirp_fft);

	//===Clean Up===//
	free(output);
	free(g);
	free(h1);

	return;
} 


//------------------------------------------------------------------------------------------------//
//==================================FFT SUPPORT FUNCTIONS=========================================//
//------------------------------------------------------------------------------------------------//

void convert_real_fft_to_cmplx_malloc( double complex* real_fft,
									   unsigned int fft_length,
									   double complex* fft )
{
	unsigned int i, j, real_length;
	double complex* imaginary_fft;

	//===Mallocs===//
	imaginary_fft = malloc(fft_length*sizeof(double complex));
	
	//===Get Imaginary===//
	real_length = fft_length/2 + 1;
	j = 0;
	for (i=1; i<real_length-1; i++){
		imaginary_fft[j++] = conj(real_fft[i]);
	}
	reverse_array_cmplx(imaginary_fft, j);

	//===Tack On===//
	copy_array_cmplx(real_fft, real_length, fft);
	for (i=0; i<j; i++){
		fft[real_length + i] = imaginary_fft[i];
	}

	//===Free===//
	free(imaginary_fft);

	return;
}

void frequencies( unsigned int length, 
				  double sampling_frequency,
				  double* frequencies )
{
	unsigned int i;
	double start;

	start = -sampling_frequency/2.0;
	for (i=0; i<length; i++){
		frequencies[i] = start;
		start += sampling_frequency/(length);
	}
	right_circular_shift_array_dbl(frequencies, length, length/2);

	return;
}

void get_nonnegative_frequencies( unsigned int num_nonnegative_nodes,
								  double sampling_rate,
								  double* nonnegative_frequencies )
{

	double temp[2*MAX_FRAME_LENGTH];
	if (num_nonnegative_nodes +1 > 2*MAX_FRAME_LENGTH){
		fprintf(stderr, "Error:: Num Nodes Is Too Large! In Function -- get_nonnegative_frequencies!\n");
		quit();
	}
	initialize_linspace_dbl(temp, num_nonnegative_nodes+1, 0, sampling_rate/2);
	copy_array_dbl(temp, num_nonnegative_nodes, nonnegative_frequencies);
	return;
}

void normalize_digital_frequencies( double* frequencies,
									unsigned int length,
									double sampling_rate )
{
	unsigned int i;
	for (i=0; i<length; i++){
		frequencies[i] = (2.0 * M_PI * frequencies[i]) / sampling_rate;
	}

	return;
}


//------------------------------------------------------------------------------------------------//
//======================================DFT FUNCTIONS=============================================//
//------------------------------------------------------------------------------------------------//

void dft_dbl( double* data,
			  unsigned int length,
			  double complex* out )
{
	int c, r;
	unsigned int num_cols, num_rows;
	double freq_nodes[MAX_SIGNAL_LENGTH];
	double complex data_cmplx[MAX_SIGNAL_LENGTH];
	double fr;
	double complex* dft;

	//===Allocate DFT===//
	dft = malloc(length*length*sizeof(double complex));
	initialize_array_cmplx(dft, length*length);

	//===Allocate Freq Nodes===//
	initialize_linspace_dbl(freq_nodes, length+1, -0.5, 0.5);
	right_circular_shift_array_dbl(freq_nodes, length, (length)/2);

	//===Create DFT Matrix===//
	num_cols = length; num_rows = length;
	for (r = 0; r<(int)num_rows; r++){
		//fr = ((double)r);
		//fr /= ((double)(num_rows));
		fr = freq_nodes[r];
		//===Compute Rth Row===//
		for (c = 0; c<(int)num_cols; c++){
			dft[c + r*num_cols] = cexp(-I * 2.0 * M_PI * fr  * ((double)c));
		}
	}	


	//===Make DFT Data===//
	combine_arrays_dbl_to_cmplx(data, NULL, length, data_cmplx);
	
	//===Perform Multiplication===//
	matrix_vector_multiply_cmplx(data_cmplx, length, dft, length, length, out);
	if (0){	
		conjugate_array_cmplx(out, length);
	}

	//===Clean Up===//
	free(dft);

	return;
}

void idft_dbl( double complex* data,
			   unsigned int length,
			   double* out )
{
	int i, c, r;
	unsigned int num_cols, num_rows;
	double freq_nodes[MAX_SIGNAL_LENGTH];
	double complex out_cmplx[MAX_SIGNAL_LENGTH];
	double fr;
	double complex* idft;

	//===Allocate DFT===//
	idft = malloc(length*length*sizeof(double complex));
	initialize_array_cmplx(idft, length*length);

	//===Allocate Freq Nodes===//
	frequencies(length, 1.0, freq_nodes);

	//===Create DFT Matrix===//
	num_cols = length; num_rows = length;
	for (r = 0; r<(int)num_rows; r++){
		//fr = ((double)r)/((double)num_rows);		
		fr = freq_nodes[r];
		//===Compute Rth Row===//
		for (c = 0; c<(int)num_cols; c++){
			idft[c + r*num_cols] = cexp(I * 2.0 * M_PI * fr * (double)c);
		}
	}	
	
	//===Perform Multiplication===//
	matrix_vector_multiply_cmplx(data, length, idft, length, length, out_cmplx);

	//===Split===//
	split_array_cmplx(out_cmplx, length, out, NULL);
	for (i=0; i<(int)length; i++) out[i] /= (double)length;
	if (0){
		reverse_array_dbl(out, length);
		right_circular_shift_array_dbl(out, length, 1);
	}

	//===Clean Up===//
	free(idft);

	return;
}

double complex run_goertzel_filter_dft( double* data,
										unsigned int length,
										double freq )
{
	unsigned int i;
	double q, q_prev, q_prev_2, ak;
	double complex dft;

	//NOTE: this algorithm suffers from severe roundoff error
	//eg: we will get real components even when only sinusoids (not cosines) are input!
	//therefore only appropriate for very small lengths (<=256 samples)

	//===Error Check===//
	if (fabs(freq) > 0.5 + 10.0*EPS){
		fprintf(stderr, "Error:: Frequency %lf Is Not Normalized! In Function -- run_goertzel_filter_dft!\n", fabs(freq));
		return 0;
	}

	//===Run Recursion===//
	q = 0; q_prev = 0; q_prev_2 = 0; ak = 0; dft = 0;
	ak = 2.0 * cos(2.0 * M_PI * freq);
	for (i=0; i<length+1; i++){
		q = ak * q_prev - q_prev_2 + (i<length)*data[i];
		q_prev_2 = q_prev;
		q_prev = q;
	}
	
	dft = q_prev - cexp(-I * 2.0 * M_PI * freq) * q_prev_2;
	//dft *= cexp(I * 2.0 * M_PI * freq * -((double)length)); //this is responsible for time direction when warping

	return dft;
}

double complex run_goertzel_filter_dft_cmplx( double complex* data,
											  unsigned int length,
											  double freq )
{
	unsigned int i;
	double complex q, q_prev, q_prev_2, ak;
	double complex dft;

	//NOTE: this algorithm suffers from severe roundoff error
	//eg: we will get real components even when only sinusoids (not cosines) are input!
	//therefore only appropriate for very small lengths (<=256 samples)

	//===Error Check===//
	if (fabs(freq) > 0.5){
		fprintf(stderr, "Error:: Frequency %lf Is Not Normalized! In Function -- run_goertzel_filter_dft_cmplx!\n", fabs(freq));
		return 0;
	}

	//===Run Recursion===//
	q = 0; q_prev = 0; q_prev_2 = 0; ak = 0; dft = 0;
	ak = 2.0 * cos(2.0 * M_PI * freq);
	for (i=0; i<length+1; i++){
		q = ak * q_prev - q_prev_2 + (i<length)*data[i];
		q_prev_2 = q_prev;
		q_prev = q;
	}
	dft = q_prev - cexp(-I * 2.0 * M_PI * freq) * q_prev_2;
	//dft *= cexp(I * 2.0 * M_PI * freq * -((double)length));
	return dft;
}


void goertzel_dft_dbl( double* data,
					   unsigned int length,
					   unsigned int num_frequency_nodes,
					   double complex* dft )
{
	unsigned int i;
	double freq_nodes[10*(MAX_POLYNOMIAL_ORDER+1)];

	//===Take Goertzel DFT===//
	initialize_array_cmplx(dft, num_frequency_nodes+1);
	frequencies(num_frequency_nodes, 1.0, freq_nodes);
	for (i=0; i<num_frequency_nodes; i++){
		dft[i] = run_goertzel_filter_dft(data, length, freq_nodes[i]);
	}

	return;
}


double run_goertzel_filter_psd( double* data,
								unsigned int length,
								double freq )
{
	unsigned int i;
	double q, q_prev, q_prev_2, ak, power;
	
	//===Error Check===//
	if (fabs(freq) > 0.5 + 10*EPS){
		fprintf(stderr, "Error:: Frequency %lf Is Not Normalized! In Function -- run_goertzel_filter_psd!\n", fabs(freq));
		return 0;
	}

	//===Run Recursion===//
	q = 0; q_prev = 0; q_prev_2 = 0; ak = 0; power = 0;
	ak = 2.0 * cos(2.0 * M_PI * freq);
	for (i=0; i<length+1; i++){
		q = ak * q_prev - q_prev_2 + (i<length)*data[i];
		q_prev_2 = q_prev;
		q_prev = q;
	}
	power = q_prev_2*q_prev_2 + q_prev*q_prev - ak*q_prev*q_prev_2;

	return power;
}

void goertzel_psd_dbl( double* data,
					   unsigned int length,
					   unsigned int num_frequency_nodes,
					   double* psd )
{
	unsigned int i;
	double freq_nodes[10*(MAX_POLYNOMIAL_ORDER+1)];

	//===Take Goertzel DFT===//
	initialize_array_dbl(psd, num_frequency_nodes+1);
	frequencies(num_frequency_nodes, 1.0, freq_nodes);
	for (i=0; i<num_frequency_nodes; i++){
		psd[i] = run_goertzel_filter_psd(data, length, freq_nodes[i]);
	}

	return;
}

void nudft_dbl( double* data,
			  	unsigned int length,
				double* freq_nodes,
				unsigned int num_freq_nodes,
			  	double complex* out )
{
	int c, r;
	unsigned int num_cols, num_rows;
	double complex* data_cmplx;
	double complex* nudft;

	//===Allocate DFT===//
	data_cmplx = malloc(length * sizeof(double complex));
	nudft = malloc(length * num_freq_nodes*sizeof(double complex));
	initialize_array_cmplx(nudft, length * num_freq_nodes);

	//===Create DFT Matrix===//
	num_cols = length; num_rows = num_freq_nodes;
	for (r = 0; r<(int)num_rows; r++){
		//===Compute Rth Row===//
		for (c = 0; c<(int)num_cols; c++){
			nudft[c + r*num_cols] = cexp(I * 2.0 * M_PI * freq_nodes[r] * -(double)c);
		}
	}	

	//===Make NUDFT Data===//
	combine_arrays_dbl_to_cmplx(data, NULL, length, data_cmplx);
	
	//===Perform Multiplication===//
	matrix_vector_multiply_cmplx(data_cmplx, length, nudft, num_freq_nodes, length, out);

	//===Clean Up===//
	free(nudft);
	free(data_cmplx);

	return;
}

void nuidft_dbl( double complex* data,
				 double* freq_nodes,
				 unsigned int num_freq_nodes,
			  	 unsigned int out_length,
			  	 double* out )
{

	int i, c, r;
	unsigned int num_cols, num_rows;
	double complex out_cmplx[MAX_SIGNAL_LENGTH];
	double fr;
	double complex* nuidft;

	//===Allocate DFT===//
	nuidft = malloc(out_length*num_freq_nodes*sizeof(double complex));
	initialize_array_cmplx(nuidft, out_length*num_freq_nodes);
	fprintf(stdout, "NUIDFT Algorithm Needs to Be Fixed -- Dimensions are not appropriate for NxM transform!\n");

	if (0){
		//===Create DFT Matrix===//
		num_cols = out_length; num_rows = num_freq_nodes;
		for (r = 0; r<(int)num_rows; r++){
			fr = freq_nodes[r];
			//===Compute Rth Row===//
			for (c = 0; c<(int)num_cols; c++){
				nuidft[c + r*num_cols] = cexp(I * 2.0 * M_PI * fr * (double)c);
			}
		}	
	}
	else{
		num_cols = out_length; num_rows = num_freq_nodes;
		for (r = 0; r<(int)num_rows; r++){
			//===Compute Rth Row===//
			for (c = 0; c<(int)num_cols; c++){
				nuidft[c + r*num_cols] = cexp(I * 2.0 * M_PI * freq_nodes[r] * -(double)c);
			}
		}	
		hermitian_transpose_matrix_cmplx(nuidft, num_rows, num_cols);
		num_cols = num_freq_nodes; num_rows = out_length;
	}

	//===Perform Multiplication===//
	matrix_vector_multiply_cmplx(data, num_freq_nodes, nuidft, num_rows, num_cols, out_cmplx);

	//===Split===//
	split_array_cmplx(out_cmplx, out_length, out, NULL);
	for (i=0; i<(int)out_length; i++) out[i] /= (double)num_freq_nodes;

	//===Clean Up===//
	free(nuidft);

	return;
}

void nudft_cmplx( double complex* data,
			  	  unsigned int length,
				  double* freq_nodes,
				  unsigned int num_freq_nodes,
			  	  double complex* out )
{
	unsigned int c, r;
	unsigned int num_cols, num_rows;
	double complex* nudft;

	//===Allocate DFT===//
	nudft = malloc(length * num_freq_nodes*sizeof(double complex));
	initialize_array_cmplx(nudft, length * num_freq_nodes);

	//===Create DFT Matrix===//
	num_cols = length; num_rows = num_freq_nodes;
	for (r = 0; r<num_rows; r++){
		//===Compute Rth Row===//
		for (c = 0; c<num_cols; c++){
			nudft[c + r*num_cols] = cpow(cexp(I * 2.0 * M_PI * freq_nodes[r] ), -c);
		}
	}	
	
	//===Perform Multiplication===//
	matrix_vector_multiply_cmplx(data, length, nudft, num_freq_nodes, length, out);

	//===Clean Up===//
	free(nudft);

	return;
}

void goertzel_nudft_dbl( double* data,
					   	 unsigned int length,
						 double* freq_nodes,
					   	 unsigned int num_frequency_nodes,
					   	 double complex* nudft )
{
	unsigned int i;

	//===Take Goertzel NUDFT===//
	initialize_array_cmplx(nudft, num_frequency_nodes+1);
	for (i=0; i<num_frequency_nodes; i++){
		nudft[i] = run_goertzel_filter_dft(data, length, freq_nodes[i]);
	}

	return;
}






//------------------------------------------------------------------------------------------------//
//====================================OTHER TRANSFORMS============================================//
//------------------------------------------------------------------------------------------------//

void hilbert_transform_dbl( double* wave,
							unsigned int length,
							double* hilbert_transform )
{
	unsigned int k;
	double Nby2, k_dbl;
	double complex fft[2*(MAX_POLYNOMIAL_ORDER+1)];

	//===Sanity Check===//
	if (length > 2*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- hilbert_transform_dbl!\n");
		quit();
	}

	//===Take FFT===//	
	fft_dbl(wave, length, fft);

	//===Compute Hilbert Transform===//
	Nby2 = ((double)length)/2.0;
	for (k=0; k<length; k++){
		k_dbl = ((double)k);
		fft[k] *= I * sign_dbl(Nby2 - k_dbl) * sign_dbl(k_dbl);
	}
	ifft_dbl(fft, length, hilbert_transform);
	
	return;
}

void compute_analytic_signal_dbl( double* data,
								  unsigned int length,
								  double complex* analytic )
{
	unsigned int i;
	double complex h[MAX_SIGNAL_LENGTH];
	double complex temp[MAX_SIGNAL_LENGTH];

	if (!ISPOW2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- compute_analytic_signal_dbl\n");
		quit();
	}									

	//===Make Unit Step===//
	initialize_array_cmplx(h, length);
	h[0] = h[length/2] = 1.0;
	for (i=1; i<length/2; i++) h[i] = 2.0;

	//===Zero Imaginary FFT===//
	fft_dbl(data, length, analytic);
	hadamard_product_cmplx(analytic, h, length, temp);

	//===Take IFFT===//
	ifft_cmplx(temp, length, analytic);
	conjugate_array_cmplx(analytic, length);
	

	return;
}

void dct_dbl( double* data,
			  unsigned int length,
			  double* dct )
{
	unsigned int k, i;
	double angle;
	for (k=0; k<length; k++){
		dct[k] = 0.0;
		for (i=0; i<length; i++){
			angle = M_PI * ((double)(k * (2*i+1)))/((double)(2*length));
			dct[k] += cos(angle) * data[i];		
		}
		dct[k] *= 2.0;
		if (k == 0) dct[k] *= sqrt(1.0/4.0*((double)length));
		else dct[k] *= sqrt(1.0/2.0*((double)length));
		dct[k] /= (double)(length);
	}
	return;
}

void dst_dbl( double* data,
			  unsigned int length,
			  double* dct )
{
	unsigned int k, i;
	double angle;
	for (k=0; k<length; k++){
		dct[k] = 0.0;
		for (i=0; i<length; i++){
			angle = M_PI * ((double)(k * (2*i+1)))/((double)(2*length));
			dct[k] += sin(angle) * data[i];		
		}
		dct[k] *= 2.0;
		if (k == 0) dct[k] *= sqrt(1.0/4.0*((double)length));
		else dct[k] *= sqrt(1.0/2.0*((double)length));
		dct[k] /= (double)(length);
	}
	return;
}

//------------------------------------------------------------------------------------------------//
//======================================Z-TRANSFORMS==============================================//
//------------------------------------------------------------------------------------------------//

void z_transform_dbl( double* data,
					  unsigned int length,
					  double radius,
					  double complex* z_transform )
{
	unsigned int n;
	double ndbl;
	double temp[4*(MAX_POLYNOMIAL_ORDER+1)];

	if (length > 4*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- z_transform_dbl!\n");
		quit();
	}

	//===Make Radial Signal===//
	for (n=0; n<length; n++){
		ndbl = n;
		temp[n] = data[n] * pow(radius, -ndbl);
	}

	//===Take FFT===//
	fft_dbl(temp, length, z_transform);

	return;
}

void z_transform_cmplx( double complex* data,
					  	unsigned int length,
					  	double radius,
					  	double complex* z_transform )
{
	unsigned int n;
	double ndbl;
	double complex temp[4*(MAX_POLYNOMIAL_ORDER+1)];

	if (length > 4*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- z_transform_dbl!\n");
		quit();
	}

	//===Make Radial Signal===//
	for (n=0; n<length; n++){
		ndbl = n;
		temp[n] = data[n] * cpow(radius, -ndbl);
	}

	//===Take FFT===//
	fft_cmplx(temp, length, z_transform);

	return;
}

void inverse_z_transform_dbl( double complex* data,
					  		  unsigned int length,
					  		  double radius,
					  		  double* inverse_z_transform )
{

	int n;
	double ndbl, N, node;
	double complex temp_cmplx[4*(MAX_POLYNOMIAL_ORDER+1)];

	if (length > 4*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- inverse_z_transform_dbl!\n");
		quit();
	}

	//===Make X(z) * e^jw Signal===//
	N = (double)length;
	for (n=0; n<(int)length; n++){
		ndbl = n;
		node = ndbl/N;
		temp_cmplx[n] = data[n] * cexp(-I * 2.0 * M_PI * node);
	}

	//===Take IFFT of X(z) * e^jw===//
	ifft_dbl(temp_cmplx, length, inverse_z_transform);

	//===Scale by (r^(n-1))/N===//
	for (n=0; n<(int)length; n++) inverse_z_transform[n] *= (pow(radius, n-1));
	right_circular_shift_array_dbl(inverse_z_transform, length, 1);

	return;
}

void z_transform_derivative_dbl( double* data,
					  			 unsigned int length,
					  			 double radius,
					  			 double complex* z_transform_derivative )
{

	int n;
	double ndbl, N;
	double temp[4*(MAX_POLYNOMIAL_ORDER+1)];
	double nodes[4*(MAX_POLYNOMIAL_ORDER+1)];

	if (length > 4*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- z_transform_derivative_dbl!\n");
		quit();
	}

	//===Make Ramped Signal n*x===//
	for (n=0; n<(int)length; n++){
		ndbl = n;
		temp[n] = -ndbl * data[n];
	}

	//===Take Z Transform of Ramped Signal n*x===//
	z_transform_dbl(temp, length, radius, z_transform_derivative);

	//===Multiply by -(z_{k})^{-1}===//
	N = (double)length;
	for (n=0; n<(int)length; n++){
		ndbl = n;
		nodes[n] = ndbl/N;
		z_transform_derivative[n] *= I * 1.0/radius * cexp(-I * 2.0 * M_PI * nodes[n]);
	}
	return;
}

void logarithmic_derivative_dbl( double* data,
					  			 unsigned int length,
					  			 double radius,
					  			 double complex* logarithmic_derivative )
{
	unsigned int i;
	double complex dXz[4*(MAX_POLYNOMIAL_ORDER+1)];
	double complex Xz[4*(MAX_POLYNOMIAL_ORDER+1)];

	if (length > 4*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- logarithmic_derivative_dbl!\n");
		quit();
	}

	//===Calculate Z Transform===//
	z_transform_dbl(data, length, radius, Xz);

	//===Calculate Derivative Of Z Transform===//
	z_transform_derivative_dbl(data,length,radius,dXz);

	//===Form Logarithmic Derivative===//
	for (i=0; i<length; i++) logarithmic_derivative[i] = dXz[i]/Xz[i];

	return;
}

void calculate_ztransform_group_delay_dbl( double* data,
									   	   unsigned int length,
										   double radius,
										   double* group_delay )
{
	unsigned int i;
	double complex logarithmic_derivative[2*MAX_FRAME_LENGTH];

	//===Calculate Logarithmic Derivative===//
	logarithmic_derivative_dbl(data, length, radius, logarithmic_derivative);

	//===Extract Group Delay===//
	for (i=0; i<length; i++){group_delay[i] = -1.0 * cimag(logarithmic_derivative[i]);}

	return;
}

//------------------------------------------------------------------------------------------------//
//======================================PHASE FUNCTIONS===========================================//
//------------------------------------------------------------------------------------------------//

void compute_phase_response_dbl( double* data,
							 	 unsigned int length,
							 	 double complex* phase )
{
	unsigned int i;
	double complex fft[MAX_SIGNAL_LENGTH];
	double angle[MAX_SIGNAL_LENGTH];
	if (!IS_POW_2(length)){
		fprintf(stderr, "Error: FFT Length Is Not A Power Of 2! In Function --- compute_phase_response_dbl\n");
		quit();
	}
	fft_dbl(data,length,fft);
	for (i=0; i<length; i++){
		angle[i] = atan2(cimag(fft[i]),creal(fft[i]));
		phase[i] = cexp(1.0*I*angle[i]);
	}
	

	return;
}

void unwrap_phase_angles( double* phase,
						  unsigned int length,
						  double cutoff )
{
	unsigned int i;
	double phase_diff[4*(MAX_POLYNOMIAL_ORDER+1)];
	double phase_diff_wrapped[4*(MAX_POLYNOMIAL_ORDER+1)];
	double cumulative_sum[4*(MAX_POLYNOMIAL_ORDER+1)];		

	//===Test For Length===//
	if (length > 4*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Phase Array Too Long! In Function -- unwrap_phase_angles!\n");
		return;
	}

	for (i=0; i<length-1; i++){

		//===Differentiate===//
		phase_diff[i] = phase[i+1] - phase[i];

		//===Map To [-pi, pi]===//
		phase_diff_wrapped[i] = phase_diff[i]+M_PI;
		phase_diff_wrapped[i] -= floor((phase_diff[i]+M_PI) / (2.0*M_PI))*(2.0*M_PI);
		phase_diff_wrapped[i] -= M_PI;
		if ( fabs(phase_diff_wrapped[i] - (-1.0*M_PI)) < 2.0*EPS && phase_diff[i] > 0){
			phase_diff_wrapped[i] = M_PI; 
		}

		//===Correct Phase Difference===//
		phase_diff[i] = phase_diff_wrapped[i] - phase_diff[i];
		if (fabs(phase_diff[i]) < cutoff){
			phase_diff[i] = 0.0;
		}

		//===Sum Phase Differences===//
		if (i==0){
			cumulative_sum[i] = phase_diff[i];
		}
		else{
			cumulative_sum[i] = cumulative_sum[i-1] + phase_diff[i];
		}

	}

	//===Correct Phases===//
	for (i=1; i<length; i++) phase[i] += cumulative_sum[i-1];

	return;
}				   

void compute_angle_response_dbl( double* data,
							 	 unsigned int length,
							 	 double* angles )
{
	unsigned int i;
	double complex fft[MAX_FRAME_LENGTH];
	double freqs[MAX_FRAME_LENGTH];

	if (length > MAX_FRAME_LENGTH){
		fprintf(stderr, "Error: Data Length Is Too Large! In Function --- compute_angle_response_dbl\n");
		quit();
	}
	
	//===Calculate FFT===//
	if (!IS_POW_2(length)){
		initialize_linspace_dbl(freqs, length+1, -0.5, 0.5);
		for (i=0; i<length; i++){
			fft[i] = run_goertzel_filter_dft(data, length, freqs[i]);
		}
		right_circular_shift_array_cmplx(fft, length, length/2);
	}
	else{
		fft_dbl(data,length,fft);
	}
	
	//===Get Angles===//
	for (i=0; i<length; i++){
		angles[i] = atan2(cimag(fft[i]),creal(fft[i]));
	}

	return;
}

void phase_shift_dbl( double* data,
					  unsigned int length,
					  double delay_seconds,
					  double sampling_rate )
{
	unsigned int i, fft_length;
	double complex fft[2*MAX_FRAME_LENGTH];
	double freqs[2*MAX_FRAME_LENGTH];
	double data_temp[2*MAX_FRAME_LENGTH];

	//===Sanity Check===//
	if (!IS_POW_2(length)){
		fprintf(stdout, "Error:: Length Is Not A Power of 2! In Function -- phase_shift_dbl!\n");
		quit();
	}
	if (2*length >= 2*MAX_FRAME_LENGTH){
		fprintf(stdout, "Error:: Length Is Too Long! In Function -- phase_shift_dbl!\n");
		quit();		
	}

	//===Copy Over===//
	fft_length = 2*length;
	initialize_array_dbl(data_temp, fft_length);
	copy_array_dbl(data, fft_length, data_temp);

	//===Calculate FFT===//
	frequencies(fft_length, 1.0, freqs);
	fft_dbl(data_temp,fft_length,fft);

	//===Gain===//
	gain_array_constant_dbl(freqs,fft_length, sampling_rate);
	
	//===Multiply===//
	for (i=0; i<fft_length; i++){
		fft[i] *= cexp(-I * 2.0 * M_PI * freqs[i] * delay_seconds); //delay is assumed positive
	}

	//===Recreate Signal===//
	ifft_dbl(fft, fft_length, data_temp);
	copy_array_dbl(data_temp, length, data);

	return;
}

void phase_shift_malloc_dbl( double* data,
						     unsigned int length,
					  		 double delay_seconds,
					  		 double sampling_rate )
{
	unsigned int i, fft_length;
	double complex* fft;
	double* freqs, *data_temp;

	//===Sanity Check===//
	if (!IS_POW_2(length)){
		fprintf(stdout, "Error:: Length Is Not A Power of 2! In Function -- phase_shift_malloc_dbl!\n");
		quit();
	}

	//===Malloc===//
	fft = malloc(2*(length+1)*sizeof(double complex));
	freqs = malloc(2*(length+1)*sizeof(double));
	data_temp = malloc(2*(length+1)*sizeof(double));

	//===Copy Over===//
	fft_length = 2*length;
	initialize_array_dbl(data_temp, fft_length);
	copy_array_dbl(data, fft_length, data_temp);

	//===Calculate FFT===//
	frequencies(fft_length, 1.0, freqs);
	fft_malloc_dbl(data_temp,fft_length,fft);

	//===Gain===//
	gain_array_constant_dbl(freqs,fft_length, sampling_rate);
	
	//===Multiply===//
	for (i=0; i<fft_length; i++){
		fft[i] *= cexp(-I * 2.0 * M_PI * freqs[i] * delay_seconds); //delay is assumed positive
	}

	//===Recreate Signal===//
	ifft_malloc_dbl(fft, fft_length, data_temp);
	copy_array_dbl(data_temp, length, data);

	//===Clean Up===//
	free(data_temp);
	free(freqs);
	free(fft);

	return;
}

//------------------------------------------------------------------------------------------------//
//=====================================MAGNITUDE FUNCTIONS========================================//
//------------------------------------------------------------------------------------------------//

void compute_magnitude_response_dbl( double* data,
							 		 unsigned int length,
							 		 double* magnitude )
{
	unsigned int i;
	double freq_nodes[MAX_SIGNAL_LENGTH];
	double complex fft[MAX_SIGNAL_LENGTH];

	//===Take FFT===//
	if (!IS_POW_2(length)){
		frequencies(length, 1.0, freq_nodes);
		for (i=0; i<length; i++){
			fft[i] = run_goertzel_filter_dft(data, length, freq_nodes[i]);
		}
	}
	else{
		fft_dbl(data,length,fft);
	}
	get_magnitude_array_cmplx(fft, length, magnitude);
	return;
}

void compute_log_magnitude_response_dbl( double* data,
								 	 	 unsigned int length,
								 	 	 double* log_magnitude )
{
	unsigned int i;
	double magnitude[MAX_SIGNAL_LENGTH];
	compute_magnitude_response_dbl(data,length,magnitude);
	for (i=0; i<length; i++){
		log_magnitude[i] = log10(pow(fabs(magnitude[i]),2.0));
	}

	return;
}

void compute_power_spectrum_dbl( double* data, 
							 	 unsigned int length,
							 	 double complex* psd_cmplx,
								 double* psd_dbl )
{
	unsigned int i;
	double complex temp;
	double freq_nodes[MAX_SIGNAL_LENGTH];
	double complex fft[MAX_SIGNAL_LENGTH];

	//===Take FFT===//
	if (!IS_POW_2(length)){
		frequencies(length, 1.0, freq_nodes);
		for (i=0; i<length; i++){
			fft[i] = run_goertzel_filter_dft(data, length, freq_nodes[i]);
		}
	}
	else{
		if (length > (MAX_SIGNAL_LENGTH-1)/2){
			fft_malloc_dbl(data,length,fft);
		}
		else{
			fft_dbl(data,length,fft);
		}
	}

	//===Compute Power Spectrum===//
	for (i=0; i<length; i++){
		if (psd_cmplx != NULL){
			psd_cmplx[i] = fft[i] * conj(fft[i]);
		}
		if (psd_dbl != NULL){
			temp = fft[i] * conj(fft[i]);
			psd_dbl[i] = cabs(temp);
		}
	}
	return;
}

void compute_power_spectrum_malloc_dbl( double* data, 
							 	 		unsigned int length,
							 	 		double complex* psd_cmplx,
								 		double* psd_dbl )
{
	unsigned int i;
	double complex temp;
	double* freq_nodes;
	double complex* fft;
	
	//===Mallocs===//
	freq_nodes = malloc(length*sizeof(double));
	fft = malloc(length*sizeof(double complex));

	//===Take FFT===//
	if (!IS_POW_2(length)){
		frequencies(length, 1.0, freq_nodes);
		for (i=0; i<length; i++){
			fft[i] = run_goertzel_filter_dft(data, length, freq_nodes[i]);
		}
	}
	else{
		if (length > (MAX_SIGNAL_LENGTH-1)/2){
			fft_malloc_dbl(data,length,fft);
		}
		else{
			fft_dbl(data,length,fft);
		}
	}

	//===Compute Power Spectrum===//
	for (i=0; i<length; i++){
		if (psd_cmplx != NULL){
			psd_cmplx[i] = fft[i] * conj(fft[i]);
		}
		if (psd_dbl != NULL){
			temp = fft[i] * conj(fft[i]);
			psd_dbl[i] = cabs(temp);
		}
	}

	//===Clean Up===//
	free(fft);
	free(freq_nodes);
	
	return;
}


void compute_power_spectrum_normalized_dbl( double* data, 
							 	 			unsigned int length,
							 	 			double complex* psd_cmplx,
								 			double* psd_dbl )
{
	unsigned int i;
	double complex temp;
	double freq_nodes[MAX_SIGNAL_LENGTH];
	double complex fft[MAX_SIGNAL_LENGTH];

	//===Sanity Check===//
	if (length > MAX_SIGNAL_LENGTH){
		fprintf(stdout, "Error:: Length Is Too Large! In Function -- compute_power_spectrum_normalized_dbl!\n");
		quit();
	}

	//===Take FFT===//
	if (!IS_POW_2(length)){
		frequencies(length, 1.0, freq_nodes);
		for (i=0; i<length; i++){
			fft[i] = run_goertzel_filter_dft(data, length, freq_nodes[i]);
		}
	}
	else{
		fft_dbl(data,length,fft);
	}

	//===Compute Power Spectrum===//
	for (i=0; i<length; i++){
		if (psd_cmplx != NULL){
			psd_cmplx[i] = fft[i] * conj(fft[i]);
			psd_cmplx[i] /= ((double)(length*length));
		}
		if (psd_dbl != NULL){
			temp = fft[i] * conj(fft[i]);
			psd_dbl[i] = cabs(temp);
			psd_dbl[i] /= ((double)(length*length));
		}
	}
	return;
}


//------------------------------------------------------------------------------------------------//
//====================================GROUP DELAY FUNCTIONS=======================================//
//------------------------------------------------------------------------------------------------//

void calculate_group_delay_dbl( double* data,
								unsigned int length,
								double sampling_frequency,
								double* group_delay )
{
	unsigned int i;
	double T, t;
	double freqs[2*MAX_FRAME_LENGTH];
	double complex ramped_data[2*MAX_FRAME_LENGTH];
	double complex fft[2*MAX_FRAME_LENGTH];
	double complex dfft[2*MAX_FRAME_LENGTH];

	//===Calculate Ramped Data===//
	T = 1.0/sampling_frequency;
	t = 0;
	for(i=0; i<length; i++){
		ramped_data[i] = -I * t * data[i];
		t += T;
	}

	//===Calculate Group Delay Via FFT===//
	if (IS_POW_2(length)){
		//===Calculate FFT===//
		fft_dbl(data, length, fft);

		//===Calculate FFT Derivative===//
		fft_cmplx(ramped_data, length, dfft);
	}
	else{
		//===Calculate FFT===//
		initialize_linspace_dbl(freqs, length+1, -0.5, 0.5);
		for (i=0; i<length; i++){
			fft[i] = run_goertzel_filter_dft(data, length, freqs[i]);
		}
		right_circular_shift_array_cmplx(fft, length, length/2);

		//===Calculate FFT Derivative===//
		for (i=0; i<length; i++){
			dfft[i] = run_goertzel_filter_dft_cmplx(ramped_data, length, freqs[i]);
		}
		right_circular_shift_array_cmplx(dfft, length, length/2);
	}

	//===Get Group Delay===//
	for(i=0; i<length; i++){
		group_delay[i] = -1.0 * cimag(dfft[i]/fft[i]);
	}

	return;
}

void calculate_warped_group_delay_dbl( double* data,
									   unsigned int length,
									   double sampling_frequency,
									   double* frequency_nodes,
									   unsigned int num_freq_nodes,
									   double* group_delay )
{

	unsigned int i;
	double T, t;
	double complex ramped_data[2*MAX_FRAME_LENGTH];
	double complex fft[2*MAX_FRAME_LENGTH];
	double complex dfft[2*MAX_FRAME_LENGTH];

	if (num_freq_nodes > 2*MAX_FRAME_LENGTH){
		fprintf(stderr, "Error:: Too Many Frequency Nodes!\n");
		quit();
	}

	//===Calculate Ramped Data===//
	T = 1.0/sampling_frequency;
	t = 0;
	for(i=0; i<length; i++){
		ramped_data[i] = -I * t * data[i];
		t += T;
	}

	//===Calculate FFT===//
	for (i=0; i<num_freq_nodes; i++){
		fft[i] = run_goertzel_filter_dft(data, length, frequency_nodes[i]);
	}

	//===Calculate FFT Derivative===//
	for (i=0; i<num_freq_nodes; i++){
		dfft[i] = run_goertzel_filter_dft_cmplx(ramped_data, length, frequency_nodes[i]);
	}

	//===Get Group Delay===//
	for(i=0; i<num_freq_nodes; i++){
		group_delay[i] = -1.0 * cimag(dfft[i]/fft[i]);
	}

	return;
}


//------------------------------------------------------------------------------------------------//
//====================================CONVOLUTION FUNCTIONS=======================================//
//------------------------------------------------------------------------------------------------//

void convolve_dbl( double *array1, 
			   	   unsigned int length1, 
			   	   double *array2, 
			   	   unsigned int length2, 
			   	   double *result)
{
	unsigned int convolution_length;
	unsigned int i, j, k;
	double temp;

	//===Length Of Convolution===//	
	convolution_length = length1+length2-1;

	//===Run Convolution===//
	for (i=0; i<convolution_length; i++){
		k = i;
		temp = 0.0;
		for (j=0; j<length2; j++){
			if((int)k>=0 && k<length1){
				temp = temp + (array1[k]*array2[j]);
			}
			k = k-1;
			result[i] = temp;
		}
	}
		
	return;
}

void convolve_cmplx( double complex *array1, 
			   		 unsigned int length1, 
			   		 double complex *array2, 
			   		 unsigned int length2, 
			   		 double complex *result )
{
	unsigned int convolution_length;
	unsigned int i, j, k;
	double complex temp;

	//===Length Of Convolution===//	
	convolution_length = length1+length2-1;

	//===Run Convolution===//
	for (i=0; i<convolution_length; i++){
		k = i;
		temp = 0.0;
		for (j=0; j<length2; j++){
			if((int)k>=0 && k<length1){
				temp = temp + (array1[k]*array2[j]);
			}
			k = k-1;
			result[i] = temp;
		}
	}
		
	return;
}

void convolve_fft_malloc_dbl( double* array1,
						 	  unsigned int length1,
						 	  double* array2,
						 	  unsigned int length2,
						 	  double* result)
{
	unsigned int fft_length;
	double complex *temp1, *temp2, *temp12;

	//===Mallocs===//
	fft_length = next_pow_2(2*(length1+length2-1));
	temp1 = malloc(fft_length*sizeof(double complex));
	temp2 = malloc(fft_length*sizeof(double complex));
	temp12 = malloc(fft_length*sizeof(double complex));

	//===Pad With Zeros===//
	pad_zeros_dbl(array1, length1, fft_length-length1); 
	pad_zeros_dbl(array2, length2, fft_length-length2); 
	
	//===Take FFTs===//
	real_fft_malloc_dbl(array1, fft_length, temp1);
	real_fft_malloc_dbl(array1, fft_length, temp2);

	//===Multiply===//
	hadamard_product_cmplx(temp1, temp2, fft_length, temp12);

	//===Take IFFT===//
	ifft_malloc_dbl(temp12, fft_length, result);

	//===Clean Up===//
	free(temp12);
	free(temp2);
	free(temp1);

	return;
}


void convolve_fft_malloc_cmplx( double complex* array1,
						 		unsigned int length1,
						 		double complex* array2,
						 		unsigned int length2,
						 		double complex* result)
{
	unsigned int fft_length;
	double complex *temp1, *temp2, *temp12;

	//===Mallocs===//
	fft_length = next_pow_2(2*(length1+length2-1));
	temp1 = malloc(fft_length*sizeof(double complex));
	temp2 = malloc(fft_length*sizeof(double complex));
	temp12 = malloc(fft_length*sizeof(double complex));

	//===Pad With Zeros===//
	pad_zeros_cmplx(array1, length1, fft_length-length1); 
	pad_zeros_cmplx(array2, length2, fft_length-length2); 
	
	//===Take FFTs===//
	fft_malloc_cmplx(array1, fft_length, temp1);
	fft_malloc_cmplx(array2, fft_length, temp2);

	//===Multiply===//
	hadamard_product_cmplx(temp1, temp2, fft_length, temp12);

	//===Take IFFT===//
	ifft_malloc_cmplx(temp12, fft_length, result);

	//===Clean Up===//
	free(temp12);
	free(temp2);
	free(temp1);

	return;
}


void normalized_correlation_dbl( double* array1,
								 unsigned int length1,
								 double* array2,
								 unsigned int length2,
								 double* result )
{
	unsigned int i, result_length;
	double maximum;
	double* copy;

	//===Copy Array===//
	copy = malloc(length2 * sizeof(double));
	copy_array_dbl(array2, length2, copy);

	//===Reverse Convolve===//
	reverse_array_dbl(copy, length2);
	convolve_dbl(array1, length1, copy, length2, result);

	//===Normalize===//
	result_length = length1+length2-1;
	maximum = result[result_length/2];
	for (i=0; i<result_length; i++) result[i] /= maximum;

	//===Clean Up===//
	free(copy);

	return;
}

void correlation_cmplx( double complex* array1,
						unsigned int length1,
						double complex* array2,
						unsigned int length2,
						double complex* result )
{
	double complex* copy;

	//===Copy Array===//
	copy = malloc((length2+1)*sizeof(double complex));
	copy_array_cmplx(array2, length2, copy);

	//===Reverse Convolve===//
	conjugate_array_cmplx(copy, length2);
	reverse_array_cmplx(copy, length2);
	convolve_cmplx(array1, length1, copy, length2, result);

	//===Clean Up===//
	free(copy);

	return;
}


//------------------------------------------------------------------------------------------------//
//====================================RESAMPLING FUNCTIONS========================================//
//------------------------------------------------------------------------------------------------//

double compute_lanczos_kernel_value_dbl( double x,
										 double a )
{
	double Lxa;

	if (fabs(x) <= 10.0*EPS){
		return 1.0;
	}
	else if (x >= -a && x <= a){
		Lxa = (a * sin(M_PI*x) * sin((M_PI*x)/a))/(M_PI*M_PI*x*x);
		return Lxa;
	}
	else{
		return 0.0;
	}

	return 0.0;
}

void resample_lanczos_dbl( double* values,
						   double* original_scale,
						   unsigned int original_length,
						   double* new_scale,
						   unsigned int new_length,
						   double kernel_size,
						   double* resampled )
{

	unsigned int i, j;
	unsigned int start_of_array, size, start, end;
	int start1, end1, size_temp;
	double min_kernel, max_kernel;
	double start_freq, end_freq;
	double kernel;
	double *original_scale_copy;
	double *new_scale_copy;

	//===Mallocs===//
	original_scale_copy = malloc(original_length*sizeof(double));
	new_scale_copy = malloc(new_length*sizeof(double));

	//===Copy Over===//
	copy_array_dbl(original_scale, original_length, original_scale_copy);
	copy_array_dbl(new_scale, new_length, new_scale_copy);

	//===Scale to [0, N]===//
	scale_array_dbl(original_scale_copy, original_length, 0, original_length);
	scale_array_dbl(new_scale_copy, new_length, 0, original_length);

	//===Run Resampling===//
	for (j=0; j<new_length; j++){

		//===Get Starting and Ending Frequencies===//
		start_freq = (MAX((floor(new_scale_copy[j]) - kernel_size + 1.0), 0));
		end_freq = MIN((floor(new_scale_copy[j]) + kernel_size), ((double)(new_scale_copy[new_length-1])));

		//===Initialize Search Bounds===//
		start1 = (int)(floor(new_scale_copy[j]) - kernel_size + 1.0);
		end1 = (int)(floor(new_scale_copy[j]) + kernel_size);

		//===Get True Search Range===//
		start_of_array = (unsigned int)MAX(start1,0);
		size_temp = 3*(end1-start1); size = size_temp;

		//===Find Starting and Ending Indices For Interpolation===//		
		start = find_closest_index_dbl(original_scale_copy + start_of_array , size, start_freq) + start_of_array;
		end = find_closest_index_dbl(original_scale_copy + start_of_array, size, end_freq) + start_of_array;

		//===Interpolate Value===//
		min_kernel = 10000000.0; max_kernel = 0.0;
		resampled[j] = 0.0; 
		for (i=start; i<end; i++){
			kernel = compute_lanczos_kernel_value_dbl(new_scale_copy[j] - original_scale_copy[i], kernel_size);
			resampled[j] += values[i] * kernel;
			if (fabs(kernel) < min_kernel){
				min_kernel = kernel;
			}
			if (fabs(kernel) > max_kernel){
				max_kernel = kernel;
			}
		}

		//fprintf(stdout, "Freq: %lf %lf %lf\n", new_scale_copy[j], min_kernel, max_kernel);

	}		

	//===Free===//
	free(new_scale_copy);
	free(original_scale_copy);

	return;
}

void resample_sinc_cmplx( double complex* values,
						  double* original_scale,
						  unsigned int original_length,
						  double* new_scale,
						  unsigned int new_length,
						  double complex* resampled )
{

	unsigned int i, j, start, end;
	unsigned int start1, end1;
	unsigned int start_of_array, end_of_array, size;
	double start_freq, end_freq;
	double kernel_size, D, omega;
	double complex kernel;
	double *original_scale_copy;
	double *new_scale_copy;

	//===Mallocs===//
	original_scale_copy = malloc(original_length*sizeof(double));
	new_scale_copy = malloc(new_length*sizeof(double));

	//===Copy Over===//
	copy_array_dbl(original_scale, original_length, original_scale_copy);
	copy_array_dbl(new_scale, new_length, new_scale_copy);

	//===Scale to [0, N]===//
	scale_array_dbl(original_scale_copy, original_length, 0, original_length);
	scale_array_dbl(new_scale_copy, new_length, 0, original_length);

	//===Set D and Kernel===//
	D = 4.0*lcm_dbl((double)original_length, (double)new_length);
	kernel_size = 10;

	//===Run Resampling===//
	for (j=0; j<new_length; j++){

		//===Get Starting and Ending Frequencies===//
		start_freq = (MAX((floor(new_scale_copy[j]) - kernel_size + 1.0), 0));
		end_freq = MIN((floor(new_scale_copy[j]) + kernel_size), ((double)(original_length)));

		//===Initialize Search Bounds===//
		start1 = (unsigned int)(MAX((floor(new_scale_copy[j]) - kernel_size + 1.0), 0));
		end1 = (unsigned int)MIN((floor(new_scale_copy[j]) + kernel_size), ((double)(original_length)));

		//===Get True Search Range===//
		start_of_array = start1 - (end1-start1);
		end_of_array = original_length - start_of_array;
		size = MIN(3*(end1-start1), end_of_array - start_of_array);

		//===Find Starting and Ending Indices For Interpolation===//		
		start = find_closest_index_dbl(original_scale_copy + start_of_array , size, start_freq) + start_of_array;
		end = find_closest_index_dbl(original_scale_copy + start_of_array, size, end_freq) + start_of_array;

		//fprintf(stdout, "Freq: %lf\n", new_scale_copy[j]);

		//===Interpolate Value===//
		resampled[j] = 0.0; 
		for (i=start; i<end; i++){
			omega = new_scale_copy[j] - original_scale_copy[i];
			kernel = cexp(I * D * omega) * (sin(D*omega*0.5)/(D*sin(omega*0.5)));
			resampled[j] += values[i] * kernel;
		}
	}		

	//===Free===//
	free(new_scale_copy);
	free(original_scale_copy);

	return;
}


//------------------------------------------------------------------------------------------------//
//======================================MIXING FUNCTIONS==========================================//
//------------------------------------------------------------------------------------------------//

void mix_signals_dbl( double* signal1,
					  unsigned int signal1_length,
					  double* signal2,
					  unsigned int signal2_length,
					  double SNR,
					  double* result )
{
	unsigned int length;

	//===Normalize===//
	normalize_vector_2norm_dbl(signal1,signal1_length);
	normalize_vector_2norm_dbl(signal2,signal2_length);

	//===Gain===//
	SNR = compute_linear_gain_constant_from_dB(SNR);
	SNR = sqrt(SNR);
	gain_array_constant_dbl(signal1, signal1_length, SNR);

	//===Add Together===//
	length = MIN(signal1_length, signal2_length);
	initialize_array_dbl(result, length);
	if (signal1_length >= signal2_length){
		copy_array_dbl(signal1, signal1_length, result);
		add_vectors_dbl(signal2, result, length, result);
	}
	else{
		copy_array_dbl(signal2, signal2_length, result);
		add_vectors_dbl(signal1, result, length, result);
	}

	return;
}

void mix_multiple_signals_dbl( double** signals,
							   unsigned int num_signals,
							   unsigned int signal_length,
							   double* signal_powers_dB,
							   double* result )
{
	unsigned int i;
	double power;
	double signal_copy[3*MAX_SIGNAL_LENGTH];

	if (signal_length > 3*MAX_SIGNAL_LENGTH){
		fprintf(stderr, "Error:: Signal Length Is Too Large! In Function -- mix_multiple_signals_dbl\n");
		quit();
	}

	//===Initialize Result===//
	initialize_array_dbl(result, signal_length);
	for (i=0; i<num_signals; i++){

		//===Copy Locally===//
		copy_array_dbl(signals[i], signal_length, signal_copy);

		//===Normalize Power===//
		normalize_vector_2norm_dbl(signal_copy, signal_length);

		//===Gain===//
		power = compute_linear_gain_constant_from_dB(signal_powers_dB[i]);
		power = sqrt(power);
		gain_array_constant_dbl(signal_copy, signal_length, power);

		//===Add Together===//
		add_vectors_dbl(signal_copy, result, signal_length, result);
	}

	return;
}


//------------------------------------------------------------------------------------------------//
//=================================INPUT GENERATION FUNCTIONS=====================================//
//------------------------------------------------------------------------------------------------//

void generate_wave_dbl( double *frequencies, 
						double *amplitudes,
						double* phase,
						unsigned int num_frequencies,
						unsigned int signal_length, 
						double sampling_rate, 
						double* wave )
{
	if (frequencies == NULL){
		fprintf(stderr, "Error: Frequencies Are NULL In Function -- generate_wave!\n");
		return;
	}
	if (!isnormal((float)num_frequencies)){
		fprintf(stderr, "Error: Num Frequencies Is Invalid In Function -- generate_wave!\n");
		return;
	}
	if (!isnormal((float)signal_length)){
		fprintf(stderr, "Error: Signal Length Is Invalid In Function -- generate_wave!\n");
		return;
	}
	if (!isnormal((float)sampling_rate)){
		fprintf(stderr, "Error: Sampling Rate Is Invalid In Function -- generate_wave!\n");
		return;
	}


	unsigned int i, j;
	double t;
	//===Create Wave===//
	for (j=0; j<signal_length; j++){
		wave[j] = 0; t = ((double)j)/sampling_rate;
		for (i=0; i<num_frequencies; i++){
			if (phase == NULL && amplitudes == NULL){
				wave[j] += 1.0* cos( (2.0 * M_PI * frequencies[i] * t) + 0.0);
			}
			else if (phase == NULL){
				wave[j] += amplitudes[i] * cos( (2.0 * M_PI * frequencies[i] * t) + 0.0);
			}
			else if (amplitudes == NULL){
				wave[j] += 1.0 * cos( (2.0 * M_PI * frequencies[i] * t) + phase[i] );
			}
			else{
				wave[j] += amplitudes[i] * cos( (2.0 * M_PI * frequencies[i] * t) + phase[i] );
			}
		}
	}

	return;
}

void generate_harmonic_wave_dbl( double fundamental,								  
								 double *amplitudes,
								 double* phase,
								 unsigned int num_harmonics,
								 unsigned int signal_length, 
								 double sampling_rate, 
								 double* wave )
{
	if (!isnormal((float)num_harmonics)){
		fprintf(stderr, "Error: Num Harmonics Is Invalid In Function -- generate_harmonic_wave!\n");
		return;
	}
	if (!isnormal((float)signal_length)){
		fprintf(stderr, "Error: Signal Length Is Invalid In Function -- generate_harmonic_wave!\n");
		return;
	}
	if (!isnormal((float)sampling_rate)){
		fprintf(stderr, "Error: Sampling Rate Is Invalid In Function -- generate_harmonic_wave!\n");
		return;
	}


	unsigned int i, j;
	double t, frequency;
	//===Create Wave===//
	for (j=0; j<signal_length; j++){
		wave[j] = 0; t = ((double)j)/sampling_rate;
		for (i=0; i<=num_harmonics; i++){
			frequency = ((double)(i+1))*fundamental;
			if (phase == NULL && amplitudes == NULL){
				wave[j] += 1.0* sin( (2.0 * M_PI * frequency * t) + 0.0);
			}
			else if (phase == NULL){
				wave[j] += amplitudes[i] * sin( (2.0 * M_PI * frequency * t) + 0.0);
			}
			else if (amplitudes == NULL){
				wave[j] += 1.0 * sin( (2.0 * M_PI * frequency * t) + phase[i] );
			}
			else{
				wave[j] += amplitudes[i] * sin( (2.0 * M_PI * frequency * t) + phase[i] );
			}
		}
	}

	return;
}


void generate_carrier_cmplx( double frequency, 
							 unsigned int signal_length, 
							 double sampling_rate, 
							 double complex* carrier )
{
	if (!isnormal((float)signal_length)){
		fprintf(stderr, "Error: Signal Length Is Invalid In Function -- generate_carrier_cmplx!\n");
		return;
	}
	if (!isnormal((float)sampling_rate)){
		fprintf(stderr, "Error: Sampling Rate Is Invalid In Function -- generate_carrier_cmplx!\n");
		return;
	}


	unsigned int j;
	double t;
	//===Create Wave===//
	for (j=0; j<signal_length; j++){
		t = ((double)j)/sampling_rate;
		carrier[j] = cexp(I * M_PI * t * frequency/(0.5*sampling_rate));
	}

	return;
}




void generate_impulse_dbl( unsigned int length,
					   	   double* impulse )
{
	unsigned int i;

	impulse[0] = 1.0;
	for (i=1; i<length; i++){
		impulse[i] = 0.0;
	}

	return;
}

void generate_unit_step_dbl( unsigned int start,
							 unsigned int length,
					   	   	 double* step )
{
	unsigned int i;

	for (i=0; i<length; i++){
		if (i>=start){
			step[i] = 1.0;
		}
		else{
			step[i] = 0.0;
		}
	}

	return;
}

void generate_unit_pulse_dbl( unsigned int start,
							  unsigned int end,
							  unsigned int length,
							  double* pulse )
{	
	unsigned int i;
	if (end > length){
		fprintf(stderr, "Error:: End Is Greater Than Length! In Function -- generate_unit_impulse_dbl!\n");
		return;
	}

	for (i=start; i<end-start; i++){
		pulse[i] = 1.0;
	}

	return;
}


//------------------------------------------------------------------------------------------------//
//===================================MEL FREQUENCY FUNCTIONS======================================//
//------------------------------------------------------------------------------------------------//

double hertz_to_mel( double hertz )
{	
	double mel;
	mel = 2595.0 * log10(1.0 + hertz/700.0);
	return mel;
}

double mel_to_hertz( double mel )
{	
	double hertz;
	hertz = 700.0*(pow(10.0,(mel/2595.0))-1.0);
	return hertz;
}

void make_mel_frequency_filterbanks( unsigned int num_frequency_nodes,
									 double sampling_rate,
						   			 unsigned int num_mel_channels,
						   			 double* mel_filterbanks )
{
	unsigned int i, center_index, idxbw;
	unsigned int mel_channel_bounds_indices[MAX_MEL_BANKS][2];
	double mel_increment;
	double freqs[MAX_FRAME_LENGTH];
	double ramp[MAX_FRAME_LENGTH];
	double mels[MAX_FRAME_LENGTH];
	double mel_centers_mel[MAX_MEL_BANKS];
	double mel_centers_hertz[MAX_MEL_BANKS];
	double mel_channel_bounds_hertz[MAX_MEL_BANKS][2];
	double* ptr;

	//===Sanity Checks===//
	if (num_frequency_nodes > MAX_FRAME_LENGTH){
		fprintf(stderr, "Error:: num_frequency_nodes Is Too Large! In Function -- make_mel_frequency_filterbanks!\n");
		quit();
	}
	if (num_mel_channels > MAX_MEL_BANKS){
		fprintf(stderr, "Error:: num_mel_channels Is Too Large! In Function -- make_mel_frequency_filterbanks!\n");
		quit();
	}
	if (mel_filterbanks == NULL){
		fprintf(stderr, "Error:: mel_filterbanks Is NULL! In Function -- make_mel_frequency_filterbanks!\n");
		quit();
	}

	//===Get Hertz Frequencies===//
	get_nonnegative_frequencies(num_frequency_nodes, sampling_rate, freqs);

	if (array_has_nans_dbl(freqs, num_frequency_nodes)){
		fprintf(stderr, "Error:: Freqs Has NaNs! In Function -- make_mel_frequency_filterbanks!\n");
		quit(); 
	}


	//===Convert To Mel===//
	for (i=0; i<num_frequency_nodes; i++) mels[i] = hertz_to_mel(freqs[i]);

	//===Get Mel Bank Centers===//
	mel_increment = find_maximum_dbl(mels, num_frequency_nodes)/(num_mel_channels + 1);
	for (i=0; i<num_mel_channels; i++){
		mel_centers_mel[i] = (i+1)*mel_increment;
	}

	//===Convert Mel Centers To Hertz===//
	for (i=0; i<num_mel_channels; i++) mel_centers_hertz[i] = mel_to_hertz(mel_centers_mel[i]);

	if (array_has_nans_dbl(mel_centers_hertz, num_mel_channels)){
		fprintf(stderr, "Error:: Mel Centers Hertz Has NaNs! In Function -- make_mel_frequency_filterbanks!\n");
		quit(); 
	}



	//===Get Hertz Frequency Bounds===//
	mel_channel_bounds_hertz[0][0] = 0;
	mel_channel_bounds_hertz[num_mel_channels-1][1] = sampling_rate/2.0;
	for (i=1; i<num_mel_channels; i++){
		mel_channel_bounds_hertz[i-1][1] = mel_centers_hertz[i];
		mel_channel_bounds_hertz[i][0] = mel_centers_hertz[i-1];
	}

	//===Get Hertz Frequency Indices===//
	for (i=0; i<num_mel_channels; i++){
		mel_channel_bounds_indices[i][0] = find_closest_index_dbl(freqs, num_frequency_nodes, mel_channel_bounds_hertz[i][0]);
		mel_channel_bounds_indices[i][1] = find_closest_index_dbl(freqs, num_frequency_nodes, mel_channel_bounds_hertz[i][1]);
	}	

	//===Create The Triangles===//
	for (i=0; i<num_mel_channels; i++){
		
		//===Initialize===//
		ptr = mel_filterbanks + i*num_frequency_nodes;
		initialize_array_dbl(ptr, num_frequency_nodes);

		//===Get Center Index===//
		center_index = find_closest_index_dbl(freqs, num_frequency_nodes, mel_centers_hertz[i]);

		//===Create 1st Ramp===//
		idxbw = center_index - mel_channel_bounds_indices[i][0] + 1;
		initialize_array_dbl(ramp, MAX_FRAME_LENGTH);
		if (idxbw > 1){
			initialize_linspace_dbl(ramp, idxbw, 0, 1);
		}

		//===Copy Over===//
		copy_array_dbl(ramp, idxbw, ptr + mel_channel_bounds_indices[i][0]);

		//===Create 2nd Ramp===//
		idxbw = mel_channel_bounds_indices[i][1] - center_index + 1;
		initialize_array_dbl(ramp, MAX_FRAME_LENGTH);
		if (idxbw > 1){
			initialize_linspace_dbl(ramp, idxbw, 1, 0);
		}

		//===Copy Over===//
		copy_array_dbl(ramp, idxbw, ptr + center_index);		
	}

	return;
}

void compute_mel_frequency_cepstral_coefficients_dbl( double* data,
													  unsigned int length,
													  double sampling_rate,
													  unsigned int num_mel_channels,
													  double* mfccs )
{
	unsigned int i, num_nonnegative_frequency_nodes;
	double L;
	double pre_emph_coeffs[2];
	double copy[MAX_FRAME_LENGTH];
	double mel_filterbanks[MAX_MEL_BANKS*MAX_FRAME_LENGTH];
	double energies[MAX_MEL_BANKS];
	double psd[MAX_FRAME_LENGTH];
	filter_t premph;

	//===Sanity Checks===//
	if (length > MAX_FRAME_LENGTH){
		fprintf(stderr, "Error:: Length Is Too Large! In Function -- compute_mel_frequency_cepstral_coefficients_dbl!\n");
		quit();
	}

	//===Create Mel Banks===//
	num_nonnegative_frequency_nodes = length/2;
	initialize_array_dbl(mel_filterbanks, MAX_MEL_BANKS*MAX_FRAME_LENGTH);
	make_mel_frequency_filterbanks(length/2, sampling_rate, num_mel_channels, mel_filterbanks);

	//===Compute Preemphasis===//
	pre_emph_coeffs[0] = 1; pre_emph_coeffs[1] = -0.97;
	initialize_filter(&(premph), pre_emph_coeffs, 2, NULL, 0); //raises frequency response with frequency (differentiator)
	run_filter(&premph, data, length, copy);
	
	//===Take Full Power Spectrum Of Data===//
	compute_power_spectrum_dbl(copy, length, NULL, psd);

	//===Sum Energy In Mel Filtered Power Spectrum===//
	for (i=0; i<num_mel_channels; i++){

		//===Inner Product With Only 'Num Nonnegative Indices' ie Positive Part Of Full PSD====//
		energies[i] = compute_inner_product_dbl(psd, mel_filterbanks + i*num_nonnegative_frequency_nodes, 
													 num_nonnegative_frequency_nodes);
		//===Take Logarithm Of Summed Energies===//
		energies[i] = log(energies[i] + EPS);

	}

	//===Take DCT Of Summed Energies===//
	dct_dbl(energies, num_mel_channels, mfccs);

	//===Lifter===//
	L = 22.0;
	for (i=0; i<num_mel_channels; i++){
		if (isnan(mfccs[i])){
			fprintf(stderr, "Error:: MFCC Is NaN! In Function -- compute_mel_frequency_cepstral_coefficients_dbl!\n");
			quit(); 
		}
		mfccs[i] *= (1.0 + (L/2.0) * sin(M_PI * ((double)i)/L));
	}

	return;
}

//------------------------------------------------------------------------------------------------//
//====================================UNIT TESTING FUNCTIONS======================================//
//------------------------------------------------------------------------------------------------//

void test_wave()
{
	unsigned int length;
	double freq, amp, phase, Fs;
	double wave[MAX_FRAME_LENGTH];
	FILE* fout;

	//===Make Wave===//
	freq = 500.0; Fs = 32000.0; amp = 1.0; phase = 0.0; length = 1024;
	generate_wave_dbl(&freq, &amp, &phase, 1, length, Fs, wave);

	//===Print Wave===//
	fout = fopen("wave.dat", "w");
	print_vector_dbl(wave, length, fout);
	fclose(fout);

	return;
}

void test_phase_unwrap()
{
	double P[100];
	double Q[100];

	//===Unwrap Phase===//
	P[0] = 0; P[1] = 7.0686; P[2] = 1.5708; P[3] = 2.3562; P[4] = 0.1963; P[5] = 0.9817; 
	P[6] = 1.7671; P[7] = 2.5525; P[8] = 6.6759; P[9] = 1.1781; P[10] = 1.9635; P[11] = 2.7489;
	P[12] = 0.5890; P[13] = 1.3744; P[14] = 2.1598; P[15] = 2.9452;
	unwrap_phase_angles(P, 16, M_PI);

	//===Compare To Truth===//
	Q[0] = 0; Q[1] = 0.78541469; Q[2] = 1.5708; Q[3] = 2.3562; Q[4] = 0.1963;
	Q[5] = 0.9817; Q[6] = 1.7671; Q[7] = 2.5525; Q[8] = 0.39271469; Q[9] = 1.1781;
	Q[10] = 1.9635; Q[11] = 2.7489; Q[12] = 0.589; Q[13] = 1.3744; Q[14] = 2.1598; 
	Q[15] = 2.9452;

	//===Print===//
	fprintf(stdout, "Computed Unwrapped Phases:\n");
	print_vector_dbl(P, 16, stdout);
	newline();
	fprintf(stdout, "True Unwrapped Phases:\n");
	print_vector_dbl(Q, 16, stdout);
	newline();

	return;
}

void test_dct()
{

	double array[MAX_FRAME_LENGTH];
	double true_dct[MAX_FRAME_LENGTH];
	double dct[MAX_FRAME_LENGTH];

	//===Make Array===//
	initialize_linspace_dbl(array, 10, 1, 10);
	
	//===Take DCT===//
	dct_dbl(array, 10, dct);

	//===Print===//
	fprintf(stdout, "Array: \n");
	print_vector_dbl(array, 10, stdout);
	newline();
	fprintf(stdout, "DCT: \n");
	print_vector_dbl(dct, 10, stdout);

	//===Print Python Answer===//
	true_dct[0] = 17.39252713;
	true_dct[1] = -9.02485113;
	true_dct[2] = 0;
	true_dct[3] = -0.9666569;
	true_dct[4] = 0;
	true_dct[5] = -0.31622777;
	true_dct[6] = 0;
	true_dct[7] = -0.12787039;
	true_dct[8] = 0;
	true_dct[9] = -0.0358573;
	newline();
	fprintf(stdout, "True DCT: \n");
	print_vector_dbl(true_dct, 10, stdout);

	return;
}

void test_mel()
{
	unsigned int num_mel_channels, num_non_nodes;
	double Fs;
	double freqs[MAX_FRAME_LENGTH];
	double* mel_filterbanks;
	FILE* fout;

	//===Mallocs===//
	mel_filterbanks = malloc(MAX_MEL_BANKS*MAX_FRAME_LENGTH*sizeof(double));

	//===Make Banks===//
	num_non_nodes = 256; num_mel_channels = 10;
	Fs = 16000;
	initialize_array_dbl(mel_filterbanks, MAX_MEL_BANKS*MAX_FRAME_LENGTH);
	make_mel_frequency_filterbanks(num_non_nodes, Fs, num_mel_channels, mel_filterbanks);

	//===Print===//
	fout = fopen("mel_banks.dat", "w");
	print_matrix_dbl(mel_filterbanks, num_mel_channels, num_non_nodes, fout);
	fclose(fout);
	fout = fopen("mel_banks_hertz.dat", "w");
	get_nonnegative_frequencies(num_non_nodes, Fs, freqs);
	print_vector_dbl(freqs, num_non_nodes, fout);
	fclose(fout);
	
	//===Clean Up===//
	free(mel_filterbanks);

	return;
}

void test_interval_chirp_transform()
{

	unsigned int signal_length, num_frequency_nodes, num_carriers;
	double start_frequency, end_frequency, omega_0, del_k, delta_w, sampling_rate;
	double freqs[10];
	double signal[MAX_SIGNAL_LENGTH];
	double complex chirp_fft[MAX_SIGNAL_LENGTH];
	FILE* fout;

	//===Create Sequence===//
	signal_length = pow(2,16); 
	num_frequency_nodes = signal_length;
	num_carriers = 4;
	freqs[0] = 500; freqs[1] = 600; freqs[2] = 700; freqs[3] = 800;
	sampling_rate = 2000.0;
	generate_wave_dbl(freqs, NULL, NULL, num_carriers, signal_length, sampling_rate, signal);

	//===Compute Chirp===//
	fprintf(stdout, "Computing Chirp FFT:\n");
	start_frequency = freqs[0]-100.0; end_frequency = freqs[num_carriers-1]+100.0;
	omega_0 = (2.0*M_PI*start_frequency)/sampling_rate; 
	del_k = (end_frequency-start_frequency)/sampling_rate;
	delta_w = (2.0*M_PI*del_k)/((double)(num_frequency_nodes-1));
	initialize_array_cmplx(chirp_fft, num_frequency_nodes);
	chirp_fft_dbl(signal, signal_length, omega_0, delta_w, num_frequency_nodes, chirp_fft);
	
	//===Print===//
	fout = fopen("wave_chirp_spectrum.dat", "w");
	print_vector_cmplx(chirp_fft, num_frequency_nodes, fout);
	fclose(fout);

	return;
}


