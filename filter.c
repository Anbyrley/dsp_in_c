#include "filter.h"

//------------------------------------------------------------------------------------------------//
//====================================ANALOG DESIGN FUNCTIONS=====================================//
//------------------------------------------------------------------------------------------------//

double prewarp_frequency_dbl( double digital_frequency )
{	

	double analog_frequency;
	analog_frequency = 2.0 * tan(digital_frequency/2.0);	
	return analog_frequency;
}

unsigned int calculate_butterworth_order( double passband_edge,
										  double stopband_edge,
										  double passband_gain,
										  double stopband_gain )
{
	unsigned int order;
	double vpnorm, vsnorm;

	//===Create Normalized Gains===//
	vpnorm = passband_gain/sqrt((1.0-pow(passband_gain,2.0)));
	vsnorm = stopband_gain/sqrt((1.0-pow(stopband_gain,2.0)));

	//===Calculate Required Order===//
	order = (unsigned int)ceil(log(vpnorm/vsnorm)/log(stopband_edge/passband_edge));

	return order;
}

void butterworth_design( unsigned int order,
						 double passband_edge,
						 double* numerator,
						 double* denominator )
{
	unsigned int v;
	double complex poles[MAX_POLYNOMIAL_ORDER];
	double complex a[MAX_POLYNOMIAL_ORDER+1];
	
	//===Calculate Poles===//
	for (v=0; v<order; v++){
		poles[v] = -sin(M_PI * (double)(2*v+1)/(double)(2*(order)));
		poles[v] += I*cos(M_PI * (double)(2*v+1)/(double)(2*(order)));
		if (fabs(cimag(poles[v])) < 1.0*EPS){
			poles[v] = creal(poles[v]) + I*0.0;
		}
	}

	//===Form Transfer Function Denominator===//
	convert_roots_to_coefficients(poles, order, a);
	split_array_cmplx(a, order+1, denominator, NULL);

	//===Scale By Passband Edge===//
	for (v=0; v<order+1; v++){
		denominator[v] /= pow(passband_edge,order-v);
	}

	//===Normalize Transfer Function===//
	numerator[0] = 1.0/denominator[0];
	make_array_monic_dbl(denominator, order+1);

	return;
}	

unsigned int calculate_chebyschev_order( double passband_ripple_dB,
										 double passband_edge,
										 double stopband_edge )
{
	unsigned int order;
	double epsilon, A;

	epsilon = sqrt(pow(10.0, passband_ripple_dB/10.0) - 1.0);
	A = pow(10.0, passband_ripple_dB);
	order = (unsigned int)ceil(acosh(sqrt((pow(A,2.0) - 1.0)/pow(epsilon, 2.0)))/acosh(stopband_edge/passband_edge)); 

	return order;
}

void chebyshev_design( unsigned int order,
					   double passband_ripple_dB,
					   double passband_edge,
					   double* numerator,
					   double* denominator )
{
	unsigned int k;
	double epsilon;
	double complex b;
	double complex poles[MAX_POLYNOMIAL_ORDER];
	double complex a[MAX_POLYNOMIAL_ORDER];  

	//===Calculate Poles===//
	epsilon = sqrt(pow(10.0, passband_ripple_dB/10.0) - 1.0);
	for (k=1; k<=order; k++){
		poles[k-1] = -1.0*sin((double)(2*k-1) * (M_PI/(2.0*(double)(order)))) * sinh((1.0/(double)(order))*asinh(1.0/epsilon));
		poles[k-1] += I*cos((double)(2*k-1) * (M_PI/(2.0*(double)(order)))) * cosh((1.0/(double)(order))*asinh(1.0/epsilon));
	}

	//===Form Transfer Function Numerator===//
	b = pow(-1.0, (double)(order));
	for (k=0; k<order; k++){
		b *= poles[k];
	}
	numerator[0] = b;	

	//===Form Transfer Function Denominator===//
	convert_roots_to_coefficients(poles, order, a);
	split_array_cmplx(a, order+1, denominator, NULL);

	//===Scale By Passband Edge===//
	for (k=0; k<order+1; k++){
		denominator[k] /= pow(passband_edge,order-k);
	}
	
	//===Scale By Filter Order===//
	if (order % 2 == 0){
		numerator[0] *= 1.0/sqrt(1.0 + pow(epsilon,2.0));
	}

	//===Normalize Transfer Function===//
	numerator[0] *= 1.0/denominator[0];
	make_array_monic_dbl(denominator, order+1);

	return;
}

void bilinear_transform( filter_t* self )
{
	int n, maxCoef, m, j, num_zeros, num_poles;
	double hC, bigT, H0;
	double complex beta, gamma, delta, eta;
	double complex work, cTwo;
	double a_z[MAX_POLYNOMIAL_ORDER], b_z[MAX_POLYNOMIAL_ORDER];
	double complex zeros[MAX_POLYNOMIAL_ORDER], poles[MAX_POLYNOMIAL_ORDER];
	double complex mu[MAX_POLYNOMIAL_ORDER];

	//===Start Up===//
	initialize_array_cmplx(mu, MAX_POLYNOMIAL_ORDER);
	initialize_array_dbl(a_z, MAX_POLYNOMIAL_ORDER);
	initialize_array_dbl(b_z, MAX_POLYNOMIAL_ORDER);

	//===Find Poles And Zeros===//
	num_zeros = self->numerator_length - 1;
	if (self->numerator_length > 1){
		initialize_array_cmplx(zeros, MAX_POLYNOMIAL_ORDER);
		find_roots_dbl(self->numerator, self->numerator_length, zeros);
	}
	num_poles = self->denominator_length - 1;
	if (self->denominator_length > 1){
		initialize_array_cmplx(poles, MAX_POLYNOMIAL_ORDER);
		find_roots_dbl(self->denominator, self->denominator_length, poles);
	}

	//===Compute Constant Gain Factor===//
	hC = 1.0; bigT = 1.0; H0 = self->numerator[0];
	work = 1.0 + I*0.0; cTwo = 2.0 + I*0.0;
	for (n=0; n<num_poles; n++){
		work = work * (cTwo - (bigT+I*0.0)*poles[n]);
		hC = hC * bigT;
	}
	hC = H0 * hC/creal(work);

	//===Run Numerator Recursion===//
	mu[0] = 1.0 + I*0.0;
	maxCoef = 0;
	self->numerator_length = num_poles+1;
	for (m=0; m<num_poles-num_zeros; m++){
		maxCoef++;
		for (j=maxCoef; j>=1; j--){
			mu[j] = mu[j] + mu[j-1];
		}
	}
	for (m=0; m<num_zeros; m++){
		maxCoef++;
		beta = (-2.0/bigT + I*0.0) - poles[n];
		for (j=maxCoef; j>=1; j--){
			mu[j] = mu[j] + beta*mu[j-1];
		}
	}
	initialize_array_dbl(self->numerator, self->numerator_length);
	for (j=0; j<=num_poles; j++){
		self->numerator[j] = creal(hC * mu[j]);
	}

	//===Run Denominator Recursion===//
	initialize_array_cmplx(mu, MAX_POLYNOMIAL_ORDER);
	self->denominator_length = num_poles+1;
	mu[0] = 1.0 + I*0.0;
	for (n=1; n<=num_poles; n++){
		gamma = (2.0/bigT + I*0.0) - poles[n-1];
		delta = (-2.0/bigT + I*0.0) - poles[n-1];
		eta = delta/gamma;
		for (j=n; j>=1; j--){
			mu[j] = mu[j] + eta * mu[j-1];
		}
	}
	initialize_array_dbl(self->denominator, self->denominator_length);
	for (j=0; j<=num_poles; j++){
		self->denominator[j] = creal(1.0 * mu[j]);
	}

	//===Normalize Transfer Function===//
	self->numerator[0] = self->numerator[0] * 1.0/self->denominator[0];
	make_array_monic_dbl(self->denominator, self->denominator_length);

	return;
}

void initialize_lowpass_ideal( double fc,
							   double Fs,
							   int order,
							   double* ideal )
{

	unsigned int mm, length;
	double alpha;
	double cutoff;
	double m[MAX_POLYNOMIAL_ORDER+1];

	//===Make Locals===//
	cutoff = fc/(0.5*Fs);
	alpha = 0.5 * (order - 1);
	
	//===Make Arrays===//
	initialize_range_dbl(m, 0, order, 1, &length);
	add_array_constant_dbl(m, length, -alpha);
	initialize_array_dbl(ideal, length);
	for (mm = 0; mm < length; mm++){
		if (cutoff * m[mm] == 0){
			ideal[mm] += cutoff;
		}
		else{
            ideal[mm] += cutoff * sin(M_PI * cutoff * m[mm])/(M_PI * cutoff * m[mm]); 
		}
	}

	return;
}

void initialize_bandpass_ideal( double f_low,
								double f_high,
								double Fs,
								int order,
								double* ideal )
{

	int v;
	double T;
	double w_low, w_high;
	double vt;

	w_low = 2.0 * M_PI * f_low;
	w_high = 2.0 * M_PI * f_high;
	T = 1.0/Fs;
	for (v=-order/2; v<order/2; v++){
		vt = ((double)v)/(1.0);
		if (v==0){
			ideal[v+order/2] = w_high - w_low;
		}
		else{
			ideal[v+order/2] = (sin(vt*w_high*T)/(vt*T)) - (sin(vt*w_low*T)/(vt*T));
		}
	}
	normalize_max_dbl(ideal, order);

	return;
}	

void initialize_highpass_ideal( double fc,
								double Fs,
								int order,
								double* ideal )
{

	unsigned int mm, length;
	double alpha;
	double cutoff;
	double m[MAX_POLYNOMIAL_ORDER+1];

	//===Make Locals===//
	cutoff = fc/(0.5*Fs);
	alpha = 0.5 * (order - 1);
	
	//===Make Arrays===//
	initialize_range_dbl(m, 0, order, 1, &length);
	add_array_constant_dbl(m, length, -alpha);
	initialize_array_dbl(ideal, length);
	for (mm = 0; mm < length; mm++){
		if (cutoff * m[mm] == 0){
			ideal[mm] += cutoff;
		}
		else{
            ideal[mm] += cutoff * sin(M_PI * cutoff * m[mm])/(M_PI * cutoff * m[mm]); 
		}
	}

	//===Spectral Inversion===//
	normalize_max_dbl(ideal, length);
	for (mm=0; mm<length; mm++){
		ideal[mm] = pow((-1),mm) * ideal[mm];
	}

	return;
}	

void initialize_bandstop_ideal( double f_low,
								double f_high,
								double Fs,
								int order,
								double* ideal )
{

	int v;
	double T;
	double w_low, w_high;
	double vt;

	w_low = 2.0 * M_PI * f_low;
	w_high = 2.0 * M_PI * f_high;
	T = 1.0/Fs;
	for (v=-order/2; v<order/2; v++){
		vt = ((double)v)/(1.0);
		if (v==0){
			ideal[v+order/2] = 1.0 - 2.0*(w_high - w_low);
		}
		else{
			ideal[v+order/2] = (sin(vt*w_low*T)/(vt*T)) - (sin(vt*w_high*T)/(vt*T));
		}
	}
	normalize_max_dbl(ideal, order);

	return;
}	

//------------------------------------------------------------------------------------------------//
//====================================WINDOW DESIGN FUNCTIONS=====================================//
//------------------------------------------------------------------------------------------------//

void initialize_cosine_window( unsigned int length,
							   double* window )
{

	unsigned int i;
	for (i=0; i<length; i++){
		window[i] = sin(M_PI/(double)(2*length)) * cos((M_PI * (i - (length/2)))/(double)length);
	}
	normalize_max_dbl(window, length);
	return;
}

void initialize_hanning_window( unsigned int length, 
								double* window )
{
	unsigned int i;
	for (i=0; i<length; i++){
		window[i] = 0.5 + 0.5*cos(2.0*M_PI*i/(length));
	}
	normalize_max_dbl(window, length);
	right_circular_shift_array_dbl(window, length, length/2);

	return;
}

void initialize_hamming_window( unsigned int length, 
								double* window )
{
	unsigned int i;
	for (i=0; i<length; i++){
		window[i] = 0.54 + 0.46*cos(2.0*M_PI*i/(length));
	}
	normalize_max_dbl(window, length);
	right_circular_shift_array_dbl(window, length, length/2);

	return;
}

void initialize_blackman_window( unsigned int length,
								 double* window )
{
	unsigned int i;
	for (i=0; i<length; i++){
		window[i] = 0.42 + 0.5*cos(2.0*M_PI*i/(length)) + 0.08*cos(4*M_PI*i/(length));
	}
	normalize_max_dbl(window, length);
	right_circular_shift_array_dbl(window, length, length/2);
	
	return;
}

void initialize_rectangular_window( unsigned int length,
									double* window )
{
	unsigned int i;
	for (i=0; i<length; i++){
		window[i] = 1.0;
	}

	return;
}

void initialize_chebyshev_window( double sidelobe_atten,
								  unsigned int length,
								  double *window )
{
    int nn, i, N;
    double M, n, sum;
    double tg, x0;

	N = length;
	sum = 0;
	tg = pow(10,sidelobe_atten/20);  // 1/r term [2], 10^gamma [2] //
    x0 = cosh((1.0/(N-1))*acosh(tg));
    M = (N-1)/2;
    if(N%2==0){
		 M = M + 0.5; // handle even length windows //
	}

	//===Run Sum===//
    for(nn=0; nn<(N/2+1); nn++){
        n = nn-M;
        sum = 0;
        for(i=1; i<=M; i++){
            sum += compute_chebyshev_polynomial(N-1,x0*cos(PI*i/N))*cos(2.0*n*PI*i/N);
        }
        window[nn] = tg + 2*sum;
        window[N-nn-1] = window[nn];
    }
	normalize_max_dbl(window, N);

    return;
}

double get_kaiser_beta( double stopband_attenuation_dB )
{
	if (stopband_attenuation_dB >= 50.0){
		  return (0.1102 * (stopband_attenuation_dB - 8.7));
	}
	else{
    	return (0.5842 * pow(stopband_attenuation_dB - 21.0, 0.4) + 0.07886 * (stopband_attenuation_dB - 21.0));
	}
	return 0;
}

unsigned int get_kaiser_filter_length( double deltaw,
									   double stopband_attenuation_dB,
									   unsigned int num_subbands )
{
	unsigned int N;
	N = (int)(((stopband_attenuation_dB - 7.95) / 14.36) * ((double)num_subbands / (deltaw)));
	return N;
}

void initialize_kaiser_window( double beta,
							   unsigned int length,
							   double* window )
{
	unsigned int i;
	double I0_input;

	for (i=0; i<length; i++){
		I0_input = sqrt(1.0 - pow(((double)(2*i)/(double)(length-1))-1.0,2.0));
		I0_input *= beta;
		window[i] = compute_bessel_I0(I0_input)/compute_bessel_I0(beta);
	}
	normalize_max_dbl(window, length);

	return;
}	

void prolate_window_design( double ws,
							unsigned int N,
							double* window )
{

	unsigned int n, m, i, num_points, M, min_index;
	double x;
	double w[MAX_NUM_QMC_POINTS], w1[MAX_NUM_QMC_POINTS];
	double X[MAX_NUM_QMC_POINTS];
	double complex eigenvalues[MAX_POLYNOMIAL_ORDER];
	double complex *P, *right_eigenvectors, *left_eigenvectors, *min_eigenvector;

	//===Allocation===//
	N += 2;
	initialize_array_dbl(window, N);
	right_eigenvectors = malloc(MAX_POLYNOMIAL_ORDER*MAX_POLYNOMIAL_ORDER*sizeof(double complex));
	left_eigenvectors = malloc(MAX_POLYNOMIAL_ORDER*MAX_POLYNOMIAL_ORDER*sizeof(double complex));

	//===Sample [0,1]===//
	num_points = MAX_NUM_QMC_POINTS;
	make_korobov_sequence(123421, num_points, w);
	for (n=0; n<num_points; n++){
		w1[n] = (M_PI - ws)*w[n] + ws;
	}

	//===Force An Odd Filter Length===//
	if (N%2==0){
		N += 1;
	}
	M = (N-1)/2;
	
	//===Form Error Integral Matrix===//
	P = malloc(M*M*sizeof(double complex));
	for (n=0; n<M; n++){
		for (m=0; m<n+1; m++){

			//===Sample Over StopBand Interval===//
			initialize_array_dbl(X, num_points);
			for (i=0; i<num_points; i++){
				X[i] += cos(n*w1[i])*cos(m*w1[i]);
			}

			//===Get Mean As Approximation Of Integral===//
			x = compute_mean_dbl(X, num_points);
			x *= 1.0/(1.0*M_PI);

			//===Set Matrix Element===//	
			P[m + n*M] = x;			
			P[n + m*M] = P[m + n*M];			
		}
	}

	//===Compute EigenSystem===//
	compute_eigensystem_cmplx(P, M, 1, eigenvalues, right_eigenvectors, left_eigenvectors);

	//===Get Minimum Eigenvector===//	
	min_index = find_minimum_magnitude_index_cmplx(eigenvalues, M);
	min_eigenvector = right_eigenvectors + min_index*M;

	//===Form Filter===//
	window[M] = min_eigenvector[0];
	for (i=1; i<M; i++){
		window[M-i] = 0.5*creal(min_eigenvector[i]);
		window[M+i] = 0.5*creal(min_eigenvector[i]);
	}
	left_shift_array_dbl(window, M+i, 1);
	normalize_sum_dbl(window,N);
	normalize_max_dbl(window,N);
	
	//===Clean Up===//
	free(P);
	free(left_eigenvectors);
	free(right_eigenvectors);
	
	return;
}

void initialize_fyrley_window( unsigned int length,
							   double* window )
{
	unsigned int i;
	double* x;

	//===Init===//
	x = malloc(length*sizeof(double));

	//===Initialize Domain===//
	initialize_linspace_dbl(x, length, 10e-6, 1-10e-6);

	//===Compute Window===//	
	for (i=0; i<length; i++){
		//e^-1/x*(1-x)^2
		window[i] = exp(-1.0/(x[i]*x[i]*((1.0-x[i])*(1.0-x[i]))));
	}

	//normalize_power_dbl(window, length);
	normalize_max_dbl(window, length);

	//===Clean Up===//
	free(x);
	
	return;
}


//------------------------------------------------------------------------------------------------//
//====================================FILTER DESIGN FUNCTIONS=====================================//
//------------------------------------------------------------------------------------------------//

void window_filter_design( filter_type type,
						   double fc1,
						   double fc2,
						   double Fs,
						   int length,
						   double* window,
						   double* filter)
{
	int i;
	double lowpass[MAX_POLYNOMIAL_ORDER+1];
	double highpass[MAX_POLYNOMIAL_ORDER+1];

	if (type == LOWPASS){
		initialize_lowpass_ideal(fc2, Fs, length, filter);
		hadamard_product_dbl(window, filter, length, filter);
	}
	else if (type == HIGHPASS){
		if (length % 2 == 0){
			fprintf(stderr, "Error:: Highpass Filters Cannot Have Even Length! In Function -- window_filter_design!\n");
			quit();
		}
		initialize_lowpass_ideal(Fs/2.0 - fc1, Fs, length, filter);
		hadamard_product_dbl(window, filter, length, filter);
		for (i=0; i<length; i++){
			filter[i] = pow((-1),i) * filter[i];
		}
	}
	else{
		if (length % 2 == 0){
			fprintf(stderr, "Error:: Bandpass Filters Cannot Have Even Length! In Function -- window_filter_design!\n");
			quit();
		}

		/*
		//===Make Lowpass===//
		initialize_lowpass_ideal(fc2, Fs, length, lowpass);
		hadamard_product_dbl(window, lowpass, length, lowpass);
		
		//===Make Highpass===//
		initialize_lowpass_ideal(Fs/2.0 - fc1, Fs, length, highpass);
		hadamard_product_dbl(window, highpass, length, highpass);
		for (i=0; i<length; i++){
			highpass[i] = pow((-1),i) * highpass[i];
		}
	
		//===Make Bandpass===//
		add_vectors_dbl(lowpass, highpass, length, filter);
		*/
		
		
		//===Make Lowpass===//
		initialize_lowpass_ideal(fc1, Fs, length, lowpass);
		hadamard_product_dbl(window, lowpass, length, lowpass);
		
		//===Make Highpass===//
		initialize_lowpass_ideal(Fs/2.0 - fc2, Fs, length, highpass);
		hadamard_product_dbl(window, highpass, length, highpass);
		for (i=0; i<length; i++){
			highpass[i] = pow((-1),i) * highpass[i];
		}
	
		//===Make Bandpass===//
		add_vectors_dbl(lowpass, highpass, length, filter);
		for (i=0; i<length; i++){
			filter[i] = pow((-1),i) * filter[i];
		}
			
	}
	return;
}

void notch_filter_design( double* notch_frequencies,
						  unsigned int num_notches,
						  double sampling_rate,
						  double notch_width,
						  unsigned int* filter_length,
						  double* numerator,
						  double* denominator )
{
	unsigned int i;
	double R;
	double complex roots[MAX_POLYNOMIAL_ORDER];
	double complex coeffs[MAX_POLYNOMIAL_ORDER+1];

	//===Make Roots===//
	for (i=0; i<num_notches; i++){
		roots[i] = cexp(I * 2.0 * M_PI * notch_frequencies[i] / sampling_rate);
	}
	for (i=num_notches; i<2*num_notches; i++){
		roots[i] = cexp(-1.0 * I * 2.0 * M_PI * notch_frequencies[i-num_notches] / sampling_rate);
	}

	//===Make Numerator Coeff===//
	num_notches *= 2;
	convert_roots_to_coefficients(roots, num_notches, coeffs);
	make_array_monic_cmplx(coeffs, num_notches+1);
	split_array_cmplx(coeffs, num_notches+1, numerator, NULL);

	//===Make Denominator Coeffs===//
	R = 1.0 - (2.0*M_PI*(notch_width)/sampling_rate)/2.0;
	for (i=0; i<num_notches+1; i++){
		denominator[i] = pow(R, num_notches-i) * numerator[i];
	}

	//===Reverse===//
	reverse_array_dbl(numerator, num_notches+1);
	reverse_array_dbl(denominator, num_notches+1);
	*filter_length = num_notches + 1;

	return;
}
						  

void eigenfilter_design( double wp,
						 double ws,
						 unsigned int N,
						 double alpha,
						 double* filter )
{
	unsigned int n, m, i, num_points, M, min_index;
	double x;
	double w[MAX_NUM_QMC_POINTS], w1[MAX_NUM_QMC_POINTS], w2[MAX_NUM_QMC_POINTS];
	double X[MAX_NUM_QMC_POINTS];
	double complex eigenvalues[MAX_POLYNOMIAL_ORDER];
	double complex *P, *right_eigenvectors, *left_eigenvectors, *min_eigenvector;

	//NOTE: For optimal results, wp - ws should be about 0.15, and no greater than 0.15.

	//===Allocation===//
	N += 2;
	initialize_array_dbl(filter, N);
	right_eigenvectors = malloc(MAX_POLYNOMIAL_ORDER*MAX_POLYNOMIAL_ORDER*sizeof(double complex));
	left_eigenvectors = malloc(MAX_POLYNOMIAL_ORDER*MAX_POLYNOMIAL_ORDER*sizeof(double complex));

	//===Sample [0,1]===//
	num_points = MAX_NUM_QMC_POINTS;
	make_korobov_sequence(123421, num_points, w);
	for (n=0; n<num_points; n++){
		w1[n] = (M_PI - ws)*w[n] + ws;
		w2[n] = wp*w[n];
	}

	//===Force An Odd Filter Length===//
	if (N%2==0){
		N += 1;
	}
	M = (N-1)/2;

	//===Form Error Integral Matrix===//
	P = malloc(M*M*sizeof(double complex));
	for (n=0; n<M; n++){
		for (m=0; m<n+1; m++){
			//===Sample Over StopBand Interval===//
			initialize_array_dbl(X, num_points);
			for (i=0; i<num_points; i++){
				X[i] += cos(n*w1[i])*cos(m*w1[i]);
			}

			//===Get Mean As Approximation Of Integral===//
			x = compute_mean_dbl(X, num_points);
			x *= alpha/(2.0*M_PI);

			//===Set Matrix Element===//	
			P[m + n*M] = x;			
			P[n + m*M] = P[m + n*M];			
			
			//===Sample Over PassBand Interval===//
			initialize_array_dbl(X, num_points);
			for (i=0; i<num_points; i++){
				X[i] += (1.0-cos(n*w2[i]))*(1.0-cos(m*w2[i]));
			}

			//===Get Mean As Approximation Of Integral===//
			x = compute_mean_dbl(X, num_points);
			x *= (1.0-alpha)/(2.0*M_PI);

			//===Set Matrix Element===//	
			P[m + n*M] += x;			
			P[n + m*M] = P[m + n*M];			
						
		}
	}

	//===Compute EigenSystem===//
	compute_eigensystem_cmplx(P, M, 1, eigenvalues, right_eigenvectors, left_eigenvectors);

	//===Get Minimum Eigenvector===//	
	min_index = find_minimum_magnitude_index_cmplx(eigenvalues, M);
	min_eigenvector = right_eigenvectors + min_index*M;

	//===Form Filter===//
	filter[M] = min_eigenvector[0];
	for (i=1; i<M; i++){
		filter[M-i] = 0.5*creal(min_eigenvector[i]);
		filter[M+i] = 0.5*creal(min_eigenvector[i]);
	}
	left_shift_array_dbl(filter, M+i, 1);
	normalize_sum_dbl(filter,N);
	
	//===Clean Up===//
	free(P);
	free(left_eigenvectors);
	free(right_eigenvectors);
	
	return;
}


void minimax_filter_design( double wp,
							double ws,
							double stopband_attn,	//dB
							double passband_ripple, //dB
							unsigned int* filter_length,
							double* filter )
{
	unsigned int i, j, k, grid_size, count, passband_refset_size, stopband_refset_size, length;
	unsigned int zero_pad_length, num_maxima, temp_uint, reference_length;
	double wc, passband_error, stopband_error, tolerance, total_error;
	double reference_grid[MAX_GRID_SIZE], a_coeffs[(MAX_POLYNOMIAL_ORDER+1)];
	double desired_response[MAX_POLYNOMIAL_ORDER+1], maxima[MAX_GRID_SIZE];
	double passband_bounds[2], stopband_bounds[2], errors[4];
	double *A, *H_real, *freq_response;
	double complex* H;

	//===Set Up===//
	length = *filter_length;
	wc = (wp+ws)/2.0;
	grid_size = MIN((unsigned int)pow(2, ceil(log2(10*length))),MAX_GRID_SIZE);
	zero_pad_length = 2*(grid_size+1-length)-3;
	tolerance = pow(10,-12);
	passband_error = 1.0 - pow(10,-passband_ripple/20.0);
	stopband_error = pow(10, -stopband_attn/20.0);
	passband_bounds[0] = 1.0 + passband_error; passband_bounds[1] = 1.0 - passband_error;
	stopband_bounds[0] = 0.0 + stopband_error; stopband_bounds[1] = 0.0 - stopband_error;
	
	//===Allocate Equation Objects===//
	A = malloc(3*(length+1)*(length+1)*sizeof(double));	
	
	//===Allocate Alternation Theorem Objects===//
	H = malloc(2*(grid_size) * sizeof(double complex));
	H_real = malloc(2*(grid_size)*sizeof(double));
	freq_response = malloc(2*(grid_size) * sizeof(double));

	//===Make Initial Passband Reference Frequency Set===//	
	passband_refset_size = MAX(1, round(length*wc/M_PI));
	initialize_linspace_dbl(reference_grid, passband_refset_size, 0, passband_refset_size-1);
	for (i=0; i<passband_refset_size; i++){
		reference_grid[i] *= wc/(double)passband_refset_size;
	}
	reference_grid[passband_refset_size] = wc;

	//===Make Initial Stopband Reference Frequency Set===//
	stopband_refset_size = length - passband_refset_size;
	initialize_linspace_dbl(reference_grid+passband_refset_size+1, stopband_refset_size, 1, stopband_refset_size);
	for (i=0; i<stopband_refset_size; i++){
		reference_grid[i+passband_refset_size+1] *= (M_PI-wc)/(double)stopband_refset_size;;
		reference_grid[i+passband_refset_size+1] += wc;
	}
	reference_length = length + 1;

	//===Run Remez Loop===//
	total_error = 1;
	count = 0;
	while(total_error > tolerance){

		//==================================//
		//===Find Response Coefficients=====//
		//==================================//

		j = 0;
		//===Form Desired Response With Room For Ripple===//
		for (i=1; i<=passband_refset_size; i++){
			desired_response[j] = passband_bounds[0] * (1.0 - pow(-1.0,i))/2.0;
			desired_response[j] += passband_bounds[1] * ((pow(-1.0,i) + 1.0)/2.0);		
			j++;
		}
		desired_response[j] = 0.5;
		j++;
		for (i=1; i<=stopband_refset_size; i++){
			desired_response[j] = stopband_bounds[1] * (1.0 - pow(-1.0,i))/2.0;
			desired_response[j] += stopband_bounds[0] * ((pow(-1.0,i) + 1.0)/2.0);
			j++;		
		} 

		//===Form Magnitude Response Cosine Matrix===//
		reference_length = passband_refset_size + stopband_refset_size + 1;
		for (i=0; i<reference_length; i++){
			for (j=0; j<reference_length; j++){
				A[j + i*(reference_length)] = cos(reference_grid[i]*j);
			}
		}

		//===Solve The Linear System===//
		solve_linear_system_dbl(A, reference_length, reference_length, desired_response, a_coeffs);

		//=======================================//
		//===Find Maxima Of Frequency Response===//
		//=======================================//

		//===Make Frequency Response===//
		zero_pad_length = next_pow_2(zero_pad_length);
		zero_pad_length -= reference_length + reference_length - 1;
		j = 0;
		for (i=0; i<reference_length; i++){
			filter[i] = a_coeffs[i];
			if (i>0){
				filter[i] *= 0.5;
			}
			j++;
		}
		for (i=j; i<j+zero_pad_length; i++){
			filter[i] = 0;
		}
		reverse_array_dbl(a_coeffs, reference_length);
		k = 0;
		for (j=i; j<i+(reference_length-1); j++){
			filter[j] = a_coeffs[k];
			k++;			
		}
		reverse_array_dbl(a_coeffs, reference_length);
		fft_dbl(filter, j, H);
		split_array_cmplx(H,j,H_real,NULL);

		//===Find Indices Of Local Maximum Maxima===//
		num_maxima = 0;
		find_local_maxima_dbl(H_real, grid_size, NULL, &num_maxima, NULL, maxima);

		//===Find Indices Of Local Minimum Maxima===//
		negate_vector_dbl(H_real, grid_size);
		temp_uint = num_maxima;
		find_local_maxima_dbl(H_real, grid_size, NULL, &num_maxima, NULL, maxima + num_maxima);
		num_maxima += temp_uint;

		//===Convert Indices Of Response Array To Radians===//
		sort_array_dbl(maxima, num_maxima);
		for (i=0; i<num_maxima; i++){
			reference_grid[i] = (double)((maxima[i]-1))*M_PI/(double)grid_size;
		}
		insert_array_dbl(reference_grid, &num_maxima, wc, 0);
		sort_array_dbl(reference_grid, num_maxima);

		//===Update Reference Set Sizes===//
		passband_refset_size = find_closest_index_dbl(reference_grid, num_maxima, wc);
		stopband_refset_size = num_maxima - passband_refset_size - 1;

		//===================//
		//===Compute Error===//
		//===================//

		//===Get Passband Frequency Response===//
		initialize_array_dbl(A, reference_length*reference_length);
		for (i=0; i<passband_refset_size; i++){
			for (j=0; j<reference_length; j++){
				A[j + i*reference_length] = cos(reference_grid[i]*j);
			}
		}
		initialize_array_dbl(freq_response, reference_length);
		matrix_vector_multiply_dbl(a_coeffs, reference_length,
							 	   A, passband_refset_size, reference_length,
							 	   freq_response);

		//===Compute Error In The Passband===//
		errors[0] = find_maximum_dbl(freq_response, passband_refset_size) - passband_bounds[0];
		errors[1] = passband_bounds[1] - find_minimum_dbl(freq_response, passband_refset_size);

		//===Get Stopband Frequency Response===//
		k = 0;
		initialize_array_dbl(A, reference_length*reference_length);
		for (i=passband_refset_size+1; i<passband_refset_size+1+stopband_refset_size; i++){
			for (j=0; j<reference_length; j++){
				A[j + k*reference_length] = cos(reference_grid[i]*j);
			}
			k++;
		}
		initialize_array_dbl(freq_response, reference_length);
		matrix_vector_multiply_dbl(a_coeffs, reference_length,
							 	   A, stopband_refset_size, reference_length,
							 	   freq_response);
		//===Compute Error In The Stopband===//
		errors[2] = find_maximum_dbl(freq_response, stopband_refset_size) - stopband_bounds[0];
		errors[3] = stopband_bounds[1] - find_minimum_dbl(freq_response, stopband_refset_size);

		//===Get Total Error===//
		total_error = find_maximum_dbl(errors, 4);

		//===Print And Break===//
		count++;
		//fprintf(stdout, "Iteration: %d	Error: %+lf\n", count, total_error);			
		if (count > 1000){
			break;
		}

	}

	//========================================//
	//===Form Filter From Found Coefficients==//
	//========================================//
	reverse_array_dbl(a_coeffs, reference_length);
	for (i=0; i<reference_length-1; i++){
		filter[i] = a_coeffs[i]/2.0;
	}
	reverse_array_dbl(a_coeffs, reference_length);
	filter[i] = a_coeffs[0];
	k = 1;
	for (j=i+1; j<i+reference_length; j++){
		filter[j] = a_coeffs[k]/2.0;
		k++;
	}
	*filter_length = j;
	normalize_sum_dbl(filter, *filter_length);

	
	//===Clean Up===//
	free(freq_response);
	free(H_real);
	free(H);
	free(A);
	
	return;
}

void least_squares_filter_design( double wp,
								  double ws,
								  double passband_weight,
								  double stopband_weight,
								  unsigned int filter_length,
								  double* filter )
{

	unsigned int i,j,n,num_points,N,M,k;
	double x;
	double w[MAX_NUM_QMC_POINTS], wP[MAX_NUM_QMC_POINTS];
	double X[MAX_NUM_QMC_POINTS], a[MAX_POLYNOMIAL_ORDER+1];
	double *Q, *b, *W;

	//===Sample [0,1]===//
	num_points = MAX_NUM_QMC_POINTS; 
	make_korobov_sequence(123421, num_points, w);
	W = malloc(num_points * sizeof(double));
	for (n=0; n<num_points; n++){
		wP[n] = M_PI * w[n];
		if (w[n] <= wp/M_PI){
			W[n] = passband_weight;
		}
		else if (w[n] >= ws/M_PI){
			W[n] = stopband_weight;
		}
		else{
			W[n] = stopband_weight;
		}
	}
	

	//===Force An Odd Filter Length===//
	N = filter_length;
	if (N%2==0){
		N += 1;
		fprintf(stderr, "Error:: Filter Length Must Be Odd! In Function -- least_squares_filter_design!\n");
		quit();
	}
	M = (N+1)/2;
	M -= 1;

	//===Allocate And Form Q===//
	Q = malloc((M+1)*(M+1)*sizeof(double));
	for (k=0; k<M+1; k++){
		for (n=0; n<M+1; n++){

			//===Sample Over [0,M_PI] Interval===//
			initialize_array_dbl(X, num_points);
			for (i=0; i<num_points; i++){
				X[i] += W[i] * cos(n*wP[i]) * cos(k*wP[i]);
			}

			//===Get Mean As Approximation Of Integral===//
			x = compute_mean_dbl(X, num_points);

			//===Set Matrix Element===//	
			Q[n + k*(M+1)] = x;				

		}
	}

	//===Form b===//
	b = malloc((M+1)*sizeof(double));
	for (k=0; k<(M+1); k++){
			//===Sample Over [0,M_PI] Interval===//
			initialize_array_dbl(X, num_points);
			for (i=0; i<num_points; i++){
				if (wP[i] <= wp){
					X[i] += W[i] * 1.0 * cos(k*wP[i]);
				}
				else{
					X[i] += W[i] * 0.0 * cos(k*wP[i]);
				}
			}

			//===Get Mean As Approximation Of Integral===//
			x = compute_mean_dbl(X, num_points);

			//===Set Matrix Element===//	
			b[k] = x;				
	}	

	//===Solve System===//
	solve_linear_system_dbl(Q, (M+1), (M+1), b, a);

	//========================================//
	//===Form Filter From Found Coefficients==//
	//========================================//
	initialize_array_dbl(filter, N);
	reverse_array_dbl(a, M+1);
	for (i=0; i<M+1-1; i++){
		filter[i] = a[i]/2.0;
	}
	reverse_array_dbl(a, M+1);
	filter[i] = a[0];
	k = 1;
	for (j=i+1; j<i+M+1; j++){
		filter[j] = a[k]/2.0;
		k++;
	}
	normalize_sum_dbl(filter, N);

	/*
	newline();
	print_vector_dbl(filter, j, stdout);
	newline();
	fprintf(stdout, "%d %d\n", j, N);
	*/

	//===Clean Up===//
	free(W); free(b); free(Q);

	return;
}

void fractional_delay_filter_design( double delay,
									 double sampling_rate,
									 unsigned int filter_length,
									 double* filter )
{
	unsigned int t;
	double x, center, samples, value;
	double window[MAX_POLYNOMIAL_ORDER+1];
	double sinc[MAX_POLYNOMIAL_ORDER+1];

	//===Convert Delay In Seconds To Samples===//
	samples = delay * sampling_rate;

	//===Create Window Function===//
	center = (double)filter_length/2; //(to make the sinc causal b/c as mathematically defined, it is not causal)
	initialize_blackman_window(filter_length, window);
	//initialize_fyrley_window(filter_length, window);
	
	//===Make Delayed Sinc Function===// 
	for (t=0; t<filter_length; t++){
		
		x = (double)t - samples;
		if (fabs(x - center) < 10e-3){
			value = 1.0;
		}
		else{
			value = sin(M_PI*(x-center));
			value /= (M_PI * (x-center));
		}
		sinc[t] = value;
	}
	
	//===Create Filter===//
	hadamard_product_dbl( window, sinc, filter_length, filter);

	return;
}

void preemphasis_filter_design( double alpha,
							   double* filter )
{
	if (filter == NULL){
		fprintf(stderr, "Error:: Filter Is NULL! In Function -- preemphasis_filter_design!\n");
		quit();
	}
	filter[0] = 1;
	filter[1] = -alpha;
	return;
}
								


//------------------------------------------------------------------------------------------------//
//====================================FILTER FUNCTIONS============================================//
//------------------------------------------------------------------------------------------------//

void initialize_filter( filter_t* self,
						double* numerator,
						unsigned int numerator_length,
						double* denominator,
						unsigned int denominator_length)
{
	if (numerator != NULL){
		copy_array_dbl(numerator, numerator_length, self->numerator);
		self->numerator_length = numerator_length;
	}
	else{
		self->numerator_length = 1;
		initialize_array_dbl(self->numerator, self->numerator_length);
		self->numerator[0] = 1;
	}	
	if (denominator != NULL){
		copy_array_dbl(denominator, denominator_length, self->denominator);
		self->denominator_length = denominator_length;
	}
	else{
		self->denominator_length = 1;
		initialize_array_dbl(self->denominator, self->denominator_length);
		self->denominator[0] = 1;
	}
	self->internal_length = MAX(denominator_length+1, numerator_length+1);
	initialize_array_dbl(self->internal, 2*(MAX_POLYNOMIAL_ORDER+1));
	return;
}

void reset_filter( filter_t* self )
{	
	initialize_array_dbl(self->internal, 2*(MAX_POLYNOMIAL_ORDER+1));
	return;
}

void copy_filter( filter_t* original,
				  filter_t* copy )
{
	copy->numerator_length = original->numerator_length;
	copy->denominator_length = original->denominator_length;
	copy->internal_length = original->internal_length;
	copy_array_dbl(original->numerator, MAX_POLYNOMIAL_ORDER+1, copy->numerator);
	copy_array_dbl(original->denominator, MAX_POLYNOMIAL_ORDER+1, copy->denominator);
	copy_array_dbl(original->internal, 2*(MAX_POLYNOMIAL_ORDER+1), copy->internal);
	return;
}

void run_filter( filter_t* self,
				 double* input,
				 unsigned int input_length,
				 double* output )
{
		unsigned int i, j, M, L, X, Y, KMAX, KMIN;
		double x, y;
		double *w, *a, *b;
		
		//===Set Locals For Easy Referencing===//
		X = input_length;
		Y = input_length;
		M = self->denominator_length;
		L = self->numerator_length;
		a = self->denominator;
		b = self->numerator;
		w = self->internal;

		//===Filter The Entire Array===//
		KMAX = MAX(M,L);		
		KMIN = MIN(X+M+L, Y);
		for (i=0; i<KMIN; i++){

			//===Grab A Single Value===//
			if (i > input_length){
				x = 0;
			}
			else{
				x = input[i];
			}

			//===Run Input Adder===//
			w[0] = x;
			for (j=1; j<M; j++){
				w[0] -= a[j] * w[j];
			}

			//===Run Output Adder===//
			y = 0;
			for (j=0; j<L; j++){
				y += b[j] * w[j];
			}
			
			//===Update The Internal State Vector===//
			for (j=KMAX-1; j>=1; j--){
				w[j] = w[j-1];
			}

			//===Save Output===//
			output[i] = y;			
		}
	
	return;
}

void run_filter_zero_phase( filter_t* self,
				 		    double* input,
				 			unsigned int input_length,
				 			double* output )
{
	unsigned int group_delay;
	double* copy;

	//NOTE: Input should have 'group_delay' number of zeros appended at the end
	//AND: the input_length have group_delay built into it

	//===Mallocs===//
	copy = malloc((MAX_POLYNOMIAL_ORDER+input_length)*sizeof(double));

	//===Filter===//
	group_delay = MAX(((self->numerator_length-1)/2), ((self->denominator_length-1)/2));
	run_filter(self, input, input_length, copy);
	left_shift_array_dbl(copy, input_length+group_delay, group_delay);

	//===Reverse And Filter And Reverse===//
	reverse_array_dbl(copy, input_length);
	run_filter(self,copy, input_length+group_delay, output);
	left_shift_array_dbl(output, input_length+group_delay, group_delay);
	reverse_array_dbl(output, input_length);

	//===Clean Up===//
	free(copy);

	return;
}

void compute_impulse_response( filter_t* self,
							   unsigned int impulse_length,
							   double* impulse_response )
{
	double impulse[MAX_FRAME_LENGTH];

	//===Create Impulse===//
	impulse[0] = 1;
	initialize_array_dbl(impulse+1, impulse_length-1);

	//===Run Filter===//
	run_filter(self, impulse, impulse_length, impulse_response);
	
	return;
}

void compute_frequency_response( filter_t* self,
								 unsigned int response_length,
								 double complex* frequency_response )
{

	unsigned int i;
	double step_size;
	double *w_axis;
	double complex *z_path;
	double complex a[MAX_POLYNOMIAL_ORDER], b[MAX_POLYNOMIAL_ORDER];

	//===Set Up===//	
	w_axis = malloc(response_length*sizeof(double));
	z_path = malloc(response_length*sizeof(double complex));
	step_size = M_PI/(double)(response_length-1);

	//===Create Complex Arrays For Evaluation===//
	combine_arrays_dbl_to_cmplx(self->numerator, NULL, self->numerator_length, b);
	combine_arrays_dbl_to_cmplx(self->denominator, NULL, self->denominator_length, a);

	//===Evaluate Rational Function Around The Unit Circle===//
	for (i=0; i<response_length; i++){
		w_axis[i] = i*step_size;
		z_path[i] = cexp(-1.0*I * w_axis[i]);
		frequency_response[i] = evaluate_polynomial_cmplx(b, self->numerator_length, z_path[i]);
		frequency_response[i] /= evaluate_polynomial_cmplx(a, self->denominator_length, z_path[i]);
	}
	
	//===Clean Up===//
	free(z_path);
	free(w_axis);

	return;
}

//------------------------------------------------------------------------------------------------//
//===================================RESAMPLING FUNCTIONS=========================================//
//------------------------------------------------------------------------------------------------//

void initialize_resampler( resampler_t* self,
						   unsigned int downsampling_factor,
						   unsigned int upsampling_factor,
						   unsigned int block_size,
						   resampler_type type )
{
	unsigned int i, filter_length;
	double wp, ws;
	double filter_coeffs[MAX_POLYNOMIAL_ORDER+1];

	//===Set Locals===//
	if (downsampling_factor == 0) downsampling_factor = 1;
	self->downsampling_factor = downsampling_factor;
	if (upsampling_factor == 0) upsampling_factor = 1;
	self->upsampling_factor = upsampling_factor;
	self->type = type;
	self->block_size = block_size;

	//===Design Filter===//
	wp = 0.875 * (1.0 / ((double)MAX(upsampling_factor,downsampling_factor))) * M_PI;
	ws = 1.025 * (1.0 / ((double)MAX(upsampling_factor,downsampling_factor))) * M_PI;
	filter_length = 255;
	minimax_filter_design(wp, ws, 120, 0.005, &filter_length, filter_coeffs);
	
	//===Initialize Filter===//
	initialize_filter(&(self->resampling_filter), filter_coeffs, filter_length, NULL, 0);

	//===Make Polyphase Filters===//
	if (type == DOWNSAMPLE){

		//===Make Filter===//
		initialize_filter(&(self->resampling_filter), filter_coeffs, filter_length, NULL, 0);

		//===Make Filter Length Polyphaseable===//
		while (filter_length % downsampling_factor != 0 ){
			pad_zeros_dbl(filter_coeffs, filter_length, 1);
			filter_length += 1;
		}

		//===Make Polyphase Filters===//
		for (i=0; i<downsampling_factor; i++){
			copy_array_dbl(self->resampling_filter.numerator, filter_length, filter_coeffs);
			left_shift_array_dbl(filter_coeffs, filter_length, i);
			downsample_array_dbl(filter_coeffs, filter_length, downsampling_factor);
			initialize_filter(&(self->polyphase_resampling_filters[i]), filter_coeffs, 
								filter_length/downsampling_factor, NULL, 0);			
		}

		//===Set Locals===//
		self->output_queue_length = block_size;
		self->input_queue_length = downsampling_factor*self->output_queue_length;

		//===Initialize Queues===//
		initialize_fifo_queue(&(self->input_queue), self->input_queue_length);
		initialize_fifo_queue(&(self->output_queue), self->output_queue_length);
		initialize_fifo_queue(&(self->tdl_queue), self->output_queue_length*self->input_queue_length);
	}
	else if (type == UPSAMPLE){

		//===Make Filter===//
		initialize_filter(&(self->resampling_filter), filter_coeffs, filter_length, NULL, 0);

		//===Make Filter Length Polyphaseable===//
		while (filter_length % upsampling_factor != 0 ){
			pad_zeros_dbl(filter_coeffs, filter_length, 1);
			filter_length += 1;
		}

		//===Make Polyphase Filters===//
		for (i=0; i<upsampling_factor; i++){
			copy_array_dbl(self->resampling_filter.numerator, filter_length, filter_coeffs);
			left_shift_array_dbl(filter_coeffs, filter_length, i);
			downsample_array_dbl(filter_coeffs, filter_length, upsampling_factor);
			initialize_filter(&(self->polyphase_resampling_filters[i]), filter_coeffs, 
								filter_length/upsampling_factor, NULL, 0);			
		}

		//===Set Locals===//
		self->input_queue_length = block_size;
		self->output_queue_length = upsampling_factor * self->input_queue_length;

		//===Initialize Queues===//
		initialize_fifo_queue(&(self->input_queue), self->input_queue_length);
		initialize_fifo_queue(&(self->output_queue), self->output_queue_length);
		initialize_fifo_queue(&(self->tdl_queue), self->output_queue_length*self->input_queue_length);
	}
	else if (type == RATIONAL){

		//===Make Filter===//
		initialize_filter(&(self->resampling_filter), filter_coeffs, filter_length, NULL, 0);

		//===Make Filter Length Polyphaseable===//
		while (filter_length % upsampling_factor != 0 ){
			pad_zeros_dbl(filter_coeffs, filter_length, 1);
			filter_length += 1;
		}

		//===Make Polyphase Filters===//
		for (i=0; i<upsampling_factor; i++){
			copy_array_dbl(self->resampling_filter.numerator, filter_length, filter_coeffs);
			left_shift_array_dbl(filter_coeffs, filter_length, i);
			downsample_array_dbl(filter_coeffs, filter_length, upsampling_factor);
			initialize_filter(&(self->polyphase_resampling_filters[i]), filter_coeffs, 
								filter_length/upsampling_factor, NULL, 0);			
		}

		//===Set Locals===//
		self->input_queue_length = downsampling_factor*block_size;
		self->output_queue_length = upsampling_factor*block_size;

		//===Initialize Queues===//
		initialize_fifo_queue(&(self->input_queue), self->input_queue_length);
		initialize_fifo_queue(&(self->output_queue), self->output_queue_length);
		initialize_fifo_queue(&(self->tdl_queue), self->input_queue_length*self->output_queue_length);
	}


	return;
}

void enqueue_resampler_dbl( resampler_t* self,
							double* data )
{
	enqueue_fifo_dbl(&(self->input_queue), data, self->input_queue_length);
	return;
}

void dequeue_resampler_dbl( resampler_t* self,
							double* data )
{
	dequeue_fifo_dbl(&(self->output_queue), self->output_queue_length, data);
	return;
}

void run_resampler( resampler_t* self )
{
	unsigned int i, j, k;
	double sum[MAX_RESAMPLING_QUEUE_LENGTH];
	double filter_data[MAX_RESAMPLING_QUEUE_LENGTH];
	double filter_input[MAX_RESAMPLING_FACTOR][MAX_RESAMPLING_QUEUE_LENGTH];
	double filter_output[MAX_RESAMPLING_FACTOR][MAX_RESAMPLING_QUEUE_LENGTH];
	filter_t* filter;

	if (self->type == DOWNSAMPLE){

		//===Get Data Filters===//
		initialize_array_dbl(sum, MAX_RESAMPLING_QUEUE_LENGTH);
		dequeue_fifo_dbl(&(self->input_queue), self->input_queue_length, filter_data);

		//===Split Up Filter Input===//
		j = 0; k = 0;
		for (j=0; j<self->output_queue_length; j++){
			for (i=0; i<self->downsampling_factor; i++){
				filter_input[i][j] = filter_data[k];
				k++;
			}
		}

		//===Run Filters===//
		for (i=0; i<self->downsampling_factor; i++){

			//===Get Filter===//
			filter = &(self->polyphase_resampling_filters[(self->downsampling_factor-1)-i]);

			//===Feed Filters===//
			run_filter(filter, filter_input[i], self->output_queue_length, filter_output[i]);
			for (j=0; j<self->output_queue_length; j++){	
				sum[j] += filter_output[i][j];
			}
		}

		//===Place In Output Queue===//
		gain_array_constant_dbl(sum, self->output_queue_length, 1.0/(double)self->downsampling_factor);
		enqueue_fifo_dbl(&(self->output_queue), sum, self->output_queue_length);
		
	}		
	else if (self->type == UPSAMPLE){

		//===Run Filters===//
		dequeue_fifo_dbl(&(self->input_queue), self->input_queue_length, filter_data);
		for (i=0; i<self->upsampling_factor; i++){

			//===Get Filter===//
			filter = &(self->polyphase_resampling_filters[i]);

			//===Feed Filters===//
			run_filter(filter, filter_data, self->input_queue_length, filter_output[i]);

		}

		//===Place Into Queue===//
		for (j=0; j<self->input_queue_length; j++){
			for (i=0; i<self->upsampling_factor; i++){
				enqueue_fifo_dbl(&(self->output_queue), &(filter_output[i][j]), 1);
			}
		}

	}
	else if (self->type == RATIONAL){

		//===Run Filters===//
		dequeue_fifo_dbl(&(self->input_queue), self->input_queue_length, filter_data);
		for (i=0; i<self->upsampling_factor; i++){

			//===Get Filter===//
			filter = &(self->polyphase_resampling_filters[i]);

			//===Feed Filters===//
			run_filter(filter, filter_data, self->input_queue_length, filter_output[i]);

		}

		//===Place Into TDL===//
		k = 0;
		for (k=0; k<self->block_size; k++){
			for (j=0; j<self->downsampling_factor; j++){
				for (i=0; i<self->upsampling_factor; i++){
					enqueue_fifo_dbl(&(self->tdl_queue), &(filter_output[i][k*self->downsampling_factor + j]), 1);
				}
			}
		}

		//===Dequeue And Downsample TDL===//
		if (fifo_queue_is_full(&(self->tdl_queue))){
			for (k=0; k<self->block_size; k++){

				//===Dequeue===//
				initialize_array_dbl(filter_data, MAX_RESAMPLING_QUEUE_LENGTH);
				dequeue_fifo_dbl(&(self->tdl_queue), self->downsampling_factor*self->upsampling_factor, filter_data);
				downsample_array_dbl(filter_data, self->downsampling_factor*self->upsampling_factor, self->downsampling_factor);

				//===Enqueue Into Output Buffer===//
				enqueue_fifo_dbl(&(self->output_queue), filter_data, self->upsampling_factor);
			}

		}
	}
	return;
}

void resample_array_dbl( resampler_t* self,
						 double* original_signal,
						 unsigned int original_signal_length,
						 unsigned int downsampling_factor,
						 unsigned int upsampling_factor,
						 resampler_type type,
						 double* resampled_signal )
{
	unsigned int i, j;
	double data[MAX_FRAME_LENGTH];

	//===Init Resampler===//
	initialize_resampler(self, downsampling_factor, upsampling_factor, 1, type);
	
	//===Resample===//
	i = 0; j = 0;
	while (i<original_signal_length){

		//===Put Data Into Input Queue===//
		while(!fifo_queue_is_full(&(self->input_queue))){

			//===Grab Data===//
			copy_array_dbl(original_signal + i, self->input_queue_length, data);

			//===Enqueue Into The Resampler===//
			enqueue_resampler_dbl(self, data);

			//===Update Input Ptr===//
			i += self->input_queue_length;	
		}

		//===Run Filters===//
		run_resampler(self);

		//===Output Data===//
		if(fifo_queue_is_full(&(self->output_queue))){

			//===Output Data From Resampler===//
			dequeue_resampler_dbl(self, data);

			//===Copy To Output Stream===//
			copy_array_dbl(data, self->output_queue_length, resampled_signal + j);
		
			//===Update Output Ptr===//
			j += self->output_queue_length;
	
		}

	}

	return;
}

void integer_upsample_dbl( double* original_signal,
						   unsigned int original_signal_length,
						   unsigned int upsample_factor,
						   double* upsampled_signal )
{

	unsigned int j, i, filter_length;
	double wp, ws, max;
	double lowpass_coeffs[MAX_POLYNOMIAL_ORDER+1];
	double* temp;
	filter_t filter;

	//===Mallocs===//
	temp = malloc(upsample_factor * original_signal_length * sizeof(double));

	//===Get Maximum===//
	max = find_maximum_dbl(original_signal, original_signal_length);

	//===Design Filter===//
	wp = 0.875 * (1.0 / ((double)upsample_factor)) * M_PI;
	ws = 1.025 * (1.0 / ((double)upsample_factor)) * M_PI;
	filter_length = 255;
	minimax_filter_design(wp, ws, 120, 0.005, &filter_length, lowpass_coeffs);
	
	//===Upsample===//
	initialize_array_dbl(temp, upsample_factor*original_signal_length);
	j = 0;
	for (i=0; i<original_signal_length; i++){
		temp[j] = original_signal[i];
		j += upsample_factor;
	}
	
	//===Run Filter===//
	initialize_filter(&filter, lowpass_coeffs, filter_length, NULL, 0);
	run_filter(&filter, temp, upsample_factor*original_signal_length, upsampled_signal);

	//===Normalize===//
	normalize_max_dbl(upsampled_signal, upsample_factor*original_signal_length);
	gain_array_constant_dbl(upsampled_signal, upsample_factor*original_signal_length, max);

	//===Clean Up===//
	free(temp);

	return;
}


//------------------------------------------------------------------------------------------------//
//===================================ALL PASS FUNCTIONS===========================================//
//------------------------------------------------------------------------------------------------//

void delay_signal_dbl( double* signal,
					   unsigned int signal_length,
					   double delay,
					   double sampling_rate,
					   double* delayed )
{
	unsigned int filter_length;
	double group_delay;
	double coeffs[2*(MAX_POLYNOMIAL_ORDER+1)];
	filter_t filter;

	//===Make Delay Filter===//
	filter_length = 256;
	fractional_delay_filter_design(delay, sampling_rate, filter_length, coeffs);
	initialize_filter(&filter, coeffs, filter_length, NULL, 0);

	//===Run Filter===//
	initialize_array_dbl(delayed, signal_length + filter_length);
	run_filter(&filter, signal, signal_length + filter_length, delayed);

	//===Take Care of Group Delay===//
	if (filter_length % 2){
		group_delay = MAX(((filter.numerator_length-1)/2), ((filter.denominator_length-1)/2));
	}
	else{
		group_delay = MAX(((filter.numerator_length)/2), ((filter.denominator_length)/2));
	}
	left_shift_array_dbl(delayed, signal_length+group_delay, group_delay);

	return;
}

void low_frequency_emphasis_dbl( double* signal,
								 unsigned int signal_length,
								 double coeff,
								 double* emphasized )
{
	double* copy;
	double coeffs[MAX_POLYNOMIAL_ORDER];
	filter_t filter;

	//===Mallocs===//
	copy = malloc(signal_length * sizeof(double));

	//===Make Filter===//	
	coeffs[0] = 1.0; coeffs[1] = coeff;
	initialize_filter(&filter, NULL, 0, coeffs, 2);

	//===Run Filter===//
	copy_array_dbl(signal, signal_length, copy);
	initialize_array_dbl(emphasized, signal_length);
	run_filter(&filter, copy, signal_length, emphasized);

	//===Clean Up===//
	free(copy);

	return;
}

//------------------------------------------------------------------------------------------------//
//====================================TEST FUNCTIONS==============================================//
//------------------------------------------------------------------------------------------------//

void test_filter()
{
	unsigned int i, impulse_length, step_length, a_length, b_length;
	double a[MAX_POLYNOMIAL_ORDER+1], b[MAX_POLYNOMIAL_ORDER+1];
	double impulse[MAX_FRAME_LENGTH], output[MAX_FRAME_LENGTH];
	double step[MAX_FRAME_LENGTH];
	filter_t filter;
	FILE* fout;

	//===Make Impulse Input===//
	impulse_length = 512;
	impulse[0] = 1;
	initialize_array_dbl(impulse+1, impulse_length-1);

	//===Make Step Input===//
	step_length = 512;
	initialize_array_constant_dbl(step, step_length, 1.0);

	//===Make Stable Roots===//
	double complex temp_coeffs[MAX_POLYNOMIAL_ORDER+1];
	double complex roots[MAX_POLYNOMIAL_ORDER];

	roots[0] = 0.75 * cexp(I * 2.0 * M_PI * 0.5);
	roots[1] = 0.35 * cexp(I * 2.0 * M_PI * 0.95);
	roots[2] = 0.55 * cexp(I * 2.0 * M_PI * -0.35);

	//===Convert Back To Process===//
	convert_roots_to_coefficients(roots, 3, temp_coeffs);
	split_array_cmplx(temp_coeffs, 4, a, NULL);

	//===Make Coefficients===//
	a_length = 4; b_length = 1;
	//a[0] = 1.0; a[1] = 0.775; a[2] = -0.425; a[3] = -0.5; a[4] = 0.1; a[5] = -0.255;
	b[0] = 1.0; b[1] = 0.5; b[2] = 0.7; b[3] = 0.9;

	//gain_array_constant_dbl(a, a_length, 2.0);

	//===Run Filter===//
	initialize_filter(&filter, b, b_length, a, a_length);
	initialize_array_dbl(output, MAX_FRAME_LENGTH);
	run_filter(&filter, impulse, impulse_length, output);

	//===Print Result===//
	fout = fopen("impulse.dat", "w");
	print_vector_dbl(output, impulse_length, fout);
	fclose(fout);
	
	//===Run Filter Real Time===//	
	reset_filter(&filter);
	initialize_filter(&filter, b, b_length, a, a_length);
	initialize_array_dbl(output, MAX_FRAME_LENGTH);
	for (i=0; i<impulse_length; i++){		
		run_filter(&filter, impulse+i, 1, output+i);
	}

	//===Print Result===//
	fout = fopen("impulse_real_time.dat", "w");
	print_vector_dbl(output, impulse_length, fout);
	fclose(fout);

	//===Run Step Filter===//
	initialize_filter(&filter, b, b_length, a, a_length);
	initialize_array_dbl(output, MAX_FRAME_LENGTH);
	run_filter(&filter, step, step_length, output);

	//===Print Result===//
	fout = fopen("step.dat", "w");
	print_vector_dbl(output, step_length, fout);
	fclose(fout);
		
	return;
}

void test_bandpass()
{
	unsigned int modulate;
	double ideal[MAX_POLYNOMIAL_ORDER+1];
	double window[MAX_POLYNOMIAL_ORDER+1];
	double filter[MAX_POLYNOMIAL_ORDER+1];
	double wave[MAX_POLYNOMIAL_ORDER+1];
	double freq, amp, phase;
	FILE* fout;

	modulate = 1;
	if (modulate){
		initialize_lowpass_ideal(100, 8000, 128, ideal);
		initialize_chebyshev_window(50, 128, window);
		hadamard_product_dbl(ideal, window, 128, filter);
		freq = 2000.0; amp = 1.0; phase = 0.0;
		generate_wave_dbl(&freq, &amp, &phase, 1, 128, 8000.0, wave);
		hadamard_product_dbl(filter, wave, 128, filter);
	}
	else{
		initialize_bandpass_ideal(1000, 2000, 8000, 128, ideal);
		initialize_chebyshev_window(50, 128, window);
		hadamard_product_dbl(ideal, window, 128, filter);
	}

	fout = fopen("bandpass.dat","w");
	print_vector_dbl(filter, 128, fout);
	fclose(fout);

	return;
}

void test_chebyshev()
{
	unsigned int order, response_length;
	double wp, Op;
	double passband_ripple_dB;
	FILE* fout;
	double a[MAX_POLYNOMIAL_ORDER];
	double b[MAX_POLYNOMIAL_ORDER];
	double impulse_response[MAX_SIGNAL_LENGTH];  
	filter_t filter;

	//===Setup Chebyshev Design===//
	passband_ripple_dB = 0.01;
	wp = 0.2*M_PI; Op = prewarp_frequency_dbl(wp);
	
	//===Run Chebyshev Design===//
	order = 2;
	chebyshev_design(order, passband_ripple_dB, Op, b, a);

	//===Print Results===//
	fout = fopen("b.dat", "w");
	print_vector_dbl(b, 1, fout);
	fclose(fout);
	fout = fopen("a.dat", "w");
	print_vector_dbl(a, order+1, fout);
	fclose(fout);

	//===Setup Filter===//
	initialize_filter(&filter, b, 1, a, order+1);

	//===Make Digital Filter===//
	bilinear_transform(&filter);
	
	//===Print Results===//
	fout = fopen("numerator.dat", "w");
	print_vector_dbl(filter.numerator, filter.numerator_length, fout);
	fclose(fout);
	fout = fopen("denominator.dat", "w");
	print_vector_dbl(filter.denominator, filter.denominator_length, fout);
	fclose(fout);

	//===Get Impulse Response===//
	response_length = 128;
	compute_impulse_response(&filter, response_length, impulse_response);

	//===Print Results===//
	fout = fopen("impulse.dat", "w");
	print_vector_dbl(impulse_response, response_length, fout);
	fclose(fout);

	return;
}

void test_butterworth()
{
	unsigned int order, response_length;
	double wp, Op;
	FILE* fout;
	double a[MAX_POLYNOMIAL_ORDER];
	double b[MAX_POLYNOMIAL_ORDER];
	double impulse_response[MAX_SIGNAL_LENGTH];  
	filter_t filter;

	//===Setup Butterworth Design===//
	wp = (1000.0/8000.0) * 2.0 * M_PI;
	Op = prewarp_frequency_dbl(wp);

	//===Run Buttworth Design===//
	order = 1;
	butterworth_design(order, Op, b, a);
	
	//===Print Results===//
	fout = fopen("b.dat", "w");
	print_vector_dbl(b, 1, fout);
	fclose(fout);
	fout = fopen("a.dat", "w");
	print_vector_dbl(a, order+1, fout);
	fclose(fout);

	//===Setup Filter===//
	initialize_filter(&filter, b, 1, a, order+1);

	//===Make Digital Filter===//
	bilinear_transform(&filter);
	
	//===Print Results===//
	fout = fopen("numerator.dat", "w");
	print_vector_dbl(filter.numerator, filter.numerator_length, fout);
	fclose(fout);
	fout = fopen("denominator.dat", "w");
	print_vector_dbl(filter.denominator, filter.denominator_length, fout);
	fclose(fout);

	//===Get Impulse Response===//
	response_length = 128;
	compute_impulse_response(&filter, response_length, impulse_response);

	//===Print Results===//
	fout = fopen("impulse.dat", "w");
	print_vector_dbl(impulse_response, response_length, fout);
	fclose(fout);

	return;
}

void view_filter_types()
{
	int order;
	double fc, Fs;
	double low_coeffs[MAX_POLYNOMIAL_ORDER+1];
	double high_coeffs[MAX_POLYNOMIAL_ORDER+1];
	double band_coeffs[MAX_POLYNOMIAL_ORDER+1];
	double window[MAX_POLYNOMIAL_ORDER+1];
	FILE* fout;

	//lowpass
	fc = 300;
	Fs = 16000;
	order = 1001;
	initialize_array_dbl(low_coeffs, MAX_POLYNOMIAL_ORDER+1);
	initialize_blackman_window(order, window);
	window_filter_design(LOWPASS, 0, fc, Fs, order, window, low_coeffs);
	fout = fopen("lowpass.dat", "w");
	print_vector_dbl(low_coeffs, order, fout);
	fclose(fout);

	//highpass
	fc = 3400;
	Fs = 16000;
	initialize_array_dbl(high_coeffs, MAX_POLYNOMIAL_ORDER+1);
	initialize_blackman_window(order, window);
	window_filter_design(HIGHPASS, fc, Fs, Fs, order, window, high_coeffs);
	fout = fopen("highpass.dat", "w");
	print_vector_dbl(high_coeffs, order, fout);
	fclose(fout);

	//bandpass
	window_filter_design(BANDPASS, 300.0, 3400.0, Fs, order, window, band_coeffs);
	fout = fopen("bandpass.dat", "w");
	print_vector_dbl(band_coeffs, order, fout);
	fclose(fout);

	//window
	fout = fopen("window.dat", "w");
	print_vector_dbl(window, order, fout);
	fclose(fout);


	return;
}

void test_notch_filter()
{
	unsigned int num_freqs, filter_length, response_length;
	double sampling_rate;
	double impulse_response[2048];
	double denominator[101], numerator[101];
	double frequencies[100];
	filter_t filter;
	FILE* fout;

	//===Make Frequencies===//
	num_freqs = 4; sampling_rate = 16000.0;
	frequencies[0] = 1000.0;
	frequencies[1] = 2750.0;
	frequencies[2] = 350.0;
	frequencies[3] = 1500.0;

	//delay caused by filter = filter_length/2

	//===Make Filter===//
	notch_filter_design(frequencies, num_freqs, sampling_rate, 10.0, &filter_length, numerator, denominator);
	initialize_filter(&filter, numerator, filter_length, denominator, filter_length);

	//===Get Impulse Response===//
	response_length = 2048;
	compute_impulse_response(&filter, response_length, impulse_response);

	//===Print Results===//
	fout = fopen("impulse.dat", "w");
	print_vector_dbl(impulse_response, response_length, fout);
	fclose(fout);

	return;
}

void test_minimax()
{
	unsigned int filter_length, response_length;
	double sampling_rate, nyquist_rate;
	double freqss[10], amps[10], phases[10];
	double modulator[1000];
	double wp_design, ws_design, wp, ws, stopband_attn, passband_ripple;
	double coeffs[1000];
	double impulse_response[4096];
	double freqs[4096];
	filter_t filter;
	FILE* fout;

	//===Design===//
	sampling_rate = 8000.0; nyquist_rate = sampling_rate/2.0;
	wp_design = 100.0; ws_design = 200.0; 
	stopband_attn = 50.0; passband_ripple = 0.1;
	filter_length = 271;

	//===Run===//
	wp = (wp_design/nyquist_rate)*M_PI;	ws = (ws_design/nyquist_rate)*M_PI; 
	minimax_filter_design(wp, ws, stopband_attn, passband_ripple, &filter_length, coeffs);
	//right_shift_array_dbl(coeffs, filter_length, 1); left_shift_array_dbl(coeffs, filter_length, 2);
	//filter_length -= 2;

	//===Modulate===//
	if (0){
		freqss[0] = 1000.0; amps[0] = 2.0; phases[0] = 0.0;
		generate_wave_dbl(freqss, amps, phases, 1, filter_length, sampling_rate, modulator);
		hadamard_product_dbl(coeffs, modulator, filter_length, coeffs);
	}

	//===Initialize===//
	initialize_filter(&filter, coeffs, filter_length, NULL, 0);

	//===Get Impulse Response===//
	response_length = 2048;
	compute_impulse_response(&filter, response_length, impulse_response);

	//===Get Frequencies===//
	frequencies(response_length, sampling_rate, freqs);

	//===Print Results===//
	fout = fopen("minimax_impulse.dat", "w");
	print_vector_dbl(impulse_response, response_length, fout);
	fclose(fout);
	fout = fopen("minimax_freqs.dat", "w");
	print_vector_dbl(freqs, response_length, fout);
	fclose(fout);
	fout = fopen("minimax_coeffs.dat", "w");
	print_vector_dbl(coeffs, filter_length, fout);
	fclose(fout);


	return;
}
