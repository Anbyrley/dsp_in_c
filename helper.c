#include "helper.h"

//================================================================================================//
//=====================================SPECIAL FUNCTIONS==========================================//
//================================================================================================//

double lcm_dbl( double a,
				double b )
{
	if (fabs(a) <= 10.0*EPS){
		return 0;
	}
	else{
		return (a*b)/gcd_dbl(a,b);
	}

	return 0;
}

double lcm_array_dbl( double* array,
					  unsigned int length )
{
	unsigned int i;
	double b;

	//===Recurse===//
	b = array[0];
	for (i=1; i<length; i++){
		b = lcm_dbl(array[i], b);
	}

	return b;
}

double gcd_dbl_2( double a, 
			      double b )
{
	//===Check===//
	if (fabs(a) < 10.0*EPS){
		return b;
	}
	while (b > 0 && fabs(b) > 10.0*EPS){
		if (a > b){
			a = a - b;
		}
		else{
			b = b - a;
		}
	}
	return a;	
}

double gcd_dbl( double a, 
			    double b )
{
	//===Recurse===//
    if (a < b){
	    return gcd_dbl(b, a);
	}
 
	//===Recurse===//
    if (fabs(b) < 10e-3){
        return a; 
	}
    else{
        return (gcd_dbl(b, a - floor(a / b) * b));
	}

	return 0;
}

double gcd_array_dbl( double* array,
					  unsigned int length )
{
	unsigned int i;
	double b;

	//===Recurse===//
	b = array[0];
	for (i=1; i<length; i++){
		b = gcd_dbl(array[i], b);
	}

	return b;
}

double degree_to_radian_dbl( double degree )
{
	double radian;
	radian = M_PI * degree / 180.0;
	return radian;
}

double diriac_dbl( double x )
{
	if (fabs(x) < 10.0*EPS){
		return 1;
	}
	else{
		return 0;
	}
	return 0;
}

double modulo_dbl( double a, 
				   double q )
{
    double b;
	b = a / q;
    return (b - floor(b)) * q;
}

double complex sinc_cmplx( double complex x )
{
	double complex sinc;
	
	sinc = csin(M_PI * x)/(M_PI * x);
	if (isnan(creal(sinc)) || isnan(cimag(sinc))){
		sinc = 1.0;
	}

	return sinc;
}

double sinc_dbl( double x )
{
	double sinc;

	sinc = sin(M_PI * x)/(M_PI*x);
	if (isnan(sinc)){
		sinc = 1.0;
	}
	return sinc;
}

double lower_incomplete_gamma_dbl( double x, 
								   double p )
{
	double a, arg, c, f;
	double value, uflo, e;

	//===Set Locals===//
	uflo = 1.0E-37;
	e = 1.0E-09;

	//===Input Check===//
	if (x<=0.0 || p<=0.0){
		fprintf(stderr, "Error: X<=0 or P<=0! In Function -- lower_incomplete_gamma_dbl!\n");
		quit();
	}		
	arg = p*log(x) - gamma_ln(p+1.0) - x;
	if (arg<log(uflo)){
		fprintf(stderr, "Error: Underflow Detected! In Function -- lower_incomplete_gamma_dbl!\n");
		quit();
	}
	f = exp(arg);
	if (fabs(f)<=10.0*EPS){
		fprintf(stderr, "Error: Underflow Detected! In Function -- lower_incomplete_gamma_dbl!\n");
		quit();
	}

	//===Begin Infinite Series===//
	c = 1.0; value = 1.0; a = p;
	while(1){
		a = a + 1.0;
		c = c * x/a;
		value = value + c;
		if (c<=e*value){
		  break;
		}
 	 }
	value *= f;
	return value;
}


double gamma_ln(double xx)
{
	int j;
	double x,y,tmp,ser;
	double coeff[6];

	if (xx < 0){
		fprintf(stdout, "Error: Input Must Be Greater Than Zero!\n");
		return 0;
	}

	coeff[0] = 76.1800917294146;
	coeff[1] = -86.50532032941677;
	coeff[2] = 24.01409824083091;
	coeff[3] = -1.231739572450155;
	coeff[4] = 0.1208650973866179e-2;
	coeff[5] = -0.5395239384953e-5;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5)*log(tmp);
	ser = 1.000000000190015;
	for (j=0; j<=5; j++) ser += coeff[j]/(++y);
	return -tmp + log(2.5066282746310005*ser/x);

}

double gamma_dbl(double xx)
{
	if (xx > 0){
		return exp(gamma_ln(xx));
	}
	else{
		return M_PI / ( sin(M_PI * xx) * (-xx) * exp(gamma_ln(-xx)));
	}
}

void get_legendre_coeffs( unsigned int order,
						  double* coeffs )
{
	unsigned int n, k;
	double coeff, tmp;

	n = order;
	for (k=0; k<=n; k++){
		tmp = (n+k-1); tmp /= 2;
		coeff = (n_choose_k_dbl(n, k) * n_choose_k_dbl(tmp , n));		
		coeffs[k] = pow(2.0, n) * coeff;
	}
	reverse_array_dbl(coeffs, order+1);

	return;
}

double compute_bessel_I0( double x )
{
	//horner's method
	double ax, ans, y;

	ax = fabs(x);
	if ( ax < 3.75 ){
		y = x/3.75;
		y *= y;
		ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.36078e-1 + y*0.45813e-2)))));
	}
	else{
		y = 3.75/ax;
		ans = (exp(ax)/sqrt(ax));
		ans *= (0.39894228 + y*(0.132859e-1 + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2
				+ y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 + y*0.392377e-2))))))));
	}

	return ans;
}

double compute_bessel_J0( double x )
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	ax = fabs(x);
	if (ax < 8.0){
		//===Direct Rational Function Fit===//
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718 +y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	}
	else{
		//===Use A Fitting===//
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}

	return ans;
}

double compute_chebyshev_polynomial( int n, 
									 double x )
{
    double res;
    if (fabs(x) <= 1){ 
		res = cos(n*acos(x));
	}
    else {
		res = cosh(n*acosh(x));
	}

    return res;
}

double sign_dbl( double x )
{
	double sign;
	if (x > 0){
		sign = 1;
	}
	if (x < 0){
		sign = -1;
	}
	if (fabs(x) <= 10e-16){
		sign = 0;
	}
	return sign;
}
				
double compute_exponential_smoothing_weight( unsigned int history )
{
	return (1.0 - 1.0/((double)(history)));
}

double ei_dbl(double x)
{
	unsigned int k;
	double sum, fact, prev, term;
	if (x <= 0){
		fprintf(stderr, "Error:: X < 0 In Function -- ei!\n");
		quit();
	}
	if (x <= 1.0e-30){
		return log(x) + 0.57721566;
	}
	if (x <= -log(6e-8)){
		sum = 0.0;
		fact = 1.0;
		for (k=1; k<=100; k++){
			fact *= x/((double)k);
			term = fact/((double)k);
			sum += term;
			if (term < 6e-8 * sum) break;
		}
		return	sum + log(x) + 0.57721566; //EULER	
	}
	else{
		sum = 0.0;
		term = 1.0;
		for (k=1; k<=100; k++){
			prev = term;
			term *= ((double)k)/x;
			if (term < 6e-8) break;
			if (term < prev){
				sum += term; //converging to add new term
			}
			else{
				sum -= prev; //diverging so subtract previous term and exit
				break;	
			}
		}
		return exp(x)*(1.0 + sum)/x;
	}
	return 0;
}

double expint_dbl( int n, double x)
{
	int i, ii, nm1;
	double a,b,c,d,del,fact,h,psi,ans;
	
	nm1 = n-1;
	if (n<0 || x < 0 || ( fabs(x) < 10.0*EPS && (n==0 || n==1))){
		fprintf(stderr, "Error:: Bad Argument! In Function -- expint!\n");
		quit();
	}
	else{
		ans = 0.0;
		if (n==0){
			ans = exp(-x)/x;
		}
		else{
			if (fabs(x) < 10.0*EPS){
				ans = 1.0/((double)nm1);
			}
			else{
				if (x > 1.0){
					b = x + (double)n;
					c = 1.0/(10e-30);
					d = 1.0/b;
					h = d;
					for (i=1; i<=100; i++){
						a = -i * (nm1 + i);
						b += 2.0;
						d = 1.0/(a*d + b);
						c = b + a/c;
						del = c*d;
						h *= del;
						if (fabs(del-1.0) < 1e-7){
							ans = h * exp(-x);
							return ans;
						}
					}
				}
				else{
					ans = (nm1 !=0 ? 1.0/((double)nm1) : -log(x) - 0.5772156649);
					fact = 1.0;
					for (i=1; i<=100; i++){
						fact *= -x/((double)i);
						if (i != nm1){
							del = -fact/((double)(i-nm1));
						}
						else{
							psi = -0.5772156649;
							for (ii=1; ii<=nm1; ii++){
								psi += 1.0/((double)ii);
							}
							del = fact * (-log(x) + psi);
						}
						ans += del;
						if (fabs(del) < fabs(ans)*1e-7){
							return ans;
						}
					}
				}

			}
		}
	}
	
	return ans;
}


//================================================================================================//
//=======================================BIT TWIDDLING============================================//
//================================================================================================//

void pause_program()
{
	fprintf(stdout,"Press ENTER key to Continue.\n");  
	getchar();  
	return;
}

void quit()
{
	exit(1);
	return;
}

unsigned int array_contains_dbl( double* array,
								 unsigned int length,
								 double value )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if (fabs(array[i]-value)<=2.0*EPS){
			return 1;
		}
	}
	return 0;
}

unsigned int array_has_nans_cmplx( double complex* array,
								   unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if ( isnan(creal(array[i])) || isnan(cimag(array[i])) ){
			return 1;
		}
	}
	return 0;
}

unsigned int array_has_nans_dbl( double* array,
								 unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if (isnan(array[i])){
			return 1;
		}
	}
	return 0;
}

unsigned int get_string_length( char* filename )
{
	unsigned int i;
	i = 0;
	while(filename[i] != '\0'){
		i++;
	}
	return i;
}

double compute_linear_gain_constant_from_dB( double dB_gain )
{
	return pow(10.0, dB_gain/10.0);
}

unsigned int next_pow_2( unsigned int length )
{
	return pow(2, ceil(log(length)/log(2)));
}

unsigned int is_real_cmplx( double complex point )
{
	unsigned int real;

	real = 0;
	if (fabs(cimag(point)) <= IS_REAL_THRESHOLD){
		real = 1;
	}
	return real;
}

unsigned int array_is_greater_than_uint( unsigned int* array,
										 unsigned int array_length,
										 unsigned int* compare,
										 unsigned int compare_length )
{
	unsigned int i, greater;

	//===Sanity Check===//
	if (compare_length != array_length && compare_length != 1){
		fprintf(stdout, "Error:: Compare Length Is Invalid! In Function -- array_is_greater_than_uint!\n");
		quit();
	}

	greater = 1;
	if (compare_length == 1){
		for (i=0; i<array_length; i++){
			if (array[i] < *compare){
				greater = 0;
			}	
		}
	}
	else{
		for (i=0; i<array_length; i++){
			if (array[i] < compare[i]){
				greater = 0;
			}
		}
	}

	return greater;
}								  


unsigned int array_is_greater_than_dbl( double* array,
										unsigned int array_length,
										double* compare,
										unsigned int compare_length )
{
	unsigned int i, greater;

	//===Sanity Check===//
	if (compare_length != array_length && compare_length != 1){
		fprintf(stdout, "Error:: Compare Length Is Invalid! In Function -- array_is_greater_than_dbl!\n");
		quit();
	}

	greater = 1;
	if (compare_length == 1){
		for (i=0; i<array_length; i++){
			if (array[i] < *compare){
				greater = 0;
			}	
		}
	}
	else{
		for (i=0; i<array_length; i++){
			if (array[i] < compare[i]){
				greater = 0;
			}
		}
	}

	return greater;
}								  

unsigned int is_array_valid_dbl( double* array,
								 unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if (!isfinite(array[i])){
			return 0;
		}
	}
	return 1;
}

unsigned int array_has_NaNs_dbl( double* array,
								 unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if (array[i]!=array[i]){
			return 1;
		}
	}
	return 0;
}

unsigned int are_equal_dbl( double val1,
						    double val2 )
{
	if ( fabs(val1 - val2) < 10.0*EPS ){
		return 1;
	}
	else{
		return 0;
	}
	return 0;
}

unsigned int are_arrays_equal_dbl( double* array1,
								   double* array2,
								   unsigned int length )
{
	unsigned int i, equal;

	//===Prove Unequality===//
	equal = 1;
	for (i=0; i<length; i++){
		equal = are_equal_dbl(array1[i], array2[i]);
	}

	return equal;
}

unsigned int are_arrays_equal_uint( unsigned int* array1,
								    unsigned int* array2,
								    unsigned int length )
{
	unsigned int i, equal;

	//===Prove Unequality===//
	equal = 1;
	for (i=0; i<length; i++){
		if (array1[i] != array2[i]){
			equal = 0;
		}
	}

	return equal;
}

unsigned int are_arrays_equal_char( char* array1,
								    char* array2,
								    unsigned int length )
{
	unsigned int i, equal;

	//===Prove Unequality===//
	equal = 1;
	for (i=0; i<length; i++){
		if (array1[i] != array2[i]){
			equal = 0;
		}
	}
	return equal;
}


unsigned int find_power_cutoff_index_dbl( double* array,
										  unsigned int length,
										  double threshold )
{
	unsigned int i, index;
	double power, threshold_power;

	//===Get Power===//
	power = compute_power_dbl(array, length);
	threshold = fabs(threshold);
	if (threshold > 1){
		threshold *= 0.01;
	}
	threshold *= power;

	//===Find Index===//
	threshold_power = 0;
	index = length-1;
	for (i=0; i<length; i++){
		threshold_power += fabs(pow(array[i],2.0));
		if (threshold_power >= threshold){
			index = i + 1;
			break;
		}
	}

	return index;
}

void find_unique_dbl( double* array,
					  unsigned int length,
					  double* unique,
					  unsigned int* new_length )
{
	unsigned int i,j, any_found;
	double unique_value;
	
	*new_length = length;
	copy_array_dbl(array, length, unique);
	sort_array_dbl(unique, length);
	for (i=0; i<*new_length; i++){
		
		//===Set Found===//
		any_found = 0;

		//===Get Value===//
		unique_value = unique[i];

		//===Eliminate Unique Values===//
		for (j=0; j<*new_length; j++){
			
			//===Value Is Equal===//
			if (j != i && are_equal_dbl(unique_value,unique[j])){

				//===Set Found===//
				any_found = 1;

				//===Left Shift By 1===//
				left_shift_array_dbl(unique + j, *new_length, 1);
				*new_length -= 1;
			}
		}

		//===Reset Search===//
		if (any_found){
			i = 0;
		}
	}

	return;
}


void find_unique_uint( unsigned int* array,
					   unsigned int length,
					   unsigned int* unique,
					   unsigned int* new_length )
{
	unsigned int i,j,any_found;
	unsigned int unique_value;
	
	*new_length = length;
	copy_array_uint(array, length, unique);
	sort_array_uint(unique, length);
	for (i=0; i<*new_length; i++){

		//===Set Found===//
		any_found = 0;

		//===Get Value===//
		unique_value = unique[i];

		//===Eliminate Unique Values===//
		for (j=0; j<*new_length; j++){

			//===Value Is Equal===//
			if (j != i && unique_value == unique[j]){

				//===Set Found===//
				any_found = 1;

				//===Left Shift By 1===//
				left_shift_array_uint(unique + j, *new_length, 1);
				*new_length -= 1;
			}
		}

		//===Reset Search===//
		if (any_found){
			i = 0;
		}

	}

	return;
}

unsigned int find_closest_index_uint( unsigned int* array,
							   	   	  unsigned int length,
							   	   	  unsigned int value )
{
	unsigned int i, index;
	unsigned int distance, min_distance;
	min_distance = UINT_MAX;
	index = UINT_MAX;

	for (i=0; i<length; i++){
		distance = abs((MAX(array[i], value)) - (MIN(array[i], value)));
		if (distance < min_distance){
			min_distance = distance;
			index = i;
		}
	}

	return index;
}

unsigned int find_closest_index_dbl( double* array,
							   	   	 unsigned int length,
							   	   	 double value )
{
	unsigned int i, index;
	double distance, min_distance;
	min_distance = DBL_MAX;
	index = UINT_MAX;

	for (i=0; i<length; i++){
		distance = fabs(array[i] - value);
		if (distance < min_distance){
			min_distance = distance;
			index = i;
		}
	}

	return index;
}

int find_insertion_index_dbl( double* array,
						  	  unsigned int length,
							  double value )
{
	int index, low_index, mid_index, high_index;
	double mid_value;
	double* copy;

	//===Sort Array===//
	copy = malloc(length * sizeof(double));
	copy_array_dbl(array, length, copy);
	sort_array_dbl(copy, length);

	//===Get Starting Indices===//
	mid_index = length/2;
	mid_value = copy[mid_index];
	if (value > mid_value){
		low_index = mid_index+1;
		high_index = length-1;
	}
	else if (value < mid_value){
		low_index = 0;
		high_index = mid_index-1;
	}
	else if ( fabs(value - mid_value) <= EPS){
		free(copy);
		return mid_index;
	}

	//===Run Binary Search===//
	index = 0;
	while (high_index >= low_index){
		mid_index = (high_index + low_index)/2;

		if (value < copy[mid_index]){
			high_index = mid_index - 1;
		}
		else if (value > copy[mid_index]){
			low_index = mid_index + 1;
		}
		else if ( fabs(value - copy[mid_index]) <= EPS){
			free(copy);
			return mid_index; 
		}
	}	
	index = -((-1 - low_index) + 1);
	
	//===Clean Up===//
	free(copy);

	return index;
}

void find_indices_where_equal_dbl( double* array,
								   unsigned int length,
								   double value,
								   unsigned int* num_indices,
								   unsigned int* indices )
{
	unsigned int i, j;
	double distance;
	j = 0;
	for (i=0; i<length; i++){
		distance = fabs(array[i] - value);
		if (distance < 2*EPS){
			if (indices != NULL){
				indices[j] = i;
			}
			j++;
		}
	}
	*num_indices = j;
	return;
}	

void find_indices_where_equal_uint( unsigned int* array,
								    unsigned int length,
								    unsigned int value,
								    unsigned int* num_indices,
								    unsigned int* indices )
{

	unsigned int i, j;
	j = 0;
	for (i=0; i<length; i++){
		if (array[i] == value){
			if (indices != NULL){
				indices[j] = i;
			}
			j++;
		}
	}
	*num_indices = j;
	return;
}



void find_indices_where_greater_dbl( double* array,
							   	   	 unsigned int length,
							   	   	 double value,
									 unsigned int* num_indices,
									 unsigned int* indices )
{
	unsigned int i, j;
	double distance;
	j = 0;
	for (i=0; i<length; i++){
		distance = fabs(array[i]) - value;
		if (distance > 0){
			if (indices != NULL) indices[j] = i;
			j++;
		}
	}
	if (num_indices != NULL) *num_indices = j;
	return;
}

void find_indices_where_nonzero_dbl( double* array,
									 unsigned int length,
									 unsigned int* num_indices,
									 unsigned int* indices )
{
	unsigned int i, j;
	j = 0;
	for (i=0; i<length; i++){
		if (fabs(array[i]) > 10e-16){
			if (indices != NULL) indices[j] = i;
			j++;
		}
	}
	if (num_indices != NULL) *num_indices = j;

	return;
}

void find_indices_where_less_dbl( double* array,
							   	  unsigned int length,
							   	  double value,
								  unsigned int* num_indices,
								  unsigned int* indices )
{
	unsigned int i, j;
	double distance;
	j = 0;
	for (i=0; i<length; i++){
		distance = fabs(array[i]) - value;
		if (distance < 0){
			if (indices != NULL) indices[j] = i;
			j++;
		}
	}
	if (num_indices != NULL) *num_indices = j;
	return;
}


									  
void insert_array_dbl( double* array,
					   unsigned int* length,
					   double value,
					   unsigned int index )
{
	unsigned int i;
	double* copy;

	//===Copy Array===//
	copy = malloc(*length * sizeof(double));
	copy_array_dbl(array, *length, copy);
	
	//===Shift Array Right at Index===//
	for (i=index; i<*(length); i++){
		array[i+1] = copy[i];
	}

	//===Place Value at Index===//
	array[index] = value; 
	*length = *length + 1;

	free(copy);

	return;
}

void remove_from_array_cmplx( double complex* array,
					   		  unsigned int* length,
					 	      unsigned int index )
{
	if (index == *length-1){
		*length -= 1;
		return;
	}
	else{
		left_shift_array_cmplx(array + index, *length - index, 1);
		*length -= 1;
	}
	return;
}


void remove_from_array_dbl( double* array,
					   		unsigned int* length,
					 	    unsigned int index )
{
	if (index == *length-1){
		*length -= 1;
		return;
	}
	else{
		left_shift_array_dbl(array + index, *length - index, 1);
		*length -= 1;
	}
	return;
}

void remove_from_array_uint( unsigned int* array,
					   		 unsigned int* length,
					 	     unsigned int index )
{
	if (index == *length-1){
		*length -= 1;
		return;
	}
	else{
		left_shift_array_uint(array + index, *length - index, 1);
		*length -= 1;
	}
	return;
}

void sort_array_by_magnitude_cmplx( double complex* array,
								  	unsigned int length )
{

    unsigned int rght, rend, left, k;
    unsigned int i,j,m;
	double complex* temp;
	temp = malloc(length*sizeof(double complex));
    for (k=1; k < length; k *= 2 ) {       
        for (left=0; left+k < length; left += k*2 ) {
            rght = left + k;        
            rend = rght + k;
            if (rend > length) rend = length; 
            m = left; i = left; j = rght; 
            while (i < rght && j < rend) { 
                if (cabs(array[i]) <= cabs(array[j])) {         
                    temp[m] = array[i]; i++;
                } else {
                    temp[m] = array[j]; j++;
                }
                m++;
            }
            while (i < rght) { 
                temp[m]=array[i]; 
                i++; m++;
            }
            while (j < rend) { 
                temp[m]=array[j]; 
                j++; m++;
            }
            for (m=left; m < rend; m++) { 
                array[m] = temp[m]; 
            }
        }
    }
	free(temp);
	return;
}

void sort_array_dbl( double* array, 
					 unsigned int length )
{
    unsigned int rght, rend, left, k;
    unsigned int i,j,m;
	double* temp;
	temp = malloc(length*sizeof(double));
    for (k=1; k < length; k *= 2 ) {       
        for (left=0; left+k < length; left += k*2 ) {
            rght = left + k;        
            rend = rght + k;
            if (rend > length) rend = length; 
            m = left; i = left; j = rght; 
            while (i < rght && j < rend) { 
                if (array[i] <= array[j]) {         
                    temp[m] = array[i]; i++;
                } else {
                    temp[m] = array[j]; j++;
                }
                m++;
            }
            while (i < rght) { 
                temp[m]=array[i]; 
                i++; m++;
            }
            while (j < rend) { 
                temp[m]=array[j]; 
                j++; m++;
            }
            for (m=left; m < rend; m++) { 
                array[m] = temp[m]; 
            }
        }
    }
	free(temp);
	return;
}

void sort_array_uint( unsigned int* array, 
					  unsigned int length )
{
    unsigned int rght, rend, left, k;
    unsigned int i,j,m;
	unsigned int* temp;
	temp = malloc(length*sizeof(unsigned int));
    for (k=1; k < length; k *= 2 ) {       
        for (left=0; left+k < length; left += k*2 ) {
            rght = left + k;        
            rend = rght + k;
            if (rend > length) rend = length; 
            m = left; i = left; j = rght; 
            while (i < rght && j < rend) { 
                if (array[i] <= array[j]) {         
                    temp[m] = array[i]; i++;
                } else {
                    temp[m] = array[j]; j++;
                }
                m++;
            }
            while (i < rght) { 
                temp[m]=array[i]; 
                i++; m++;
            }
            while (j < rend) { 
                temp[m]=array[j]; 
                j++; m++;
            }
            for (m=left; m < rend; m++) { 
                array[m] = temp[m]; 
            }
        }
    }
	free(temp);
	return;
}


void sort_array_indices_by_magnitude_cmplx( double complex* array,
											unsigned int length,
											unsigned int* indices )
{
	double* magnitudes;

	//===Extract Magnitudes===//
	magnitudes = malloc(length * sizeof(double));
	get_magnitude_array_cmplx(array, length, magnitudes);

	//===Sort===//	
	sort_array_indices_dbl(magnitudes, length, indices);

	//===Clean Up===//
	free(magnitudes);

	return;
}

void sort_array_indices_dbl( double* array,
					   		 unsigned int length,
					   		 unsigned int* indices )
{
	unsigned int i, n, M, indxt, ir, swap_temp, j, k, l, m;
	int jstack;
	double a;
	unsigned int* indx;
	int* istack;
	double* arr;

	//===Allocates===//
	indx = calloc(length+1, sizeof(unsigned int));
	arr = calloc(length+1, sizeof(double));
	istack = calloc(length+1, sizeof(int));

	//===Copy Array===//
	for (m=0; m<length; m++){
		arr[m+1] = array[m];
	}

	//===Make Index Array===//
	jstack = 0; n = length; M = length+1;
	l=1; ir=length;
	for (j=1; j<=n; j++){
		indx[j] = j;
	}

	//===Run Sort===//
	for (;;){
		if (ir-l < M){
			for (j=l+1; j<=ir; j++){
				indxt = indx[j];
				a = arr[indxt];
				for (i=j-1; i>=1; i--){
					if (arr[indx[i]] <= a){
						break;
					}
					indx[i+1] = indx[i];
				}
				indx[i+1] = indxt;
			}
			if (jstack == 0){
				break;
			}
			ir = istack[jstack--];
			l = istack[jstack--];
		}
		else{
			k = (l+ir) >> 1;
			SWAP(indx[k], indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]){
				SWAP(indx[l], indx[ir]);
			}
			if (arr[indx[l+1]] > arr[indx[ir]]){
				SWAP(indx[l+1], indx[ir]);
			}
			if (arr[indx[l]] > arr[indx[l+1]]){
				SWAP(indx[l], indx[l+1]);
			}
			i = l + 1;
			j = ir;
			indxt = indx[l+1];
			a = arr[indxt];
			for (;;){
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i], indx[j]);
			}
			indx[l+1] = indx[j];
			indx[j] = indxt;
			jstack += 2;
			if (ir-i+1 > j-1){
				istack[jstack] = ir;
				istack[jstack-1] = i;
				ir = j-1;
			}
			else{
				istack[jstack] = j-1;
				istack[jstack-1] = l;
				l=i;
			}
		}
	}

	//===Copy Over===//
	for (m=0; m<length; m++){
		indices[m] = indx[m+1] - 1;
	}

	//===Free===//
	free(istack);
	free(arr);
	free(indx);	
	
	return;
}

//REALLY SHOULD BE DECIMATE ARRAY
void downsample_array_dbl( double* signal,
				   		   unsigned int N,
				   		   unsigned int order )
{
	unsigned int i, M;
	double signal_temp[4*(MAX_POLYNOMIAL_ORDER+1)];

	if (N > 4*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Input Array Is Too Large! In Function -- downsample_array_dbl!\n");
		quit();
	} 

	//===Set Length Of New Signal===//
	if (N % order == 0){
		M = N/order;
	}
	else{
		M = N/order + 1;
	}

	//===Perform Downsampling===//
	copy_array_dbl(signal, N, signal_temp);
	initialize_array_dbl(signal, N);
	for (i=0; i<M; i++){
		signal[i] = signal_temp[order*i];
	}

	return;
}	

void downsample_array_uint( unsigned int* signal,
				   		    unsigned int N,
				   		    unsigned int order )
{
	unsigned int i, M;
	unsigned int signal_temp[4*(MAX_POLYNOMIAL_ORDER+1)];

	if (N > 4*(MAX_POLYNOMIAL_ORDER+1)){
		fprintf(stderr, "Error:: Input Array Is Too Large! In Function -- downsample_array_uint!\n");
		quit();
	} 

	//===Set Length Of New Signal===//
	if (N % order == 0){
		M = N/order;
	}
	else{
		M = N/order + 1;
	}

	//===Perform Downsampling===//
	copy_array_uint(signal, N, signal_temp);
	initialize_array_uint(signal, N);
	for (i=0; i<M; i++){
		signal[i] = signal_temp[order*i];
	}

	return;
}	


void upsample_array_dbl( double* signal,
			     		 unsigned int N,
				 		 unsigned int order )
{

	unsigned int i, M;
	double signal_temp[MAX_POLYNOMIAL_ORDER+1];

	if (N > MAX_POLYNOMIAL_ORDER+1){
		fprintf(stderr, "Error:: Input Array Is Too Large! In Function -- upsample_array_dbl!\n");
		quit();
	} 

	//===Perform Upsampling===//
	M = order*N;
	copy_array_dbl(signal, N, signal_temp);
	initialize_array_dbl(signal, M);
	for (i=0; i<N; i++){
		signal[order*i] = signal_temp[i];
	}

	return;
}


//================================================================================================//
//=====================================ARRAY MANIPULATION=========================================//
//================================================================================================//

void initialize_array_cmplx( double complex* array,
					   	   	 unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = 0;
	}

	return;
}

void initialize_array_dbl( double* array,
					   	   unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = 0;
	}

	return;
}

void initialize_array_int( int* array,
					   	   unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = 0;
	}

	return;
}


void initialize_array_uint( unsigned int* array,
					   	   	unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = 0;
	}

	return;
}

void initialize_array_char( char* array,
					   	   	unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = '\0';
	}

	return;
}


void initialize_array_uchar( unsigned char* array,
					   	   	 unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = 0;
	}

	return;
}

void initialize_array_constant_cmplx( double complex* array,
					   	   			  unsigned int length,
									  double complex constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = constant;
	}

	return;
}


void initialize_array_constant_dbl( double* array,
					   	   			unsigned int length,
									double constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = constant;
	}

	return;
}

void initialize_array_constant_uint( unsigned int* array,
					   	   			 unsigned int length,
									 unsigned int constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = constant;
	}

	return;
}


void copy_array_char( char* data,
				 	  unsigned int length,
				 	  char* copy )
{
	unsigned int i;
	for (i=0; i<length; i++){
		copy[i] = data[i];
	}

	return;
}

void copy_array_int( int* data,
				 	 unsigned int length,
				 	 int* copy )
{
	unsigned int i;
	for (i=0; i<length; i++){
		copy[i] = data[i];
	}
	return;
}

void copy_array_uint( unsigned int* data,
				 	  unsigned int length,
				 	  unsigned int* copy )
{
	unsigned int i;
	for (i=0; i<length; i++){
		copy[i] = data[i];
	}

	return;
}

void copy_array_dbl( double* data,
				 	 unsigned int length,
				 	 double* copy )
{
	unsigned int i;
	for (i=0; i<length; i++){
		copy[i] = data[i];
	}

	return;
}

void copy_array_cmplx( double complex* data,
				 	   unsigned int length,
				 	   double complex* copy )
{
	unsigned int i;
	for (i=0; i<length; i++){
		copy[i] = data[i];
	}

	return;
}

void combine_arrays_dbl_to_cmplx( double* real,
								  double* imag,
								  unsigned int length,
								  double complex* output )
{
	unsigned int i;

	if (real != NULL && imag != NULL){
		for (i=0; i<length; i++){
			output[i] = real[i] + I*imag[i];
		}
	}
	if (real == NULL){
		for (i=0; i<length; i++){
			output[i] = 0.0 + I*imag[i];
		}
	}
	if (imag == NULL){
		for (i=0; i<length; i++){
			output[i] = real[i] + I*0.0;
		}
	}

	return;
}

void split_array_cmplx( double complex* array,
				 	    unsigned int length,
				 	    double* real,
						double* imag )
{
	unsigned int i;

	if (real != NULL){
		for (i=0; i<length; i++){
			real[i] = creal(array[i]);
		}
	}

	if (imag != NULL){
		for (i=0; i<length; i++){
			imag[i] = cimag(array[i]);
		}
	}

	return;
}

void swap_real_imag_cmplx( double complex* array,
						   unsigned int length )
{
	unsigned int i;
	double a, b;
	double complex swap_temp;

	for (i=0; i<length; i++){
		a = creal(array[i]); b = cimag(array[i]);
		SWAP(a,b);
		array[i] = a + I*b;
	}

	return;
}

void zero_array_real_cmplx( double complex* array,
				 	        unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = 0.0 + cimag(array[i])*I;
	}

	return;
}

void zero_array_imag_cmplx( double complex* array,
				 	        unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = creal(array[i]) + 0.0*I;
	}

	return;
}


void abs_imaginary_cmplx( double complex* array,
						  unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = creal(array[i]) + I*fabs(cimag(array[i]));
	}
	
	return;
}

void abs_array_cmplx( double complex* array,
				      unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = cabs(array[i]);
	}
	
	return;
}


void abs_array_dbl( double* array,
				    unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = fabs(array[i]);
	}
	
	return;
}

void pow_array_cmplx( double complex* array,
					  unsigned int length,
				      double complex power )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = cpow(array[i], power);
	}

	return;
}


void pow_array_dbl( double* array,
					unsigned int length,
				    double power )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = pow(array[i],power);
	}

	return;
}

void copy_matrix_cmplx( double complex* matrix,
				  		unsigned int num_rows,
				  		unsigned int num_cols,
				  		double complex* copy )
{

	unsigned int i, j;
	for (i=0; i<num_rows; i++){
		for (j=0; j<num_cols; j++){
			copy[j +i*num_cols] = matrix[j +i*num_cols];
		}
	}


	return;
}

void copy_matrix_dbl( double* matrix,
				  	  unsigned int num_rows,
				  	  unsigned int num_cols,
				  	  double* copy )
{

	unsigned int i, j;
	for (i=0; i<num_rows; i++){
		for (j=0; j<num_cols; j++){
			copy[j +i*num_cols] = matrix[j +i*num_cols];
		}
	}


	return;
}

void copy_matrix_uint( unsigned int* matrix,
				  	   unsigned int num_rows,
				  	   unsigned int num_cols,
				  	   unsigned int* copy )
{

	unsigned int i, j;
	for (i=0; i<num_rows; i++){
		for (j=0; j<num_cols; j++){
			copy[j +i*num_cols] = matrix[j +i*num_cols];
		}
	}


	return;
}


void normalize_angle_cmplx( double complex* data,
							unsigned int length )
{
	unsigned int i;
	double magnitude, max_angle, angle;
	double complex phase;

	max_angle = find_maximum_angle_cmplx(data, length);
	for (i=0; i<length; i++){
		magnitude = cabs(data[i]);
		angle = atan2(cimag(data[i]), creal(data[i])) / max_angle;
		phase = cexp(1.0*I*angle);
		data[i] = magnitude * phase;
	}

	return;
}

void normalize_vector_2norm_cmplx( double complex* vector,
					  			   unsigned int length )
{
	unsigned int i;
	double norm;
	norm = compute_vector_2norm_cmplx(vector, length);
	if (fabs(norm) > 10e-300){
		for (i=0; i<length; i++){
			vector[i] *= (1.0/norm);
		}
	}		
	return;
}

void normalize_vector_2norm_dbl( double* vector,
					  			 unsigned int length )
{
	unsigned int i;
	double norm;
	norm = compute_vector_2norm_dbl(vector, length);
	if (fabs(norm) > 10e-300){
		for (i=0; i<length; i++){
			vector[i] *= (1.0/norm);
		}
	}
		
	return;
}

void normalize_max_dbl( double* data,
						unsigned int length )
{
	unsigned int i;
	double maximum;

	maximum = find_maximum_dbl(data, length);
	for (i=0; i<length; i++){
		data[i] /= maximum;
	}

	return;
} 

void normalize_range_dbl( double* data,
						  unsigned int length )
{
	unsigned int i;
	double max, min;
	
	//===Normalize To [0,1]===//
	max = find_maximum_dbl(data, length);
	min = find_minimum_dbl(data, length);
	for (i=0; i<length; i++){
		data[i] = (data[i]-min)/(max-min);
	}

	return;
}
	

void normalize_abs_max_dbl( double* data,
							unsigned int length )
{
	unsigned int i;
	double maximum;

	maximum = find_abs_maximum_dbl(data, length);
	for (i=0; i<length; i++){
		data[i] /= maximum;
	}

	return;
} 


void normalize_magnitude_cmplx( double complex* data,
								unsigned int length )
{
	unsigned int i;
	double magnitude;

	magnitude = find_maximum_magnitude_cmplx(data,length);
	for (i=0; i<length; i++){
		data[i] = ((1.0/magnitude) * cabs(data[i])) * cexp(1.0 * I * atan2(cimag(data[i]),creal(data[i])) );
	}

	return;
}

void normalize_variance_dbl( double* data,
							 unsigned int length )
{
	unsigned int i;
	double std;

	std = sqrt(compute_variance_dbl(data, length));
	for (i=0; i<length; i++){
		data[i] /= std;
	}	

	return;
}


void normalize_sum_dbl( double* data,
						unsigned int length )
{
	unsigned int i;
	double sum;

	sum = compute_sum_dbl(data,length);
	for (i=0; i<length; i++){
		data[i] /= sum;
	}	

	return;
}

void normalize_power_dbl( double* data,
						  unsigned int length )
{
	unsigned int i;
	double power;

	power = compute_power_dbl(data,length);
	for (i=0; i<length; i++){
		data[i] /= sqrt(power);
	}	

	return;
}


void make_array_monic_dbl( double* array, 
					   	   unsigned int length )
{
	//===Locals===//
	unsigned int i;
	double lead;

	//===Monicize===//
	lead = array[0];
	for (i=0; i<length; i++){
		array[i] /= lead;
	}
	
	return;
}

void make_array_monic_cmplx( double complex* array, 
					   	   	 unsigned int length )
{
	//===Locals===//
	unsigned int i;
	double complex lead;

	//===Monicize===//
	lead = array[0];
	for (i=0; i<length; i++){
		array[i] /= lead;
	}
	
	return;

}

void reverse_array_uint( unsigned int* array, 
						 unsigned int length )
{
	//===Locals===//
	unsigned int c, end;
	unsigned int t;

	//===Reverse Array===//
	end = length - 1; 
	for (c=0; c<length/2; c++){
		t = array[c];
		array[c] = array[end];
		array[end] = t;
		end--;
	}

	return;
}


void reverse_array_dbl( double* array, 
						unsigned int length )
{
	//===Locals===//
	unsigned int c, end;
	double t;

	//===Reverse Array===//
	end = length - 1; 
	for (c=0; c<length/2; c++){
		t = array[c];
		array[c] = array[end];
		array[end] = t;
		end--;
	}

	return;
}

void reverse_array_cmplx( double complex* array, 
						  unsigned int length )
{
	//===Locals===//
	unsigned int c, end;
	double complex t;

	//===Reverse Array===//
	end = length - 1; 
	for (c=0; c<length/2; c++){
		t = array[c];
		array[c] = array[end];
		array[end] = t;
		end--;
	}

	return;
}

void reorder_array_uint( unsigned int* array,
						 unsigned int* new_indices,
						 unsigned int length )
{
	unsigned int i;
	unsigned int* copy;

	//===Malloc===//
	copy = malloc(length*sizeof(unsigned int));

	//===Copy===//
	for (i=0; i<length; i++){
		copy[i] = array[new_indices[i]];
	}
	copy_array_uint(copy, length, array);

	//===Clean Up===//
	free(copy);

	return;
}

void reorder_array_dbl( double* array,
						unsigned int* new_indices,
						unsigned int length )
{
	unsigned int i;
	double* copy;

	//===Malloc===//
	copy = malloc(length*sizeof(double));

	//===Copy===//
	for (i=0; i<length; i++){
		copy[i] = array[new_indices[i]];
	}
	copy_array_dbl(copy, length, array);

	//===Clean Up===//
	free(copy);

	return;
}

void reorder_array_cmplx( double complex* array,
						  unsigned int* new_indices,
						  unsigned int length )
{
	unsigned int i;
	double complex* copy;

	//===Malloc===//
	copy = malloc(length*sizeof(double complex));

	//===Copy===//
	for (i=0; i<length; i++){
		copy[i] = array[new_indices[i]];
	}
	copy_array_cmplx(copy, length, array);

	//===Clean Up===//
	free(copy);

	return;
}


void extract_array_coordinates_cmplx( double complex* array,
									  unsigned int length,
									  double* real,
									  double* imag )
{
	unsigned int i;	
	for (i=0; i<length; i++){
		if (real != NULL){
			real[i] = creal(array[i]);
		}
		if (imag != NULL){
			imag[i] = cimag(array[i]);
		}
	}

	return;
}

void initialize_linspace_dbl( double* linspace,
							  unsigned int length,
							  double start,
							  double end )
{
	unsigned int s;
	double step_size;
	step_size = (end - start)/(length-1);
	for (s=0; s<length; s++){	
		linspace[s] = start + s*step_size; 
	}
	//fprintf(stdout, "Last Linspace: %lf	End: %lf\n", linspace[length-1], end);
	return;
}

void initialize_linspace_uint( unsigned int* linspace,
							   unsigned int length,
							   unsigned int start,
							   unsigned int end )
{
	unsigned int s;
	double step_size;
	step_size = (end - start)/(length-1);
	for (s=0; s<length; s++){	
		linspace[s] = start + (unsigned int)(s*step_size);
	}
	return;
}

void initialize_linspace_int( int* linspace,
							  int length,
							  int start,
							  int end )
{
	int s;
	double step_size;
	step_size = (end - start)/(length-1);
	for (s=0; s<length; s++){	
		linspace[s] = start + (int)(s*step_size);
	}
	return;
}


void initialize_logspace_dbl( double* logspace,
							  unsigned int length,
							  double start,
							  double end )
{
	unsigned int s;
	double step_size;
	step_size = 1.0/((double)length-1.0);
	for (s=0; s<length; s++){	
		logspace[s] = pow(end, ((double)s)*step_size) * pow(start, 1.0-(((double)s)*step_size));
	}
	return;
}


void initialize_range_uint( unsigned int* range,
							unsigned int start,
							unsigned int end )
{
	unsigned int s;
	for (s=0; s<end-start; s++){	
		range[s] = start + s;
	}
	return;
}

void initialize_range_dbl( double* array,
						   double start,
						   double end,
						   double step,
						   unsigned int* length )
{
	unsigned int new_length;
	double new_x;

	array[0] = start;
	new_x = array[0];
	new_length = 1;
	while (1){
		new_x = array[new_length-1] + step;
		if (new_x > (end + EPS)){
			break;
		}
		else{
			array[new_length] = new_x;
			new_length += 1;
		}
	}
	if (length != NULL){
		*length = new_length;
	}
	return;
}   

void pad_zeros_cmplx( double complex* array,
					  unsigned int length,
					  unsigned int num_zeros )
{
	unsigned int i;
	for (i=length; i<num_zeros+length; i++){
		array[i] = 0;
	}

}


void pad_zeros_dbl( double* array,
					unsigned int length,
					unsigned int num_zeros )
{
	unsigned int i;
	for (i=length; i<num_zeros+length; i++){
		array[i] = 0;
	}

}

void gain_array_constant_cmplx( double complex* array,
							  	unsigned int length,
							  	double complex constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] *= constant;
	}
	return;
}		


void gain_array_constant_dbl( double* array,
							  unsigned int length,
							  double constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] *= constant;
	}
	return;
}		

void gain_array_constant_uint( unsigned int* array,
							   unsigned int length,
							   unsigned int constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] *= constant;
	}
	return;
}		


void add_array_constant_uint( unsigned int * array,
							  unsigned int length,
							  unsigned int constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] += constant;
	}
	return;
}

void add_array_constant_int( int * array,
							 int length,
							 int constant )
{
	int i;
	for (i=0; i<length; i++){
		array[i] += constant;
	}
	return;
}

void add_array_constant_dbl( double* array,
							 unsigned int length,
							 double constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] += constant;
	}
	return;
}


void scale_array_dbl( double* array,
					  unsigned int length,
					  double start,
					  double end )
{
	unsigned int i;
	double minn, maxx;
	double* copy;

	//===Mallocs===//
	copy = malloc(length*sizeof(double));
	copy_array_dbl(array, length, copy);

	//===Scale To [0, max+min]===//
	minn = find_minimum_dbl(copy, length);
	if (minn < 0) add_array_constant_dbl(copy, length, fabs(minn));
	else add_array_constant_dbl(copy, length, -1.0*minn);

	//===Normalize to [0, 1]===//
	maxx = fabs(find_maximum_dbl(copy, length));
	gain_array_constant_dbl(copy, length, 1.0/(maxx));

	//===Scale to [start, end]===//
	for (i=0; i<length; i++){
		array[i] = start + copy[i]*(end-start);
	}
	
	//===Clean Up===//
	free(copy);

	return;
}


void append_array_unique_uint( unsigned int* array,
							   unsigned int* array_length,
							   unsigned int value )
{
	unsigned int i, inside;

	//===See If Already Inside===//
	inside = 0;
	for (i=0; i<*array_length; i++){
		if (array[i] == value){
			inside = 1;
		}
	}

	//===Append===//
	if (!inside){
		array[*array_length] = value;
		*array_length = *array_length + 1;		
	}
	
	return;
}


void append_array_dbl( double* array,
					   unsigned int* array_length,
					   double* appendix,
					   unsigned int appendix_length )
{
	unsigned int i, j;
	j = 0;
	for (i=*array_length; i<*array_length+appendix_length; i++){
		array[i] = appendix[j++];
	}
	*array_length = *array_length + appendix_length;
	return;
}


void right_circular_shift_array_dbl( double* array,
									 unsigned int array_length,
									 unsigned int shift_length )
{
	double* saved; 
	saved = malloc(shift_length * sizeof(double));
	initialize_array_dbl(saved, shift_length);

	//===Save Last Shift Length Values===//
	copy_array_dbl(array + array_length - shift_length, shift_length, saved);

	//===Shift Right By Shift Length===//
	right_shift_array_dbl(array, array_length, shift_length);

	//===Place Saved Values Up Front===//
	copy_array_dbl(saved, shift_length, array);

	//===Clean Up===//
	free(saved);

	return;
}

void left_circular_shift_array_cmplx( double complex* array,
									  unsigned int array_length,
									  unsigned int shift_length )
{
	double complex* saved; 
	saved = malloc(shift_length * sizeof(*saved));
	initialize_array_cmplx(saved, shift_length);

	//===Save Last Shift Length Values===//
	copy_array_cmplx(array, shift_length, saved);

	//===Shift Left By Shift Length===//
	left_shift_array_cmplx(array, array_length, shift_length);

	//===Place Saved Values Up Front===//
	copy_array_cmplx(saved, shift_length, array + array_length - shift_length);

	//===Clean Up===//
	free(saved);

	return;
}								

void left_circular_shift_array_dbl( double* array,
									unsigned int array_length,
									unsigned int shift_length )
{
	double* saved; 
	saved = malloc(shift_length * sizeof(*saved));
	initialize_array_dbl(saved, shift_length);

	//===Save Last Shift Length Values===//
	copy_array_dbl(array, shift_length, saved);

	//===Shift Left By Shift Length===//
	left_shift_array_dbl(array, array_length, shift_length);

	//===Place Saved Values Up Front===//
	copy_array_dbl(saved, shift_length, array + array_length - shift_length);

	//===Clean Up===//
	free(saved);

	return;
}				

void left_circular_shift_array_uint( unsigned int* array,
									 unsigned int array_length,
									 unsigned int shift_length )
{
	unsigned int* saved; 
	saved = malloc(shift_length * sizeof(*saved));
	initialize_array_uint(saved, shift_length);

	//===Save Last Shift Length Values===//
	copy_array_uint(array, shift_length, saved);

	//===Shift Left By Shift Length===//
	left_shift_array_uint(array, array_length, shift_length);

	//===Place Saved Values Up Front===//
	copy_array_uint(saved, shift_length, array + array_length - shift_length);

	//===Clean Up===//
	free(saved);

	return;
}				

void right_shift_array_uint( unsigned int* array,
							 unsigned int array_length,
							 unsigned int delay_length )
{
	unsigned int i;
	unsigned int* copy;

	//===Delay===//
	if (delay_length >= array_length){
		initialize_array_uint(array, array_length);
	}
	else{
		copy = malloc((array_length+delay_length)*sizeof(unsigned int));
		initialize_array_uint(copy, array_length+delay_length);
		for (i=0; i<array_length; i++){
			copy[i+delay_length] = array[i];
		}
		copy_array_uint(copy, array_length, array);
		free(copy);
	}
	return;
}


void right_shift_array_dbl( double* array,
							unsigned int array_length,
							unsigned int delay_length )
{
	unsigned int i;
	double* copy;

	//===Delay===//
	if (delay_length >= array_length){
		initialize_array_dbl(array, array_length);
	}
	else{
		copy = malloc((array_length+delay_length)*sizeof(double));
		initialize_array_dbl(copy, array_length+delay_length);
		for (i=0; i<array_length; i++){
			copy[i+delay_length] = array[i];
		}
		copy_array_dbl(copy, array_length, array);
		free(copy);
	}
	return;
}

void left_shift_array_uint( unsigned int* array,
						    unsigned int array_length,
						    unsigned int shift_length )
{
	unsigned int i;
	unsigned int* copy;

	//===Shift===//
	if (shift_length >= array_length){
		initialize_array_uint(array, array_length);
	}
	else{
		copy = malloc((array_length+shift_length)*sizeof(unsigned int));
		initialize_array_uint(copy, array_length+shift_length);
		for (i=0; i<array_length-shift_length; i++){
			copy[i] = array[i+shift_length];
		}
		copy_array_uint(copy, array_length, array);
		free(copy);
	}
	return;
}


void left_shift_array_dbl( double* array,
						   unsigned int array_length,
						   unsigned int shift_length )
{
	unsigned int i;
	double* copy;

	//===Shift===//
	if (shift_length >= array_length){
		initialize_array_dbl(array, array_length);
	}
	else{
		copy = malloc((array_length+shift_length)*sizeof(double));
		initialize_array_dbl(copy, array_length+shift_length);
		for (i=0; i<array_length-shift_length; i++){
			copy[i] = array[i+shift_length];
		}
		copy_array_dbl(copy, array_length, array);
		free(copy);
	}
	return;
}

void left_shift_array_cmplx( double complex* array,
						   	 unsigned int array_length,
						   	 unsigned int shift_length )
{
	unsigned int i;
	double complex* copy;

	//===Shift===//
	if (shift_length >= array_length){
		initialize_array_cmplx(array, array_length);
	}
	else{
		copy = malloc((array_length+shift_length)*sizeof(double complex));
		initialize_array_cmplx(copy, array_length+shift_length);
		for (i=0; i<array_length-shift_length; i++){
			copy[i] = array[i+shift_length];
		}
		copy_array_cmplx(copy, array_length, array);
		free(copy);
	}
	return;
}

void remove_array_element_uint( unsigned int* array,
								unsigned int* array_length,
								unsigned int element )
{
	unsigned int i, num_indices, num_removed;
	unsigned int indices[MAX_POLYNOMIAL_ORDER+1];
	
	//===Find All Indicies Where Element Exists===//
	num_indices = 0;
	for (i=0; i<*array_length; i++){
		if (array[i] == element){
			indices[num_indices++] = i;
		}
	}

	num_removed = 0;
	while(num_removed != num_indices){

		//===Left Shift Over This Index===//
		left_shift_array_uint(array+indices[num_removed],*array_length-indices[num_removed],1);

		//===Fix Length===//
		*array_length = *array_length - 1;

		//===Update Removed===//
		num_removed++;
	}

	return;
}

void repeat_array_dbl( double* array,
					   unsigned int length,
					   unsigned int num_repeats )
{
	unsigned int i;
	double *copy;

	//===Malloc===//
	copy = malloc(length*sizeof(double));

	//===Repeat===//
	copy_array_dbl(array, length, copy);
	for (i=0; i<num_repeats; i++){
		copy_array_dbl(copy, length, array + i*length);
	}
	
	//===Free===//
	free(copy);

	return;
}	

void left_shift_matrix_columns_cmplx( double complex* matrix,
									  unsigned int num_rows,
									  unsigned int num_cols,
									  unsigned int num_shifts )
{
	unsigned int i, n;
	double complex column[MAX_POLYNOMIAL_ORDER+1];	

	for (n=0; n<num_shifts; n++){
		for (i=1; i<num_cols-1; i++){

			//===Copy Current Column===//
			copy_matrix_column_cmplx(matrix, num_rows, num_cols, i, column);

			//===Replace Previous Column===//		
			replace_matrix_column_cmplx(matrix, num_rows, num_cols, column, i-1);
		}

		//===Replace Final Column With Zeros===//
		initialize_array_cmplx(column, num_rows);
		replace_matrix_column_cmplx(matrix, num_rows, num_cols, column, i);
	}
	
	return;
}

void remove_matrix_column_cmplx( double complex* matrix,
								 unsigned int num_rows,
								 unsigned int* num_cols,
								 unsigned int column_number )
{
	unsigned int r;
	double complex row[MAX_POLYNOMIAL_ORDER+1];
	double complex* copy;
	
	//===Copy Matrix===//
	copy = malloc(num_rows * (*num_cols) * sizeof(double complex));
	copy_array_cmplx(matrix, num_rows*(*num_cols), copy);	

	//===Zero Original===//
	initialize_array_cmplx(matrix, num_rows*(*num_cols));

	//===Remove Column===//	
	for (r=0; r<num_rows; r++){

		//===Copy Row===//
		copy_matrix_row_cmplx(copy, *num_cols, r, row);
		
		//===Replace Column Index With Zero===//
		row[column_number] = 0.0;

		//===Left Shift Zero Out===//
		left_shift_array_cmplx(row + column_number, *num_cols-column_number, 1);

		//===Place Row Back In===//
		replace_matrix_row_cmplx(matrix, *num_cols-1, row, r);
	}
	
	//===Update Number Of Columns===//
	*num_cols -= 1;

	//===Clean Up===//
	free(copy);


	return;
}

void remove_matrix_column_dbl( double* matrix,
							   unsigned int num_rows,
							   unsigned int* num_cols,
							   unsigned int column_number )
{
	unsigned int r;
	double row[MAX_POLYNOMIAL_ORDER+1];
	double* copy;
	
	//===Copy Matrix===//
	copy = malloc(num_rows * (*num_cols) * sizeof(double));
	copy_array_dbl(matrix, num_rows*(*num_cols), copy);	

	//===Zero Original===//
	initialize_array_dbl(matrix, num_rows*(*num_cols));

	//===Remove Column===//	
	for (r=0; r<num_rows; r++){

		//===Copy Row===//
		copy_matrix_row_dbl(copy, *num_cols, r, row);
		
		//===Replace Column Index With Zero===//
		row[column_number] = 0.0;

		//===Left Shift Zero Out===//
		left_shift_array_dbl(row + column_number, *num_cols-column_number, 1);

		//===Place Row Back In===//
		replace_matrix_row_dbl(matrix, *num_cols-1, row, r);
	}
	
	//===Update Number Of Columns===//
	*num_cols -= 1;

	//===Clean Up===//
	free(copy);

	return;
}

void remove_matrix_row_dbl( double* matrix,
							unsigned int* num_rows,
							unsigned int num_cols,
							unsigned int row_number )
{
	unsigned int c;
	double column[MAX_POLYNOMIAL_ORDER+1];
	double* copy;
	
	//===Copy Matrix===//
	copy = malloc(*(num_rows) * num_cols * sizeof(double));
	copy_array_dbl(matrix, (*num_rows)*num_cols, copy);	

	//===Zero Original===//
	initialize_array_dbl(matrix, (*num_rows)*num_cols);

	//===Remove Row===//	
	for (c=0; c<num_cols; c++){

		//===Copy Column===//
		copy_matrix_column_dbl(copy, *num_rows, num_cols, c, column);
		
		//===Replace Row Index With Zero===//
		column[row_number] = 0.0;

		//===Left Shift Zero Out===//
		left_shift_array_dbl(column + row_number, num_cols-row_number, 1);

		//===Place Column Back In===//
		replace_matrix_column_dbl(matrix, *num_rows, num_cols, column, c);
	}
	
	//===Update Number Of Rows===//
	*num_rows -= 1;

	//===Clean Up===//
	free(copy);

	return;
}



void right_circular_shift_array_cmplx( double complex* array,
									   unsigned int array_length,
									   unsigned int shift_length )
{
	double complex* saved; 
	saved = malloc(shift_length * sizeof(*saved));
	initialize_array_cmplx(saved, shift_length);

	//===Save Last Shift Length Values===//
	copy_array_cmplx(array + array_length - shift_length, shift_length, saved);

	//===Shift Right By Shift Length===//
	right_shift_array_cmplx(array, array_length, shift_length);

	//===Place Saved Values Up Front===//
	copy_array_cmplx(saved, shift_length, array);

	//===Clean Up===//
	free(saved);

	return;
}				

void right_circular_shift_array_uint( unsigned int* array,
									  unsigned int array_length,
									  unsigned int shift_length )
{
	unsigned int* saved; 
	saved = malloc(shift_length * sizeof(*saved));
	initialize_array_uint(saved, shift_length);

	//===Save Last Shift Length Values===//
	copy_array_uint(array + array_length - shift_length, shift_length, saved);

	//===Shift Right By Shift Length===//
	right_shift_array_uint(array, array_length, shift_length);

	//===Place Saved Values Up Front===//
	copy_array_uint(saved, shift_length, array);

	//===Clean Up===//
	free(saved);

	return;
}				


void right_shift_array_cmplx( double complex* array,
							  unsigned int array_length,
							  unsigned int delay_length )
{
	unsigned int i;
	double complex* copy;

	//===Delay===//
	if (delay_length >= array_length){
		initialize_array_cmplx(array, array_length);
	}
	else{
		copy = malloc((array_length+delay_length)*sizeof(*copy));
		initialize_array_cmplx(copy, array_length+delay_length);
		for (i=0; i<array_length; i++){
			copy[i+delay_length] = array[i];
		}
		copy_array_cmplx(copy, array_length, array);
		free(copy);
	}
	return;
}

void eliminate_array_zeros_dbl( double* array,
								unsigned int array_length,
								unsigned int* new_length )
{
	unsigned int i, j;
	j = array_length;
	for (i=0; i<j; i++){
		if (fabs(array[i]) < 10e-16){
			left_shift_array_dbl(array + i, array_length-i, 1);
			j--; i--;
		}
	}
	if (new_length != NULL){
		*new_length = j;
	}
	return;
}


void eliminate_array_zeros_uint( unsigned int* array,
								 unsigned int array_length,
								 unsigned int* new_length  )
{
	unsigned int i, j;
	j = array_length;
	for (i=0; i<j; i++){
		if (array[i] == 0){
			left_shift_array_uint(array + i, array_length-i, 1);
			j--; i--;
		}
	}
	if (new_length != NULL){
		*new_length = j;
	}
	return;
}

void eliminate_array_infs_dbl( double* array,
							   unsigned int array_length )
{
	unsigned int i, j;
	j = array_length;
	for (i=0; i<j; i++){
		if (isinf(array[i])){
			left_shift_array_dbl(array + i, array_length-i, 1);
			j--; i--;
		}
	}
	return;
}

void append_matrix_column_dbl( double* matrix,
							   unsigned int num_rows,
							   unsigned int* num_cols,
							   double* column )
{
	unsigned int c, num_old_cols, num_new_cols;
	double temp_column[MAX_POLYNOMIAL_ORDER+1];
	double* matrix_temp;

	//===Malloc===//	
	matrix_temp = malloc(num_rows*((*num_cols)+1)*sizeof(double));

	//===Save===//	
	num_old_cols = *num_cols;  num_new_cols = num_old_cols + 1;
	copy_matrix_dbl(matrix, num_rows, num_old_cols, matrix_temp);

	//===Copy Over===//	
	initialize_matrix_dbl(matrix, num_rows, num_old_cols);
	initialize_matrix_dbl(matrix, num_rows, num_new_cols);
	for (c=0; c<num_old_cols; c++){
		copy_matrix_column_dbl(matrix_temp, num_rows, num_old_cols, c, temp_column);
		replace_matrix_column_dbl(matrix, num_rows, num_new_cols, temp_column, c);
	}
	replace_matrix_column_dbl(matrix, num_rows, num_new_cols, column, c);
	*num_cols = num_new_cols;
	

	//===Clean Up===//
	free(matrix_temp);

	return;
}

void insert_matrix_column_dbl( double* matrix,
							   unsigned int num_rows,
							   unsigned int* num_cols,
							   double* column, 
							   unsigned int column_number )
{

	unsigned int c, num_old_cols, num_new_cols;
	double temp_column[MAX_POLYNOMIAL_ORDER+1];
	double* matrix_temp;

	//===Malloc===//	
	matrix_temp = malloc(num_rows*((*num_cols))*sizeof(double));

	//===Save===//	
	num_old_cols = *num_cols;  num_new_cols = num_old_cols + 1;
	copy_matrix_dbl(matrix, num_rows, num_old_cols, matrix_temp);

	//===Copy Over===//	
	initialize_matrix_dbl(matrix, num_rows, num_old_cols);
	initialize_matrix_dbl(matrix, num_rows, num_new_cols);
	for (c=0; c<num_new_cols; c++){
		if (c == column_number){
			replace_matrix_column_dbl(matrix, num_rows, num_new_cols, column, c);
		}
		else{
			if (c < column_number){
				copy_matrix_column_dbl(matrix_temp, num_rows, num_old_cols, c, temp_column);
				replace_matrix_column_dbl(matrix, num_rows, num_new_cols, temp_column, c);
			}
			else{
				copy_matrix_column_dbl(matrix_temp, num_rows, num_old_cols, c-1, temp_column);
				replace_matrix_column_dbl(matrix, num_rows, num_new_cols, temp_column, c);
			}
		}
	}
	*num_cols = num_new_cols;
	

	return;
}


//================================================================================================//
//=====================================MAXIMUM AND MINIMUM========================================//
//================================================================================================//

double find_maximum_pairwise_difference_dbl( double* data,
								   			 unsigned int length )
{
	double maximum_difference;
	double* sorted_data;

	//===Mallocs===//
	sorted_data = malloc(length*sizeof(double));

	//===Sort Array===//
	copy_array_dbl(data, length, sorted_data);
	sort_array_dbl(sorted_data, length);

	//===Find Maximum Difference===//
	maximum_difference = sorted_data[length-1] - sorted_data[0];

	//===Clean Up===//
	free(sorted_data);

	return maximum_difference;
}

double find_minimum_abs_pairwise_difference_dbl( double* data,
								   	   			 unsigned int length )
{
	unsigned int i;
	double difference, minimum_difference;
	double* sorted_data;

	//===Mallocs===//
	sorted_data = malloc(length*sizeof(double));

	//===Sort Array===//
	copy_array_dbl(data, length, sorted_data);
	sort_array_dbl(sorted_data, length);

	//===Find Maximum Increment===//
	minimum_difference = DBL_MAX;
	for (i=1; i<length; i++){
		difference = sorted_data[i] - sorted_data[i-1];
		if (fabs(difference) < fabs(minimum_difference)){
			minimum_difference = difference;
		}
	}

	//===Clean Up===//
	free(sorted_data);

	return minimum_difference;
}

double find_maximum_magnitude_cmplx( double complex* data, 
					 	 			 unsigned int length )
{
	unsigned int i;
	double maximum;

	maximum = -10e300;
	for (i=0; i<length; i++){
		if (cabs(data[i]) > maximum){
			maximum = cabs(data[i]);
		}
	}

	return maximum;
}

double find_maximum_angle_cmplx( double complex* data, 
					 	 	     unsigned int length )
{
	unsigned int i;
	double angle, maximum;

	maximum = -10e300;
	for (i=0; i<length; i++){
		angle = fabs(atan2(cimag(data[i]), creal(data[i])));
		if (angle > maximum){
			maximum = angle;
		}
	}

	return maximum;
}

double complex compute_sum_cmplx( double complex* data,
								  unsigned int length )
{
	unsigned int i;
	double complex sum;

	sum = 0;
	for (i=0; i<length; i++){
		sum += data[i];
	}

	return sum;
}

double compute_sum_dbl( double* data,
					 	unsigned int length )
{
	unsigned int i;
	double sum;

	sum = 0;
	for (i=0; i<length; i++){
		sum += data[i];
	}

	return sum;
}

double find_maximum_dbl( double* data, 
					 	 unsigned int length )
{
	unsigned int i;
	double maximum;

	maximum = -10.0*EPS;
	for (i=0; i<length; i++){
		if (data[i] > maximum){
			maximum = data[i];
		}
	}

	return maximum;
}

double find_abs_maximum_dbl( double* data, 
					 	 	 unsigned int length )
{
	unsigned int i;
	double maximum;
	maximum = -(10.0*EPS);
	for (i=0; i<length; i++){
		if (fabs(data[i]) > fabs(maximum)){
			maximum = data[i];
		}
	}
	return fabs(maximum);
}


unsigned int find_maximum_uint( unsigned int* data, 
					 	 		unsigned int length )
{
	unsigned int i;
	unsigned int maximum;

	maximum = 0;
	for (i=0; i<length; i++){
		if (data[i] > maximum){
			maximum = data[i];
		}
	}

	return maximum;
}


double find_median_dbl( double* data,
						unsigned int length )
{
	double copy[MAX_SIGNAL_LENGTH];

	//===Copy Sort And Return Middle Number===//	
	copy_array_dbl(data, length, copy);
	sort_array_dbl(copy, length);

	return copy[length/2];
}

double find_median_malloc_dbl( double* data,
							   unsigned int length )
{
	double value;
	double* copy;

	//===Mallocs===//
	copy = malloc(length*sizeof(double));

	//===Copy Sort And Return Middle Number===//	
	copy_array_dbl(data, length, copy);
	sort_array_dbl(copy, length);

	//===Get Median===//
	value = copy[length/2];

	//===Free===//
	free(copy);

	return value;
}


//NOTE: This function needs not to threshold for minimax design!
void find_local_maxima_dbl( double* data,
							unsigned int length,
							double* threshold,
							unsigned int* num_maxima,
							unsigned int* positions,
							double* maxima )
{

	unsigned int i, j;
	double thresh;

	//===Find Maxima===//
	j = 0; thresh = 0; if (threshold != NULL) thresh = *threshold;
	for (i=1; i<length-1; i++){
		if (data[i] >= data[i+1] + thresh && data[i] > data[i-1] + thresh){
			if (positions != NULL){
				positions[j] = i;
			}
			if (maxima != NULL){
				maxima[j] = data[i];
			}
			j++;
		}
	}
	if (num_maxima != NULL){
		*num_maxima = j;
	}
	return;
}


unsigned int find_maximum_index_dbl( double* data,
						 			 unsigned int length)
{
	unsigned int i, index;
	double max;

	max = -DBL_MAX;
	index = 0;
	for (i=0; i<length; i++){
		if (data[i] > max){
			max = fabs(data[i]);
			index = i;
		}
	}
	return index;
}


double find_minimum_dbl( double* data, 
					 	 unsigned int length )
{
	unsigned int i;
	double minimum;

	minimum = 10e300;
	for (i=0; i<length; i++){
		if (data[i] < minimum){
			minimum = data[i];
		}
	}

	return minimum;
}

unsigned int find_minimum_uint( unsigned int* data, 
					 	 		unsigned int length )
{
	unsigned int i;
	unsigned int minimum;

	minimum = UINT_MAX;
	for (i=0; i<length; i++){
		if (data[i] < minimum){
			minimum = data[i];
		}
	}

	return minimum;
}

unsigned int find_minimum_index_dbl( double* data,
						 			 unsigned int length)
{
	unsigned int i, index;
	double min;

	//===Find Minimum===//
	min = find_minimum_dbl(data, length);
	
	//===Find First Instance Of Minimum===//
	index = 0;
	for (i=0; i<length; i++){
		if (fabs(data[i] - min) < 2.0*EPS){
			index = i;
			break;
		}
	}
	return index;
}

unsigned int find_maximum_magnitude_index_cmplx( double complex* data, 
		  						 				 unsigned int length )
{
	unsigned int i, max_index;
	double maximum;

	max_index = 0;
	maximum = -DBL_MAX;
	for (i=0; i<length; i++){
		if (cabs(data[i]) > cabs(maximum)){
			maximum = data[i];
			max_index = i;
		}
	}

	return max_index;
}


unsigned int find_minimum_magnitude_index_cmplx( double complex* data, 
		  						 				 unsigned int length )
{
	unsigned int i, min_index;
	double minimum;

	min_index = 0;
	minimum = DBL_MAX;
	for (i=0; i<length; i++){
		if (cabs(data[i]) < cabs(minimum)){
			minimum = data[i];
			min_index = i;
		}
	}

	return min_index;
}

void initialize_matrix_uchar( unsigned char* matrix,
							  unsigned int num_rows,
							  unsigned int num_cols )
{
	unsigned int r,c;

	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			matrix[c + r*num_cols] = 0;
		}
	}

	return;
}

void initialize_matrix_dbl( double* matrix,
							unsigned int num_rows,
							unsigned int num_cols )
{
	unsigned int r,c;

	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			matrix[c + r*num_cols] = 0;
		}
	}

	return;
}

void initialize_matrix_cmplx( double complex* matrix,
							  unsigned int num_rows,
							  unsigned int num_cols )
{
	unsigned int r,c;

	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			matrix[c + r*num_cols] = 0;
		}
	}

	return;
}


unsigned int find_matrix_maximum_uint( unsigned int* matrix,
								  	   unsigned int num_rows,
								  	   unsigned int num_cols )
{
	unsigned int r, c, maximum;

	//===Find Maximum===//
	maximum = 0;
	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			if (matrix[c + r*num_cols] > maximum){
				maximum = matrix[c+ r*num_cols];
			}
		}
	}

	return maximum;
}

double find_matrix_maximum_dbl( double* matrix,
								unsigned int num_rows,
								unsigned int num_cols )
{
	unsigned int r, c;
	double maximum;

	//===Find Maximum===//
	maximum = 0;
	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			if (matrix[c + r*num_cols] > maximum){
				maximum = matrix[c+ r*num_cols];
			}
		}
	}

	return maximum;
}



//================================================================================================//
//=====================================CONVERSION ROUTINES========================================//
//================================================================================================//

double q115_to_dbl( int16_t x )
{
	return ( ((double)x) / (1<<15) );
}

int16_t dbl_to_q115( double x )
{
	int16_t y;
	y = (int16_t) (x * (1<<15) + (x>=0 ? 0.5 : -0.5));
	if (y >= 32767){
		y = 32767;
	}
	else if (y <= -32767){
		y = -32767;
	}
	return y;
}

//NOTE: Was previously pow to 2 -- 7-14-16
void get_magnitude_array_cmplx( double complex* array,
						  		unsigned int length,
						  		double* magnitude )
{
	unsigned int i;
	for (i=0; i<length; i++){
		//magnitude[i] = pow(cabs(array[i]),2.0);
		magnitude[i] = cabs(array[i]);
	}

	return;
}

void get_angle_array_cmplx( double complex* data,
							unsigned int length,
							double* angles )
{
	unsigned int i;
	for (i=0; i<length; i++){
		angles[i] = atan2(cimag(data[i]), creal(data[i]));
	}

	return;
}

void convert_vector_uint_to_dbl( unsigned int* vector,
								 unsigned int length,
								 double* dbl_vector )
{
	unsigned int i;
	for (i=0; i<length; i++){
		dbl_vector[i] = (double)(vector[i]);
	}
	return;
}


void convert_matrix_uint_to_uchar( unsigned int* matrix,
					   			   unsigned int num_rows,
						    	   unsigned int num_cols,
								   unsigned char* char_matrix)
{
	unsigned int r, c, maximum;
	double* matrix_temp;
	matrix_temp = malloc(num_rows*num_cols*sizeof(double));


	maximum = find_matrix_maximum_uint(matrix, num_rows, num_cols);
	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			matrix_temp[c + r*num_cols] = (matrix[c + r*num_cols]);
			matrix_temp[c + r*num_cols] /= (double)maximum;
			matrix[c + r*num_cols] = (unsigned int)(matrix_temp[c + r*num_cols] * 255.0);
		}
	}
	free(matrix_temp);

	if (char_matrix != NULL){
		for (r=0; r<num_rows; r++){
			for (c=0; c<num_cols; c++){
				char_matrix[c+r*num_cols] = matrix[c+r*num_cols];
			}
		}
	}

	return;
}

void conjugate_array_cmplx( double complex* array,
							unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] = conj(array[i]);
	}
	return;
}



//================================================================================================//
//==================================COMBINATORICS ROUTINES========================================//
//================================================================================================//

unsigned int n_choose_k_uint( unsigned int n,
							  unsigned int k )
{
	unsigned int i, out;
	double n_choose_k;

	n_choose_k = 1;
	for (i=1; i<=k; i++){
		n_choose_k *= (((double)(n+1-i))/((double)i));
	}
	out = (unsigned int)n_choose_k;
	return out;
}

double n_choose_k_dbl( double n,
					   double k )
{
	double num, denom;

	//===Check Numerator===//
	if ( fabs(n+1) < 10e-8){
		return 0;
	}

	//===Check Denominator 1===//
	if ( fabs(k+1) < 10e-8){
		return 0;
	}
	if ( k+1 < 0 && fabs( fabs(k+1) - floor( fabs(k+1) )) < 10e-8 ){
		return 0;
	} 
	if ( k+1 < 0 && fabs( fabs(k+1) - ceil(fabs(k+1))) < 10e-8 ){
		return 0;
	} 

	//===Check Denominator 2===//
	if ( fabs(n-k+1) < 10e-8){
		return 0;
	}
	if ( n-k+1 < 0 && fabs( fabs(n-k+1) - floor( fabs(n-k+1) )) < 10e-8 ){
		return 0;
	} 
	if ( n-k+1 < 0 && fabs( fabs(n-k+1) - ceil(fabs(n-k+1))) < 10e-8 ){
		return 0;
	} 

	//===Calculate===//
	num = gamma_dbl(n+1);
	denom = gamma_dbl(k+1) * gamma_dbl(n-k+1);

	//fprintf(stdout, "N: %lf K: %lf NK1: %lf\n", n, k, n-k+1);
	//fprintf(stdout, "N: %lf D1: %5.5e D2: %5.5e\n", num, gamma_dbl(k+1), gamma_dbl(n-k+1));

	return num/denom;
}



//================================================================================================//
//=======================================STATISTICS ROUTINES======================================//
//================================================================================================//

double complex compute_mean_cmplx( double complex* array,
					 	   		   unsigned int length )
{
	unsigned int i;
	double complex mean;

	mean = 0;
	for (i=0; i<length; i++){
		mean += array[i];
	}
	mean /= (double complex)length;

	return mean;
}

double compute_mean_dbl( double* array,
					 	 unsigned int length )
{
	unsigned int i;
	double mean;

	mean = 0;
	for (i=0; i<length; i++){
		mean += array[i];
	}
	mean /= (double)length;

	return mean;
}

unsigned int compute_mean_uint( unsigned int* array,
					 	 		unsigned int length )
{
	unsigned int i;
	unsigned int mean;

	mean = 0;
	for (i=0; i<length; i++){
		mean += array[i];
	}
	mean /= length;

	return mean;
}

double compute_geometric_mean_dbl( double* array,
								   unsigned int length )
{
	unsigned int i;
	double mean;

	mean = 1;
	for (i=0; i<length; i++){
		mean *= array[i];
	}
	mean = pow(mean, 1.0/((double)length));

	return mean;
}

void recurse_mean_dbl( double* array,
					   unsigned int length,
					   unsigned int iteration,
					   double* mean )
{
	unsigned int i;
	for (i=0; i<length; i++){
		mean[i] = mean[i] + (1.0/(double)(iteration+1)) * (array[i] - mean[i]);
	}

	return;
} 

void remove_mean_cmplx( double complex* array,
					  	unsigned int length )
{
	unsigned int i;
	double complex mean;
	
	mean = compute_mean_cmplx(array, length);
	for (i=0; i<length; i++){
		array[i] -= mean;
	}	

	return;
}

void remove_mean_dbl( double* array,
					  unsigned int length )
{
	unsigned int i;
	double mean;
	
	mean = compute_mean_dbl(array, length);
	for (i=0; i<length; i++){
		array[i] -= mean;
	}	

	return;
}

void add_constant_dbl( double* array,
					   unsigned int length,
					   double constant )
{
	unsigned int i;
	for (i=0; i<length; i++){
		array[i] += constant;
	}	

	return;
}

double compute_variance_dbl( double* array,
						 	 unsigned int length )
{
	unsigned int i;
	double mean, variance;

	mean = compute_mean_dbl(array, length);
	variance = 0;
	for (i=0; i<length; i++){
		variance += pow((array[i] - mean),2.0);
	}
	variance /= (((double)length) - 1.0);

	return variance;
}

double compute_sample_skewness_dbl( double* array,
									unsigned int length )
{
	unsigned int i;
	double m2, m3, g1;
	double mean;

	//===Compute M2 and M3===//
	mean = compute_mean_dbl(array, length);
	g1 = m2 = m3 = 0;
	for (i=0; i<length; i++){
		m2 += pow(array[i] - mean, 2.0);
		m3 += pow(array[i] - mean, 3.0);
	}
	m2 /= (double)length;
	m3 /= (double)length;

	//===Make Sample Skewness==//
	g1 = m3/pow(m2, 1.5);
	g1 *= sqrt( (double)(length * (length-1)))/((double)(length-2)+EPS);

	return g1;
}

double compute_sample_kurtosis_dbl( double* array,
									unsigned int length )
{
	unsigned int i;
	double m2, m4, kurtosis, excess_kurtosis;
	double mean;

	//===Compute M2 and M3===//
	mean = compute_mean_dbl(array, length);
	kurtosis = excess_kurtosis = m2 = m4 = 0;
	for (i=0; i<length; i++){
		m2 += pow(array[i] - mean, 2.0);
		m4 += pow(array[i] - mean, 4.0);
	}
	m2 /= (double)length;
	m4 /= (double)length;

	//===Make Kurtosis==//
	kurtosis = m4/pow(m2, 2.0);
	excess_kurtosis = kurtosis - 3.0;

	//===Make Sample Excess Kurtosis===//
	excess_kurtosis *= (double)(length+1);
	excess_kurtosis += 6;
	excess_kurtosis *= ((double)(length-1)+EPS);
	excess_kurtosis /= (((double)(length-2)+EPS)*((double)(length-3)+EPS));

	return excess_kurtosis;
}

void recurse_mean_cmplx( double complex* array,
					   	 unsigned int length,
					   	 unsigned int iteration,
					   	 double complex* mean )
{
	unsigned int i;
	for (i=0; i<length; i++){
		mean[i] = mean[i] + (1.0/(double)(iteration+1)) * (array[i] - mean[i]);
	}

	return;
} 


void normalize_standard_dbl( double* array,
							 unsigned int length )
{
	unsigned int i;
	double mean, std;

	mean = compute_mean_dbl(array, length);
	std = sqrt(compute_variance_dbl(array,length));

	for (i=0; i<length; i++){
		array[i] -= mean;
		array[i] /= std;
	}

	return;
}

void threshold_array_dbl( double* array,
						  unsigned int length,
						  double threshold )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if (array[i] <= threshold){
			array[i] = 0;
		}
	}

	return;
}

void threshold_array_abs_dbl( double* array,
						  	  unsigned int length,
						  	  double threshold )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if (fabs(array[i]) <= threshold){
			array[i] = 0;
		}
	}

	return;
}


void clip_array_dbl( double* array,
					 unsigned int length,
					 double minimum,
					 double maximum )
{
	unsigned int i;
	if (maximum < minimum){
		fprintf(stderr, "Error:: Maximum Is Less Than Minimum! In Function -- clip_array_dbl\n");
		return;
	}	
	for (i=0; i<length; i++){
		if (array[i] < minimum){
			array[i] = minimum;
		}
		if (array[i] > maximum){
			array[i] = maximum;
		}
	}
	return;
}

int compute_mean_error_int( int* vector1,
						    int* vector2,
						    unsigned int length )
{	
	unsigned int i;
	int error;

	error = 0;
	for (i=0; i<length; i++){
		error += abs(((int)(vector1[i] - vector2[i])));
	}
	error /= (int)length;

	return error;
}

double complex compute_mean_squared_error_cmplx( double complex* vector1,
									     		 double complex* vector2,
									   	 		 unsigned int length )
{	
	unsigned int i;
	double complex error;

	error = 0;
	for (i=0; i<length; i++){
		error += cpow(vector1[i]-vector2[i],2.0);
	}
	error /= (double complex)length;

	return error;
}


double compute_mean_squared_error_dbl( double* vector1,
									   double* vector2,
									   unsigned int length )
{	
	unsigned int i;
	double error;

	error = 0;
	for (i=0; i<length; i++){
		error += pow(vector1[i]-vector2[i],2.0);
	}
	error /= (double)length;

	return error;
}

double compute_manhattan_distance_dbl( double* vector1,
									   double* vector2,
									   unsigned int length )
{
	unsigned int i;
	double distance;

	distance = 0;
	for (i=0; i<length; i++){
		distance += fabs(vector1[i] - vector2[i]);
	}

	return distance;
}

double compute_minkowski_distance_dbl( double* vector1,
									   double* vector2,
									   unsigned int length,
									   double power )
{
	
	unsigned int i;
	double distance;

	distance = 0;
	for (i=0; i<length; i++){
		distance += pow(fabs(vector1[i] - vector2[i]),power);
	}
	distance = pow(distance, 1.0/power);

	return distance;

}

double compute_chebyshev_distance_dbl( double* vector1,
									   double* vector2,
									   unsigned int length )
{	
	unsigned int i;
	double tmp_distance, distance;
	
	distance = 0;
	for (i=0; i<length; i++){
		tmp_distance = fabs(vector1[i] - vector2[i]);
		if (tmp_distance > distance) distance = tmp_distance;
	}	

	return distance;
}

double compute_dice_distance_dbl( double* vector1,
								  double* vector2,
								  unsigned int length )
{
	unsigned int i;
	double numer, denom1, denom2, distance;

	numer = denom1 = denom2 = 0;
	for (i=0; i<length; i++){
		numer += pow(vector1[i] - vector2[i], 2.0);
		denom1 += vector1[i] * vector1[i];
		denom2 += vector2[i] * vector2[i];
	}
	distance = numer/(denom1 + denom2);

	return distance;
}


double compute_maximum_error_cmplx( double complex* vector1,
								  	double complex* vector2,
								  	unsigned int length )
{	
	unsigned int i;
	double error, max_error;

	max_error = error = 0;
	for (i=0; i<length; i++){
		error = cabs(cpow(vector1[i]-vector2[i],2.0));
		//for debug
		//fprintf(stdout, "(%lf, %lf) (%lf, %lf) Error: %lf\n", creal(vector1[i]), cimag(vector1[i]), creal(vector2[i]), cimag(vector2[i]), error);
		if (error > max_error){
			max_error = error;
		}
	}

	return max_error;
}


double compute_maximum_error_dbl( double* vector1,
								  double* vector2,
								  unsigned int length )
{	
	unsigned int i;
	double error, max_error;

	max_error = error = 0;
	for (i=0; i<length; i++){
		error = pow(vector1[i]-vector2[i],2.0);
		if (error > max_error){
			max_error = error;
		}
	}

	return max_error;
}


double compute_power_cmplx( double complex* array,
					  		unsigned int length )
{
	unsigned int i;
	double power;

	power = 0;
	for (i=0; i<length; i++){
		power += cabs(cpow(array[i], 2.0));
	}

	return power;
}

double compute_power_dbl( double* array,
					  	  unsigned int length )
{
	unsigned int i;
	double power;

	power = 0;
	for (i=0; i<length; i++){
		power += fabs(pow(array[i], 2.0));
	}

	return power;
}


double compute_vector_2norm_cmplx( double complex* vector,
					  			   unsigned int length )
{
	return sqrt(compute_power_cmplx(vector, length));
}


double compute_vector_2norm_dbl( double* vector,
					  			 unsigned int length )
{
	return sqrt(compute_power_dbl(vector, length));
}

double compute_vector_1norm_dbl( double* vector,
								 unsigned int length )
{
	unsigned int i;
	double norm;
	norm = 0;
	for (i=0; i<length; i++){
		norm += vector[i];
	}
	norm /= (double)length;
	return norm;
}

double compute_matrix_1norm_cmplx( double complex* matrix,
								   unsigned int num_rows,
								   unsigned int num_cols )
{
	unsigned int c;
	double norm;
	double complex column[MAX_POLYNOMIAL_ORDER+1];
	double column_sum[MAX_POLYNOMIAL_ORDER+1];
	
	for (c=0; c<num_cols; c++){
		//===Get Column===//
		copy_matrix_column_cmplx(matrix, num_rows, num_cols, c, column);

		//===Abs Column===//
		abs_array_cmplx(column, num_rows);
	
		//===Compute Sum===//
		column_sum[c] = cabs(compute_sum_cmplx(column,num_rows));
	}
	norm = find_maximum_dbl(column_sum, num_cols); 
	return norm;
}				


double compute_matrix_1norm_dbl( double* matrix,
								 unsigned int num_rows,
								 unsigned int num_cols )
{
	unsigned int c;
	double norm;
	double column[MAX_POLYNOMIAL_ORDER+1];
	double column_sum[MAX_POLYNOMIAL_ORDER+1];
	
	for (c=0; c<num_cols; c++){
		//===Get Column===//
		copy_matrix_column_dbl(matrix, num_rows, num_cols, c, column);

		//===Abs Column===//
		abs_array_dbl(column, num_rows);
	
		//===Compute Sum===//
		column_sum[c] = compute_sum_dbl(column,num_rows);
	}
	norm = find_maximum_dbl(column_sum, num_cols); 
	return norm;
}				
	
double compute_matrix_infnorm_dbl( double* matrix,
								   unsigned int num_rows,
								   unsigned int num_cols )
{
	unsigned int r;
	double norm;
	double row[MAX_POLYNOMIAL_ORDER+1];
	double row_sum[MAX_POLYNOMIAL_ORDER+1];
	
	for (r=0; r<num_rows; r++){
		//===Get Column===//
		copy_matrix_row_dbl(matrix, num_cols, r, row);

		//===Abs Column===//
		abs_array_dbl(row, num_cols);
	
		//===Compute Sum===//
		row_sum[r] = compute_sum_dbl(row,num_cols);
	}
	norm = find_maximum_dbl(row_sum, num_rows); 
	return norm;
}				

double compute_matrix_fnorm_dbl( double* matrix,
								 unsigned int num_rows,
								 unsigned int num_cols )
{
	unsigned int r,c;
	double norm;
	norm = 0;
	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			norm += matrix[c + r*num_cols] * matrix[c + r*num_cols];
		}
	}	
	norm = pow(norm, 0.5);
	return norm;
}

void exponentially_smooth_array_dbl( double* array,
									 double* update,
									 unsigned int length,
									 double alpha,
									 double* result )
{
	unsigned int i;
	if (alpha > 1 || alpha < 0 ){
		fprintf(stderr, "Error:: Alpha Is Not [0,1]! In Function -- exponentially_smooth_array_dbl!\n");
		quit();
	}
	
	//===Smooth===//
	for (i=0; i<length; i++){
		result[i] = alpha*array[i] + (1.0-alpha)*update[i];
	}

	return;
}


//================================================================================================//
//=====================================LINEAR ALGEBRA ROUTINES====================================//
//================================================================================================//

void hadamard_product_cmplx( double complex* vector1,
					   	   	 double complex* vector2,
					   	   	 unsigned int length,
					   	     double complex* result )
{
	unsigned int i;
	for (i=0; i<length; i++){
		result[i] = vector1[i] * vector2[i];
	}
	return;
}

void hadamard_product_dbl( double* vector1,
					   	   double* vector2,
					   	   unsigned int length,
					   	   double* result )
{
	unsigned int i;
	for (i=0; i<length; i++){
		result[i] = vector1[i] * vector2[i];
	}
	return;
}

double complex compute_inner_product_cmplx( double complex* vector1,
					   	   		  			double complex* vector2,
					   	   		  			unsigned int length )
{
	unsigned int i;
	double complex product;

	//===Conjugate And Multiply===//
	product = 0;
	conjugate_array_cmplx(vector1, length);	
	for (i=0; i<length; i++){
		product += vector1[i] * vector2[i];
	}

	//===Fix Conjugation===//
	conjugate_array_cmplx(vector1, length);	

	return product;
}

double compute_inner_product_dbl( double* vector1,
					   	   		  double* vector2,
					   	   		  unsigned int length )
{
	unsigned int i;
	double product;

	product = 0;
	for (i=0; i<length; i++){
		product += vector1[i] * vector2[i];
	}
	return product;
}

double compute_vector_distance_cmplx( double complex* vector1,
									  double complex* vector2,
									  unsigned int length,
									  double p )
{
	unsigned int i;
	double distance;

	distance = 0;
	for (i=0; i<length; i++){
		distance += pow(cabs(vector1[i] - vector2[i]), p);
	}
	distance = pow(distance, 1.0/p);

	return distance;
}	


double compute_vector_distance_dbl( double* vector1,
									double* vector2,
									unsigned int length,
									double p )
{
	unsigned int i;
	double distance;

	distance = 0;
	for (i=0; i<length; i++){
		distance += pow((vector1[i] - vector2[i]), p);
	}
	distance = pow(distance, 1.0/p);

	return distance;
}	

void subtract_vectors_cmplx( double complex* vector1,
				  	    	 double complex* vector2,
				  	    	 unsigned int length,
				  	    	 double complex* result )
{
	unsigned int i;
	for (i=0; i<length; i++){
		result[i] = vector1[i] - vector2[i];
	}
	return;
}

void subtract_vectors_dbl( double* vector1,
				  	       double* vector2,
				  	       unsigned int length,
				  	       double* result )
{
	unsigned int i;
	for (i=0; i<length; i++){
		result[i] = vector1[i] - vector2[i];
	}
	return;
}


void add_vectors_cmplx( double complex* vector1,
				  	    double complex* vector2,
				  	    unsigned int length,
				  	    double complex* result )
{
	unsigned int i;
	for (i=0; i<length; i++){
		result[i] = vector1[i] + vector2[i];
	}
	return;
}

void add_vectors_dbl( double* vector1,
				  	  double* vector2,
				  	  unsigned int length,
				  	  double* result )
{
	unsigned int i;
	for (i=0; i<length; i++){
		result[i] = vector1[i] + vector2[i];
	}
	return;
}

void add_matrices_cmplx( double complex* matrix1,
					     double complex* matrix2,
					     unsigned int num_rows,
					     unsigned int num_cols,
					     double complex* result )
{
	unsigned int i;
	for (i=0; i<num_rows*num_cols; i++){
		result[i] = matrix1[i] + matrix2[i];
	}
	return;
}

void subtract_matrices_cmplx( double complex* matrix1,
					          double complex* matrix2,
					          unsigned int num_rows,
					          unsigned int num_cols,
					          double complex* result )
{
	unsigned int i;
	for (i=0; i<num_rows*num_cols; i++){
		result[i] = matrix1[i] - matrix2[i];
	}
	return;
}

void add_matrices_dbl( double* matrix1,
					   double* matrix2,
					   unsigned int num_rows,
					   unsigned int num_cols,
					   double* result )
{
	unsigned int i;
	for (i=0; i<num_rows*num_cols; i++){
		result[i] = matrix1[i] + matrix2[i];
	}
	return;
}

void negate_vector_dbl( double* vector,
						unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		vector[i] *= -1.0;
	}
	return;
}

void invert_vector_cmplx( double complex* vector,
						  unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if (cabs(vector[i]) <= 10e-16){
			vector[i] = 0.0;
		}
		else{
			vector[i] = 1.0/vector[i];
		}
	}
	return;
}

void invert_vector_dbl( double* vector,
						unsigned int length )
{
	unsigned int i;
	for (i=0; i<length; i++){
		if (fabs(vector[i]) <= 10e-16){
			vector[i] = 0.0;
		}
		else{
			vector[i] = 1.0/vector[i];
		}
	}
	return;
}


void divide_vectors_dbl( double* vector1,
						 double* vector2,
						 unsigned int length,
						 double* result )
{
	unsigned int i,j;
	unsigned int indices[MAX_SIGNAL_LENGTH];
	j = 0;
	for (i=0; i<length; i++){
		if (fabs(vector2[i]) <= 10e-16){
			result[i] = 0.0;
			indices[j] = i;
			j++;
		}
		else{
			result[i] = vector1[i]/vector2[i];
		}
	}
	for (i=0; i<j; i++){
		if (indices[i]-1 > 0 && indices[i] + 1 < length){
			result[indices[i]] = 0.5*(result[indices[i]-1] + result[indices[i]+1]); 
		}
	}


	return;
}  

//================================================================================================//
//===============================MATRIX MANIPULATION ROUTINES=====================================//
//================================================================================================//

void set_matrix_diagonal_cmplx( double complex* matrix,
							  	unsigned int num_cols,
							  	double complex* elements,
							  	unsigned int num_elements )
{
	unsigned int i;
	for (i=0; i<num_elements; i++){
		matrix[i + i*num_cols] = elements[i];		
	}

	return;
}

void set_matrix_diagonal_dbl( double* matrix,
							  unsigned int num_cols,
							  double* elements,
							  unsigned int num_elements )
{
	unsigned int i;
	for (i=0; i<num_elements; i++){
		matrix[i + i*num_cols] = elements[i];		
	}

	return;
}

void set_matrix_offdiagonal_dbl( double* matrix,
								 unsigned int num_rows,
							  	 unsigned int num_cols,
							  	 double* elements,
							  	 unsigned int num_elements,
								 int off_diagonal_index )
{
	unsigned int i;
	int i1;
	i1 = off_diagonal_index;
	for (i=0; i<num_elements; i++){
		if (i + (i+i1)*num_cols > 0 && i + (i+i1)*num_cols < num_rows*num_cols){
			matrix[i + (i+i1)*num_cols] = elements[i];		
		}
	}

	return;
}


void copy_matrix_diagonal_cmplx( double complex* matrix,
							  	 unsigned int num_rows,
								 unsigned int num_cols,
							  	 double complex* diagonal )
{
	unsigned int i;
	for (i=0; i<num_rows; i++){
		diagonal[i] = matrix[i + i*num_cols];		
	}

	return;
}

void copy_matrix_diagonal_dbl( double* matrix,
							   unsigned int num_rows,
							   unsigned int num_cols,
							   double* diagonal )
{
	unsigned int i;
	for (i=0; i<num_rows; i++){
		diagonal[i] = matrix[i + i*num_cols];		
	}

	return;
}

void copy_matrix_offdiagonal_dbl( double* matrix,
							   	  unsigned int num_rows,
							   	  unsigned int num_cols,
								  int off_diagonal_index,
							   	  double* off_diagonal )
{

	unsigned int i, k;
	int i1;
	i1 = off_diagonal_index; k = 0;
	for (i=0; i<num_rows; i++){
		if (i + (i+i1)*num_cols > 0 && i + (i+i1)*num_cols < num_rows*num_cols){
			off_diagonal[k++] = matrix[i + (i+i1)*num_cols];		
		}
	}

	return;
}

void copy_matrix_row_cmplx( double complex* matrix,
							unsigned int num_cols,
							unsigned int row_number,
							double complex* row )
{	
	unsigned int j;
	for (j=0; j<num_cols; j++){
		row[j] = matrix[j + row_number*num_cols];
	}
	return;
}

void copy_matrix_row_dbl( double* matrix,
					      unsigned int num_cols,
						  unsigned int row_number,
						  double* row )
{	
	unsigned int j;
	for (j=0; j<num_cols; j++){
		row[j] = matrix[j + row_number*num_cols];
	}
	return;
}

void copy_matrix_row_uint( unsigned int* matrix,
					       unsigned int num_cols,
						   unsigned int row_number,
						   unsigned int* row )
{	
	unsigned int j;
	for (j=0; j<num_cols; j++){
		row[j] = matrix[j + row_number*num_cols];
	}
	return;
}

void copy_matrix_column_cmplx( double complex* matrix,
							   unsigned int num_rows,
							   unsigned int num_cols,
							   unsigned int column_number,
							   double complex* column )
{	
	unsigned int j;
	for (j=0; j<num_rows; j++){
		column[j] = matrix[column_number + j*num_cols];
	}
	return;
}

void copy_matrix_column_dbl( double* matrix,
							 unsigned int num_rows,
						     unsigned int num_cols,
							 unsigned int column_number,
							 double* column )
{	
	unsigned int j;
	for (j=0; j<num_rows; j++){
		column[j] = matrix[column_number + j*num_cols];
	}
	return;
}

void copy_matrix_column_uint( unsigned int* matrix,
							  unsigned int num_rows,
						      unsigned int num_cols,
							  unsigned int column_number,
							  unsigned int* column )
{	
	unsigned int j;
	for (j=0; j<num_rows; j++){
		column[j] = matrix[column_number + j*num_cols];
	}
	return;
}


void replace_matrix_column_cmplx( double complex* matrix,
								  unsigned int num_rows,
								  unsigned int num_cols,
								  double complex* column,
								  unsigned int column_number )
{
	unsigned int i;
	for (i=0; i<num_rows; i++){
		matrix[column_number + num_cols*i] = column[i];
	}
	return;
}

void replace_matrix_column_dbl( double* matrix,
								unsigned int num_rows,
								unsigned int num_cols,
								double* column,
								unsigned int column_number )
{
	unsigned int i;
	for (i=0; i<num_rows; i++){
		matrix[column_number + num_cols*i] = column[i];
	}
	return;
}

void replace_matrix_row_cmplx( double complex* matrix,
							   unsigned int num_cols,
							   double complex* row,
							   unsigned int row_number )
{
	unsigned int i;
	for (i=0; i<num_cols; i++){
		matrix[i + num_cols*row_number] = row[i];
	}
	return;
}

void replace_matrix_row_dbl( double* matrix,
							 unsigned int num_cols,
							 double* row,
							 unsigned int row_number )
{
	unsigned int i;
	for (i=0; i<num_cols; i++){
		matrix[i + num_cols*row_number] = row[i];
	}
	return;
}

void replace_matrix_row_uint( unsigned int* matrix,
							  unsigned int num_cols,
							  unsigned int* row,
							  unsigned int row_number )
{
	unsigned int i;
	for (i=0; i<num_cols; i++){
		matrix[i + num_cols*row_number] = row[i];
	}
	return;
}



void partition_matrix_columns_cmplx( double complex* matrix,
								     unsigned int num_rows,
								     unsigned int num_cols,
								     unsigned int column_break,
								   	 double complex* block1,
									 double complex* block2 )
{
	unsigned int i, c;
	double complex column[MAX_POLYNOMIAL_ORDER+1];

	if (block1 != NULL){
		for (c=0; c<column_break; c++){
			copy_matrix_column_cmplx(matrix, num_rows, num_cols, c, column);
			replace_matrix_column_cmplx(block1, num_rows, column_break, column, c);
		}
	}
	if (block2 != NULL){
		i = 0;
		for (c=column_break; c<num_cols; c++){
			copy_matrix_column_cmplx(matrix, num_rows, num_cols, c, column);
			replace_matrix_column_cmplx(block2, num_rows, num_cols - column_break, column, i);
			i++;
		}
	}

	return;
}

void partition_matrix_rows_cmplx( double complex* matrix,
								  unsigned int num_rows,
								  unsigned int num_cols,
								  unsigned int row_break,
								  double complex* block1,
								  double complex* block2 )
{
	unsigned int i, r;
	double complex row[MAX_POLYNOMIAL_ORDER+1];

	if (block1 != NULL){
		for (r=0; r<row_break; r++){
			copy_matrix_row_cmplx(matrix, num_cols, r, row);
			replace_matrix_row_cmplx(block1, num_cols, row, r);
		}
	}
	if (block2 != NULL){
		i = 0;
		for (r=row_break; r<num_rows; r++){
			copy_matrix_row_cmplx(matrix, num_cols, r, row);
			replace_matrix_row_cmplx(block2, num_cols, row, i);
			i++;
		}
	}

	return;
}


void transpose_matrix_cmplx( double complex* matrix,
							 unsigned int num_rows,
							 unsigned int num_cols )
{
	unsigned int i;
	double complex* copy;
	double complex* row;

	copy = malloc(num_rows * num_cols * sizeof(double complex));
	row = malloc(num_cols * sizeof(double complex));

	//===Initialize Copy===//
	initialize_array_cmplx(copy, num_rows * num_cols);

	//===Copy Rows As Columns===//
	for (i=0; i<num_rows; i++){

		//===Zero Row===//
		initialize_array_cmplx(row, num_cols);

		//===Get Row===//
		copy_matrix_row_cmplx(matrix, num_cols, i, row);

		//===Save As Column===//
		replace_matrix_column_cmplx(copy, num_cols, num_rows, row, i);

	}

	//===Copy Over===//
	copy_array_cmplx(copy, num_rows * num_cols, matrix);

	//===Clean Up===//
	free(row);
	free(copy);

	return;
}

void transpose_matrix_dbl( double* matrix,
						   unsigned int num_rows,
						   unsigned int num_cols )
{
	unsigned int i;
	double* copy;
	double* row;

	copy = malloc(num_rows * num_cols * sizeof(double));
	row = malloc(num_rows * num_cols * sizeof(double));

	//===Initialize Copy===//
	initialize_array_dbl(copy, num_rows * num_cols);

	//===Copy Rows As Columns===//
	for (i=0; i<num_rows; i++){

		//===Zero Row===//
		initialize_array_dbl(row, num_cols);

		//===Get Row===//
		copy_matrix_row_dbl(matrix, num_cols, i, row);

		//===Save As Column===//
		replace_matrix_column_dbl(copy, num_cols, num_rows, row, i);

	}


	//===Copy Over===//
	copy_array_dbl(copy, num_rows * num_cols, matrix);

	//===Clean Up===//
	free(row);
	free(copy);

	return;
}

void hermitian_transpose_matrix_cmplx( double complex* matrix,
							 		   unsigned int num_rows,
							 		   unsigned int num_cols )
{
	transpose_matrix_cmplx(matrix, num_rows, num_cols);
	conjugate_array_cmplx(matrix, num_rows*num_cols);	
	return;
}

void append_zeros_matrix_cmplx( double complex* matrix,
							  	unsigned int num_rows,
							  	unsigned int num_cols,
							  	unsigned int extra_rows,
							  	unsigned int extra_cols )
{
	unsigned int r, c;
	double complex* copy;
	copy = malloc((num_rows + extra_rows) * (num_cols + extra_cols) * sizeof(double complex));

	//===Copy Matrix===//
	initialize_array_cmplx(copy, (num_rows + extra_rows) * (num_cols + extra_cols));
	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			copy[c + r*(num_cols+extra_cols)] = matrix[c + r*(num_cols)];
		}	
	}

	//===Copy Back===//
	copy_array_cmplx(copy, (num_rows + extra_rows) * (num_cols + extra_cols), matrix);	

	//===Clean Up===//
	free(copy);

	return;
}

void append_zeros_matrix_dbl( double* matrix,
							  unsigned int num_rows,
							  unsigned int num_cols,
							  unsigned int extra_rows,
							  unsigned int extra_cols )
{
	unsigned int r, c;
	double* copy;
	copy = malloc((num_rows + extra_rows) * (num_cols + extra_cols) * sizeof(double));

	//===Copy Matrix===//
	initialize_array_dbl(copy, (num_rows + extra_rows) * (num_cols + extra_cols));
	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			copy[c + r*(num_cols+extra_cols)] = matrix[c + r*(num_cols)];
		}	
	}

	//===Copy Back===//
	copy_array_dbl(copy, (num_rows + extra_rows) * (num_cols + extra_cols), matrix);	

	//===Clean Up===//
	free(copy);

	return;
}

void normalize_matrix_columns_dbl( double* matrix,
								   unsigned int num_rows,
								   unsigned int num_cols )
{
	unsigned int c;
	double column[MAX_POLYNOMIAL_ORDER+1];
	for (c=0; c<num_cols; c++){
		copy_matrix_column_dbl(matrix, num_rows, num_cols, c, column);
		normalize_max_dbl(column, num_rows);
		replace_matrix_column_dbl(matrix, num_rows, num_cols, column, c);
	}
	return;
}

void normalize_matrix_rows_dbl( double* matrix,
								unsigned int num_rows,
								unsigned int num_cols )
{
	unsigned int r;
	double row[MAX_POLYNOMIAL_ORDER+1];
	for (r=0; r<num_rows; r++){
		copy_matrix_row_dbl(matrix, num_cols, r, row);
		normalize_max_dbl(row, num_cols);
		replace_matrix_row_dbl(matrix, num_cols, row, r);
	}
	return;
}

void normalize_matrix_rows_sum_cmplx( double complex* matrix,
									  unsigned int num_rows,
									  unsigned int num_cols )
{
	unsigned int r, i;
	double complex row[MAX_POLYNOMIAL_ORDER+1];
	double complex sum;
	for (r=0; r<num_rows; r++){
		copy_matrix_row_cmplx(matrix, num_cols, r, row);
		sum = compute_sum_cmplx(row, num_cols);
		for (i=0; i<num_cols; i++) row[i] /= sum;
		replace_matrix_row_cmplx(matrix, num_cols, row, r);
	}
	return;
}

void normalize_matrix_rows_sum_dbl( double* matrix,
									unsigned int num_rows,
									unsigned int num_cols )
{
	unsigned int r, i;
	double row[MAX_POLYNOMIAL_ORDER+1];
	double sum;
	for (r=0; r<num_rows; r++){
		copy_matrix_row_dbl(matrix, num_cols, r, row);
		sum = compute_sum_dbl(row, num_cols);
		for (i=0; i<num_cols; i++) row[i] /= sum;
		replace_matrix_row_dbl(matrix, num_cols, row, r);
	}
	return;
}

void normalize_matrix_frobenius_dbl( double* matrix,
									 unsigned int num_rows,
									 unsigned int num_cols )
{
	unsigned int r, c;
	double frobenius_norm;

	//===Calculate Frobenius Norm===//
	frobenius_norm = 0;
	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			frobenius_norm += (fabs(matrix[c + r*num_cols])*fabs(matrix[c + r*num_cols]));
		}
	}
	frobenius_norm = sqrt(frobenius_norm);

	//===Normalize Frobenius Norm===//
	for (r=0; r<num_rows; r++){
		for (c=0; c<num_cols; c++){
			matrix[c + r*num_cols] = matrix[c + r*num_cols]/frobenius_norm;
		}
	}

	return;
}

void normalize_matrix_1norm_cmplx( double complex* matrix,
								   unsigned int num_rows,
								   unsigned int num_cols )
{
	unsigned int c;
	double complex column[MAX_POLYNOMIAL_ORDER+1];
	double complex sum_cmplx;
	double sum, max_sum;

	//===Figure Out Maximum Column Sum===//
	max_sum = -DBL_MAX;
	for (c=0; c<num_cols; c++){
		copy_matrix_column_cmplx(matrix, num_rows, num_cols, c, column);
		sum_cmplx = (compute_sum_cmplx(column, num_rows));
		sum = cabs(sum_cmplx);
		if (sum > max_sum){
			max_sum = sum;
		}
	}

	//===Gain Matrix===//
	sum_cmplx = max_sum + I*0.0;
	gain_array_constant_cmplx(matrix, num_rows*num_cols, 1.0/sum_cmplx);

	return;
}


void normalize_matrix_rows_2norm_cmplx( double complex* matrix,
									    unsigned int num_rows,
									    unsigned int num_cols )
{
	unsigned int r;
	double complex row[MAX_POLYNOMIAL_ORDER+1];
	for (r=0; r<num_rows; r++){
		copy_matrix_row_cmplx(matrix, num_cols, r, row);
		normalize_vector_2norm_cmplx(row,num_cols);
		replace_matrix_row_cmplx(matrix, num_cols, row, r);
	}
	return;
}



void initialize_identity_matrix_cmplx( double complex* matrix,
									   unsigned int num_rows,
									   unsigned int num_cols )
{
	unsigned int i,j;
	for (i=0; i<num_rows; i++){
		for (j=0; j<num_cols; j++){
			if (i == j){
				matrix[j + i*num_cols] = 1.0;
			}
			else{
				matrix[j + i*num_cols] = 0.0;
			}
		}
	}
	return;
}

void initialize_exchange_matrix_cmplx( double complex* matrix,
									   unsigned int num_rows,
									   unsigned int num_cols )
{
	unsigned int i,j;
	for (i=0; i<num_rows; i++){
		for (j=0; j<num_cols; j++){
			if (i == (num_cols-1)-j){
				matrix[j + i*num_cols] = 1.0;
			}
			else{
				matrix[j + i*num_cols] = 0.0;
			}
		}
	}
	return;
}

void initialize_identity_matrix_dbl( double* matrix,
									 unsigned int num_rows,
									 unsigned int num_cols )
{
	unsigned int i,j;
	for (i=0; i<num_rows; i++){
		for (j=0; j<num_cols; j++){
			if (i == j){
				matrix[j + i*num_cols] = 1.0;
			}
			else{
				matrix[j + i*num_cols] = 0.0;
			}
		}
	}
	return;
}

void initialize_exchange_matrix_dbl( double* matrix,
									 unsigned int num_rows,
									 unsigned int num_cols )
{
	unsigned int i,j;
	for (i=0; i<num_rows; i++){
		for (j=0; j<num_cols; j++){
			if (i == (num_cols-1)-j){
				matrix[j + i*num_cols] = 1.0;
			}
			else{
				matrix[j + i*num_cols] = 0.0;
			}
		}
	}
	return;
}


void form_circulant_matrix_dbl( double* data,
							    unsigned int length,
								double* matrix )
{
	unsigned int i;
	double temp_data[MAX_FRAME_LENGTH];

	//===Make Matrix===//
	copy_array_dbl(data, length, temp_data);
	for (i=0; i<length; i++){
		if (i!=0){

			//===Shift Right By 1===//
			right_shift_array_dbl(temp_data, length, 1);

			//===Append Ith sample to the 0th Position===//
			temp_data[0] = data[i];
		}

		//===Replace Ith Column Circulant Matrix===//
		replace_matrix_column_dbl(matrix, length, length, temp_data, i);

	}
	return;
}

void initialize_unit_vector_dbl( double* array,
								 unsigned int length,
								 unsigned int coordinate )
{
	initialize_array_dbl(array, length);
	array[coordinate] = 1;
	return;
}


//================================================================================================//
//=======================================I/O ROUTINES=============================================//
//================================================================================================//

void newline()
{
	fprintf(stdout, "\n");
	return;
}

void hello_world()
{
	fprintf(stdout, "Hello World!\n");
	return;
}

void print_row_vector_uint( unsigned int* vector,
						    unsigned int vector_size,
						    FILE* file )
{

	unsigned int i;

	//===Print===//
	for (i=0; i<vector_size; i++){
		fprintf(file, "%d ", vector[i]);
	}
	fprintf(file, "\n");

	return;
}


void print_vector_uint( unsigned int* vector,
				   	    unsigned int vector_size,
				   	    FILE* file )
{
	unsigned int i;

	for (i=0; i<vector_size; i++){
		fprintf(file, "%+d\n", vector[i]);
	}
	
	return;
}

void print_row_vector_int( int* vector,
						   int vector_size,
						   FILE* file )
{

	int i;

	//===Print===//
	for (i=0; i<vector_size; i++){
		fprintf(file, "%d ", vector[i]);
	}
	fprintf(file, "\n");

	return;
}

void print_vector_int( int* vector,
				   	   unsigned int vector_size,
				   	   FILE* file )
{
	unsigned int i;

	for (i=0; i<vector_size; i++){
		fprintf(file, "%+d\n", vector[i]);
	}
	
	return;
}

void print_row_vector_dbl( double* vector,
						   unsigned int vector_size,
						   FILE* file )
{

	unsigned int i;

	if (isatty(fileno(file))){
		//===Print To Stdout===//
		for (i=0; i<vector_size; i++){
			fprintf(file, "%+6.6f ", vector[i]);
			//fprintf(file, "%+8.8e ", vector[i]);
		}
	}
	else{
		//===Print To File===//
		for (i=0; i<vector_size; i++){
			fprintf(file, "%+32.32e ", vector[i]);
		}
	}
	fprintf(file, "\n");

	return;
}

void print_vector_dbl( double* vector,
				   	   unsigned int vector_size,
				   	   FILE* file )
{
	unsigned int i;

	if (isatty(fileno(file))){
		//===Print To Stdout===//
		for (i=0; i<vector_size; i++){
			fprintf(file, "%+lf\n", vector[i]);
			//fprintf(file, "%+3.3f\n", vector[i]);
			//fprintf(file, "%+16.16e\n", vector[i]);
			//fprintf(file, "%+4.4e\n", vector[i]);
			//fprintf(file, "%+8.8e\n", vector[i]);
		}
	}
	else{
		//===Print To File===//
		for (i=0; i<vector_size; i++){
			//fprintf(file, "%+16.16e\n", vector[i]);
			fprintf(file, "%+32.32e\n", vector[i]);
			//fprintf(file, "%+lf\n", vector[i]);
			//fprintf(file, "%+3.3f\n", vector[i]);
			//fprintf(file, "%+4.4e\n", vector[i]);
		}

	}
	
	return;
}

void print_row_vector_cmplx( double complex* vector,
						   	 unsigned int vector_size,
						   	 FILE* file )
{

	unsigned int i;

	if (isatty(fileno(file))){
		//===Print To Stdout===//
		for (i=0; i<vector_size; i++){
			fprintf(file, "%+3.3e%+3.3ej ", creal(vector[i]),cimag(vector[i]));
			//fprintf(file, "%+8.8e ", vector[i]);
			//fprintf(file, "%+lf%+lfj ", creal(vector[i]),cimag(vector[i]));
		}
	}
	else{
		//===Print To File===//
		for (i=0; i<vector_size; i++){
			fprintf(file, "%+32.32e%+3.32ej ", creal(vector[i]),cimag(vector[i]));
		}
	}
	fprintf(file, "\n");

	return;
}


void print_vector_cmplx( double complex* vector,
				   	     unsigned int vector_size,
				   	     FILE* file )
{
	unsigned int i;

	if (isatty(fileno(file))){
		//===Print To Stdout===//	
		for (i=0; i<vector_size; i++){
			fprintf(file, "%+lf%+lfj\n", creal(vector[i]),cimag(vector[i]));
			//fprintf(file, "%+16.16e %+16.16e\n", creal(vector[i]),cimag(vector[i]));
			//fprintf(file, "%+3.3e%+3.3ej\n", creal(vector[i]),cimag(vector[i]));
			//fprintf(file, "%+3.3e%+3.3ej\n", creal(vector[i]),cimag(vector[i]));
			//fprintf(file, "%+6.6f%+6.6fj\n", creal(vector[i]),cimag(vector[i]));
		}
	}
	else{
		//===Print To File===//
		for (i=0; i<vector_size; i++){
			//fprintf(file, "%+16.16e %+16.16e\n", creal(vector[i]),cimag(vector[i]));
			fprintf(file, "%+6.6e %+6.6e\n", creal(vector[i]),cimag(vector[i]));
		}
	}

	
	return;
}

void print_matrix_uint( unsigned int* matrix,
				   	    int num_rows,
				   	    int num_cols,
				   	    FILE* file )
{
	int i, j;
	if (isatty(fileno(file))){
		//===Print To Stdout===//	
		for (i=0; i<num_rows; i++){
			for (j=0; j<num_cols; j++){
				fprintf(file, "%+d ", matrix[j + i*num_cols]);
			}
			fprintf(file, "\n");
		}
	}
	else{
		//===Print To File===//
		for (i=0; i<num_rows; i++){
			for (j=0; j<num_cols; j++){
				fprintf(file, "%+d ", matrix[j + i*num_cols]);
			}
			fprintf(file, "\n");
		}
	}
	return;
}

void print_matrix_dbl( double* matrix,
				   	   int num_rows,
				   	   int num_cols,
				   	   FILE* file )
{
	int i, j;
	if (isatty(fileno(file))){
		//===Print To Stdout===//	
		for (i=0; i<num_rows; i++){
			for (j=0; j<num_cols; j++){
				//fprintf(file, "%+16.16e ", matrix[j + i*num_cols]);
				fprintf(file, "%+lf ", matrix[j + i*num_cols]);
				//fprintf(file, "%+6.6f ", matrix[j + i*num_cols]);
				//fprintf(file, "%+8.8e ", matrix[j + i*num_cols]);
			}
			fprintf(file, "\n");
		}
	}
	else{
		//===Print To File===//
		for (i=0; i<num_rows; i++){
			for (j=0; j<num_cols; j++){
				//fprintf(file, "%+16.16e ", matrix[j + i*num_cols]);
				fprintf(file, "%+3.3f ", matrix[j + i*num_cols]);
				//fprintf(file, "%+lf ", matrix[j + i*num_cols]);
			}
			fprintf(file, "\n");
		}
	}


	return;
}

void print_matrix_cmplx( double complex* matrix,
				   		 int num_rows,
				   		 int num_cols,
				   		 FILE* file )
{
	int i, j;
	if (isatty(fileno(file))){
		//===Print To Stdout===//	
		for (i=0; i<num_rows; i++){
			for (j=0; j<num_cols; j++){
				//fprintf(file, "%+3.3e%+3.3ej ", creal(matrix[j + i*num_cols]), cimag(matrix[j + i*num_cols]));
				fprintf(file, "%+3.3f%+3.3fj ", creal(matrix[j + i*num_cols]), cimag(matrix[j + i*num_cols]));
				//fprintf(file, "%+6.6f%+6.6fj ", creal(matrix[j + i*num_cols]), cimag(matrix[j + i*num_cols]));
				//fprintf(file, "%+lf%+lfj ", creal(matrix[j + i*num_cols]), cimag(matrix[j + i*num_cols]));
			}
			fprintf(file, "\n");
		}
	}
	else{
		//===Print To File===//
		for (i=0; i<num_rows; i++){
			for (j=0; j<num_cols; j++){
				//fprintf(file, "%+16.16e %+16.16e ", creal(matrix[j + i*num_cols]), cimag(matrix[j + i*num_cols]));
				fprintf(file, "%+6.6lf %+6.6e ", creal(matrix[j + i*num_cols]), cimag(matrix[j + i*num_cols]));
			}
			fprintf(file, "\n");
		}
	}

	return;
}

void read_dat_file_cmplx( char* filename,
						  unsigned int* data_length,
						  double complex* data_out )
{
	FILE* fin;
	int item;
	double real, imag;
	unsigned int i;

	//===Open And Read===//
	fin = fopen(filename, "r");
	i = 0; item = 666;
	while (1){ 
		
		//===Read Real===//
		item = fscanf(fin,"%lf",&real);
		if (item == EOF){
			break;
		}

		//===Read Imaginary===//
		item = fscanf(fin,"%lf",&imag);
		if (item == EOF){
			break;
		}
		
		//===Set Data===//
		data_out[i] = real + I*imag;
		i++;
	}
	if (data_length != NULL){
		*data_length = i;
	}
	
	//===Clean Up===//
	fclose(fin);

	return;
}

unsigned int get_dat_file_length_dbl( char* filename )
{
	FILE* fin;
	unsigned int length;
	int item;
	double real;
	unsigned int i;

	//===Open And Read===//
	fin = fopen(filename, "r");
	i = 0; item = 666;
	while (1){ 
		
		//===Read Real===//
		item = fscanf(fin,"%lf",&real);
		if (item == EOF){
			break;
		}
		
		i++;
	}
	length = i;
	
	//===Clean Up===//
	fclose(fin);

	return length;
}

void read_dat_file_dbl( char* filename,
						unsigned int* data_length,
						double* data_out )
{
	FILE* fin;
	int item;
	double real;
	unsigned int i;

	//===Open And Read===//
	fin = fopen(filename, "r");
	i = 0; item = 666;
	while (1){ 
		
		//===Read Real===//
		item = fscanf(fin,"%lf",&real);
		if (item == EOF){
			break;
		}
		
		//===Set Data===//
		data_out[i] = real;
		i++;
	}
	if (data_length != NULL){
		*data_length = i;
	}

	//===Clean Up===//
	fclose(fin);

	return;
}

void read_dat_file_uint( char* filename,
						 unsigned int* data_length,
						 unsigned int* data_out )
{
	FILE* fin;
	int item;
	unsigned int data;
	unsigned int i;

	//===Open And Read===//
	fin = fopen(filename, "r");
	i = 0; item = 666;
	while (1){ 
		
		//===Read Real===//
		item = fscanf(fin,"%d",&data);
		if (item == EOF){
			break;
		}
		
		//===Set Data===//
		data_out[i] = data;
		i++;
	}
	if (data_length != NULL){
		*data_length = i;
	}

	//===Clean Up===//
	fclose(fin);

	return;
}

void print_column_major_matrix_dbl( double* matrix,
									int num_rows,
									int num_cols,
									FILE* file )
{
	int i, j;
	if (isatty(fileno(file))){
		//===Print To Stdout===//	
		for (i=0; i<num_rows; i++){
			for (j=0; j<num_cols; j++){
				//fprintf(file, "%+16.16e ", matrix[i + j*(num_rows)]);
				fprintf(file, "%+lf ", matrix[i + j*(num_rows)]);
			}
			fprintf(file, "\n");
		}
	}
	else{
		//===Print To File===//
		for (i=0; i<num_rows; i++){
			for (j=0; j<num_cols; j++){
				fprintf(file, "%+16.16e ", matrix[i + j*(num_rows)]);
				//fprintf(file, "%+lf ", matrix[i + j*(num_rows)]);
			}
			fprintf(file, "\n");
		}
	}

	return;
}


//should probably have a function that reads multicolumn dat files!!

void row_to_column_major_matrix_cmplx( double complex* matrix,
							  		   unsigned int num_rows,
							  		   unsigned int num_cols )
{
	unsigned int i;
	double complex* copy;
	double complex column[MAX_POLYNOMIAL_ORDER+1];

	//===First Copy The Matrix===//
	copy = malloc(num_rows * num_cols * sizeof(double complex));
	copy_array_cmplx(matrix, num_rows*num_cols, copy);

	//So basically, I just need to grab each column and place them next to each other
	for (i=0; i<num_cols; i++){
		copy_matrix_column_cmplx(copy, num_rows, num_cols, i, column);
		copy_array_cmplx(column, num_rows, matrix + i*num_rows);
	}

	//===Clean Up===//
	free(copy);

	return;
}

void row_to_column_major_matrix_dbl( double* matrix,
							  		 unsigned int num_rows,
							  		 unsigned int num_cols )
{
	unsigned int i;
	double* copy;
	double column[MAX_POLYNOMIAL_ORDER+1];

	//===First Copy The Matrix===//
	copy = malloc(num_rows * num_cols * sizeof(double));
	copy_array_dbl(matrix, num_rows*num_cols, copy);

	//So basically, I just need to grab each column and place them next to each other
	for (i=0; i<num_cols; i++){
		copy_matrix_column_dbl(copy, num_rows, num_cols, i, column);
		copy_array_dbl(column, num_rows, matrix + i*num_rows);
	}

	//===Clean Up===//
	free(copy);

	return;
}

void column_to_row_major_matrix_cmplx( double complex* matrix,
							  		   unsigned int num_rows,
							  		   unsigned int num_cols )
{
	unsigned int i;
	double complex* copy;
	double complex row[MAX_POLYNOMIAL_ORDER+1];

	//===First Copy The Matrix===//
	copy = malloc(num_rows * num_cols  * sizeof(double complex));
	copy_array_cmplx(matrix, num_rows*num_cols, copy);

	//So basically, I just need to grab each column (now a row) and place them in the columns of the matrix
	for (i=0; i<num_cols; i++){
		copy_matrix_row_cmplx(copy, num_rows, i, row);
		replace_matrix_column_cmplx(matrix, num_rows, num_cols, row, i);
	}

	//===Clean Up===//
	free(copy);

	return;
}


void column_to_row_major_matrix_dbl( double* matrix,
							  		 unsigned int num_rows,
							  		 unsigned int num_cols )
{
	unsigned int i;
	double* copy;
	double row[MAX_POLYNOMIAL_ORDER+1];

	//===First Copy The Matrix===//
	copy = malloc(num_rows * num_cols  * sizeof(double));
	copy_array_dbl(matrix, num_rows*num_cols, copy);

	//So basically, I just need to grab each column (now a row) and place them in the columns of the matrix
	for (i=0; i<num_cols; i++){
		copy_matrix_row_dbl(copy, num_rows, i, row);
		replace_matrix_column_dbl(matrix, num_rows, num_cols, row, i);
	}

	//===Clean Up===//
	free(copy);

	return;
}

//================================================================================================//
//======================================OTHER ROUTINES============================================//
//================================================================================================//

void make_korobov_sequence( unsigned int base,
						    unsigned int num_samples, 
							double* sequence )
{
	unsigned int i;

	sequence[0] = fmod(base * 1.0/(double)num_samples, 1.0);
	for (i=1; i<num_samples; i++){
		sequence[i] = fmod(sequence[i-1] + 1.0/(double)num_samples, 1.0);
	}	
	for (i=0; i<num_samples; i++){
		if (sequence[i] < 0.5){
			sequence[i] *= 2.0;
		}
		else{
			sequence[i] = 2.0 * (1.0 - sequence[i]);
		}
	}

	return;
}

