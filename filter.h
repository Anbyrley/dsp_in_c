/** @file filter.h
*   @brief Contains functions for filtering.
*
*
*  @author Alex N. Byrley (anbyrley)
*  @date November 2015
*  @bug No known bugs
*/

#ifndef FILTER_H
#define FILTER_H

//================================================================================================//
//===================================STANDARD INCLUDES============================================//
//================================================================================================//

#include "macros.h"
#include "helper.h"
#include "ring.h"
#include "dsp.h"
#include "polynomial.h"
#include "linear_algebra.h"

//================================================================================================//
//======================================DATA STRUCTURES===========================================//
//================================================================================================//

//================================================================================================//
/** @enum resampler_type
*   @brief This enum tells whether the resampler does: downsample, upsample, rational or arbitrary.
*/
//================================================================================================//
typedef enum{
	DOWNSAMPLE, 
	UPSAMPLE,
	RATIONAL,
} resampler_type;


//================================================================================================//
/** @enum filter_type
*   @brief This enum tells whether the filter performs: lowpass, highpass, bandpass or bandstop.
*/
//================================================================================================//
typedef enum{
	LOWPASS, 
	HIGHPASS,
	BANDPASS,
	BANDSTOP,
} filter_type;


//================================================================================================//
/** @struct filter_t
*   @brief This structure is the typedef for the filter_t object.
*/
//================================================================================================//
typedef struct filter_s filter_t;
typedef struct filter_s{
	unsigned int numerator_length;
	unsigned int denominator_length;
	unsigned int internal_length;
	double numerator[MAX_POLYNOMIAL_ORDER+1];
	double denominator[MAX_POLYNOMIAL_ORDER+1];
	double internal[2*(MAX_POLYNOMIAL_ORDER+1)];
} filter_t;

//================================================================================================//
/** @struct resampler_initializer_t
*   @brief This structure is the typedef for the resampler_initializer_t object.
*/
//================================================================================================//
typedef struct resampler_initializer_s resampler_initializer_t;
typedef struct resampler_initializer_s{
	resampler_type type;
	unsigned int block_size;
	unsigned int downsampling_factor;
	unsigned int upsampling_factor;
} resampler_initializer_t;

//================================================================================================//
/** @struct resampler_t
*   @brief This structure is the typedef for the resampler_t object.
*/
//================================================================================================//
typedef struct resampler_s resampler_t;
typedef struct resampler_s{
	//unsigned int output_enqueue_count;	//when % downsampling_factor == 0, enqueue into the output queue
	resampler_type type;
	unsigned int block_size;
	unsigned int downsampling_factor;
	unsigned int upsampling_factor;
	unsigned int input_queue_length;
	unsigned int output_queue_length;
	fifo_queue_t input_queue;
	fifo_queue_t output_queue;
	fifo_queue_t tdl_queue;
	filter_t resampling_filter;
	filter_t polyphase_resampling_filters[MAX_NUM_CMFB_FILTERS];
} resampler_t;



//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//

//------------------------------------------------------------------------------------------------//
//==================================ANALOG DESIGN FUNCTIONS=======================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function converts a digital frequency from [0,M_PI] its analog equivalent.
*
* Note: Digital frequencies must be pre-warped before passing them as parameters to analog designs.
*	
* @param[in] double digital_frequency
*
* @return double analog_frequency
*/
//================================================================================================//
double prewarp_frequency_dbl(double);


//================================================================================================//
/**
* @brief This function calculates the required order for a butterworth filter given the specs.
*
* @param[in] double passband_edge
* @param[in] double stopband_edge
* @param[in] double passband_gain
* @param[in] double stopband_gain
*
* @return unsigned int order
*/
//================================================================================================//
unsigned int calculate_butterworth_order(double,double,double,double);


//================================================================================================//
/**
* @brief This function designs an IIR Analog Butterworth Filter.
*
* Note: Digital frequencies [0,M_PI] must be pre-warped before use.
*
* @param[in] unsigned int order
* @param[in] double passband_edge
* @param[out] double* numerator
* @param[out] double* denominator
*
* @return NONE
*/
//================================================================================================//
void butterworth_design(unsigned int,double,double*,double*);


//================================================================================================//
/**
* @brief This function calculates the required order for an chebyshev typeI filter given the specs.
*
* @param[in] double passband_ripple_dB
* @param[in] double passband_edge
* @param[in] double stopband_edge
*
* @return unsigned int order
*/
//================================================================================================//
unsigned int calculate_chebyschev_order(double,double,double);


//================================================================================================//
/**
* @brief This function designs an IIR Analog Chebyshev Type I Filter.
*
* Note: Digital frequencies [0,M_PI] must be prewarped before use. There is a tradeoff between
* the amount of passband ripple and the rolloff of the filter. High allowed ripple means a 
* steep rolloff, while tight ripple requirements result in a filter with a shallow roll off.
*
* @param[in] unsigned int order
* @param[in] double passband_ripple_dB
* @param[in] double passband_edge
* @param[out] double* numerator
* @param[out] double* denominator
*
* @return NONE
*/
//================================================================================================//
void chebyshev_design(unsigned int,double,double,double*,double*);


//================================================================================================//
/**
* @brief This function converts an analog filter to a digital filter in place.
*
* Algorithm Taken From: Digital Filter Designer's Handbook by C. Britton Rorabaugh pg. 297 
*
* @param[in,out] filter_t* self
*
* @return NONE
*/
//================================================================================================//
void bilinear_transform( filter_t* self );


//================================================================================================//
/**
* @brief This function initializes an ideal lowpass filter.
*
* @param[in] double fc
* @param[in] double Fs
* @param[in] int order
* @param[in,out] double* ideal
*
* @return NONE
*/
//================================================================================================//
void initialize_lowpass_ideal(double,double,int,double*);


//================================================================================================//
/**
* @brief This function initializes an ideal bandpass filter.
*
* @param[in] double f_low
* @param[in] double f_high
* @param[in] double Fs
* @param[in] int order
* @param[in,out] double* ideal
*
* @return NONE
*/
//================================================================================================//
void initialize_bandpass_ideal(double,double,double,int,double*);


//================================================================================================//
/**
* @brief This function initializes an ideal highpass filter.
*
* @param[in] double fc
* @param[in] double Fs
* @param[in] int order
* @param[in,out] double* ideal
*
* @return NONE
*/
//================================================================================================//
void initialize_highpass_ideal(double,double,int,double*);


//================================================================================================//
/**
* @brief This function initializes an ideal bandstop filter.
*
* @param[in] double f_low
* @param[in] double f_high
* @param[in] double Fs
* @param[in] int order
* @param[in,out] double* ideal
*
* @return NONE
*/
//================================================================================================//
void initialize_bandstop_ideal(double,double,double,int,double*);


//------------------------------------------------------------------------------------------------//
//==================================WINDOW DESIGN FUNCTIONS=======================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function initializes a cosine window.
*
* @param[in] unsigned int length
* @param[out] double* window
*
* @return NONE
*/
//================================================================================================//
void initialize_cosine_window(unsigned int, double*);


//================================================================================================//
/**
* @brief This function initializes a hanning window.
*
* @param[in] unsigned int length
* @param[out] double* window
*
* @return NONE
*/
//================================================================================================//
void initialize_hanning_window(unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a hamming window.
*
* @param[in] unsigned int length
* @param[out] double* window
*
* @return NONE
*/
//================================================================================================//
void initialize_hamming_window(unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a blackman window.
*
* @param[in] unsigned int length
* @param[out] double* window
*
* @return NONE
*/
//================================================================================================//
void initialize_blackman_window(unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a rectangular window.
*
* @param[in] unsigned int length
* @param[out] double* window
*
* @return NONE
*/
//================================================================================================//
void initialize_rectangular_window(unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a chebyshev window in the time domain.
*
* Time domain implementation was found in: Antoniou, A., "Digital Filters", McGraw-Hill, 2000.
*
* @param[in] double sidelobe_atten (dB)
* @param[in] unsigned int length
* @param[out] double* window
*
* @return NONE
*/
//================================================================================================//
void initialize_chebyshev_window(double,unsigned int,double*);


//================================================================================================//
/**
* @brief This function returns the required beta parameter given filter specifications.
*
* @param[in] double stopband_attenuation_dB
*
* @return double beta
*/
//================================================================================================//
double get_kaiser_beta(double);


//================================================================================================//
/**
* @brief This function returns the required order of a kaiser window given filter specifications.
*
* @param[in] double deltaw
* @param[in] double stopband_attenuation_dB
* @param[in] unsigned int num_subbands;
*
* @return unsigned int N
*/
//================================================================================================//
unsigned int get_kaiser_filter_length(double,double,unsigned int);


//================================================================================================//
/**
* @brief This function initializes a kaiser window in the time domain.
*
* @param[in] double beta
* @param[in] unsigned int length
* @param[out] double* window
*
* @return NONE
*/
//================================================================================================//
void initialize_kaiser_window(double,unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a prolate window design.
*
* @param[in] double ws
* @param[in] unsigned int N
* @param[in,out] double* window
*
* @return NONE
*/
//================================================================================================//
void prolate_window_design(double,unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a Fyrley window.
*
* @param[in] unsigned int length
* @param[out] double* window
*
* @return NONE
*/
//================================================================================================//
void initialize_fyrley_window(unsigned int, double*);


//------------------------------------------------------------------------------------------------//
//===================================FILTER DESIGN FUNCTIONS=======================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function initializes a filter via the window method.
*
* The window to be used must be passed in as a parameter.
*
* @param[in] filter_type
* @param[in] double fc1
* @param[in] double fc2
* @param[in] double Fs
* @param[in] unsigned int length
* @param[in] double* window
* @param[in,out] double* filter
*
* @return NONE
*/
//================================================================================================//
void window_filter_design(filter_type,double,double,double,int,double*,double*);


//================================================================================================//
/**
* @brief This function initializes a notch filter.
*
* @param[in] double* notch_frequencies
* @param[in] unsigned int num_notches
* @param[in] double sampling_rate
* @param[in] double notch_width
* @param[out] unsigned int *filter_length
* @param[out] double* numerator
* @param[out] double* denominator
*
* @return NONE
*/
//================================================================================================//
void notch_filter_design(double*,unsigned int,double,double,unsigned int*,double*,double*);


//================================================================================================//
/**
* @brief This function initializes an eigenfilter.
*
*  NOTE: wp and ws must be in the range [0, M_PI]. Since freq axis is assumed [-M_PI, M_PI].
*  Also select alpha (0, 1.0) since 0 or 1 weights (stopband,passband) with nothing causing failure.
*
*  NOTE: Taken from Eigenfilters: A New Approach to Least-Squares FIR Filter Design and Applications
*  Including Nyquist Filters by P. P. Vaidyanathan and T.Q. Nguyen.
*
* @param[in] double wp
* @param[in] double ws
* @param[in] unsigned int N
* @param[in] double alpha
* @param[in,out] double* filter
*
* @return NONE
*/
//================================================================================================//
void eigenfilter_design(double,double,unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function initializes a filter using the minimax error criteria.
*
* NOTE: wp and ws must be in the range [0, M_PI]. Since freq axis is assumed [-M_PI, M_PI].
*
* NOTE: Taken from EL 713 Lecture Notes by I. Selesnick
*
*
* @param[in] double wp
* @param[in] double ws
* @param[in] double stopband_attn (dB)
* @param[in] double passband_ripple (dB)
* @param[in,out] unsigned int* filter_length
* @param[in,out] double* filter
*
* @return NONE
*/
//================================================================================================//
void minimax_filter_design(double,double,double,double,unsigned int*,double*);


//================================================================================================//
/**
* @brief This function initializes a filter using the L2 error criteria.
*
* NOTE: wp and ws must be in the range [0, M_PI]. Since freq axis is assumed [-M_PI, M_PI].
*
* NOTE: Taken from EL 713 Lecture Notes by I. Selesnick
*
* @param[in] double wp
* @param[in] double ws
* @param[in] double passband_weight
* @param[in] double stopband_weight
* @param[in] unsigned int filter_length
* @param[in,out] double* filter
*
* @return NONE
*/
//================================================================================================//
void least_squares_filter_design(double,double,double,double,unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a fractional delay filter.
*
* @param[in] double delay (seconds)
* @param[in] double sampling_rate
* @param[in] unsigned int filter_length
* @param[out] double* filter
*
* @return NONE
*/
//================================================================================================//
void fractional_delay_filter_design(double,double,unsigned int,double*);


//================================================================================================//
/**
* @brief This function initializes a preemphasis delay filter. (1 - alpha)x
*
* @param[in] double alpha
* @param[out] double* filter
*
* @return NONE
*/
//================================================================================================//
void preemphasis_filter_design(double,double*);


//------------------------------------------------------------------------------------------------//
//======================================FILTER FUNCTIONS==========================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function initializes the filter_t object.
*
* @param[in,out] filter_t* self
* @param[in] double* numerator
* @param[in] unsigned int numerator_length
* @param[in] double* denominator
* @param[in] unsigned int denominator_length
*
* @return NONE
*/
//================================================================================================//
void initialize_filter(filter_t*,double*,unsigned int,double*,unsigned int);


//================================================================================================//
/**
* @brief This function resets the filter_t object.
*
* @param[in,out] filter_t* self
*
* @return NONE
*/
//================================================================================================//
void reset_filter(filter_t*);


//================================================================================================//
/**
* @brief This function copies the filter_t object.
*
* @param[in] filter_t* original
* @param[out] filter_t* copy
*
* @return NONE
*/
//================================================================================================//
void copy_filter(filter_t*,filter_t*);


//================================================================================================//
/**
* @brief This function runs input through the filter_t object.
*
* @param[in,out] filter_t* self
* @param[in] double* input
* @param[in] unsigned int input_length
* @param[in,out] double* output
*
* @return NONE
*/
//================================================================================================//
void run_filter(filter_t*,double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function runs input through the filter_t object, reverses, filters, and reverses.
*
* NOTE: Input Length should include zeros for filter's group delay.
*
* @param[in,out] filter_t* self
* @param[in] double* input
* @param[in] unsigned int input_length
* @param[in,out] double* output
*
* @return NONE
*/
//================================================================================================//
void run_filter_zero_phase(filter_t*,double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes the impulse response of the filter_t object.
*
* @param[in,out] filter_t* self
* @param[in] unsigned int impulse_length
* @param[out] double* impulse_response
*
* @return NONE
*/
//================================================================================================//
void compute_impulse_response(filter_t*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes the frequency response of a digital filter from [0, M_PI].
*
* This function evaluates a rational function aka the transfer function around the unit circle.
*
* @param[in] filter_t* self
* @param[in] unsigned int response_length
* @param[out] double complex* frequency_response
*
* @return NONE
*/
//================================================================================================//
void compute_frequency_response(filter_t*,unsigned int,double complex*);


//------------------------------------------------------------------------------------------------//
//===================================RESAMPLING FUNCTIONS=========================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function initializes the resampler_t object.
*
* @param[in,out] resampler_t* self
* @param[in] unsigned int downsampling_factor
* @param[in] unsigned int upsampling_factor
* @param[in] unsigned int block_size
* @param[in] resampler_type type
*
* @return NONE
*/
//================================================================================================//
void initialize_resampler(resampler_t*,unsigned int,unsigned int,unsigned int,resampler_type);


//================================================================================================//
/**
* @brief This function enqueues data into the resampler_t object.
*
* NOTE: The data must the length of the input_queue.
*
* @param[in,out] resampler_t* self
* @param[in] double* data
*
* @return NONE
*/
//================================================================================================//
void enqueue_resampler_dbl(resampler_t*,double*);


//================================================================================================//
/**
* @brief This function dequeues data out of the resampler_t object.
*
* NOTE: The data must the length of the output_queue.
*
* @param[in,out] resampler_t* self
* @param[out] double* data
*
* @return NONE
*/
//================================================================================================//
void dequeue_resampler_dbl(resampler_t*,double*);


//================================================================================================//
/**
* @brief This function runs the resampler.
*
* @param[in,out] resampler_t* self
*
* @return NONE
*/
//================================================================================================//
void run_resampler(resampler_t*);


//================================================================================================//
/**
* @brief This function runs the resampler on an array.
*
* @param[in,out] resampler_t* self
* @param[in] double* original_signal
* @param[in] unsigned int original_signal_length
* @param[in] unsigned int downsampling_factor
* @param[in] unsigned int upsampling_factor
* @param[in] resampler_type type
* @param[out] double* resampled_signal
*
* @return NONE
*/
//================================================================================================//
void resample_array_dbl(resampler_t*,double*,unsigned int,unsigned int,
						unsigned int,resampler_type,double*);


//------------------------------------------------------------------------------------------------//
//===================================ALL PASS FUNCTIONS===========================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function delays a signal via a fractional delay filter.
*
* @param[in] double* signal
* @param[in] unsigned int signal_length
* @param[in] double delay
* @param[in] double sampling_rate
* @param[out] double* delayed
*
* @return NONE
*/
//================================================================================================//
void delay_signal_dbl(double*, unsigned int, double, double, double*);


//================================================================================================//
/**
* @brief This function boosts the low frequencies of a signal.
*
* @param[in] double* signal
* @param[in] unsigned int signal_length
* @param[in] double coeff
* @param[out] double* emphasized
*
* @return NONE
*/
//================================================================================================//
void low_frequency_emphasis_dbl(double*,unsigned int,double,double*);


//------------------------------------------------------------------------------------------------//
//====================================TEST FUNCTIONS==============================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function tests the impulse response of a filter and prints it to a file.
*
* @return NONE
*/
//================================================================================================//
void test_filter();


//================================================================================================//
/**
* @brief This function generates an AR process and prints it to a file.
*
* @return NONE
*/
//================================================================================================//
void generate_ar_process();


//================================================================================================//
/**
* @brief This function tests the fractional delay filter and prints results to files.
*
* @return NONE
*/
//================================================================================================//
void test_fractional_delay();


//================================================================================================//
/**
* @brief This function tests the impulse response of a chebyshev filter and prints it to a file.
*
* @return NONE
*/
//================================================================================================//
void test_chebyshev();


//================================================================================================//
/**
* @brief This function tests the impulse response of a buttworth filter and prints it to a file.
*
* @return NONE
*/
//================================================================================================//
void test_butterworth();


//================================================================================================//
/**
* @brief This function tests the resampler and shows the program flow for real-time use.
*
* @return NONE
*/
//================================================================================================//
void test_resampler();


//================================================================================================//
/**
* @brief This function prints the various filter types (low, high, band, notch) for viewing.
*
* @return NONE
*/
//================================================================================================//
void view_filter_types();


//================================================================================================//
/**
* @brief This function prints the impulse response of a notch filter.
*
* @return NONE
*/
//================================================================================================//
void test_notch_filter();


#endif //FILTER_H//
