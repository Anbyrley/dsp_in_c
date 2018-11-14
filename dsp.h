/** @file dsp.h
*   @brief Contains functions for working with FFT.
*
*
*  @author Alex N. Byrley (anbyrley)
*  @date September 2015
*  @bug No known bugs
*/

#ifndef DSP_H
#define DSP_H

//================================================================================================//
//===================================STANDARD INCLUDES============================================//
//================================================================================================//

#include <math.h>
#include <complex.h>


//================================================================================================//
//===================================LOCAL INCLUDES============================================//
//================================================================================================//

#include "macros.h"
#include "helper.h"
#include "linear_algebra.h"
#include "filter.h"



//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//


//------------------------------------------------------------------------------------------------//
//=====================================FFT FUNCTIONS==============================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function runs an FFT for a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void fft_dbl(double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function runs an IFFT to a double precision array.
*
* @param[in] double complx* data
* @param[in] unsigned int length
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void ifft_dbl(double complex*,int,double*);


//================================================================================================//
/**
* @brief This function runs an FFT for a double complex array.
*
* @param[in] double complex* data
* @param[in] int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void fft_cmplx(double complex*,int,double complex*);


//================================================================================================//
/**
* @brief This function runs an IFFT for a double complex array.
*
* @param[in] double complex* data
* @param[in] int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void ifft_cmplx(double complex*,int,double complex*);


//================================================================================================//
/**
* @brief This function runs a normalized FFT for a double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void fft_normalized_dbl(double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function runs a normalized FFT for a large double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void fft_normalized_malloc_dbl(double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function runs an FFT for a double precision array using malloc.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void fft_malloc_dbl(double*, unsigned int, double complex*);


//================================================================================================//
/**
* @brief This function runs an FFT for a double complex array using a malloc operation.
*
* @param[in] double complex* data
* @param[in] int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void fft_malloc_cmplx(double complex*,unsigned int, double complex*);


//================================================================================================//
/**
* @brief This function runs an IFFT for a double complex array using a malloc operation.
*
* @param[in] double complex* data
* @param[in] int length
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void ifft_malloc_cmplx(double complex*, int, double complex*);


//================================================================================================//
/**
* @brief This function computes a Chirp IFFT for a double precision array using a malloc operation.
*
* @param[in] double* data
* @param[in] unsigned int signal_length
* @param[in] double reference_frequency
* @param[in] double frequency_spacing
* @param[in] unsigned int num_frequency_nodes
* @param[out] double complex* chirp_fft
*
* @return NONE
*/
//================================================================================================//
void chirp_fft_dbl(double*,unsigned int,double,double,unsigned int,double complex*);


//------------------------------------------------------------------------------------------------//
//==================================FFT SUPPORT FUNCTIONS=========================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function converts the output of a real fft to that of a complex fft.
*
* @param[in] double complex* real_fft
* @param[in] unsigned int fft_length
* @param[out] double complex* complex_fft
*
* @return NONE
*/
//================================================================================================//
void convert_real_fft_to_cmplx_malloc(double complex*, unsigned int, double complex*);


//================================================================================================//
/**
* @brief This function returns the frequencies associated with the fft.
*
* This function returns an array [ [0, Fs/2), [-Fs/2, 0) ].
*
* @param[in] unsigned int length
* @param[in] double sampling_frequency
* @param[out] double* frequencies
*
* @return NONE
*/
//================================================================================================//
void frequencies(unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function returns the nonnegative frequencies associated with the fft.
*
* This function returns an array [0, Fs/2).
*
* @param[in] unsigned int num_nonnegative_nodes
* @param[in] double sampling_rate
* @param[out] double* nonnegative_frequencies
*
* @return NONE
*/
//================================================================================================//
void get_nonnegative_frequencies(unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function normalizes an array of digital frequencies to the range [-pi, pi)
*
* @param[in,out] double* frequencies
* @param[in] unsigned int length
* @param[in] double sampling_rate
*
* @return NONE
*/
//================================================================================================//
void normalize_digital_frequencies(double*,unsigned int,double);


//------------------------------------------------------------------------------------------------//
//=====================================DFT FUNCTIONS==============================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function computes a DFT given double precision data.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* data_dft
*
* @return NONE
*/
//================================================================================================//
void dft_dbl(double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function computes a double precision IDFT.
*
* @param[in] double complex* data_dft
* @param[in] unsigned int length
* @param[out] double* data
*
* @return NONE
*/
//================================================================================================//
void idft_dbl(double complex*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes the dft value of a given frequency via goertzel recursion.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[in] double freq
*
* @return double dft
*/
//================================================================================================//
double complex run_goertzel_filter_dft(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function computes the dft value of a given frequency via goertzel recursion for cmplx.
*
* @param[in] double complex* data
* @param[in] unsigned int length
* @param[in] double freq
*
* @return double dft
*/
//================================================================================================//
double complex run_goertzel_filter_dft_cmplx(double complex*,unsigned int,double);


//================================================================================================//
/**
* @brief This function computes the dft of a double precision signal via goertzel recursion.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[in] unsigned int num_frequency_nodes
* @param[out] double complex* dft
*
* @return NONE
*/
//================================================================================================//
void goertzel_dft_dbl(double*,unsigned int,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function computes the power of a given frequency via goertzel recursion.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[in] double freq
*
* @return double power
*/
//================================================================================================//
double run_goertzel_filter_psd(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function computes the power spectrum via goertzel recursion.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[in] unsigned int num_frequency_nodes
* @param[in] double freq
* @param[out] double* psd
*
* @return NONE
*/
//================================================================================================//
void goertzel_psd_dbl(double*, unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes a NUDFT given double precision data and nodes [-0.5,0.5).
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[in] double* freq_nodes
* @param[in] unsigned int num_freq_nodes
* @param[out] double complex* data_dft
*
* @return NONE
*/
//================================================================================================//
void nudft_dbl(double*,unsigned int,double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function computes an Inverse NUDFT given double precision data and nodes [-0.5,0.5).
*
* @param[in] double complex* data
* @param[in] double* freq_nodes
* @param[in] unsigned int num_freq_nodes
* @param[in] unsigned int out_length
* @param[out] double* out
*
* @return NONE
*/
//================================================================================================//
void nuidft_dbl(double complex*,double*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes a NUDFT given double complex data and nodes [-0.5,0.5).
*
* @param[in] double complex* data
* @param[in] unsigned int length
* @param[in] double* freq_nodes
* @param[in] unsigned int num_freq_nodes
* @param[out] double complex* data_dft
*
* @return NONE
*/
//================================================================================================//
void nudft_cmplx(double complex*,unsigned int,double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function computes the nudft of a double precision signal via goertzel recursion.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[in] double* freq_nodes
* @param[in] unsigned int num_frequency_nodes
* @param[out] double complex* dft
*
* @return NONE
*/
//================================================================================================//
void goertzel_nudft_dbl(double*,unsigned int,double*,unsigned int,double complex*);





//------------------------------------------------------------------------------------------------//
//====================================OTHER TRANSFORMS============================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function runs calculates the Hilbert Transform of a double precision array.
*
* @param[in] double* data
* @param[in] int length
* @param[out] double* analytic
*
* @return NONE
*/
//================================================================================================//
void hilbert_transform_dbl(double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes the analytic signal of passed in data
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* analytic
*
* @return NONE
*/
//================================================================================================//
void compute_analytic_signal_dbl(double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function runs calculates the Form II DCT directly for a double precision array.
*
* @param[in] double* data
* @param[in] int length
* @param[out] double* dct
*
* @return NONE
*/
//================================================================================================//
void dct_dbl(double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function runs calculates the Form I DST directly for a double precision array.
*
* @param[in] double* data
* @param[in] int length
* @param[out] double* dct
*
* @return NONE
*/
//================================================================================================//
void dst_dbl(double*,unsigned int,double*);


//------------------------------------------------------------------------------------------------//
//======================================Z-TRANSFORMS==============================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function runs a z-transform for a double precision array.
*
* @param[in] double* data
* @param[in] int length
* @param[in] double radius
* @param[out] double complex* z_transform
*
* @return NONE
*/
//================================================================================================//
void z_transform_dbl(double*,unsigned int,double,double complex*);


//================================================================================================//
/**
* @brief This function runs a z-transform for a double complex array.
*
* @param[in] double complex* data
* @param[in] int length
* @param[in] double radius
* @param[out] double complex* z_transform
*
* @return NONE
*/
//================================================================================================//
void z_transform_cmplx(double complex*,unsigned int,double,double complex*);


//================================================================================================//
/**
* @brief This function runs an inverse z-transform for a double complex array.
*
* @param[in] double complex* data
* @param[in] int length
* @param[in] double radius
* @param[out] double* z_transform
*
* @return NONE
*/
//================================================================================================//
void inverse_z_transform_dbl(double complex*,unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function takes the derivative of the z-transform for a double precision array.
*
* @param[in] double* data
* @param[in] int length
* @param[in] double radius
* @param[out] double complex* z_transform_derivative
*
* @return NONE
*/
//================================================================================================//
void z_transform_derivative_dbl(double*,unsigned int,double,double complex*);


//================================================================================================//
/**
* @brief This function takes the logarithmic derivative for a double precision array.
*
* @param[in] double* data
* @param[in] int length
* @param[in] double radius
* @param[out] double complex* logarithmic_derivative
*
* @return NONE
*/
//================================================================================================//
void logarithmic_derivative_dbl(double*,unsigned int,double,double complex*);


//================================================================================================//
/**
* @brief This calculates the group delay around a circle of radius r.
*
* @param[in] double* data
* @param[in] int length
* @param[in] double radius
* @param[out] double* group_delay
*
* @return NONE
*/
//================================================================================================//
void calculate_ztransform_group_delay_dbl(double*,unsigned int, double, double*);


//------------------------------------------------------------------------------------------------//
//====================================PHASE FUNCTIONS=============================================//
//------------------------------------------------------------------------------------------------//


//================================================================================================//
/**
* @brief This function computes the phase response of the data.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* phase
*
* @return NONE
*/
//================================================================================================//
void compute_phase_response_dbl(double*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function unwraps angles by changing absolute jumps greater than cutoff to their 
* TWO_PI complement.
*
* @param[in,out] double* phase
* @param[in] unsigned int length
* @param[in] double cutoff
*
* @return NONE
*/
//================================================================================================//
void unwrap_phase_angles(double*,unsigned int,double);


//================================================================================================//
/**
* @brief This function computes the angle response of the data.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double* angles
*
* @return NONE
*/
//================================================================================================//
void compute_angle_response_dbl(double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function delays the signal via a phase shift.
*
* @param[in,out] double* data
* @param[in] unsigned int length
* @param[in] double delay_seconds
* @param[in] double sampling_rate
*
* @return NONE
*/
//================================================================================================//
void phase_shift_dbl(double*,unsigned int,double,double);


//================================================================================================//
/**
* @brief This function delays the signal via a phase shift for long signals using malloc.
*
* @param[in,out] double* data
* @param[in] unsigned int length
* @param[in] double delay_seconds
* @param[in] double sampling_rate
*
* @return NONE
*/
//================================================================================================//
void phase_shift_malloc_dbl(double*,unsigned int,double,double);


//------------------------------------------------------------------------------------------------//
//=====================================MAGNITUDE FUNCTIONS========================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function computes the magnitude response of the data.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double* magnitude
*
* @return NONE
*/
//================================================================================================//
void compute_magnitude_response_dbl(double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function returns the log magnitude spectrum of the fourier array.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double* log_magnitude
*
* @return NONE
*/
//================================================================================================//
void compute_log_magnitude_response_dbl(double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function returns the sample power spectral density.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* psd_cmplx
* @param[out] double* psd_dbl
*
* @return NONE
*/
//================================================================================================//
void compute_power_spectrum_dbl(double*,unsigned int,double complex*,double*);


//================================================================================================//
/**
* @brief This function returns the sample power spectral density for a large double precision array.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[out] double complex* psd_cmplx
* @param[out] double* psd_dbl
*
* @return NONE
*/
//================================================================================================//
void compute_power_spectrum_malloc_dbl(double*,unsigned int,double complex*,double*);


//------------------------------------------------------------------------------------------------//
//====================================GROUP DELAY FUNCTIONS=======================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function computes the group delay of sampled data.
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[in] double sampling_frequency
* @param[out] double* group_delay
*
* @return NONE
*/
//================================================================================================//
void calculate_group_delay_dbl(double*,unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function computes the group delay of sampled data at given frequency nodes in [-.5,.5]
*
* @param[in] double* data
* @param[in] unsigned int length
* @param[in] double sampling_frequency
* @param[in] double* frequency_nodes
* @param[in] unsigned int num_freq_nodes
* @param[out] double* group_delay
*
* @return NONE
*/
//================================================================================================//
void calculate_warped_group_delay_dbl(double*,unsigned int,double,double*,unsigned int,double*);


//------------------------------------------------------------------------------------------------//
//=================================CONVOLUTION FUNCTIONS==========================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function convolves two arrays of type double.
*
* @param[in] double* array1
* @param[in] unsigned int length1
* @param[in] double* array2
* @param[in] unsigned int length2
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void convolve_dbl(double*, unsigned int, double*, unsigned int, double*);


//================================================================================================//
/**
* @brief This function convolves two arrays of type double.
*
* @param[in] double complex* array1
* @param[in] unsigned int length1
* @param[in] double complex* array2
* @param[in] unsigned int length2
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void convolve_cmplx(double complex*, unsigned int, double complex*, unsigned int, double complex*);


//================================================================================================//
/**
* @brief This function convolves two large arrays of type double using mallocs.
*
* @param[in] double* array1
* @param[in] unsigned int length1
* @param[in] double* array2
* @param[in] unsigned int length2
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void convolve_fft_malloc_dbl(double*,unsigned int,double*,unsigned int, double*);


//================================================================================================//
/**
* @brief This function convolves two large arrays of type double complex using mallocs.
*
* @param[in] double complex* array1
* @param[in] unsigned int length1
* @param[in] double complex* array2
* @param[in] unsigned int length2
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void convolve_fft_malloc_cmplx(double complex*,unsigned int,double complex*,unsigned int,double complex*);


//================================================================================================//
/**
* @brief This function returns the normalized cross correlation of two arrays.
*
* @param[in] double* array1
* @param[in] unsigned int length1
* @param[in] double* array2
* @param[in] unsigned int length2
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void normalized_correlation_dbl(double*,unsigned int,double*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function returns the cross correlation of two double complex arrays.
*
* @param[in] double complex* array1
* @param[in] unsigned int length1
* @param[in] double complex* array2
* @param[in] unsigned int length2
* @param[out] double complex* result
*
* @return NONE
*/
//================================================================================================//
void correlation_cmplx(double complex*,unsigned int,double complex*,unsigned int,double complex*);


//------------------------------------------------------------------------------------------------//
//====================================RESAMPLING FUNCTIONS========================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function computes the lanczos kernel for a given value.
*
* @param[in] double x
* @param[in] double a
*
* @return L(x, a)
*/
//================================================================================================//
double compute_lanczos_kernel_value_dbl(double,double);


//================================================================================================//
/**
* @brief This function resamples a double precision signal to a new scale using the lanczos kernel.
*
* @param[in] double* values
* @param[in] double* original_scale
* @param[in] unsigned int original_length
* @param[in] double* new_scale
* @param[in] unsigned int new_length
* @param[in] double kernel_size
* @param[out] double* resampled
*
* @return NONE
*/
//================================================================================================//
void resample_lanczos_dbl(double*,double*,unsigned int,double*,unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function resamples a double complex signal to a new scale using the sinc kernel.
*
* @param[in] double complex* values
* @param[in] double* original_scale
* @param[in] unsigned int original_length
* @param[in] double* resampled_scale
* @param[in] unsigned int resampled_length
* @param[out] double complex* resampled
*
* @return NONE
*/
//================================================================================================//
void resample_sinc_cmplx(double complex*,double*,unsigned int,double*,unsigned int,double complex*);


//------------------------------------------------------------------------------------------------//
//======================================MIXING FUNCTIONS==========================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function mixes two signals with a certain SNR.
*
* NOTE: SNR gains up first signal with respect to the other.
*
* @param[in] double* signal1
* @param[in] unsigned int signal1_length
* @param[in] double* signal2
* @param[in] unsigned int signal2_length
* @param[in] double SNR
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void mix_signals_dbl(double*,unsigned int,double*,unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function mixes multiple signals according to their given power.
*
* NOTE: Each signal must have the same length, and powers must be in dB.
*
* @param[in] double** signals
* @param[in] unsigned int num_signals
* @param[in] unsigned int signal_length
* @param[in] double* signal_powers_dB
* @param[out] double* result
*
* @return NONE
*/
//================================================================================================//
void mix_multiple_signals_dbl(double**,unsigned int,unsigned int,double*,double*);


//------------------------------------------------------------------------------------------------//
//=================================INPUT GENERATION FUNCTIONS=====================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function generates a bandlimited signal with specified frequencies and amplitudes.
*
* @param[in] double* frequencies
* @param[in] double* amplitudes
* @param[in] double* phase
* @param[in] unsigned int num_frequencies
* @param[in] unsigned int signal_length
* @param[in] double sampling_rate
* @param[out] double* wave
*
* @return NONE
*/
//================================================================================================//
void generate_wave_dbl(double*,double*,double*,unsigned int,unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function generates a bandlimited harmonic signal with specified f0 and amplitudes.
*
* @param[in] double fundamental
* @param[in] double* amplitudes
* @param[in] double* phase
* @param[in] unsigned int num_harmonics
* @param[in] unsigned int signal_length
* @param[in] double sampling_rate
* @param[out] double* wave
*
* @return NONE
*/
//================================================================================================//
void generate_harmonic_wave_dbl(double,double*,double*,unsigned int,unsigned int,double,double*);


//================================================================================================//
/**
* @brief This function generates a bandlimited carrier at a specified frequency.
*
* @param[in] double frequency
* @param[in] unsigned int signal_length
* @param[in] double sampling_rate
* @param[out] double complex* carrier
*
* @return NONE
*/
//================================================================================================//
void generate_carrier_cmplx(double,unsigned int,double,double complex*);


//================================================================================================//
/**
* @brief This function generates an impulse.
*
* @param[in] unsigned int length
* @param[out] double* impulse
*
* @return NONE
*/
//================================================================================================//
void generate_impulse_dbl(unsigned int,double*);


//================================================================================================//
/**
* @brief This function generates a unit step from start.
*
* @param[in] unsigned int start
* @param[in] unsigned int length
* @param[out] double* step
*
* @return NONE
*/
//================================================================================================//
void generate_unit_step_dbl(unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function generates a unit pulse from start to end.
*
* @param[in] unsigned int start
* @param[in] unsigned int end
* @param[in] unsigned int length
* @param[out] double* pulse
*
* @return NONE
*/
//================================================================================================//
void generate_unit_pulse_dbl(unsigned int,unsigned int,unsigned int,double*);


//------------------------------------------------------------------------------------------------//
//===================================MEL FREQUENCY FUNCTIONS======================================//
//------------------------------------------------------------------------------------------------//

//================================================================================================//
/**
* @brief This function converts a value in hertz to that in mel frequency.
*
* @param[in] double hertz
*
* @return double mel
*/
//================================================================================================//
double hertz_to_mel(double);


//================================================================================================//
/**
* @brief This function converts a value in mel freqwuency to that in hertz.
*
* @param[in] double mel
*
* @return double hertz
*/
//================================================================================================//
double mel_to_hertz(double);


//================================================================================================//
/**
* @brief This function creates magnitude spetrum mel filterbanks.
*
* @param[in] unsigned int num_frequency_nodes
* @param[in] double sampling_rate
* @param[in] unsigned int num_mel_channels
* @param[out] double* mel_filterbanks
*
* @return NONE
*/
//================================================================================================//
void make_mel_frequency_filterbanks(unsigned int,double,unsigned int,double*);


//================================================================================================//
/**
* @brief This function computes the MFCCs from a given data array.
*
* @param[in] double* data 
* @param[in] unsigned int length
* @param[in] double sampling_rate
* @param[in] unsigned int num_mel_channels
* @param[out] double* mfccs
*
* @return NONE
*/
//================================================================================================//
void compute_mel_frequency_cepstral_coefficients_dbl(double*,unsigned int,double,unsigned int,double*);


#endif //DSP_H//
