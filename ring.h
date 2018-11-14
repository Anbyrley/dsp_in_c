/** @file ring.h
*   @brief Contains functions for the ring buffer.
*
*
*  @author Alex N. Byrley (anbyrley)
*  @date August 2016
*  @bug No known bugs
*/

#ifndef RING_H
#define RING_H
 


//================================================================================================//
//=====================================LOCAL INCLUDES=============================================//
//================================================================================================//

#include "macros.h"
#include "helper.h"

//================================================================================================//
//======================================DATA STRUCTURES===========================================//
//================================================================================================//

//================================================================================================//
/** @union fifo_data
*   @brief This union contains all needed data types for a fifo object.
*/
//================================================================================================//
union buffer_data{
	unsigned int uint[MAX_QUEUE_LENGTH];
	int sint[MAX_QUEUE_LENGTH];
	double dbl[MAX_QUEUE_LENGTH];
	double complex cmplx[MAX_QUEUE_LENGTH];
} buffer_data;


//================================================================================================//
/** @struct ring_buffer_t
*   @brief This structure is the typedef for the ring_buffer_t object.
*/
//================================================================================================//
typedef struct ring_buffer_s ring_buffer_t;
typedef struct ring_buffer_s{
	union buffer_data data;
	unsigned int full;
	unsigned int tail;
	unsigned int queue_length;
} ring_buffer_t;


//================================================================================================//
/** @struct fifo_queue_t
*   @brief This structure is the typedef for the fifo_queue_t object.
*/
//================================================================================================//
typedef struct fifo_queue_s fifo_queue_t;
typedef struct fifo_queue_s{
	union buffer_data data;
	unsigned int head;
	unsigned int tail;
	unsigned int queue_length;
	unsigned int head_wrap : 1;
	unsigned int tail_wrap : 1;
} fifo_queue_t;


//================================================================================================//
//==================================FUNCTION DECLARATIONS=========================================//
//================================================================================================//

//================================================================================================//
/**
* @brief This function initializes the ring buffer.
*
* @param[in,out] ring_buffer_t* self 
* @param[in] unsigned int queue_length
*
* @return NONE
*/
//================================================================================================//
void initialize_ring_buffer(ring_buffer_t*,unsigned int);


//================================================================================================//
/**
* @brief This function copies a ring buffer.
*
* @param[in] ring_buffer_t* original
* @param[out] ring_buffer_t* copy
*
* @return NONE
*/
//================================================================================================//
void copy_ring_buffer(ring_buffer_t*,ring_buffer_t*);


//================================================================================================//
/**
* @brief This function writes all the data in the buffer to an array and then zeros its own data.
*
* @param[in,out] ring_buffer_t* self 
* @param[out] double* out_data
*
* @return NONE
*/
//================================================================================================//
void flush_ring_buffer(ring_buffer_t*, double*);


//================================================================================================//
/**
* @brief This function adds to the ring buffer.
*
* @param[in,out] ring_buffer_t* self 
* @param[in] double* data
* @param[in] unsigned int data_length
*
* @return NONE
*/
//================================================================================================//
void enqueue_ring_buffer_dbl(ring_buffer_t*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function reads from the ring buffer.
*
* @param[in,out] ring_buffer_t* self 
* @param[in] unsigned int data_length
* @param[out] double* out_data
*
* @return NONE
*/
//================================================================================================//
void read_ring_buffer_dbl(ring_buffer_t*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function calculates the mean of the data in the ring buffer.
*
* @param[in] ring_buffer_t* self 
*
* @return double mean
*/
//================================================================================================//
double calculate_ring_buffer_mean_dbl(ring_buffer_t*);


//================================================================================================//
/**
* @brief This function calculates the variance of the data in the ring buffer.
*
* @param[in] ring_buffer_t* self 
*
* @return double variance
*/
//================================================================================================//
double calculate_ring_buffer_variance_dbl(ring_buffer_t*);


//================================================================================================//
/**
* @brief This function finds the median of the data in the ring buffer.
*
* @param[in] ring_buffer_t* self 
*
* @return double median
*/
//================================================================================================//
double find_ring_buffer_median_dbl(ring_buffer_t* self);


//================================================================================================//
/**
* @brief This function performs unit tests on the ring buffer.
*
* @return NONE
*/
//================================================================================================//
void test_ring_buffer();


//================================================================================================//
/**
* @brief This function initializes the fifo queue.
*
* @param[in,out] fifo_queue_t* self
* @param[in] unsigned int queue_length
*
* @return NONE
*/
//================================================================================================//
void initialize_fifo_queue(fifo_queue_t*,unsigned int);


//================================================================================================//
/**
* @brief This function resets the fifo queue.
*
* @param[in,out] fifo_queue_t* self
*
* @return NONE
*/
//================================================================================================//
void reset_fifo_queue(fifo_queue_t*);


//================================================================================================//
/**
* @brief This function copies the fifo queue.
*
* @param[in] fifo_queue_t* original
* @param[out] fifo_queue_t* copy
*
* @return NONE
*/
//================================================================================================//
void copy_fifo_queue(fifo_queue_t*,fifo_queue_t*);


//================================================================================================//
/**
* @brief This function enters data into the fifo queue.
*
* @param[in,out] fifo_queue_t* self
* @param[in] double* in_data
* @param[in] unsigned int num_frames
*
* @return NONE
*/
//================================================================================================//
void enqueue_fifo_dbl(fifo_queue_t*,double*,unsigned int);


//================================================================================================//
/**
* @brief This function takes data from the fifo queue.
*
* @param[in,out] fifo_queue_t* self
* @param[in] unsigned int num_frames
* @param[out] double* out_data
*
* @return NONE
*/
//================================================================================================//
void dequeue_fifo_dbl(fifo_queue_t*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function copies data from the fifo queue without dequeueing it.
*
* @param[in,out] fifo_queue_t* self
* @param[in] unsigned int num_frames
* @param[out] double* out_data
*
* @return NONE
*/
//================================================================================================//
void copy_fifo_data_dbl(fifo_queue_t*,unsigned int,double*);


//================================================================================================//
/**
* @brief This function tests if the fifo queue is full or not.
*
* @param[in] fifo_queue_t* self
*
* @return unsigned int is_full
*/
//================================================================================================//
unsigned int fifo_queue_is_full(fifo_queue_t*);


//================================================================================================//
/**
* @brief This function tests if the fifo queue can accept the given data length.
*
* @param[in] fifo_queue_t* self
* @param[in] unsigned int data_length
*
* @return unsigned int can_enqueue
*/
//================================================================================================//
unsigned int fifo_queue_has_room(fifo_queue_t*,unsigned int);


//================================================================================================//
/**
* @brief This function tests if the fifo queue is empty or not.
*
* @param[in] fifo_queue_t* self
*
* @return unsigned int is_empty
*/
//================================================================================================//
unsigned int fifo_queue_is_empty(fifo_queue_t*);


//================================================================================================//
/**
* @brief This function averages the frames in the fifo_queue.
*
* @param[in] fifo_queue_t* self
* @param[in] unsigned int num_frames
* @param[in] unsigned int frame_length
* @param[out] double* average_frame
*
* @return NONE
*/
//================================================================================================//
void average_fifo_frames(fifo_queue_t*,unsigned int,unsigned int,double*);


//================================================================================================//
/**
* @brief This function unit tests the fifo_queue.
*
* @return NONE
*/
//================================================================================================//
void test_fifo_queue();


#endif //RING_H//
