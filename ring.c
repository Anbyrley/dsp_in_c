#include "ring.h"

void initialize_ring_buffer( ring_buffer_t* self,
							 unsigned int queue_length)
{

	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Ring Buffer Is NULL! In Function -- initialize_ring_buffer!\n");
		quit();
	}
	if (queue_length > MAX_QUEUE_LENGTH){
		fprintf(stderr, "Error:: Qeueue Length Is Too Large! In Function -- initialize_ring_buffer!\n");
		quit();
	}

	//===Set Data Fields===//
	self->queue_length = queue_length;

	//===Set Head and Tail===//
	self->tail = 0;
	self->full = 0;

	//===Zero Array===//
	initialize_array_cmplx(self->data.cmplx, MAX_QUEUE_LENGTH);

	return;
}

void copy_ring_buffer( ring_buffer_t* original,
					   ring_buffer_t* copy )
{
	copy->queue_length = original->queue_length;
	copy->tail = original->tail;
	copy->full = original->full;
	copy_array_cmplx(original->data.cmplx, MAX_QUEUE_LENGTH, copy->data.cmplx);
	return;
}

void flush_ring_buffer( ring_buffer_t* self,
						double* out_data )
{
	if (out_data != NULL){
		read_ring_buffer_dbl(self, self->queue_length, out_data);
	}
	self->tail = 0;	
	self->full = 0;
	initialize_array_cmplx(self->data.cmplx, MAX_QUEUE_LENGTH);
	return;
}

void enqueue_ring_buffer_dbl( ring_buffer_t* self,
						  	  double* in_data,
						  	  unsigned int data_length )
{
	unsigned int i, k;

	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Ring Buffer Is NULL! In Function -- enqueue_ring_buffer_dbl!\n");
		quit();
	}
	if (in_data == NULL){
		fprintf(stderr, "Error:: In Data Is NULL! In Function -- enqueue_ring_buffer_dbl!\n");
		quit();
	}
	if (!isnormal((double)data_length)){
		fprintf(stderr, "Error:: Data Length In Invalid! In Function -- enqueue_ring_buffer_dbl!\n");
		quit();
	}

	//===Run Enqueue===//
	k = 0;
	for (i=0; i<data_length; i++){			

		//===Add Data===//
		*(self->data.dbl + self->tail) = *(in_data + k++);			

		//===Increment Tail===//
		self->tail += 1;

		//===Check For Wrap Around===//	
		if (self->tail == self->queue_length){
			self->tail = 0;
			self->full = 1;
		}
	}

	return;
}

void read_ring_buffer_dbl( ring_buffer_t* self,
						   unsigned int data_length,
						   double* out_data )
{
	int j;
	unsigned int i, k;

	k = 0;
	j = self->tail; 
	if (j==0){
		j = self->queue_length;
	}
	for (i=0; i<data_length; i++){			

		//===Get Data===//
		*(out_data + k++) = *(self->data.dbl + j - 1); 			
		j--;

		//===Check For Wrap Around===//	
		if (j <= 0){
			j = self->queue_length;
		}
	}

	return;
}

void print_ring_buffer_data_dbl( ring_buffer_t* self,
								 FILE* file )
{

	int k;
	unsigned int i;
	k = self->tail;
	for (i=0; i<self->queue_length; i++){
		fprintf(file, "%lf\n", *(self->data.dbl + k - 1));
		k--;
		if (k <= 0) k = self->queue_length;		
	}

}

double calculate_ring_buffer_mean_dbl(ring_buffer_t* self)
{
	return compute_mean_dbl(self->data.dbl, self->queue_length);
}

double calculate_ring_buffer_variance_dbl(ring_buffer_t* self)
{
	return compute_variance_dbl(self->data.dbl, self->queue_length);
}

double find_ring_buffer_median_dbl(ring_buffer_t* self)
{

	if (self->full){
		return find_median_dbl(self->data.dbl, self->queue_length);
	}
	else{
		return find_median_dbl(self->data.dbl, self->tail);
	}
}

void test_ring_buffer()
{
	unsigned int i, out_length, num_frames, frame_length;
	double mean, variance, in_values[2*MAX_QUEUE_LENGTH], out_values[2*MAX_QUEUE_LENGTH];
	ring_buffer_t* ring_buffer;

	//===Mallocs===//
	ring_buffer = malloc(sizeof(ring_buffer_t));

	//===Create Values To Enter===//
	for (i=0; i<2*MAX_QUEUE_LENGTH; i++){
		in_values[i] = i;
		out_values[i] = 0;
	}

	//===Initialize Buffer===//
	frame_length = 3; num_frames = 1;
	initialize_ring_buffer(ring_buffer, frame_length*num_frames);

	//===Add To Buffer===//
	enqueue_ring_buffer_dbl(ring_buffer, in_values, 5);

	//===Print Buffer===//
	fprintf(stdout, "Buffer Values\n");
	print_ring_buffer_data_dbl(ring_buffer, stdout);
	newline();

	//===Calculate Mean===//
	mean = calculate_ring_buffer_mean_dbl(ring_buffer);
	fprintf(stdout, "Mean: %lf\n", mean);
	newline();

	//===Calculate Variance===//
	variance = calculate_ring_buffer_variance_dbl(ring_buffer);
	fprintf(stdout, "Variance: %lf\n", variance);
	newline();

	//===Get Data====//
	out_length = 8;
	read_ring_buffer_dbl(ring_buffer, out_length, out_values);

	//===Print Results===//
	newline();
	fprintf(stdout, "Output:\n");
	print_vector_dbl(out_values, out_length, stdout);
	newline();

	//===Print Buffer===//
	fprintf(stdout, "Buffer Values\n");
	print_ring_buffer_data_dbl(ring_buffer, stdout);
	newline();
	
	//===Free===//
	free(ring_buffer);

	return;
}

//================================================================================================//
//=========================================Queue Methods==========================================//
//================================================================================================//

void initialize_fifo_queue( fifo_queue_t* self,
							unsigned int queue_length )
{

	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Fifo Queue Is NULL! In Function -- initialize_fifo_queue!\n");
		quit();
	}
	if (!isnormal((double)queue_length)){
		fprintf(stderr, "Error:: Queue Length Is Invalid! In Function -- initialize_fifo_queue!\n");
		quit();
	}
	if (queue_length > MAX_QUEUE_LENGTH){
		fprintf(stderr, "Error:: Queue Length Is Too Large! In Function -- initialize_fifo_queue!\n");
		quit();
	}

	//===Set Locals===//
	self->head = 0;
	self->tail = 0;
	self->head_wrap = 0;
 	self->tail_wrap = 0;
	self->queue_length = queue_length;	

	//===Initialize===//
	initialize_array_cmplx(self->data.cmplx, MAX_QUEUE_LENGTH);

	return;
}

void reset_fifo_queue( fifo_queue_t* self )
{

	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Fifo Queue Is NULL! In Function -- reset_fifo_queue!\n");
		quit();
	}

	//===Reset Head And Tail===//
	self->head = self->tail = 0;
	self->head_wrap = self->tail_wrap = 0;

	//===Reset Data===//	
	initialize_array_cmplx(self->data.cmplx, MAX_QUEUE_LENGTH);

	return;
}

void copy_fifo_queue( fifo_queue_t* original,
					  fifo_queue_t* copy )
{
	//===Error Pre-Check===//
	if (original == NULL){
		fprintf(stderr, "Error:: Original Is NULL! In Function -- copy_fifo_queue!\n");
		quit();
	}
	if (copy == NULL){
		fprintf(stderr, "Error:: Copy Is NULL! In Function -- copy_fifo_queue!\n");
		quit();
	}

	//===Copy Locals===//
	copy->head = original->head;
	copy->tail = original->tail;
	copy->head_wrap = original->head_wrap;
	copy->tail_wrap = original->tail_wrap;
	copy->queue_length = original->queue_length;

	//===Copy Data===//
	copy_array_cmplx(original->data.cmplx, original->queue_length, copy->data.cmplx);

	return;
}

void enqueue_fifo_dbl( fifo_queue_t* self, 
				   	   double* in_data, 
				   	   unsigned int data_length )
{
	unsigned int i;
	double temp[MAX_POLYNOMIAL_ORDER+1];

	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Fifo Queue Is NULL! In Function -- enqueue_fifo_dbl!\n");
		quit();
	}
	if (in_data == NULL){
		fprintf(stderr, "Error:: In Data Is NULL! In Function -- enqueue_fifo_dbl!\n");
		quit();
	}
	if (!isnormal((double)data_length)){
		fprintf(stderr, "Error:: Data Length Is Invalid! In Function -- enqueue_fifo_dbl!\n");
		quit();
	}

	//===Run Enqueue===//
	for (i=0; i<data_length; i++){

		//===Check if Full===//
		if (self->tail == self->head &&
			self->tail_wrap != self->head_wrap) {
			fprintf(stderr, "Error:: Fifo Queue Is Full! In Function -- enqueue_fifo_dbl!\n");
			dequeue_fifo_dbl(self, 1, temp);
		}
				
		//===Add Data===//
		*(self->data.dbl + self->tail) = *(in_data + i);			

		//===Increment Tail===//
		self->tail += 1;
		
		//===Check For Wrap Around===//	
		if (self->tail == self->queue_length){
			self->tail = 0;
			(self->tail_wrap)++;
		}

	}

	return;
}

void dequeue_fifo_dbl( fifo_queue_t* self,
				   	   unsigned int data_length, 
				   	   double* out_data)
{
	unsigned int i;

	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Fifo Queue Is NULL! In Function -- dequeue_fifo_dbl!\n");
		quit();
	}
	if (out_data == NULL){
		fprintf(stderr, "Error:: Out Data Is NULL! In Function -- dequeue_fifo_dbl!\n");
		quit();
	}
	if (!isnormal((double)data_length)){
		fprintf(stderr, "Error:: Data Length Is Invalid! In Function -- dequeue_fifo_dbl!\n");
		quit();
	}


	//===Run Dequeue===//
	for (i=0; i<data_length; i++){
			
		//===Check For Empty===//	
		if (self->head == self->tail &&
			self->head_wrap == self->tail_wrap){
			fprintf(stderr, "Error:: Fifo Queue Is Empty! In Function -- dequeue_fifo_dbl!\n");
			return;
		}	
		
		//===Save Data===//				
		*(out_data + i) = *(self->data.dbl + self->head);
	
		//===Incremement Head===//
		self->head += 1;

		//===Check For Wrap Around===//	
		if (self->head == self->queue_length){
			self->head = 0;
			(self->head_wrap)++;
		}	

	}

	return;
}

void copy_fifo_data_dbl( fifo_queue_t* self,
						 unsigned int data_length,
						 double* out_data )
{
	unsigned int i, saved_head, saved_head_wrap;

	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Fifo Queue Is NULL! In Function -- copy_fifo_data_dbl!\n");
		quit();
	}
	if (out_data == NULL){
		fprintf(stderr, "Error:: Out Data Is NULL! In Function -- copy_fifo_data_dbl!\n");
		quit();
	}
	if (!isnormal((double)data_length)){
		fprintf(stderr, "Error:: Data Length Is Invalid! In Function -- copy_fifo_data_dbl!\n");
		quit();
	}

	//===Run Dequeue===//
	saved_head = self->head; saved_head_wrap = self->head_wrap;
	for (i=0; i<data_length; i++){
			
		//===Check For Empty===//	
		if (self->head == self->tail &&
			self->head_wrap == self->tail_wrap){
			fprintf(stderr, "Error:: Fifo Queue Is Empty! In Function -- copy_fifo_data_dbl!\n");
			return;
		}	
		
		//===Save Data===//				
		*(out_data + i) = *(self->data.dbl + self->head);
	
		//===Incremement Head===//
		self->head += 1;

		//===Check For Wrap Around===//	
		if (self->head == self->queue_length){
			self->head = 0;
			(self->head_wrap)++;
		}	
	}

	//===Replace===//
	self->head = saved_head;
	self->head_wrap = saved_head_wrap;

	return;
}

unsigned int fifo_queue_is_full( fifo_queue_t* self )
{
	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Fifo Queue Is NULL! In Function -- fifo_queue_is_full!\n");
		quit();
	}

	//===Check if Full===//
	if ((self->tail == self->head) &&
		(self->tail_wrap != self->head_wrap)) {
		return 1;
	}
	else{
		return 0;
	}
	return 0;
}

unsigned int fifo_queue_has_room( fifo_queue_t* self,
								  unsigned int data_length )
{
	unsigned int i;
	unsigned int faux_tail, faux_head;
	unsigned int faux_tail_wrap;
	unsigned int faux_head_wrap;

	//===Save===//
	faux_tail = self->tail;
	faux_head = self->head;
	faux_tail_wrap = self->tail_wrap;
	faux_head_wrap = self->head_wrap;

	//===Run Through Length===//
	for (i=0; i<data_length; i++){
		//===Check if Full===//
		if (faux_tail == faux_head &&
			faux_tail_wrap != faux_head_wrap) {
			return 0;
		}

		//===Increment Tail===//
		faux_tail += 1;
		
		//===Check For Wrap Around===//	
		if (faux_tail == self->queue_length){
			faux_tail = 0;
			faux_tail_wrap = (faux_tail_wrap + 1)%2;
		}

	}
	return 1;
}

unsigned int fifo_queue_is_empty( fifo_queue_t* self )
{
	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Fifo Queue Is NULL! In Function -- fifo_queue_is_empty!\n");
		quit();
	}
	
	//===Check For Empty===//	
	if (self->head == self->tail &&
		self->head_wrap == self->tail_wrap){
		return 1;
	}	
	else{
		return 0;
	}
	return 0;
}

void average_fifo_frames_dbl( fifo_queue_t* self,  
						  	  unsigned int num_frames,
						  	  unsigned int frame_length,
						  	  double* average_frame )
{

	unsigned int i,f;

	//===Error Pre-Check===//
	if (self == NULL){
		fprintf(stderr, "Error:: Fifo Queue Is NULL! In Function -- average_fifo_frames_dbl!\n");
		quit();
	}
	if (average_frame == NULL){
		fprintf(stderr, "Error:: Average Frame Is NULL! In Function -- average_fifo_frames_dbl!\n");
		quit();
	}
	if (!isnormal((double)num_frames)){
		fprintf(stderr, "Error:: Num Frames Is Invalid! In Function -- average_fifo_frames_dbl!\n");
		quit();
	}
	if (!isnormal((double)frame_length)){
		fprintf(stderr, "Error:: Frame Length Is Invalid! In Function -- average_fifo_frames_dbl!\n");
		quit();
	}

	//===Initialize Average Frame===//
	initialize_array_dbl(average_frame, frame_length);

	//===Sum Data===//
	for (f=0; f<num_frames*frame_length; f++){
		average_frame[f%frame_length] += *(self->data.dbl + f);
	}
	
	//===Normalize By Frames Added To This Frame===//
	for(i=0; i<frame_length; i++){
		average_frame[i] /= (double)num_frames;
	}

	return;
}

void test_fifo_queue()
{
	double frame[MAX_POLYNOMIAL_ORDER+1];
	fifo_queue_t* queue;

	//===Mallocs===//
	queue = malloc(sizeof(fifo_queue_t));

	//===Initialize===//
	initialize_fifo_queue(queue, 2*3);

	//===Put Two In, Get 1st Out===//
	frame[0] = 1; frame[1] = 2; frame[2] = 3;
	newline();
	fprintf(stdout, "Entered:\n");
	print_vector_dbl(frame, 3, stdout);
	enqueue_fifo_dbl(queue, frame, 3);
	frame[0] = 4; frame[1] = 5; frame[2] = 6;
	fprintf(stdout, "Entered:\n");
	print_vector_dbl(frame, 3, stdout);
	newline();
	enqueue_fifo_dbl(queue, frame, 3);
	dequeue_fifo_dbl(queue, 3, frame);
	fprintf(stdout, "Retreived:\n");
	print_vector_dbl(frame, 3, stdout);
	newline();

	//===Reset===//
	reset_fifo_queue(queue);

	//===Put Three In, Get 2nd Out===//
	frame[0] = 1; frame[1] = 2; frame[2] = 3;
	newline();
	fprintf(stdout, "Entered:\n");
	print_vector_dbl(frame, 3, stdout);
	enqueue_fifo_dbl(queue, frame, 3);
	frame[0] = 4; frame[1] = 5; frame[2] = 6;
	fprintf(stdout, "Entered:\n");
	print_vector_dbl(frame, 3, stdout);
	enqueue_fifo_dbl(queue, frame, 3);
	frame[0] = 7; frame[1] = 8; frame[2] = 9;
	fprintf(stdout, "Entered:\n");
	print_vector_dbl(frame, 3, stdout);
	newline();
	enqueue_fifo_dbl(queue, frame, 3);
	dequeue_fifo_dbl(queue, 3, frame);
	fprintf(stdout, "Retreived:\n");
	print_vector_dbl(frame, 3, stdout);
	newline();

	
	return;
}
