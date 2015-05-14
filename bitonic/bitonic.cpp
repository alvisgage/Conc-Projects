#include <stdio.h> 
#include <time.h>
#include <math.h> 
#include <stdlib.h> 
#include "mpi.h" 


//global local_size to be accessed without sending as arg
int local_size;
int compareMyType(const void * a, const void * b) {
    return ( * (int *)a - * (int *)b );
}

int getMax(const int * arr, int size){
	int max = 0;
	for (int i = 1; i < size; i++){
		if (arr[i] > max){
			max = arr[i];
		}
	}
	return max;
}
int getMin(const int * arr, int size){
	int min = local_size + 1;
	for (int i = 1; i < size; i++){
		if (arr[i] < min){
			min = arr[i];
		}
	}
	return min;
}
void HalfCleaner(bool up, int j, int partner, int*local_array){
	int i, buf, *send;
	//printf("me: %d, partner: %d\n", j, b);

	if(up){ //if sending to higher proc, send the lowest number in array
		send = &local_array[0];
	}else { //else send largest number in array
	   send = &local_array[local_size - 1];
	}

    int send_counter = 0;
    int * buffer_send = new int[local_size+1];

	//send either the smalles number or largest number in array
    MPI_Send(send, 1,  MPI_INT,  partner, 0, MPI_COMM_WORLD);

    // Receive new min of sorted numbers
    int recv_counter;
    int * buffer_recieve = new int[local_size+1];
	//get number from partner
    MPI_Recv(&buf, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //get all numbers in array that are larger or smaller than number recvd
    for (i = 0; i < local_size; i++) {
        if ((up && local_array[i] > buf) || (!up && local_array[i] < buf)) {
			send_counter++;
            buffer_send[send_counter] = local_array[i];
        } else { break; }
    }

    buffer_send[0] = send_counter; //how many we are sending

    // send numbers to partner
    MPI_Send(buffer_send,send_counter, MPI_INT, partner,0, MPI_COMM_WORLD);

    // get numbers from partner
    MPI_Recv(buffer_recieve, local_size, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	recv_counter = buffer_recieve[0]; //how many we are recving

	if (up){
		local_array[0] = getMax(buffer_recieve, recv_counter);
	} else{
		local_array[local_size - 1] = getMin(buffer_recieve, buffer_recieve[0]);
	}

    // resort the array
    qsort(local_array, local_size, sizeof(int), compareMyType);

    return;	
}

int main(int argc, char * argv[]) {

	double start;
	double stop;
	int my_rank;
	int p;
	int tag = 0;
	int * local_array;
	int num_elements;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    num_elements = atoi(argv[1]);
    local_size = num_elements / p;
    local_array = new int[local_size];

	//each proc generates its own local array so proc 0 doesn't waste time scattering one
    srand(time(0)); 
    for (int i = 0; i < local_size; i++) {
        local_array[i] = rand() % num_elements;
    }

    // Cube Dimension
    int dimensions = (int)(log2(p)); 
   
    if (my_rank == 0) {
        start = MPI_Wtime();
    }	   
    // local sort
    qsort(local_array, local_size, sizeof(int), compareMyType);

    // half cleaner calls
    for (int i = 0; i < dimensions; i++) {
        for (int j = i; j >= 0; j--) { 
			int partner = my_rank ^ ((j << 1) + (j+1));
			//int partner = my_rank ^ (1 << j);
			if (my_rank > partner) {
				HalfCleaner(false, j, partner, local_array);	//lower proc
			} else {
				HalfCleaner(true, j, partner, local_array); //higher proc
			}
        }
    } 

    MPI_Barrier(MPI_COMM_WORLD);

    if (my_rank == 0) {
		stop = MPI_Wtime();

        /*for (int i = 0; i < local_size; i++) {
            printf("%d ",local_array[i]);           
        }
		for (int j = 1; j < p; j++){
			int * tmp;
			tmp = new int[local_size];
			MPI_Recv(tmp, local_size, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int i = 0; i < local_size; i++){
				printf("%d ", tmp[i]);
			}
		}
        printf("\n\n");*/

        printf("time: %f\n", stop - start);
		printf("p: %d\n", p);
		printf("n: %d\n", num_elements);
    } else{
		//MPI_Send(local_array, local_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

   
    MPI_Finalize();
    return 0;
}