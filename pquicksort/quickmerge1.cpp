#include <iostream>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include "mpi.h"
#include <fstream>
#include <queue>

using namespace std;

int compareMyType (const void * a, const void * b);
struct Node
{
    int value;
    int index;
};		
bool operator<(const Node& a, const Node& b)
{
  return a.value > b.value;
}
main(int argc, char* argv[])
{
	int my_rank;
	int p;
	int source;
	int dest;
	int tag=0;
	int start, end;
	int n = 0;
	int local_n = 0;
	int i = 0;
	int local_i = 0;
	int index = 0; //for various index counts in for loops
	double time_start, time_stop;
	int *all_numbers;
	//local vars:
	int loc_n, *loc_arr, *loc_splitters, *loc_count;
	//global vars:
	int *glo_splitters, *glo_count, *all_loc_splitters;
	//int size = sizeof(arr)/sizeof(*arr);	//number of numbers
	
	
	MPI_Init(&argc, &argv);	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 	
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	//generate some numbers
	n = 2048;
	all_numbers = new int[n];
	srand(time(0));
	if(my_rank == 0){
		for (i = 0; i < n; i++){
			all_numbers[i] = rand()% 100; 
		}
	}

	//1. local sort -> all proc
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	loc_n = n/p;
	loc_arr = new int[loc_n];
	//scatter large array
	MPI_Scatter(all_numbers, loc_n, MPI_INT, loc_arr, loc_n, MPI_INT, 0, MPI_COMM_WORLD);
	qsort(loc_arr, loc_n, sizeof(int), compareMyType);
	//this part is working


	//2. find local splitters -> all proc
	loc_splitters = new int[p];
	index = 0;
	for(int i = 0; i < loc_n; i++){
		if(i%(loc_n/p)==0){//if i splits evenly by the number of elements in local array
			loc_splitters[index] = loc_arr[i];
			index++;  
		}
	}
	//splitters found -- this part working
	
	//3. find global splitters -> proc 0
	MPI_Gather(loc_splitters, p, MPI_INT, all_loc_splitters, p, MPI_INT, 0, MPI_COMM_WORLD);
   	if(my_rank==0){
		qsort(all_loc_splitters, p*p, sizeof(int), compareMyType);
	}
	//global splitters are gathered and sorted -- this part working
	
	//proc 0 picks splitters by:
	//	bcast global splitters
	MPI_Bcast(all_loc_splitters, p * p, MPI_INT, 0, MPI_COMM_WORLD);
	//	each proc generates a count that is <= each global splitter
	loc_count = new int[p*p];
	index = 0;
	int sum = 0;
	while(index < p*p){
		for(int i = 0; i < loc_n; i++){
			if (loc_arr[i] <= all_loc_splitters[index]){
				sum++;
			}
		}
		loc_count[index] = sum;
		sum = 0;
		index++;
	}
	//	proc 0 calls mpi reduce
	glo_count = new int[p*p];
	MPI_Reduce(loc_count, glo_count, p * p, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	//this part is working

	//	find count that is closest to rank * n/p, use that splitter
	if(my_rank==0){
		int processor = 1;//start with 1, because all counts are greater than 0 * n/p
		i = 0;
		for(i = 0; i < p*p; i++){
			if(glo_count[i] > processor * n/p){
				if(glo_count[i] - processor >= glo_count[i-1] - processor){
					glo_splitters[processor-1] = all_loc_splitters[i-1];
				}
				else{
					glo_splitters[processor-1] = all_loc_splitters[i];
				}
				processor++;				
			}
			glo_splitters[p-1] = all_loc_splitters[p*p-1];
		}
	}//not sure if this is working
	MPI_Bcast(glo_splitters, p, MPI_INT, 0, MPI_COMM_WORLD);

	//4. distribute data -> all to all, barrier
	//alltoallv will redistribute the local arrays if each proc has a sent_num and recv_num
	int * sendbuf = local_array;
    int * sent_num = new int[p];
    int * sdispls = new int[p];
    int * recvbuf; //should be large 
    int * recv_num = new int[p];
    int * recv_disp = new int[p];
        
    //find the sendcounts and sdispls arrays for each process
    int num;                    //total of elements in between splitters
    int displs;                 //displacement
    j = 0;                      //goes through local array
    for(i = 0; i < p; i++)      //goes through global splitters
    {
        num = 0;
        displs = j;        
        while(local_array[j] <= global_splitters[i] && j < local_n)
        {
            num++;
            j++;
        }        
        sent_num[i] = num;
        sdispls[i] = displs;
    }
        
    //give the recvcounts array to each process
    MPI_Alltoall(sent_num, 1, MPI_INT, recv_num, 1, MPI_INT, MPI_COMM_WORLD);

    //get the rdispls array for each process
    recv_disp[0] = 0;
    for(i = 1; i < p; i++)
        recv_disp[i] = recv_disp[i - 1] + recv_num[i - 1];

    //set up recvbuf array
    int arraysize = recv_disp[p - 1] + recv_num[p - 1];
    recvbuf = new int[arraysize];
   
	//sendbuf - starting address
	//sent_num - array specifying number of elements sent to each processor
	//sent_disp - array specifying displacement relative to sendbuf to process j
	//recv_num - array specifying max number that can be received
	//recv_disp - array specifying displacement relative to sendbuf FROM process i
	MPI_Alltoallv(sendbuf, sent_num, sent_disp, MPI_INT, recvbuf, recv_num, recv_disp, MPI_INT, MPI_COMM_WORLD);
	
	//5. local merge -> all proc
	// take first elements of every array into a new array
	//		call build heap
	//		take root (it will be smalles)
	//		move pointer in array that the root came from
	//		insert new array, call heapify
	//6. proc 0 merge 

	

	MPI_Finalize();
}

int compareMyType (const void * a, const void * b) //from cplusplus.com
{
  return ( *(int*)a - *(int*)b );
}