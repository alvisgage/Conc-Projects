#include <fstream>
#include <iostream>
#include <cstdlib>      //for random function and qsort
#include <ctime>
#include <queue>        //for min heap
#include "mpi.h"

using namespace std;

//Node to store value and index for the priority queue
struct Node
{
    int value;
    int index;
};

//override < opperator for min heap priority queue
bool operator<(const Node& a, const Node& b)
{
  return a.value > b.value;
}

int compareMyType(const void* a, const void* b)
{
    return(*(int*)a - *(int*)b);
}

int main(int argc, char** argv)
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
	double time_start, time_stop;
	int *splitter_arr, *total_splitters, *splitter_count, *total_splitter_count;
	int *arr, *local_arr;
	//int size = sizeof(arr)/sizeof(*arr);	//number of numbers
	//set spliter arrays
	splitter_arr = new int[p];
	total_splitters = new int[p*p];//p^2 splitters total
	splitter_count = new int[p*p]; //summing subarray numbers <= each splitter
	total_splitter_count = new int[p*p]; //proc 0 does mpi_reduce on the splitter numbers
	//time_start = MPI_Wtime();
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    n = 2048;
	arr = new int[n];
	srand(time(0));
	
	if(my_rank == 0){
		for (i = 0; i < n; i++){
			arr[i] = rand()% 100; 
		}
	}

    //send number of elements
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    local_n = n/p;
	local_arr = new int[local_n];
	//scatter array... man i hope im doing this right
	MPI_Scatter(arr, local_n, MPI_INT, local_arr, local_n, MPI_INT, 0, MPI_COMM_WORLD);
   	
	qsort(local_arr, local_n, sizeof(int), compareMyType);

    //each proc has p splitters, total has p^2 splitters
	//go through local array, if i % local_n / p == 0, add that number to splitter array
	int index = 0;
	for (i = 0; i < local_n; i++){
		if (i%(local_n/p)==0){
			splitter_arr[index] = local_arr[i];
			index++;
		}
	}

    //send splitters to 0
	//mpi gather should get all local splitter arrays
	MPI_Gather(splitter_arr, p, MPI_INT, total_splitters, p, MPI_INT, 0, MPI_COMM_WORLD);
   	if (my_rank == 0)
	{				
		//splitters should be in total_splitters
		qsort(total_splitters, p*p, sizeof(int), compareMyType);
	}



   //send allllll the splitters back to everyone
	MPI_Bcast(total_splitters, p*p, MPI_INT, 0, MPI_COMM_WORLD);
	
    //for each splitter, count how many numbers in subarray are <= that splitter
	for (i = 0; i < p*p; i++){
		splitter_count[i] = splitter_count[i-1];//set next count to previous count
		for(int j = 0; j < local_n && local_arr[j] <= total_splitters[i]; j++){
			splitter_count[i] = splitter_count[i] + 1;//increment current count	
		}
	}


    //proc 0???
    MPI_Reduce(splitter_count, total_splitter_count, p*p, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
    //process 0 finds the global splitters
    if(my_rank == 0)
    {
        int splitSpot = 1;      //find the counts that are closest to splitSpot * n/p
        for(i = 0; i < p * p; i++)      //goes through allCounts array
        {
            if(allCounts[i] > splitSpot * n/p)
            {
                if(allCounts[i] - splitSpot > allCounts[i-1] - splitSpot)
                    global_splitters[splitSpot-1] = all_splitters[i-1];
                else
                    global_splitters[splitSpot-1] = all_splitters[i];
                    
                splitSpot++;
            }
        }
        global_splitters[p - 1] = all_splitters[p * p - 1];
    }
 
    //broadcast the global splitters
    MPI_Bcast(global_splitters, p, MPI_INT, 0, MPI_COMM_WORLD);

    //4. DISTRIBUTE DATA
    //variables for MPI_Alltoallv
    int * sendbuf = local_array;
    int * sendcounts = new int[p];
    int * sdispls = new int[p];
    int * recvbuf;
    int * recvcounts = new int[p];
    int * rdispls = new int[p];
        
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
        
        sendcounts[i] = num;
        sdispls[i] = displs;
    }
        
    //give the recvcounts array to each process
    MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

    //get the rdispls array for each process
    rdispls[0] = 0;
    for(i = 1; i < p; i++)
        rdispls[i] = rdispls[i - 1] + recvcounts[i - 1];

    //set up recvbuf array
    int arraysize = rdispls[p - 1] + recvcounts[p - 1];
    recvbuf = new int[arraysize];

    //all to all call to redistribute the arrays
    MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_INT, recvbuf, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);

    //5. LOCAL MERGE
    local_array = new int[arraysize];   //clear local_array

    //make a copy of rdispls for the priority queue
    int qdispls[p];
    int q2displs[p];
    j = 0;
    for(i = 0; i < p; i++)
    {
        if(rdispls[i] < arraysize)
        {
            if(i == 0 || i > 0 && rdispls[i] != rdispls[i-1])
            {
                qdispls[j] = q2displs[j] = rdispls[i];
                j++;
            }
        }
    }
    
    //create a min heap with nodes that store value and index
    priority_queue<Node> min_heap;
    Node z[p];
    for(i = 0; i < j; i++)
    {
        z[i].index = i;
        z[i].value = recvbuf[q2displs[i]];
        min_heap.push(z[i]);
    }

    //merge arrays by using the queue
    Node tmp;
    for(i = 0; i < arraysize; i++)
    {
        tmp = min_heap.top();
        local_array[i] = tmp.value;
        q2displs[tmp.index]++;

        //pop out top value
        if(!min_heap.empty())
            min_heap.pop();

        //push a new value in
        if(tmp.index == j - 1)      //prevent going past array size
        {
            if(q2displs[tmp.index] < arraysize)
            {
                tmp.value = recvbuf[q2displs[tmp.index]];
                min_heap.push(tmp);
            }
        }
        else    //prevent going to next block
        {
            if(q2displs[tmp.index] < qdispls[tmp.index+1])
            {
                tmp.value = recvbuf[q2displs[tmp.index]];
                min_heap.push(tmp);
            }
        }
    }

    // process 0 collects the sorted data in array, and writes the data to a file
    int * acounts = new int[p];
    int * adispls = new int[p];
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //find acounts and adispls
    MPI_Gather(&arraysize, 1, MPI_INT, acounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(my_rank == 0)
    {    
        adispls[0] = 0;
        for(i = 1; i < p; i++)
            adispls[i] = adispls[i - 1] + acounts[i - 1];
    }

    //gather all local_arrays into array
    MPI_Gatherv(local_array, arraysize, MPI_INT, array, acounts, adispls, MPI_INT, 0, MPI_COMM_WORLD);

    stop = MPI_Wtime();
    
    //check if array is sorted
    if(my_rank == 0)
    {
        bool isSorted = true;
        for(i = 1; i < n; i++)
        {
            if(array[i] < array[i-1])
                isSorted = false;
        }
        
        if(isSorted == true)
        {
            cout << "The array is sorted.\n";
            cout << "First number: " << array[0] << endl;
            cout << "Last numner: " << array[n-1] << endl;
        }
        else
            cout << "The array is not sorted.\n";
    }
  
    MPI_Finalize();

 }  // main
