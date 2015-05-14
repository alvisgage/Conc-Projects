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

int compare(const void* a, const void* b)
{
    return(*(int*)a - *(int*)b);
}

int main(int argc, char** argv)
{
    int* array;                 //overall array containing n numbers
    int* local_array;           //local processor array containing n/p numbers
    int* local_splitters;       //local processor array containing local splitters
    int* all_splitters;         //process 0 array containing all the local splitters
    int* global_splitters;      //array containing the global splitters
    int* dist_array;            //local array containing the distributed numbers
    int n;              // number of  elements to sort
    int local_n;        // = n/p
    int p;              // number of processes
    int my_rank;
    int i, j;
    int tag = 0;
    double start, stop; // for timing
    MPI_Status status;

    // start MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0)
    {        
        //create an array of size n taken from runtime arguments
        n = atoi(argv[1]);
        array = new int[n];
        if(array == NULL)
        {
             cout << "array allocation failed!!\n";
             MPI_Abort(MPI_COMM_WORLD, 4);
             exit(0);
        }
        
        //fill array with random numbers
        srand(time(0));
        for(i = 0; i < n; i++)
        {
            array[i] = rand()%n;
        }
        
        start = MPI_Wtime();
    }

    //1. LOCAL SORT
    // broadcast the value of n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //create local arrays
    local_n = n/p;
    local_array = new int[local_n];
    if(local_array == NULL)
    {
         cout << "local_array allocation failed!!\n";
         MPI_Abort(MPI_COMM_WORLD, 4);
         exit(0);
    }
    
    // scatter the array over the processes
    MPI_Scatter(array, local_n, MPI_INT, local_array, local_n, MPI_INT, 0, MPI_COMM_WORLD);
    
    //each process use quicksort to sort the local array
    qsort(local_array, local_n, sizeof(int), compare);

    //2. FIND LOCAL SPLITTERS
    //each process finds local splitters (p splitters per process)
    local_splitters = new int[p];
    j = 0;
    for(i = local_n/p; i < local_n; i = i + (local_n/p))
    {
        local_splitters[j] = local_array[i-1];
        j++;
    }
    local_splitters[j] = local_array[i-1];

    //3. FIND GLOBAL SPLITTERS
    //process 0 gathers all the local splitters
    global_splitters = new int[p];
    all_splitters = new int[p * p];
    MPI_Gather(local_splitters, p, MPI_INT, all_splitters, p, MPI_INT, 0, MPI_COMM_WORLD);
    
    //process 0 sorts all the splitters
    if(my_rank == 0)
        qsort(all_splitters, p * p, sizeof(int), compare);
/*   
    //find global splitters without communication
    if(my_rank == 0)
    {        
        j = 0;
        for(i = p; i < p * p; i += p)
        {
            global_splitters[j] = all_splitters[i-1];
            j++;
        }
        global_splitters[j] = all_splitters[p * p - 1];
    }
*/
    //process 0 broadcasts all the splitters to every process
    MPI_Bcast(all_splitters, p * p, MPI_INT, 0, MPI_COMM_WORLD);
    
    //every process generates a count that is less than/equal to the splitter
    int * counts = new int[p * p];      //stores local counts
    j = 0;
    for(i = 0; i < p * p; i++)  //goes through all_splitters array
    {
        while(j < local_n && local_array[j] <= all_splitters[i])
            j++;
        counts[i] = j;
    }
    
    //Sum up the counts to process 0
    int * allCounts = new int[p * p];   //stores global counts
    MPI_Reduce(counts, allCounts, p * p, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
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
    if (recvbuf == NULL)
    {
         cout << "recvbuf allocation failed!!\n";
         MPI_Abort(MPI_COMM_WORLD, 4);
         exit(0);
    }

    //all to all call to redistribute the arrays
    MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_INT, recvbuf, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);

    //5. LOCAL MERGE
    local_array = new int[arraysize];   //clear local_array
    if (local_array == NULL)
    {
         cout << "local_array allocation failed!!\n";
         MPI_Abort(MPI_COMM_WORLD, 4);
         exit(0);
    }

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
/*
    //output sorted array to file
    if(my_rank == 0)
    {
        ofstream outfile;
        outfile.open(argv[2]);
        if(!outfile)
            cout << "outfile failed to open \n";
        for(i = 0; i < n; i++)
            outfile << array[i] << " ";
        outfile.close();
    }
*/    
    MPI_Finalize();

 }  // main
