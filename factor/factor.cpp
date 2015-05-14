#include <iostream>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include "mpi.h"

using namespace std;
bool isSquare(long int sum);
main(int argc, char* argv[])
{
	int my_rank; 						/*processor rank*/
	int p;								/* processors*/
	int source;							/* sending processor*/
	int dest;							/*receiving processor */
	int tag=0;
	long int i;							/* starting index for each proc*/
	long int limit;						/* ending index*/
	long long arr[2];					/* (j+i)(j-i)*/
	long long n = atol(argv[1]);			/* product of two primes*/
	double start,stop;					/* timers */
	double tmpP;


	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	//MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* broadcast product of primes*/
	
	start = MPI_Wtime(); //start time
	
	limit = pow(n,2)/2;	  
	i = (my_rank * n)/(3*p);				/* start index*/
	limit = (((my_rank + 1) * n)/(3*p)) -1;	/* limit index */

	bool square = false;

	while(!(square) && i < limit)
	{		    												 
		i++;
		if (isSquare(n + (i * i))){  /* if n + i^2 is a square*/
			square = true; 
		}	 
	}	 

	if (square){/* found square */
		//set j+i, j-1
		arr[0] = (sqrt(n+(i*i)) + i);
		arr[1] = (sqrt(n+(i*i)) - i);
		MPI_Send(&arr, 2, MPI_LONG_LONG, 0, tag, MPI_COMM_WORLD);
	}  	

	//MPI_Send(&arr, 2, MPI_LONG, 0, tag, MPI_COMM_WORLD);
	//MPI_Bcast(&arr, 1, MPI_LONG, 0, MPI_COMM_WORLD);

	if (my_rank == 0)
	{	  
	   //receive 
	    MPI_Recv(&arr, 2, MPI_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		cout << arr[0] << ", " << arr[1] << endl;
		stop = MPI_Wtime();
		tmpP = stop - start;
		cout << "time for " << p << " processors: " << tmpP << endl;;   
	}

	//run again with 1 processor, tmp1
	//Sp = tmp1/tmpP
	//Ep = Sp/p
		
	MPI_Finalize();
}

bool isSquare(long int sum)
{
	bool result = false;
	int mod = sum % 16;
	if (mod == 0 || mod == 1 || mod == 4 || mod == 16)
	{
		if(pow(int(sqrt(sum) + .1),2) ==  sum){
			result = true;
		}
	}
	return result;
}