#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
void print(int * m, int size, const char* msg, int rank, bool whole);
void initializeMatrix(int * m, int size, bool empty);
void multiply(int * a, int * b, int * c, int size);
void swap(int * m, int * recvbuf, int * tmp);
int main(int argc, char *argv[])
{
	int rank, p;
	int n, local_n;
	int dimensions[2];
	int wrap_arounds[2];
	int left,right,up,down;
	int *a,*b,*c;
	int *recvbuf,*tmp;
	double start, stop;
	MPI_Comm comm;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	
	srand(time(0));
	
	dimensions[0] = dimensions[1] = 0;
	wrap_arounds[0] = wrap_arounds[1] = 1;
	MPI_Dims_create(p,2,dimensions);
	n = atoi(argv[1]);
	
	local_n=n/dimensions[0];	
	
	a = new int[local_n * local_n];
	b = new int[local_n * local_n];
	recvbuf = new int[local_n * local_n];
	c = new int[local_n * local_n];
	
	initializeMatrix(a, local_n, false);
	initializeMatrix(b, local_n, false);
	initializeMatrix(c, local_n, true);
	
	MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,wrap_arounds,1,&comm);
	MPI_Cart_shift(comm,0,1,&left,&right);
	MPI_Cart_shift(comm,1,1,&up,&down);
	start=MPI_Wtime();
	
	for(int i=0;i<dimensions[0];i++) {
		if(i==dimensions[0]-1)
			break;
		multiply(a, b, c, local_n);
		
		MPI_Sendrecv(a,local_n*local_n,MPI_INT,left,1,recvbuf,local_n*local_n,MPI_INT,right,1,comm,MPI_STATUS_IGNORE);
		swap(a, recvbuf, tmp);
		MPI_Sendrecv(b,local_n*local_n,MPI_INT,up,2,recvbuf,local_n*local_n,MPI_INT,down,2,comm,MPI_STATUS_IGNORE);
		swap(b, recvbuf, tmp);		
	}
	MPI_Barrier(comm);
	stop=MPI_Wtime();
	
	if(rank==0){
		printf("p: %d\tn: %d\ttime: %f\n",p, local_n, stop-start);
		//printf("%d\n",local_n);
		//print(a, local_n, "matrix a\n", rank, true);
		//print(b, local_n, "matrix b\n", rank, true);
		//print(c, local_n, "matrix c\n", rank, true);
	}
	
	free(a); 
	free(b); 
	free(recvbuf); 
	free(c);
	MPI_Finalize();
	return 0;
}
void swap(int * m, int * recvbuf, int * tmp){
	tmp = recvbuf;
	recvbuf = m;
	m = tmp;
}
void multiply(int * a, int * b, int * c, int size){
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			for(int k=0;k<size;k++)
				c[i*size+k]+=a[i*size+j]*b[j*size+k];
}
void initializeMatrix(int * m, int size, bool empty){
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++) {
			int r = 0;
			if (!empty)
				r = rand() % RAND_MAX + 1;
			m[i*size+j]=r;
		}
}
void print(int * m, int size, const char* msg, int rank, bool whole){
	printf(msg);
	printf(" rank: %d\n", rank);
	if(whole){/*print entire matrix*/
		for (int i = 0; i < size * size; i++){
			printf("%d ", m[i]);
			if ( (i+1) % size == 0){
				printf("\n");
			}
		}
	}
	else{
		for(int i=0;i<size;i++)
			for(int j=0;j<size;j++) {
				printf("%d ", m[i*size+j]);
			}
	}
	printf("\n");
}