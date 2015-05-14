#include <stdlib.h>
#include "mpi.h"
using namespace std; 

void print( double *M, int x, int y ) 
{
  int i, j;

  for( i=0; i<y; i++ ) {
    for( j=0; j<x; j++ )  
      cout << M[ i*x+j ] << "  ";
    cout << endl;
  }
}
void print_solution( double *A, double *y, int n )
{
  int i, j;
  double *x = new double[n];

  for( i=n-1; i>=0; i-- ) {
    x[i] = y[i];
    for( j=n-1; j>i; j-- )
      x[i] -= x[j] * A[i*n+j];
  }
  for( i=0; i<n; i+=4 ) {
    for( j=i; j<i+4 && j<n; j++ )
      cout << "x[" << j << "] = " << x[j] << "  ";
    cout << endl;
  }
}

main( int argc, char *argv[] )
{
  double *A, *b, *y;
  double start, stop;
  int n, r, p, my_rank, local_n;
double *local_A, *local_b;
  n = atoi(argv[1]);
	cout << "n: " << n << endl;


	if(my_rank==0){
		A = new double[n*n];	
		b = new double[n];
		y = new double[n];
}
	/*values were resulting in nan, inf, -inf*/
	/*this should set solution to be integer*/

	for (int i=0; i<n; i++){
		b[i] = 0.0;
		for (int j=0; j<n; j++){
			r = rand() % n;
			A[i*n+j] = r;
			b[i] += j*r;
		}
	} 
  MPI_Init (&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&p);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); 
	local_n = n/p;
  
	if(my_rank==0){
		cout << "mat A" << endl;
		print(A, n, n);
		cout << "b: " << endl;
		print(b, n, 1);
		
	}	
	local_A = new double[n*local_n];
	local_b = new double[local_n];

  start = MPI_Wtime();
  /*distribute A as local matrix*/
	if( my_rank == 0 ) {
		for( int h=p-1; h>=0; h-- ) {
			for( int i=h; i<n; i+=p )
				for( int j=0; j<n; j++ ) {
					local_A[ (i/p)*n+j ] = A[ i*n+j ];
				}
			if( h != 0 ) 
				MPI_Send( local_A, local_n*n, MPI_DOUBLE, h, 10, MPI_COMM_WORLD );
		}
	}
	else{
		MPI_Recv( local_A, local_n*n, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	/*end A distribution*/	
  /*distribute b as local b*/
	if( my_rank == 0 ) {
		for( int h=p-1; h>=0; h-- ) {
			for( int i=h; i<n; i+=p )
				for( int j=0; j<1; j++ ) {
					local_b[ (i/p)+j ] = b[ i+j ];
				}
			if( h != 0 ) 
				MPI_Send( local_b, (local_n), MPI_DOUBLE, h, 10, MPI_COMM_WORLD );
		}
	}
	else{
		MPI_Recv( local_b, (local_n), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	}
	/*end b distribution*/  
	
	/*gaussian elimination*/
	double *tmp = new double[local_n];
	double *row = new double[n+1];
	int receiver, sender, pos;
	receiver = (p+(my_rank-1)) % p;
	sender = (my_rank+1) % p;
	for(int k=0; k<n; k++ ) {
		if( (k%p) == my_rank ) {
			for(int j=k+1; j<n; j++ ){					
				local_A[ ((k/p)*n)+j ] /= local_A[ ((k/p)*n)+k ];
			}
			tmp[ k/p ] = local_b[ k/p ] / local_A[ ((k/p)*n)+k ];
			local_A[ ((k/p)*n)+k ] = 1.0;

			for(int j=0; j<n; j++ ){
				row[j] = local_A[ ((k/p)*n)+j ];
			}
			row[n] = tmp[ k/p ];
			MPI_Send( row, n+1, MPI_DOUBLE, sender, 20, MPI_COMM_WORLD );	
		}
		else {			
			MPI_Recv( row, n+1, MPI_DOUBLE, receiver, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			if( sender != k%p ) {	
				MPI_Send( row, n+1, MPI_DOUBLE, sender, 20, MPI_COMM_WORLD );
			}
		}
		/*if (my_rank <= k%p){
			pos = (int)k/p+1;
		}
		else{
			pos = (int)k/p;
		}*/
		pos = (my_rank <= k&p) ? (int) k/p+1 : (int) k/p;
		for(int i=pos; i<local_n; i++ ) {			
			for(int j=k+1; j<n; j++ ) {
				local_A[ (i*n)+j ] -= local_A[ (i*n)+k ] * row[j];
			}
			local_b[i] -= local_A[ (i*n)+k ] * row[n];
			local_A[ (i*n)+k ] = 0.0;
		}
	} 
	/*end division and elimination*/ 
	
		/*send local_a back to 0*/
	if( my_rank == 0 ) {
		for( int h=0; h<p; h++ ) {
			if(h != 0){
				MPI_Recv( local_A, (local_n)*n, MPI_DOUBLE, h, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			}
		for( int i=h; i<n; i+=p )
			for( int j=0; j<n; j++ )
				A[ i*n+j ] = local_A[ (i/p)*n+j ];
		}
	}
	else
		MPI_Send( local_A, (local_n)*n, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
	/*end sending local_a */
	
	/*send tmp back to 0*/
	if( my_rank == 0 ) {
		for( int h=0; h<p; h++ ) {
			if(h != 0){
				MPI_Recv( tmp, (local_n), MPI_DOUBLE, h, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			}
		for( int i=h; i<n; i+=p )
			for( int j=0; j<1; j++ )
				y[ i+j ] = tmp[ (i/p)+j ];
		}
	}
	else
		MPI_Send( tmp, (local_n), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );

	/*end sending tmp*/
	MPI_Barrier( MPI_COMM_WORLD );    
	stop = MPI_Wtime();			


  if ( my_rank==0 )
  {   
	
	 /*for (int i = 0; i < n; i++){
		cout << y[i] << " ";
	 }*/
	 print(A, n, n);
     print_solution( A, y, n );
	 cout << endl;
     cout << "n = " << n << "  p = " << p;
     cout << "  x[" << n-1 << "] = " << y[n-1] << endl;
     cout << "  Total Time = " << stop-start << endl;
  }
  MPI_Finalize();
}
