#include <stdlib.h>
#include "mpi.h"
#include <math.h>
using namespace std;

int p,my_rank, n; 
void subst(int my_rank,int col,double *local_a, double *local_b, double *row){
	int pos = col/p+1;
	double tmp = 0.0;
	if (my_rank > col%p){
		pos--;
	}
    for(int i=pos; i<n/p; i++ ) {	
      for(int j=col+1; j<n; j++ ) {
		tmp = local_a[ i*n+col ] * row[j];
        local_a[ i*n+j ] -= tmp;
      }
      local_b[i] -= local_a[ i*n+col ] * row[n];
    }
}
main( int argc, char *argv[] )
{
  double *A, *b, *local_a, *local_b, *row;
  double start, stop;
  int r, local_n;
    int recv, send, pos;
  
  n = atoi(argv[1]);
	if( my_rank==0 ) {
		A = new double[n*n];	
		b = new double[n];
	}
	
	/*values were being solved as nan, inf, -inf*/
	/*this should set up values to be integers*/
	for( int i=0; i<n; i++ ) {	
		b[i] = 0.0;
		for( int j=0; j<n; j++ ) {	
			r = rand() % n+1;		
			A[i*n+j] = r;
			b[i] += j*r;
		}
	}

  MPI_Init (&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&p);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	local_n = n/p;
	/*print matrix and b array for testing*/
	/*if(my_rank==0){
		cout << "mat A" << endl;
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				cout << A[i*n+j] << "\t";
			}
			cout << endl;
		}
		cout << "b: " << endl;
		for (int i = 0; i < n; i++){
			cout << b[i] << "\t";			
		}
		cout << endl;		
	}*/

  start = MPI_Wtime();	
  /*init local arraya*/
  ////////////////////////////////////////////////////////////////////////////
  local_a = new double[ n*local_n ];		
  local_b = new double[ local_n ];
  row = new double[n+1];
  ////////////////////////////////////////////////////////////////////////////

  /*distribute A to proc*/
  ///////////////////////////////////////////////////////////////////////
  if( my_rank == 0 ) {
    for( int h=p-1; h>=0; h-- ) {
      for( int i=h; i<n; i+=p )
        for( int j=0; j<n; j++ ) {
          local_a[ (i/p)*n+j ] = A[ i*n+j ];
        }
      if( h > 0 ) 
        MPI_Send( local_a, local_n*n, MPI_DOUBLE, h, 1, MPI_COMM_WORLD );
    }
  }
  else{
    MPI_Recv( local_a, local_n*n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	}
////////////////////////////////////////////////////////////////////////////////
  /*distribute b to proc*/
  ///////////////////////////////////////////////////////////////////////////////
  if( my_rank == 0 ) {
    for( int h=p-1; h>=0; h-- ) {
      for( int i=h; i<n; i+=p )
        for( int j=0; j<1; j++ ) {
          local_b[ (i/p)+j ] = b[ i+j ];
        }
      if( h > 0 ) 
        MPI_Send( local_b, local_n, MPI_DOUBLE, h, 1, MPI_COMM_WORLD );
    }
  }
  else{
    MPI_Recv( local_b, local_n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
  ///////////////////////////////////////////////////////////////////////////////
  /*elimination*/
  ///////////////////////////////////////////////////////////////////////////////
	/*predecessor and successor ranks to send row*/
  recv = (p+(my_rank-1)) % p;
  send = (my_rank+1) % p;
  for(int col=0; col<n; col++ ) {
    if(my_rank != col%p) {			
      MPI_Recv( row, n+1, MPI_DOUBLE, recv, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      if( send != col%p ) {	
        MPI_Send( row, n+1, MPI_DOUBLE, send, 2, MPI_COMM_WORLD );
      }
    }
	else{
		for(int j=col+1; j<n; j++ )	{				
			local_a[ (col/p)*n+j ] = local_a[ (col/p)*n+j ] / local_a[ (col/p)*n+col ];
		}
		local_b[ col/p ] = local_b[ col/p ] / local_a[ (col/p)*n+col ];
		local_a[ (col/p)*n+col ] = 1.0;
		for(int j=0; j<n; j++ ){
			row[j] = local_a[ (col/p)*n+j ];
		}
		row[n] = local_b[ col/p ];
		MPI_Send( row, n+1, MPI_DOUBLE, send, 2, MPI_COMM_WORLD );
	}
	///////////////////////////////////////////////////////////////////////////
	/*substitution*/
	///////////////////////////////////////////////////////////////////////////	
	subst(my_rank, col, local_a, local_b, row);
  } 
  //////////////////////////////////////////////////////////////////
  /*sending local_a back to 0*/
  /////////////////////////////////////////////////////////////////////////
  if( my_rank == 0 ) {
    for( int h=0; h<p; h++ ) {
      if( h > 0 ) 
        MPI_Recv( local_a, local_n*n, MPI_DOUBLE, h, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      for( int i=h; i<n; i+=p )
        for( int j=0; j<n; j++ )
          A[ i*n+j ] = local_a[ (i/p)*n+j ];
    }
  }
  else{
    MPI_Send( local_a, local_n*n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
	}
  ////////////////////////////////////////////////////////////////////////
  /*sending local_y back to proc 0*/
  ///////////////////////////////////////////////////////////////////////////
  if( my_rank == 0 ) {
    for( int h=0; h<p; h++ ) {
      if( h > 0 ) 
        MPI_Recv( local_b, (local_n), MPI_DOUBLE, h, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      for(int i=h; i<n; i+=p )
        for(int j=0; j<1; j++ )
          b[ i+j ] = local_b[ (i/p)+j ];
    }
  }
  else{
    MPI_Send( local_b, (local_n), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
	}
  ///////////////////////////////////////////////////////////////////////////
  
  
  MPI_Barrier( MPI_COMM_WORLD );       
  stop = MPI_Wtime();			


  if ( my_rank==0 )
  {
	/*back substitution once array is solved*/
    for(int i=n-1; i>=0; i-- ) {
		for(int j=n-1; j>i; j-- ){
			b[i] -= b[j] * A[i*n+j];
		}
	}
	for( int i=0; i<n; i+=4 ) {
		for(int j=i; j<i+4 && j<n; j++ )
			cout << "b[" << j << "] = " << round(b[j]) << "  ";
		cout << endl;
	}	
     cout << "n = " << n << "  p = " << p << endl;;
     cout << stop-start << endl;
  }
  MPI_Finalize();
}
