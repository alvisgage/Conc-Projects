#include <iostream>
#include <string.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <cstdlib>
#include "mpi.h"

using namespace std;

int parent(int i){
	return i / 2;}
	
int left(int i){
	return 2 * i;}
	
int right(int i){
	return 2 * i + 1;}
	
struct heapObj{
	int val;
	int source;};
	
void min_heapify(heapObj * h, int size, int i){
	int l = left(i);
	int r = right(i);
	int smallest;
	if (l < size && h[l].val < h[i].val){
		smallest = l;}
	else{
		smallest = i;}
	if (r < size && h[r].val < h[smallest].val){
		smallest = r;}
	if (smallest != i){
		heapObj temp = h[i];
		h[i] = h[smallest];
		h[smallest] = temp;
		min_heapify(h, size, smallest);}}
		
void printHeap(heapObj * h, int size){
	if (size <= 0){
		cout << "Heap empty" << endl;}
	else{
		for (int i = 0; i < size; i++){
			cout << "Val = " << h[i].val << " Source = " << h[i].source << endl;}}}
		
void printArr(int * arr, int size){
	if (size <= 0){
		cout << "Array empty" << endl;}
	else{
		for (int i = 0; i < size; i++){
			cout << arr[i] << " ";}
		cout << endl;}}
		
void print2dArr(int ** arr, int size, int * arrSizes){
	for (int i = 0; i < size; i++){
		for (int j = 0; j < arrSizes[i]; j++){
			cout << arr[i][j] << " ";}}
	cout << endl;}
	
void print2dArr2(int ** arr, int size, int * arrSizes){
	for (int i = 0; i < size; i++){
		cout << i << " = ";
		for (int j = 0; j < arrSizes[i]; j++){
			cout << arr[i][j] << " ";}
		cout << endl;}}

int comparator(const void * left, const void * right){
	return ( *(int*)left - *(int*)right );}
	
void errorCheck(int error_code){
	char error_string[BUFSIZ];
	int length_of_error_string, error_class;
	MPI_Error_class(error_code, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);
	cout << error_string << endl;}

main(int argc, char* argv[]) {
	int my_rank;
	int p;
	int globalSize = atoi(argv[1]);
	int localSize;
	double start, end;
	MPI_Status status;
	int * localArr;
	int * medians = new int [p];
	srand(time(0));
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
	
	localSize = globalSize / p;
	localArr = new int[localSize];
	
	if (my_rank == 0){
		start = MPI_Wtime();
		for (int i = 0; i < localSize; i++){
			localArr[i] = rand() % 1000000;}
		if (p > 1){
			for (int j = 1; j < p; j++){
				int * sendRands = new int[localSize];
				for (int i = 0; i < localSize; i++){
					sendRands[i] = rand() % 1000000;}
				MPI_Send(sendRands, localSize, MPI_INT, j, 0, MPI_COMM_WORLD);
				delete[] sendRands;}}}
	else{
		int err = MPI_Recv(localArr, localSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);}
		
	MPI_Barrier(MPI_COMM_WORLD);
	
	qsort(localArr, localSize, sizeof(int), comparator);
	
	cout << my_rank << " localArr size " << localSize << endl;
	
	if (p > 1){
		if (my_rank == 0){
			int mediansIndex = 0;
			int * splitVals = new int[p*p];
			for (int i = 0; i < p; i++){
				splitVals[mediansIndex] = localArr[(i+1) * (localSize / p) - 1];
				mediansIndex++;}
			int * recvVals = new int [p];
			for (int i = 1; i < p; i++){
				MPI_Recv(recvVals, p, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				for (int j = 0; j < p; j++){
					splitVals[mediansIndex] = recvVals[j];
					mediansIndex++;}}
			qsort(splitVals, p * p, sizeof(int), comparator);
			for (int i = 0; i < p; i++){
				medians[i] = splitVals[(i * p) + p - 1];}
			delete[] splitVals;
			delete[] recvVals;}
		else{
			int * sendVals =  new int [p];
			for (int i = 0; i < p; i++){
				sendVals[i] = localArr[(i+1) * (localSize / p) - 1];}
			MPI_Send(sendVals, p, MPI_INT, 0, 0, MPI_COMM_WORLD);
			delete[] sendVals;}

		MPI_Bcast(medians, p, MPI_INT, 0, MPI_COMM_WORLD);
		
		if (my_rank == 0){
			cout << "Medians" << endl;
			printArr(medians, p);}

		int ** sendArrs = new int* [p];
		int ** recvArrs = new int* [p];
		int * sendCounts = new int [p];
		int * recvCounts = new int [p];
		for (int i = 0; i < p; i++){
			sendCounts[i] = 0;}
		int currProc = 0;
		
		//Find sendCounts
		for (int i = 0; i < localSize; i++){
			if (localArr[i] <= medians[currProc]){
				sendCounts[currProc]++;}
			else if (currProc < p - 1){
				currProc++;
				while (localArr[i] > medians[currProc] && currProc < p){
					currProc++;}
				sendCounts[currProc]++;}
			else{
				sendCounts[p-1]++;}}
				
		delete[] medians;
		
		//Allocate sendArrs
		for (int i = 0; i < p; i++){
			cout << "Alloc start sendArrs " << i << " by " << my_rank << " of size " << sendCounts[i] << endl;
			sendArrs[i] = new int[sendCounts[i]];
			cout << "Alloc finish sendArrs " << i << " by " << my_rank << endl;}
			
		cout << my_rank << " allocated sendArrs" << endl;
		
		//Set sendArrs values
		int localArrIndex = 0;
		for (int i = 0; i < p; i++){
			for (int j = 0; j < sendCounts[i]; j++){
				sendArrs[i][j] = localArr[localArrIndex];
				localArrIndex++;}}
				
		cout << my_rank << " set sendArrs vals" << endl;
						
		delete[] localArr;
		
		//Send and receive counts for arrays to be exchanged
		for (int i = 0; i < p; i++){
			for (int j = 0; j < p; j++){
				if (my_rank == i){
					if (i == j){
						recvCounts[my_rank] = sendCounts[my_rank];
						recvArrs[my_rank] =  new int [sendCounts[my_rank]];
						for (int x = 0; x < sendCounts[my_rank]; x++){
							recvArrs[my_rank][x] = sendArrs[my_rank][x];}}
					else{
						MPI_Send(&sendCounts[j], 1, MPI_INT, j, 0, MPI_COMM_WORLD);}}
				else{
					if (my_rank == j){
						MPI_Recv(&recvCounts[i], 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);}}}}
		
		//Send and receive arrays
		for (int i = 0; i < p; i++){
			for (int j = 0; j < p; j++){
				if (my_rank == i){
					if (my_rank != j){
						MPI_Send(sendArrs[j], sendCounts[j], MPI_INT, j, 0, MPI_COMM_WORLD);}}
				else if (my_rank == j){
					recvArrs[i] =  new int [recvCounts[i]];
					MPI_Recv(recvArrs[i], recvCounts[i], MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);}}}
		
		delete[] sendCounts;
		delete[] sendArrs;
			
		int resultArrSize = 0;
		for (int i = 0; i < p; i++){
			resultArrSize += recvCounts[i];}
		int * resultArr =  new int [resultArrSize];
		int * recvIndexes =  new int [p];
		for (int i = 0; i < p; i++){
			recvIndexes[i] = 0;}
		
		if (my_rank == 0){
			end = MPI_Wtime();
			cout << "Final array for " << p << " after " << end - start << endl;}
		cout << my_rank << " final array count = " << resultArrSize << endl;
		
		heapObj * h = new heapObj [p];
		int hCurSize = 0; //max is p
		
		//Insert the first value from each array into heap
		for (int i = 0; i < p; i++){
			if (recvCounts[i] > 0){
				heapObj x;
				x.val = recvArrs[i][0];
				x.source = i;
				h[hCurSize] = x;
				recvIndexes[i]++;
				hCurSize++;}}
			
		//Build the heap
		for (int i = hCurSize / 2; i >= 0; i--){
			min_heapify(h, hCurSize, i);}
		
		for (int i = 0; i < resultArrSize; i++){
			heapObj smallest = h[0];
			if (recvIndexes[smallest.source] < recvCounts[smallest.source]){
				heapObj newest;
				newest.val = recvArrs[smallest.source][recvIndexes[smallest.source]];
				newest.source = smallest.source;
				h[0] = newest;
				recvIndexes[smallest.source]++;
				min_heapify(h, hCurSize, 0);}
			else{
				h[0] = h[hCurSize - 1];
				hCurSize--;
				min_heapify(h, hCurSize, 0);}
			resultArr[i] = smallest.val;}
			
		delete[] recvCounts;
		delete[] recvArrs;
		delete[] recvIndexes;
		delete[] h;
		
		int ** finalArrs = new int*[p-1];
		int * finalArrsCounts = new int[p-1];
		
		if (my_rank == 0){
			for (int i = 1; i < p; i++){
				MPI_Recv(&finalArrsCounts[i-1], 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);}}
		else{
			MPI_Send(&resultArrSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);}
			
		for (int i = 1; i < p; i++){
			if (my_rank == 0){
				finalArrs[i-1] = new int[finalArrsCounts[i-1]];
				MPI_Recv(finalArrs[i-1], finalArrsCounts[i-1], MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);}
			else if (my_rank == i){
				MPI_Send(resultArr, resultArrSize, MPI_INT, 0, 0, MPI_COMM_WORLD);}}}
		
	if (my_rank == 0){
		end = MPI_Wtime();
		cout << "Run time for " << p << " processors = " << end - start << endl;}
	
	MPI_Finalize();}