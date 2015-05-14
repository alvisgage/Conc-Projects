#include <iostream>
#include "mpi.h"
using namespace std;
void Par_dfs(STACK_T local_stack, MPI_Comm comm);
void Service_requests(STACK_T local_stack, MPI_Comm comm);
int main(int argc, char* argv[]){
	NODE_T root;
	MPI_Comm comm;
	STACK_T local_stack;
	NODE_T* node_list;
	NODE_T node;
	int p;
	int my_rank;
	
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &my_rank);
	
	if(my_rank==0){
		Generate(root, &node_list, p);
	}
	
	Scatter(node_list, &node, p);
	Initialize(node, &local_stack);
	
	do{
		//search for a while
		Par_dfs(local_stack, com);
		
		//service request for work
		Service_requests(local_stack, comm);
		
		//if local_stack isn't empty, return
		//if local_stack is empty, send requests until we get work or recv msg terminating program
	}while (Work_remains(local_stack, comm));
	
	Update_solution(comm);
	Free(local_stack);
	Print_solution(comm);
	
}

void Par_dfs(STACK_T local_stack, MPI_Comm comm){
	int count;
	NODE_T node;
	float temp_solution;
	
	//search local subtree for a while
	count = 0;
	while(!Empty(local_stack) && (count < MAX_WORK)){
		node = Pop(local_stack);
		if (Solution(node)){
			temp_solution = Evaluate(node);
			if(temp_solution < Best_solution(comm)){
				Local_solution_update(temp_solution, node);
				Bcast_solution(comm);
			}
		} else if (Feasible(node, comm)){
			Expand(node, local_stack);
		}
		count++;
	}	
}

void Service_requests(STACK_T local_stack, MPI_Comm comm){
	STACK_T send_stack;
	int destination;
	
	while(Work_requests_pending(comm)){
		destination = Get_dest(comm);
		if(Nodes_available(local_stack)){
			Split(local_stack, &send_stack);
			Send_work(destination, send_stack, comm);
		} else {
			Send_reject(desination, comm);
		}
	}
}

int Work_remains(STACK_T local_stack, MPI_Comm comm){
	int work_available;
	int request_sent = FALSE;
	int work_request_process;
	
	if (!Empty(local_stack)){
		return TRUE;
	} else{
		Out_of_work(comm);
		while(TRUE){
			
		}
	}
}