#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
using namespace std;
typedef struct {
        int v1;
        int v2;
        int weight;
} Edge;

struct list_edge {
        Edge edge;
        struct list_edge * next;
};

list_edge * curr, *head, *mstTree;
int number_edges = 6;
int number_vertexs;
int* root;
int* size;
void printData(list_edge * list) {
        list_edge *a = list;
         do{
                printf("%d-------------%d----------------%f\n", a->edge.v1, a->edge.v2,
                                a->edge.weight);
                a = a->next;
        }while (a->next!=NULL);
}

void init_graph() {
        FILE * pFile;
        int i;
        head = NULL;
        //read
	
        pFile = fopen("graph.txt", "r");
        fscanf(pFile, "%d", &number_vertexs);
        fscanf(pFile, "%d", &number_edges);
        root = (int*) malloc(number_vertexs * sizeof(int));
        size = (int*) malloc(number_vertexs * sizeof(int));
        head = (list_edge*) malloc(sizeof(list_edge));
        curr=head;
        for (i = 0; i < number_edges; i++) {
                list_edge * newEdge= (list_edge*) malloc(sizeof(list_edge));
                fscanf(pFile, "%d", &curr->edge.v1);
                fscanf(pFile, "%d", &curr->edge.v2);
                fscanf(pFile, "%d", &curr->edge.weight);
                newEdge->next=NULL;
                curr->next = newEdge;
                curr=newEdge;
        }

        //print data
		
        //printData(head);

        //init root

        for (i = 0; i < number_vertexs; i++) {
                root[i] = i;
                size[i]=1;
        }
		//printf("done reading file\n");
}
int findRoot(int i) {
        //find root
        while (i != root[i]) {
                i = root[i];
        }
        return i;
}
int less(Edge e, Edge f) {
        return (e.weight < f.weight);

}
void unionEdge(int i, int j) {

        if (i == j) {
                return;
        }

        // make smaller root point to larger one
        if (size[i] < size[j]) {
                root[i] = j;
                size[j] += size[i];
        } else {
                root[j] = i;
                size[i] += size[j];
        }
}
void boruvkaMst(int rank) {
        list_edge *tree;
        Edge* closest = (Edge *) malloc(number_vertexs * sizeof(Edge));
		Edge* tmp = (Edge *) malloc(number_vertexs * sizeof(Edge));
        mstTree = (list_edge*) malloc(sizeof(list_edge));
        tree=mstTree;
        int mstTreeSize = 0;
        int countEdge=0;
        int countLoop=0;
		
        do {
                for (int i = 0; i < number_vertexs; i++) {
                        closest[i].weight= -1;
                }
                curr = head;
                //find closest edge
                while (curr->next!=NULL) {
                        int i = findRoot(curr->edge.v1);
                        int j = findRoot(curr->edge.v2);
                        if (i == j) {
                                curr = curr->next;
                                continue;
                        }
                        if (closest[i].weight == -1
                                        ||curr->edge.weight < closest[i].weight) {
                                closest[i] = curr->edge;
                        }
                        if (closest[j].weight == -1
                                        || curr->edge.weight < closest[j].weight) {
                                closest[j] = curr->edge;
                        }
                        curr = curr->next;
                }
                countEdge=0;
                for (int i = 0; i < number_vertexs; i++) {
                        Edge e = closest[i];
                        if (e.weight != -1) {
                                int u = findRoot(e.v1);
                                int v = findRoot(e.v2);
                                if (u!= v) {
                                        countEdge++;
                                        list_edge* newEdge = (list_edge*) malloc(sizeof(list_edge));
                                        tree->edge=e;
                                    tree->next = newEdge;
                                        tree=newEdge;
                                        mstTreeSize++;
                                        unionEdge(u, v);
                                }
                        }
                }
                countLoop++;
        }
        while(mstTreeSize<number_vertexs-1);
   printf("number loop :%d\n",countLoop);
}
int main(int argc, char *argv[]) {
		int p, rank;
			MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	//create a MPI Struct Datatype for Edge
    const int edgeItems = 3;
    int blocklengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Datatype edgeT;
    MPI_Aint offsets[3];
    offsets[0] = offsetof(Edge, v1);
    offsets[1] = offsetof(Edge, v2);
    offsets[2] = offsetof(Edge, weight);
    MPI_Type_create_struct(edgeItems, blocklengths, offsets, types, &edgeT);
    MPI_Type_commit(&edgeT);
	MPI_Status status;
		if (rank==0){
			init_graph();
			//printf("proc 0 sending\n");
			//send list
			curr = head;
			while(curr->next != NULL){
				printf("rank: %d\tv1: %d\tv2: %d\tweight: %d\n", rank, curr->edge.v1, curr->edge.v2, curr->edge.weight);
				curr = curr->next;
			}
			for (int i = 1; i < p; i++){
				MPI_Send(head, 2, edgeT, i, 1, MPI_COMM_WORLD);
			}	
		}
		else{
			head = (list_edge*) malloc(sizeof(list_edge));
			printf("proc %d receiving\n", rank);
			//recieve list
			Edge e;
			MPI_Recv(head, 2, edgeT, 0, 1, MPI_COMM_WORLD, &status);
			curr = head;
			while(curr->next != NULL){
				printf("rank: %d\tv1: %d\tv2: %d\tweight: %d\n", rank, curr->edge.v1, curr->edge.v2, curr->edge.weight);
				curr = curr->next;
			}
			
		}
		//printf("-------------------------------------------------------\n");
       // boruvkaMst(rank);
        //printData(mstTree);

		
			MPI_Finalize();
			
        return 0;
}

