#include <iostream>
#include <vector>
#include <fstream>
#include "mpi.h"
#include <vector>
#include <algorithm>

using namespace std;
#define foreach BOOST_FOREACH
struct Edge
{
    int v1;
    int v2;
    int weight;
};
//check if edge is in tree
bool edgeInTree(vector<Edge> Tree, Edge e){
	bool result = false;
	for (int i = 0; i < Tree.size(); i++){
		if (Tree[i].v1 == e.v1)
			if(Tree[i].v2 == e.v2)
				if(Tree[i].weight == e.weight)
					result = true;
	}
	return result;
}
//check if vertex is in tree
bool v2InTree(vector<Edge> Tree, int v2){
	bool result = false;
	for (int i = 0; i < Tree.size(); i++)
		if (Tree[i].v2 == v2)
			result = true;
	return result;
}
//if item is in vector
bool vectorContains(vector<int> v, int item){
	bool result = false;
	for(int i = 0; i < v.size(); i++){
		if (v[i] == item)
			result = true;
	}
	return result;
}
//delete value from vector
void vectorRemove(vector<int> &v, int value){
	for (int i =0; i < v.size(); i++){
		if (v[i] == value)
			v.erase(v.begin()+i);
	}
}
//function for sort
struct myfunction{
	bool operator()(const Edge& left, const Edge& right) const{
		return left.v1 < right.v1;
	}
};
//sort tree by v1
void orderTree(vector<Edge> &Tree){
	std::sort(Tree.begin(), Tree.end(), myfunction());
}
//check if all vertices in tree
void errCheckTree(vector<Edge> &Tree){
	bool errors = false;
	for (int i = 1; i < Tree.size(); i+=2){
		if (Tree[i-1].v1 != i || Tree[i-1].v2 != i+1)
			errors = true;
	}
	if (errors)
		printf("INVALID PATH IN TREE\n");
	else
		printf("tree path valid\n");
}
//prints out entire edge list
void printData(Edge * a, int size, int rank)
{
    for(int i = 0; i < size; i++)
        cout << rank << " " << a[i].v1 << " to " << a[i].v2 << " weight: " << a[i].weight << endl;
}
//print tree
void printTree(vector<Edge> Tree){
	for(int x = 0; x < Tree.size(); x++){
		printf("from\t%d\tto\t%d\tweight\t%d\n", Tree[x].v1, Tree[x].v2, Tree[x].weight);
	}
}

//initialize graph from input file
Edge * init_graph(int & numVertices, int & numEdges)
{
    //open "graph.txt"
    fstream graphFile;
    graphFile.open("graph.txt", ios::in);

    //find number of edges and vertices
    graphFile >> numVertices;
    graphFile >> numEdges;

    //create a vector of edges
    Edge * edges = new Edge[numEdges];
    for(int i = 0; i < numEdges; i++)
    {
        graphFile >> edges[i].v1;
        graphFile >> edges[i].v2;
        graphFile >> edges[i].weight;
    }
	//printData(edges, numEdges, 0);
    return edges;
}

//initialize root array
int * initRoot(int numVertices)
{
    int * root = new int[numVertices];
    for(int i = 0; i < numVertices; i++)
        root[i] = i;
    return root;
}

//initialize size array
int * initSize(int numVertices)
{
    int * size = new int[numVertices];
    for(int i = 0; i < numVertices; i++)
        size[i] = 1;
    return size;
}
//return root of vertex
int findRoot(int a, int * root)
{
    return root[a];
}

//merge the smaller component into the larger component
void unionVertices(int a, int b, int * size, int * root, int numVertices)
{
    if(size[a] < size[b])
    {
        int tmp = root[a];
        for(int i = 0; i < numVertices; i++)
        {
            if(root[i] == tmp)
                root[i] = b;
        }
        size[b] += size[a];
    }
    else
    {
        int tmp = root[b];
        for(int i = 0; i < numVertices; i++)
        {
            if(root[i] == tmp)
                root[i] = a;
        }
        size[a] += size[b];
    }
}

//perform Boruvka to find mst
Edge * boruvkaMST(Edge * edges, int * size, int * root, int numVertices, int numEdges)
{
	printf("beginning of boruvka\n");
    //array of edges to store the minimum weight edges
    Edge * minWeightEdges = new Edge[numVertices];

    //create mstTree to store mst edges
    Edge * mstTree = new Edge[numVertices-1];

    int mstTreeSize = 0;
    do{
        //initialize every minimum weight edge with -1
        for(int i = 0; i < numVertices; i++)
            minWeightEdges[i].weight = -1;

        //go through every edge to find minimum weight edge
        for(int i = 0; i < numEdges; i++)
        {
            //find the roots of the vertices of the edges
            int k = findRoot(edges[i].v1, root);
            int j = findRoot(edges[i].v2, root);

            //store new minimum weight edge of vertex i and j if weight is smaller than current edge weight
            if(k != j)
            {
                if(minWeightEdges[k].weight == -1 || edges[i].weight < minWeightEdges[k].weight)
                    minWeightEdges[k] = edges[i];

                if(minWeightEdges[j].weight == -1 || edges[i].weight < minWeightEdges[j].weight)
                    minWeightEdges[j] = edges[i];
            }
        }
        
        //go through all the minWeightEdges
        for(int i = 0; i < numVertices; i++)
        {
            Edge e = minWeightEdges[i];
            if(e.weight != -1)
            {
                int u = findRoot(e.v1, root);
                int v = findRoot(e.v2, root);

                if(u != v)
                {
                    //add edge to MST tree
                    mstTree[mstTreeSize] = e;
                    
                    //merge vertices u and v together
                    unionVertices(u, v, size, root, numVertices);
                    
                    mstTreeSize++;
                }
            }
        }
    }while(mstTreeSize < numVertices-1);
    printf("end of boruvka\n");
    return mstTree;
}

//assign which processsor each vertex should go
int * assignVerticesToProc(int numVertices, int p)
{
    int * vertexProcLocation = new int[numVertices];
    int localNumVertices;
    int counter = 0;
    for(int i = 0; i < p; i++)
    {
        localNumVertices = numVertices/p;
            
        for(int j = 0; j < localNumVertices; j++)
        {
            vertexProcLocation[counter] = i;
            counter++;
        }
    }
    
    return vertexProcLocation;
}
//assign number of edges each processor receives
int * assignNumEdgesPerProc(Edge * edges, int numVertices, int numEdges, int p)
{
    //find how high of a vertex number each proc will get
    int * numVerticesPerProc = new int[p];
    int localNumVertices = 0;
    for(int i = 0; i < p; i++)
    {
        localNumVertices += numVertices/p;
        numVerticesPerProc[i] = localNumVertices;
    }
    
    int * numEdgesPerProc = new int[p];
    int pCounter = 0;
    int edgeCounter = 0;
    for(int i = 0; i < numEdges; i++)
    {
        if(edges[i].v1 <= numVerticesPerProc[pCounter])
            edgeCounter++;
        else
        {
            numEdgesPerProc[pCounter] = edgeCounter;
            pCounter++;
            edgeCounter = 1;
        }
    }
    numEdgesPerProc[pCounter] = edgeCounter;
    
    return numEdgesPerProc;
}
int main(int argc, char* argv[])
{
    //MPI initializations
    int rank;
    int p;
	double start, stop;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
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
   
    
    Edge * edges;
    Edge * localEdges;
    int numVertices;
    int numEdges;
    int localNumEdges;
    int * root;
    int * size;
    
    int * vertexProcLocation;
    int * numEdgesPerProc;

    //1. INITIALIZE GRAPH AND DISTRIBUTE BETWEEN PROCESSORS
    if(rank == 0)
    {
		start = MPI_Wtime();
        //initialize edge list
        edges = init_graph(numVertices, numEdges);
        root = initRoot(numVertices);
        size = initSize(numVertices);

        //assign which processor each vertex should go
        //e.g. vertexProcLocation[1] = 0 means vertex 1 goes to process 0
        vertexProcLocation = assignVerticesToProc(numVertices, p);

        //assign number of edges each process should get
        numEdgesPerProc = assignNumEdgesPerProc(edges, numVertices, numEdges, p);
    }
    
    //send the number of edges each processor is going to get
    MPI_Scatter(numEdgesPerProc, 1, MPI_INT, &localNumEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);   
    
    //divide edges between processors
    if(rank == 0)
    {
        int counter = 0;
        for(int i = 0; i < p; i++)
        {
            MPI_Send(&edges[counter], numEdgesPerProc[i], edgeT, i, 0, MPI_COMM_WORLD);
            counter += numEdgesPerProc[i];
        }
    }
    localEdges = new Edge[localNumEdges];
    MPI_Recv(localEdges, localNumEdges, edgeT, 0, 0, MPI_COMM_WORLD, &status);
    //printf("%d\tnum edges:\t%d\n", rank, localNumEdges);

    
	//each process gets min edge for its assigned vertices
		//create list of each vertex in that edge
		//for each in list, get min edge
	vector<int> vertices;
	for (int i=0; i<localNumEdges;i++){
		int v = localEdges[i].v1;
		if (!vectorContains(vertices, v)){
			//v not found
			vertices.insert(vertices.begin(), 1, v);
		}
	}	
	

	
	Edge * minEdge = new Edge[vertices.size()];
	int minEdgeCounter = 0;
	for(int j = 0; j < vertices.size(); j++){
		int v = vertices[j];
		Edge tempEdge;
		tempEdge.v1 = -1;
		tempEdge.v2 = -1;
		tempEdge.weight = -1;
		minEdge[minEdgeCounter] = tempEdge;
		for (int i=0; i<localNumEdges; i++){
			if(localEdges[i].v1==v){
				if(minEdge[minEdgeCounter].weight == -1 || localEdges[i].weight < minEdge[minEdgeCounter].weight)
					{
						minEdge[minEdgeCounter].v1 = localEdges[i].v1;
						minEdge[minEdgeCounter].v2 = localEdges[i].v2;
						minEdge[minEdgeCounter].weight = localEdges[i].weight;
					}
			}
		}
		minEdgeCounter++;
	}

	//send these edges to 0
		//assume every proc has the same number of vertices
	MPI_Barrier(MPI_COMM_WORLD);
	//printData(minEdge, vertices.size(), rank);
	MPI_Send(&minEdge[0], vertices.size(), edgeT, 0, 1, MPI_COMM_WORLD);

	if(rank==0){
		vector<Edge> Tree;		
		Edge * otherProcEdges = new Edge[vertices.size()];
		for (int i=0; i < p; i++){
			//printf("0 receiving...\n");
			MPI_Recv(otherProcEdges, vertices.size(), edgeT, i, 1, MPI_COMM_WORLD, &status);
			//printData(otherProcEdges, vertices.size(), i);
			//printf("0 recvd %d...\n", i);
			//0 goes through proc in order and adds edges to tree IF 1. is not alread in tree and v2 is not already in tree
			for (int j = 0; j < vertices.size(); j++){
					//write these functions
				if (!edgeInTree(Tree, otherProcEdges[j]) && !v2InTree(Tree, otherProcEdges[j].v2)){
					Tree.push_back(otherProcEdges[j]);
					//printf("\nadding edge...\n");
					//printTree(Tree);
				}
			}
		}		
		//printTree(Tree);
		//printf("0 added edges to tree...\n");
	//get list of v1's that are not in edge list (ignoring last vertex)
	////start with list of all vertices except last, for each edge in tree, if v1 is in vertex list, remove it
	//get list of v2's that are not in edge list (ignoring first vertex)
	////start with list of all vertices except first, for each edge in tree, if v2 is in vertex list, remove it
		vector<int> v1sNotInList;
		vector<int> v2sNotInList;
		for (int i=0; i< numVertices; i++){
			if(i != numVertices -1)
				v1sNotInList.insert(v1sNotInList.begin(), 1, i);
			if(i > 0)
				v2sNotInList.insert(v2sNotInList.begin(), 1, i);
		}
		for (int i = 0; i < Tree.size(); i++){
			if(vectorContains(v1sNotInList, Tree[i].v1))
				vectorRemove(v1sNotInList, Tree[i].v1);
			if(vectorContains(v2sNotInList, Tree[i].v2))
				vectorRemove(v2sNotInList, Tree[i].v2);
		}
		//printf("got v1 and v2 list...\n");
		for(int j = 0; j < v2sNotInList.size(); j++){
			int v = v2sNotInList[j];
			for(int i=0; i < numEdges; i++){
				if(edges[i].v2 == v && vectorContains(v1sNotInList, edges[i].v1)){
					Tree.insert(Tree.begin(), sizeof(Edge), edges[i]);
				}
			}
		}
		stop = MPI_Wtime();
		orderTree(Tree);
		errCheckTree(Tree);
		printf("printing tree...\n");
		printTree(Tree);
		printf("n:\t%d\tp:\t%d\ttime:\t%d\n", numVertices, p, stop - start);
	}	
     
	//for each v in v2 list, for each edge in entire edge list, if edge.v2 = v and edge.v1 is in v1 list, add edge to tree
	//this should connect remaining components
    //printData(localEdges, localNumEdges, rank);
    cout << endl;
    
    MPI_Finalize();
}
