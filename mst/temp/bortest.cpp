#include <iostream>
#include <vector>
#include <fstream>
#include "mpi.h"

using namespace std;

struct Edge
{
    int v1;
    int v2;
    int weight;
};

//prints out entire edge list
void printData(Edge * a, int size, int rank)
{
    for(int i = 0; i < size; i++)
        cout << rank << " " << a[i].v1 << " to " << a[i].v2 << " weight: " << a[i].weight << endl;
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
        if(i < numVertices%p)
            localNumVertices = numVertices/p + 1;
        else
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
        if(i < numVertices%p)
            localNumVertices += numVertices/p + 1;
        else
            localNumVertices += numVertices/p;
        numVerticesPerProc[i] = localNumVertices;
    }
    
    int * numEdgesPerProc = new int[p];
    int pCounter = 0;
    int edgeCounter = 0;
    for(int i = 0; i < numEdges; i++)
    {
        if(edges[i].v1 < numVerticesPerProc[pCounter])
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
    
    //timing of runtime
    double start, stop, runtime;
    
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
    
    //2. IMPLEMENT BARUVKA
    

    printData(localEdges, localNumEdges, rank);
    cout << endl;
    
    MPI_Finalize();
}
