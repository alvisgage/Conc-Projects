#include <iostream>
#include <fstream>
#include "mpi.h"

using namespace std;

//each edge contains its root, v1, v2, and weight
struct Edge
{
    int v1;
    int v2;
    int weight;
};

//prints out entire edge list
void printData(Edge * a, int size)
{
    for(int i = 0; i < size; i++)
        cout << a[i].v1 << " " << a[i].v2 << " " << a[i].weight << endl;
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

//stores number of edges per vertex
int * assignNumEdgesPerVertex(Edge * edges, int numVertices, int numEdges)
{
	int * numEdgesPerVertex = new int[numVertices];
	int counter = 0;
	for(int i = 0; i < numEdges; i++)
	{
		if(i == 0 || edges[i].v1 == edges[i-1].v1)
			counter++;
		else
		{
			numEdgesPerVertex[edges[i-1].v1] = counter;
			counter = 1;
		}
	}
	numEdgesPerVertex[edges[numEdges-1].v1] = counter;
	
	return numEdgesPerVertex;
}

//stores number of vertices each processor has
int * assignNumVerticesPerProc(int numVertices, int p)
{
    int * numVerticesPerProc = new int[p];
    for(int i = 0; i < p; i++)
    {
        if(i < numVertices%p)
            numVerticesPerProc[i] = numVertices/p + 1;
        else
            numVerticesPerProc[i] = numVertices/p;
    }
    return numVerticesPerProc;    
}

//assign which processor each vertex is located
int * assignVertexProcLocation(int * numVerticesPerProc, int numVertices, int p)
{
    int * vertexProcLocation = new int[numVertices];
    int counter = 0;
    for(int i = 0; i < p; i++)
    {
        for(int j = 0; j < numVerticesPerProc[i]; j++)
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

//return root of vertex
int findRoot(int a, int * root)
{
    return root[a];
}

//return processor location of root of vertex
int findLocation(int a, int * vertexProcLocation, int * root)
{
	return vertexProcLocation[root[a]];
}

void sortEdgesBasedOnLocation(Edge * edges, int numEdges, int * vertexProcLocation, int * root)
{
	for(int i = 0; i < numEdges; i++)
	{
		int x = i;
		for(int j = i; j < numEdges; j++)
		{
			int a = findLocation(edges[x].v2, vertexProcLocation, root);
			int b = findLocation(edges[j].v2, vertexProcLocation, root);
			
			if(a > b)
				x = j;
		}
		Edge temp = edges[i];
		edges[i] = edges[x];
		edges[x] = temp;
	}
}
//find the minimum weight edge for each vertex
Edge * findMinWeightEdges(Edge * localEdges, int localNumVertices, int localNumEdges, int startVertex, int * root)
{
    //initialize edge weights to -1 (meaning no minimum weight edge was found)
    Edge * minWeightEdges = new Edge[localNumVertices];
    for(int i = 0; i < localNumVertices; i++)
        minWeightEdges[i].weight = -1;

    //go through every edge to find the minimum weight edge
    for(int i = 0; i < localNumEdges; i++)
    {
        //find the roots of the vertices of the edges
        int v1 = findRoot(localEdges[i].v1, root);
        int v2 = findRoot(localEdges[i].v2, root);
        
        //store minimum weight edge into minWeightEdges (v1-startVertex)
        if(v1 != v2) //make sure edge is non-cycling
        {
            if(minWeightEdges[v1-startVertex].weight == -1 || localEdges[i].weight < minWeightEdges[v1-startVertex].weight)
                minWeightEdges[v1-startVertex] = localEdges[i];
        }
    }
    
    return minWeightEdges;
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
	
	//process 0 owns these variables and these don't change (except adding to the mstTree)
	Edge * edges;
	Edge * mstTree;
	int numVertices;
	int numEdges;
	int mstTreeSize;
	
	//everyone should have these (only the root vertices are changed)
	int * size;
	int * numEdgesPerVertex;
	int * vertexProcLocation;
	
	//everyone should have these (every vertex is changed)
	int * root;
	
	//information relevant to local processor
	Edge * localEdges;
    Edge * localMinWeightEdges;
    int localNumVertices;
    int localNumEdges;
    int startVertex;            //each processor starts with this vertex
	
	//information for sending joining vertices together
	Edge * tmpLocalMinWeightEdges;
	int * numSendVertices;
	int * numSendEdges;
	int * numReceiveVertices;
	int * numReceiveEdges;
	int totalNumReceive;
	int * sdispls;
	int * rdispls;
	int totalSend;
	int totalReceive;
	int * sendTheseVertex;
	Edge * sendTheseEdges;
	int counter;
	int counter2;
	
	//things that are only used for distributing vertices through the processes (gets deleted)
	int * numEdgesPerProc;
	int * numVerticesPerProc;
	int * info;
	int * localInfo;
	
	if(rank == 0)
	{
		//initialize graph
		edges = init_graph(numVertices, numEdges);
		
		//arrays for distributing the edges to processors
		numVerticesPerProc = assignNumVerticesPerProc(numVertices, p);
		numEdgesPerProc = assignNumEdgesPerProc(edges, numVertices, numEdges, p);
		
		//arrays that each process should have
		root = initRoot(numVertices);
		size = initSize(numVertices);
		numEdgesPerVertex = assignNumEdgesPerVertex(edges, numVertices, numEdges);
		vertexProcLocation = assignVertexProcLocation(numVerticesPerProc, numVertices, p);
		
		//mstTree
		mstTree = new Edge[numVertices-1];
		mstTreeSize = 0;
	}
	
	//send information to processors
	MPI_Bcast(&numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		info = new int[3*p];
		int startVertex = 0;
		for(int i = 0; i < p; i++)
		{
			info[3*i] = numEdgesPerProc[i];
			info[3*i+1] = numVerticesPerProc[i];
			info[3*i+2] = startVertex;
			startVertex += numVerticesPerProc[i];
		}
	}
	localInfo = new int[3];
    MPI_Scatter(info, 3, MPI_INT, localInfo, 3, MPI_INT, 0, MPI_COMM_WORLD);   
    localNumEdges = localInfo[0];
    localNumVertices = localInfo[1];
    startVertex = localInfo[2];
	
	//distribute edges between processors
	if(rank == 0)
    {
        counter = 0;
        for(int i = 0; i < p; i++)
        {
            MPI_Send(&edges[counter], numEdgesPerProc[i], edgeT, i, 0, MPI_COMM_WORLD);
            counter += numEdgesPerProc[i];
        }
    }
    localEdges = new Edge[localNumEdges];
    MPI_Recv(localEdges, localNumEdges, edgeT, 0, 0, MPI_COMM_WORLD, &status);
		
	//delete useless information
	if(rank == 0)
	{
		delete [] numEdgesPerProc;
		delete [] numVerticesPerProc;
		delete [] info;
	}
	delete [] localInfo;
	
	//broadcast information about all vertices to processors
	if(rank != 0)
	{
		root = new int[numVertices];
		size = new int[numVertices];
		numEdgesPerVertex = new int[numVertices];
		vertexProcLocation = new int[numVertices];
	}
	MPI_Bcast(root, numVertices, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(size, numVertices, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(numEdgesPerVertex, numVertices, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vertexProcLocation, numVertices, MPI_INT, 0, MPI_COMM_WORLD);
	
	//find minimum weight edges
	localMinWeightEdges = findMinWeightEdges(localEdges, localNumVertices, localNumEdges, startVertex, root);
	
	sortEdgesBasedOnLocation(localMinWeightEdges, localNumVertices, vertexProcLocation, root);
	
	//find out how many vertices each process is sending/receiving to each process
	numSendVertices = new int[p];
	numSendEdges = new int[p];
	numReceiveVertices = new int[p];
	numReceiveEdges = new int[p];
	for(int i = 0; i < p; i++)
	{
		numSendVertices[i] = 0;
		numSendEdges[i] = 0;
	}
	for(int i = 0; i < localNumVertices; i++)
	{
		int x = findLocation(localMinWeightEdges[i].v2, vertexProcLocation, root);
		numSendVertices[x]++;
	}	
	MPI_Alltoall(numSendVertices, 1, MPI_INT, numReceiveVertices, 1, MPI_INT, MPI_COMM_WORLD);
	
	//alltoall send of vertices
	totalNumReceive = 0;
	sdispls = new int[p];
	rdispls = new int[p];
	for(int i = 0; i < p; i++)
	{
		totalNumReceive += numReceiveVertices[i];
		if(i == 0)
		{
			sdispls[i] = 0;
			rdispls[i] = 0;
		}
		else
		{
			sdispls[i] = sdispls[i-1] + numSendVertices[i-1];
			rdispls[i] = rdispls[i-1] + numReceiveVertices[i-1];
		}
	}
	tmpLocalMinWeightEdges = new Edge[totalNumReceive];
	MPI_Alltoallv(localMinWeightEdges, numSendVertices, sdispls, edgeT, tmpLocalMinWeightEdges, numReceiveVertices, rdispls, edgeT, MPI_COMM_WORLD);
	
	//mark which vertices to send and find out how many edges receiving
	totalReceive = 0;
	totalSend = 0;
	sendTheseVertex = new int[localNumVertices];
	counter = 0;
	for(int i = 0; i < localNumVertices; i++)
	{
		int tmp1v1 = findRoot(localMinWeightEdges[i].v1, root);
		int tmp1v2 = findRoot(localMinWeightEdges[i].v2, root);
		for(int j = 0; j < totalNumReceive; j++)
		{
			int tmp2v1 = findRoot(tmpLocalMinWeightEdges[j].v1, root);
			int tmp2v2 = findRoot(tmpLocalMinWeightEdges[j].v2, root);
			if(tmp1v1 == tmp2v2 && tmp1v2 == tmp2v1)	//make sure edges match
			{
				if(size[tmp1v1] > size[tmp1v2] || (size[tmp1v1] == size[tmp1v2] && tmp1v1 < tmp1v2))	//v1 absorbs v2 (receive edges from v2)
				{
					int x = findLocation(tmp1v2, vertexProcLocation, root);
					numReceiveEdges[x] += numEdgesPerVertex[tmp1v2];
					totalReceive += numEdgesPerVertex[tmp1v2];
				}
				else	//v2 absorbs v1 (send edges to v2)
				{
					int x = findLocation(tmp1v2, vertexProcLocation, root);
					numSendEdges[x] += numEdgesPerVertex[tmp1v1];
					totalSend += numEdgesPerVertex[tmp1v1];
					sendTheseVertex[counter] = tmp1v1;
					counter++;
				}
			}
		}
	}
	
	//create Edge array to store all edges being sent
	sendTheseEdges = new Edge[totalSend];
	counter2 = 0;
	for(int i = 0; i < localNumEdges; i++)
	{
		for(int j = 0; j < counter; j++)
		{
			if(findRoot(localEdges[i].v1, root) == sendTheseVertex[j])
			{
				sendTheseEdges[counter2] = localEdges[i];
				counter2++;
			}
		}
	}
	
	sortEdgesBasedOnLocation(sendTheseEdges, totalSend, vertexProcLocation, root);
	printData(sendTheseEdges, totalSend);
	
	//new tmp array based on edges receiving and edges sending
	//send/receive arrays
	//cleanup phase???
	//update size of vertex, where each vertex is located, root of each vertex, edges per vertex
	cout << endl << endl;
	MPI_Finalize();
}