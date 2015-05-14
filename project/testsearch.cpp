#include <iostream>
#include <list>
#include "mpi.h"
#include <cstdlib>
using namespace std;
 
int rank;
int p;
long n;
long local_n;

struct Edge
{
    long v;
    long w;
};

class Graph
{
    long V;    // No. of vertices    
	int *names; //pointer to names of vertices
    void DFSUtil(long v, bool visited[]);  // A function used by DFS
public:
    Graph(long V, long l);   // Constructor
    void addEdge(long v, long w);   // function to add an edge to graph
    void DFS(long v, int stop);    // prints DFS traversal of the complete graph
	void OrderList();
	list<long> *adj;    // Pointer to an array containing adjacency lists
	list<long> path;
	long *pathTaken;
	
};
 
Graph::Graph(long V, long l)
{
    this->V = V;
    adj = new list<long>[V];
	pathTaken = new long[l + 1];
	//names = new int[V];
	/*for (int i = 0; i < V; i++){
		names[i] = i + rank * (n/p);
	}*/
}
 
void Graph::addEdge(long v, long w)
{
	adj[v].push_back(w);
	
	/*int newV = v%(n/p);
	int newW = w%(n/p);
	if (v >= this->lowerLimit && v < this->upperLimit)
		if (w >= this->lowerLimit && w < this->upperLimit)*/
			//adj[newV].push_back(newW); // Add w to v’s list
	

}
 
void Graph::DFSUtil(long v, bool visited[])
{
    // Mark the current node as visited and print it
    visited[v] = true;
	path.push_back(v);
	pathTaken[path.size()-1] = v;
    //cout << rank << "\t" << v << "\n";
 
    // Recur for all the vertices adjacent to this vertex
    list<long>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i)
        if(!visited[*i])
            DFSUtil(*i, visited);
}
 
// The function to do DFS traversal. It uses recursive DFSUtil()
void Graph::DFS(long v, int stop = -1)
{
	//printf("rank:\t%d\tin DFS\n", rank);
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for(long i = 0; i < V; i++)
		if (stop != -1 && i > stop)
			visited[i] = true;
		else
			visited[i] = false;
 
    // Call the recursive helper function to print DFS traversal
    // starting from all vertices one by one
    /*for(int i = 0; i < V; i++)
        if(visited[i] == false)
			DFSUtil(i, visited);*/
		//while(path.size() < V){
           // DFSUtil(v, visited);
			//v += 1;
		//}
		DFSUtil(v, visited);
}
 //order adj list
 void Graph::OrderList(){	
	for(long i =0; i < V; i++){
		adj[i].sort();
	}
 }
 //print path taken
 void PrintPath(list<int> path){
	//cout << rank << endl;
	list<int>::const_iterator i;
    for(i = path.begin(); i != path.end(); i++)
		cout << *i << " ";
	cout << endl;
 }
 void ComparePaths(list<int> sequential, list<int> parallel){
	bool result = true;
	list<int>::const_iterator seq_it;
	list<int>::const_iterator par_it;
	for(seq_it = sequential.begin(), par_it = parallel.begin(); seq_it != sequential.end(); seq_it++, par_it++){
		if (*seq_it != *par_it)
			result = false;
	}
	if (result)
		printf("path successful");
	else
		printf("FAILURE!");
 }
 //if (v, w) belongs to root
 bool BelongsToRoot(int root, long v){
	bool result = false;
	while(v > root){
		if (v%2==0)
			v = (v/2) - 1;
		else
			v = ((v+1)/2) - 1;
	}
	if (v == root)
		result = true;
		
	return result;
 }
int main(int argc, char* argv[])
{
    //MPI initializations
	n = atoi(argv[1]);
	
	double start, stop;
	int *procRoots;
	long edgesPerProc;
	//long *pathTaken;
	list<int> finalPath;
	Edge * localEdges;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	//create a MPI Struct Datatype for Edge
    const int edgeItems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_LONG, MPI_LONG};
    MPI_Datatype edgeT;
    MPI_Aint offsets[2];
    offsets[0] = offsetof(Edge, v);
    offsets[1] = offsetof(Edge, w);
    MPI_Type_create_struct(edgeItems, blocklengths, offsets, types, &edgeT);
    MPI_Type_commit(&edgeT);	
	
	if (rank==0)
		local_n = p-2;
	else
		local_n = (n-((p-1)*2))/(p-1);
		
	Graph g(n, local_n);
	
	if (rank == 0){
		start = MPI_Wtime();
		// Create a graph given in the above diagram
		Edge ** edges;
		edges = new Edge*[p-1];
		//edgesPerProc = (n-(p-2))/(p-1);
		edgesPerProc = (n-((p-1)*2))/(p-1);
		procRoots = new int[p-1];
		long counter = 1;
		long tmpCounter = 0;
		long *edgeCounter = new long[p-1];
		int root = -1;
		for (int i 0; i < p-1; i++){
			edges[i-1] = new Edge[edgesPerProc];
			edgeCounter[i] = 0;
		}
		for (long i=0; i < n/2 - 1; i++, counter++){
			g.addEdge(i, i + counter);
			g.addEdge(i, i + counter + 1);
			for (int i = 1; i < p; i++){
				root = i + (p-3);
				if(BelongsToRoot(root, i+counter)){
					Edge e;
					e.v = i;
					e.w = i+counter;
					tmpCounter = edgeCounter[i-1];
					edges[i-1][tmpCounter] = e;
					tmpCounter++;
					edgeCounter[i-1] = tmpCounter;
				}
			}
		}
		for (int i = 1; i < p; i++){
			/*edges[i-1] = new Edge[edgesPerProc];
			root = i + (p-3);
			procRoots[i-1] = root;
			counter = 0;
			for (long j = root; j < n; j++){
				list<long>::iterator k;
				for(k = g.adj[j].begin(); k != g.adj[j].end(); k++){
					//edge (v, w) = (j, k)
					if (BelongsToRoot(root, j)){
						Edge e;
						e.v = j;
						e.w = *k;
						edges[i-1][counter] = e;
						counter++;
					}
				}
			}*/
			/*for (int j = 0; j < edgesPerProc; j++){
				printf("proc:\t%d\tgets\t(%d, %d)\n",i,edges[i-1][j].v, edges[i-1][j].w);
			}*/
			//printf("\nsending edges to proc:\t%d\n", i);
			MPI_Send(&edges[i-1][0], edgesPerProc, edgeT, i, 0, MPI_COMM_WORLD);
		}
		if (p>1){
			g.DFS(0, procRoots[p-2]);
		}
		else{
			g.DFS(0);
		}
	}
	else{
		localEdges = new Edge[local_n];
		//printf("proc:\t%d\treceiving from 0\n", rank);
		MPI_Recv(localEdges, local_n, edgeT, 0, 0, MPI_COMM_WORLD, &status);

		for (long i = 0; i < local_n; i++){
			//printf("rank:\t%d\tedge.v\t%d\tedge.w\t%d\n",rank,localEdges[i].v, localEdges[i].w); 
			g.addEdge(localEdges[i].v, localEdges[i].w);
		}
		g.DFS(rank + (p-3));
		/*pathTaken = new long[local_n + 1];
		list<long>::iterator i;
		long count = 0;
		for(i = g.path.begin(); i != g.path.end(); i++, count++){
			pathTaken[count] = *i;
		}*/
		//printf("rank:\t%d\tcount:\t%d\n",rank,count);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		//g.PrintPath();
		long ** paths;
		paths = new long*[p-1];
		for (int i = 1; i < p; i++){
			paths[i-1] = new long[edgesPerProc + 1];
			//receive edges, print
			//printf("receiving %d edges from:\t%d\n",edgesPerProc + 1, i);
			MPI_Recv(paths[i-1], edgesPerProc + 1, MPI_LONG, i, 2, MPI_COMM_WORLD, &status);
			//printf("done receiving from:\t%d\n", i);
		}
		//print proc0 path until root is reach, print proc list that has root
		//the following is only needed to display/verify data
		//comment out this section to reduce time
		/*
		list<int>::const_iterator it;
		int rootIndex = 0;
		for(it = g.path.begin(); it != g.path.end(); it++){
			if (*it == procRoots[rootIndex]){
				for (int i = 0; i < edgesPerProc + 1; i++){
					//cout << paths[rootIndex][i] << "[not0]\t";
					finalPath.push_back(paths[rootIndex][i]);
				}
				rootIndex++;
			} else {
				//cout << *it << "\t";
				finalPath.push_back(*it);
			}
		}
		while(rootIndex < (p-1)){
			for (int i = 0; i < edgesPerProc + 1; i++){
					//cout << paths[rootIndex][i] << "\t";
					finalPath.push_back(paths[rootIndex][i]);
			}
			rootIndex++;
		}
		*/
		stop = MPI_Wtime();
		printf("n:\t%d\tp:\t%d\ttime:\t%f\n", n, p, stop - start);
		//PrintPath(finalPath);
	} else {
		//send
		//printf("rank\t%d\tsending %d edges to proc 0\n", rank, g.path.size());
		//printf("rank\t%d\tpath\t%d\n", rank, g.pathTaken[0]);
		MPI_Send(&g.pathTaken[0], local_n + 1, MPI_LONG, 0, 2, MPI_COMM_WORLD);
		//printf("rank\t%d\tdone sending to 0\n", rank);
	}

	/*if (rank==0 && p > 1)
		g.DFS(0, p-3);
	else if (rank == 0)
		g.DFS(0);
	else
		g.DFS(rank + (p-3));*/

	//only add edge if v1 is in assigned vertex list
	//rank * n/p to rank + 1 * n/p
    //cout << "Following is Depth First Traversal (starting from vertex 2) \n";
	//g.OrderList();
	//g.DFS(0);
	//g.PrintPath();
	//run full  search again and compare paths to make sure it ran successfully - TAKES TOO LONG
	/*if (rank==0 && p > 1){
		local_n = n;
		g.path.erase(g.path.begin(), g.path.end());
		g.DFS(0);
		//PrintPath(g.path);
		//PrintPath(finalPath);
		ComparePaths(g.path, finalPath);
	}*/
	
	MPI_Finalize();
	
	//to do:
	//DONE split edges up and send them to correct proc 
	////we know what root the proc will call, just need to find out if other edges are dependent on that root
	////to find root: v = v/2 - 1; recurs until v = root or v < root
	//each proc should add edges to graph (will have to keep track of names some how, perhaps when adding edges, add the name to an array, then resend this array in order)
	////to keep track of names, use adjacency list array, adj[newV].pushback(name)
	//send array back to 0
	//0 should just output in order
	////if theres time, figure out way to print partial 0, then 1, then 0, then 2, etc
	///////0 keeps an array of all roots, print until you get to a root, then print that proc's array
 
    return 0;
}
