#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;
int main(int argc, char* argv[])
{	
	int n = atoi(argv[1]);
	ofstream graphFile;
    graphFile.open("newgraph.txt", ios::trunc);
	
	graphFile << n << endl;
	graphFile << ((n+1)*n)/2 << endl;
	
	for(int i =1; i <= n; i++)
		for (int j = i+1; j <= n; j++){
			graphFile << i << "\t";
			graphFile << j << "\t";
			graphFile << i+j * 2 << endl;
		}
}