#include <iostream>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;
bool isSquare(long int sum);
main(int argc, char* argv[])
{
	ofstream outputFile("numbers.txt");

	for (int i = 0; i < 1600; i++){
		outputFile << rand() % RAND_MAX << endl;	
	}
}
