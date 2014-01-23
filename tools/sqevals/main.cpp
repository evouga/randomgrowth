#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
	if(argc != 2)
	{
		cerr << "Usage: sqevals (spectrum file)" << endl;
		return -1;
	}
	
	ifstream ifs(argv[1]);
	if(!ifs)
	{
		cerr << "Cannot open spectrum file " << argv[1] << endl;
		return -1;
	}
	int nummodes;
	ifs >> nummodes;
	for(int i=0; i<nummodes && ifs; i++)
	{
		double sqeval;
		ifs >> sqeval;
		if(i != 0)
			cout << " ";
		cout << sqeval;
	}
	cout << endl;
	if(!ifs)
	{
		cerr << "Error reading spectrum file " << argv[1] << endl;
		return -1;
	}
	return 0;
}
		
