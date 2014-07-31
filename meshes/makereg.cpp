#include <iostream>
#include <cstdlib>

using namespace std;

int main()
{
	for(int i=0; i<10; i++)
	{
		for(int j=0; j<10; j++)
		{
			double x = -1.0 + 2.0*double(i)/double(9);
			double y = -1.0 + 2.0*double(j)/double(9);
			cout << "v " << x << " " << y << " 0" << endl;
		}
	}

	for(int i=0; i<9; i++)
	{
		for(int j=0; j<9; j++)
		{
			if(rand() % 1)
			{
				int idx1 = 1+i+10*j;
				cout << "f " << idx1 << " " << idx1+1 << " " << idx1+10 << endl;
				cout << "f " << idx1+10 << " " << idx1+1 << " " << idx1+11 << endl;
			}
			else
			{
				int idx1 = 1+i+10*j;
				cout << "f " << idx1 << " " << idx1+11 << " " << idx1+10 << endl;
				cout << "f " << idx1+11 << " " << idx1 << " " << idx1+1 << endl;
			}
		}
	}
}
