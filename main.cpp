#include "Image.h"
#include "utilities.h"


#include "matchdist.h"
#include "entropyminimize.h"
#include "qpprogram.h"

int main(int argc, char * argv[])
{
	clock_t start = clock();

	//entropy minimization
	EntropyMinimize(argc,argv);

	//matching distribution - KL or Bhattachayya
	//MatchDist(argc,argv);

	//Quadratic Programming
	//QPprogram(argc, argv);
	
	clock_t finish = clock();
	cout<<"main function time "<<(double)(finish-start)/CLOCKS_PER_SEC<<" seconds!"<<endl;
	return -1;
}