#ifndef param_h
#define param_h

#include "TTree.h"
#include "TFile.h"

//5,15,35,70,126,210,330,495,715,1001
/*const int nParam=495;
const int nMat=8;
const int nXf=8;
const int nXpf=5;
const int nYf=5;
const int nYpf=4;*/

//const int nParam=210;
//const int nMat=6;
//const int nXf=6;
//const int nXpf=5;
//const int nYf=4;
//const int nYpf=4;

// --- For path length ---- //
//const int nParam=35; 
//const int nMat=3;
//const int nXf=3;
//const int nXpf=3;
//const int nYf=3;
//const int nYpf=3;

const int nParam=70;
const int nMat=4;
const int nXf=4;
const int nXpf=4;
const int nYf=4;
const int nYpf=4;

/*const int nParam=70;
const int nMat=4;
const int nXf=4;
const int nXpf=4;
const int nYf=4;
const int nYpf=4;*/
//const int Maxevent=2000;
//const int Maxevent=3000;
//const int Maxevent=4000;
//const int Maxevent=5000;
//const int Maxevent=6000;
//const int Maxevent=7000;
const int Maxevent=10000;

//4,10,20,35,56,84,120,165,220,286
//const int nfParam = 84;
//const int nfMat = 6;

void calcMatrix(int nTrig, double* yVal, double* eVal, 
		double* xf, double* xpf,
		double* yf, double* ypf,
		char* Mname);
int readvalue(TFile *f, int nMax, double fVal[][4], double tVal[][8]);
bool DCana(int N, TTree* tree);

#endif
