#ifndef recon_h
#define recon_h

#include <TFile.h>

const int nHoleX=9;
const int nHoleY=6;
//5,15,35,70,126,210,330,495,715,1001
/*const int nParamT=495;
const int nMatT=8;
const int nXf=8;
const int nXpf=8;
const int nYf=8;
const int nYpf=8;*/

/*
const int nParamT=210;
const int nMatT=6;
const int nXf=6;
const int nXpf=5;
const int nYf=4;
const int nYpf=4;
*/

// --- For path length ---- //
const int nParamT=35; 
const int nMatT=3;
const int nXf=3;
const int nXpf=3;
const int nYf=3;
const int nYpf=3;

/*const int nParamT=70;
const int nMatT=4;
const int nXf=4;
const int nXpf=4;
const int nYf=4;
const int nYpf=4;*/

//4,10,20,35,56,84,120,165,220,286
const Int_t nParamS = 35;
const Int_t nMatS = 4;
const Int_t nXpt = 4;
const Int_t nYpt = 4;
const Int_t nMom = 4;

int readvalue(TFile *f, int nMax, double fVal[][4], 
					double tVal[][6], double sVal[][4]);
double calcf2t(double* P, double xf, double xpf, 
                 double yf, double ypf);
double calct2s(double* P, double xpt, double ypt, double mom);

#endif
