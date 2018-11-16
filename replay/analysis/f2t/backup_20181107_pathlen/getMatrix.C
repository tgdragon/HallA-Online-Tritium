////////////////////////////////////////////
// HES Backword Matrix making Program w/ SVD 
//
// D.Kawama
////////////////////////////////////////////

#include "getMatrix.h"
#include <fstream>
#include <iostream>
#include <TMatrixD.h>
#include <TFile.h>
#include <TTree.h>
using namespace std;


int main()
{
  double tVal[Maxevent][8];
  double fVal[Maxevent][4];
  double len[Maxevent];
  double elen[Maxevent];
  double mom[Maxevent];
  double emom[Maxevent];
  double xpt[Maxevent];
  double expt[Maxevent];
  double ypt[Maxevent];
  double eypt[Maxevent];
  double xf[Maxevent];
  double xpf[Maxevent];
  double yf[Maxevent];
  double ypf[Maxevent];
  double zt[Maxevent];
  int nTrig=0;

 char name_Mmom[100];
 char name_Mxpt[100];
 char name_Mypt[100];
 char name_Mlen[100];
 char name_Mzt[100];
 char name_ifile[100];
  TFile* file;
  TTree* tree;
 
  ifstream input("parMatrix.dat");
  input >> name_Mmom;
  input >> name_Mxpt;
  input >> name_Mypt;
  input >> name_Mlen;
  input >> name_Mzt;
  input >> name_ifile;
  input.close();
  
  file = new TFile(name_ifile);
  tree = (TTree*)file->Get("tree");

  cout << "reading file..."<<endl;
  nTrig=readvalue(file,Maxevent,fVal,tVal);
  cout << "read "<<nTrig<<" events"<<endl;
  for (int i=0;i<Maxevent;i++){
    xpt[i]  = tVal[i][0]; 
    ypt[i]  = tVal[i][1]; 
    mom[i]  = tVal[i][2]; 
    len[i]  = tVal[i][3]; 
    zt[i]   = tVal[i][4]; 
    
    xf[i]   = fVal[i][0]; 
    xpf[i]  = fVal[i][1]; 
    yf[i]   = fVal[i][2]; 
    ypf[i]  = fVal[i][3]; 
  }
  //cout << "Momentum-reconstruction matrix is being calculated..." <<endl;
  //calcMatrix(nTrig, mom, emom, xf, xpf, yf, ypf, name_Mmom);
  //cout << "x'-reconstruction matrix is being calculated..." <<endl;
  //calcMatrix(nTrig, xpt, expt, xf, xpf, yf, ypf, name_Mxpt);
  //cout << "y'-reconstruction matrix is being calculated..." <<endl;
  //calcMatrix(nTrig, ypt, eypt, xf, xpf, yf, ypf, name_Mypt);
  cout << "Path-length reconstruction matrix is being calculated..." <<endl;
  calcMatrix(nTrig, len, elen, xf, xpf, yf, ypf, name_Mlen);
  //cout << "Zt reconstruction matrix is being calculated..." <<endl;
  //calcMatrix(nTrig, zt, elen, xf, xpf, yf, ypf, name_Mzt);
  return 0;
}

