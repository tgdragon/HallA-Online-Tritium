//#include "pathl_opt.hh"
#include <fstream>
#include <iostream>
//#include <TMatrixD.h>
#include <math.h>

using namespace std;
//const int nMatT=6;
//const int nXf=6;
//const int nXpf=5;
//const int nYf=4;
//const int nYpf=4;

// -----3rd order -----
const int nMatT=3; 
const int nXf=3;
const int nXpf=3;
const int nYf=3;
const int nYpf=3;

//////////////////////////////////////////////////
double calcf2t(double* P, double xf, double xpf, 
                 double yf, double ypf)
//////////////////////////////////////////////////
{

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0;
  
  for (int n=0;n<nMatT+1;n++){
  	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){
		
		if (a+b+c+d==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf){
		    x = pow(xf,double(a))*pow(xpf,double(b))*
		      pow(yf,double(c))*pow(ypf,double(d));
		  }
		  else{
		    x = 0.;
		  }
		  Y += x*P[npar];
		  npar++;
		}
		
	      }
	    }
	  }
  	}
  }
  
  return Y;
}

/*
//////////////////////////////////////////////////
double calct2s(double* P, double xpt, double ypt, double mom)
//////////////////////////////////////////////////
{
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0;
  
  for (int n=0;n<nMatS+1;n++){
    for (c=0;c<n+1;c++){
      for (b=0;b<n+1;b++){
  	for (a=0;a<n+1;a++){
	  
	  if (a+b+c==n){
	    if (a<=nXpt && b<=nYpt && c<=nMom ){
	      x = pow(xpt,double(a))*pow(ypt,double(b))*
		pow(mom,double(c));
	    }
	    else{
	      x = 0.;
	    }
	    Y += x*P[npar];
	    npar++;
	  }
	  
  	}
      }
    }
  }
  
  return Y;
}
*/
