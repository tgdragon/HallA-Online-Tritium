#include "recon.h"
#include <fstream>
#include <iostream>
#include <TMatrixD.h>
#include <math.h>
using namespace std;


//////////////////////////////////////////////////
double calcf2t(double* P, double xf, double xpf, 
                 double yf, double ypf)
//////////////////////////////////////////////////
{

  double Y=0.;
  double x=1.; 
  Int_t npar=0;
  Int_t a=0,b=0,c=0,d=0;
  
  for (Int_t n=0;n<nMatT+1;n++){
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

//////////////////////////////////////////////////
double calcf2t_3rd(double* P, double xf, double xpf, 
		   double yf, double ypf)
//////////////////////////////////////////////////
{

  double Y=0.;
  double x=1.; 
  Int_t npar=0;
  Int_t a=0,b=0,c=0,d=0;
  
  for (Int_t n=0;n<nMatT+1;n++){
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

//////////////////////////////////////////////////
double calcf2t_4th(double* P, double xf, double xpf, 
		   double yf, double ypf)
//////////////////////////////////////////////////
{

  double Y=0.;
  double x=1.; 
  Int_t npar=0;
  Int_t a=0,b=0,c=0,d=0;
  
  for (Int_t n=0;n<nMatT_4+1;n++){
  	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){
		
		if (a+b+c+d==n){
		  if (a<=nXf_4 && b<=nXpf_4 && c<=nYf_4 && d<=nYpf_4){
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

//////////////////////////////////////////////////
double calct2s(double* P, double xpt, double ypt, double mom)
//////////////////////////////////////////////////
{

  double Y=0.;
  double x=1.; 
  Int_t npar=0;
  Int_t a=0,b=0,c=0;
  
  for (Int_t n=0;n<nMatS+1;n++){
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

