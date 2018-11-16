#ifndef pathlopt_h
#define pathlopt_h

//#include <TFile.h>


const double  XFPm=-0.7, XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18; 
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74;
const double  XFPr=1.3, XpFPr=0.27; 
const double  YFPr=0.1, YpFPr=0.10; 
const double  Xptr=0.15, Yptr=0.08, Momr=0.18;
const double  PLm = 25.4, PLr=0.7;
const double  Ztm = -0.15, Ztr=0.35;

//const int nParamT=210; // 6th order
const int nParamT=35; // 3rd order
//const int nMatT=6;
//const int nXf=6;
//const int nXpf=5;
//const int nYf=4;
//const int nYpf=4;
double Plen_opt[nParamT];

double calcf2t(double* P, double xf, double xpf, double yf, double ypf);
double tune(double* pa, int j);
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
double tune2(double* pa, int j);
void fcn2(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
double tune3(double* pa, int j);
void fcn3(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
double tune4(double* pa, int j);
void fcn4(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);

double Calc_FPcor(double* val, double* par);
//double Calc_FPcor();

//const int tunemax = 100;
const int tunemax = 10000;
//const int tunemax = 8000;
//const int tunemax = 5000;
//const int tunemax = 2000;
double tval0[tunemax]; // x 
double tval1[tunemax]; // xp 
double tval2[tunemax]; // y 
double tval3[tunemax]; // yp
double tval4[tunemax]; // len
double betabeta[tunemax]; // beta
double timet[tunemax];  // time at target
double s2plen[tunemax]; // time at target

double tval0_L[tunemax]; // x 
double tval1_L[tunemax]; // xp 
double tval2_L[tunemax]; // y 
double tval3_L[tunemax]; // yp
double tval4_L[tunemax]; // len
double betabeta_L[tunemax]; // beta
double timet_L[tunemax];  // time at target
double s2plen_L[tunemax]; // time at target

int ntune;

const int npar_rtime_ycor = 2;
double par_rtime_ycor[npar_rtime_ycor];

//const int npar_rtime_ycor_L = 2;
const int npar_rtime_ycor_L = 2;
double par_rtime_ycor_L[npar_rtime_ycor_L];


const int npar_pathl_L_cor = 2;
double par_pathl_L_cor[npar_pathl_L_cor];

//double tval4[tunemax]; // path length

#endif
