/////////////////////////////
// HES Backword Matrix making Program w/ SVD 
// readvalue.c
//
// D.Kawama
//////////////////////////////

#include <fstream>
#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TCut.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <math.h>
//#include "recon.h"
using namespace std;

////////////////////////////////////////////////////////////
int readvalue(TFile *f, int nMax, double fVal[][4], 
	      double tVal[][7], double sVal[][4])
////////////////////////////////////////////////////////////
{
  int nentries=0;
  int nTrig=0;
  
 TTree *tree = (TTree*)f->Get("tree");
 nentries =(int)tree->GetEntries();
 //nentries = 10000;
 TRandom3* rnd = new TRandom3();
 
 int jj=0;
 double XFP=0.,XpFP=0.,YFP=0.,YpFP=0.,Mom=0.;
 double XSPLP=0.,YSPLP=0.;
 double Xpt=0.,Ypt=0.,LFP=0.,LEDC=0.,Xt=0.,Yt=0., Zt=0.0;
 double XSS=0.,YSS=0.;
 double PL=0.0;
 double XpSS=0.,YpSS=0.;
 double XQ1i=0.,YQ1i=0.;
 double XQ2i=0.,YQ2i=0.;
 double XDi=0.,YDi=0.;
 double XDe=0.,YDe=0.;
 double XFP1=0.,YFP1=0.;
 double XpFP1=0.,YpFP1=0.;
 char SPLTrig=0,KTrig=0;
// const double  XFPm=0., XpFPm=0.; 
// const double  YFPm=0., YpFPm=0.; 
// const double  Xptm=0., Yptm=0., Momm=0.;
// const double  XFPr=25., XpFPr=0.1; 
// const double  YFPr=3., YpFPr=0.01; 
// const double  Xptr=1., Yptr=1., Momr=1.;
 TCut trig="ESPLTrig==1";
 
 //for DCana
 const double sigma=0.021;
 const double FPpos=0.0;
 const int nLayerX=16;
 const int nLayerY=8;
 const int nFit=2;
 double XEDC[nLayerX];
 double YEDC[nLayerY];
 int max = 100;
 double xdet[max], ydet[max], xpdet[max], ypdet[max];
 double xptar[max], yptar[max];
 //double ztar[max], momentum[max], pathl[max];
 double momentum[max], pathl[max];
 double ztarR[max], ztarL[max];
 double s2nhit;
 int LR_flag = 2;
 
 //tree->SetBranchAddress("EXSPLP",&XSPLP);
 //tree->SetBranchAddress("EYSPLP",&YSPLP);
 //tree->SetBranchAddress("EXSS",&XSS);
 //tree->SetBranchAddress("EXpSS",&XpSS);
 //tree->SetBranchAddress("EYSS",&YSS);
 //tree->SetBranchAddress("EYpSS",&YpSS);
 //tree->SetBranchAddress("EXQ1i",&XQ1i);
 //tree->SetBranchAddress("EYQ1i",&YQ1i);
 //tree->SetBranchAddress("EXDi",&XDi);
 //tree->SetBranchAddress("EYDi",&YDi);
 //tree->SetBranchAddress("EXDe",&XDe);
 //tree->SetBranchAddress("EYDe",&YDe);
 //tree->SetBranchAddress("EXFP1",&XFP1);
 //tree->SetBranchAddress("EYFP1",&YFP1);
 //tree->SetBranchAddress("EXpFP1",&XpFP1);
 //tree->SetBranchAddress("EYpFP1",&YpFP1);
 //tree->SetBranchAddress("ELFP3",&LFP);
 /*tree->SetBranchAddress("XEDCFP",&XFP);
 tree->SetBranchAddress("XpEDCFP",&XpFP);
 tree->SetBranchAddress("YEDCFP",&YFP);
 tree->SetBranchAddress("YpEDCFP",&YpFP);*/
 tree->SetBranchAddress("R.tr.vz",&ztarR);
 tree->SetBranchAddress("L.tr.vz",&ztarL);
 if(LR_flag==2){
   tree->SetBranchAddress("R.tr.x" ,&xdet); // dispersive plane
   tree->SetBranchAddress("R.tr.th",&xpdet);
   tree->SetBranchAddress("R.tr.y", &ydet); // non-dispersive plane
   tree->SetBranchAddress("R.tr.ph",&ypdet);
   tree->SetBranchAddress("R.tr.tg_th",&xptar);
   tree->SetBranchAddress("R.tr.tg_ph",&yptar);
   tree->SetBranchAddress("R.tr.p",&momentum);
   tree->SetBranchAddress("R.tr.pathl",&pathl);
 }
 else{
   tree->SetBranchAddress("L.tr.x" ,&xdet); // dispersive plane
   tree->SetBranchAddress("L.tr.th",&xpdet);
   tree->SetBranchAddress("L.tr.y", &ydet); // non-dispersive plane
   tree->SetBranchAddress("L.tr.ph",&ypdet);
   tree->SetBranchAddress("L.tr.tg_th",&xptar);
   tree->SetBranchAddress("L.tr.tg_ph",&yptar);
   tree->SetBranchAddress("L.tr.p",&momentum);
   tree->SetBranchAddress("L.tr.pathl",&pathl);
 }
 
 
 //tree->SetBranchAddress("R.s2.nthit",&s2nhit);
 //tree->SetBranchAddress("KXFP",&XFP);
 //tree->SetBranchAddress("KXpFP",&XpFP);
 //tree->SetBranchAddress("KYFP",&YFP);
 //tree->SetBranchAddress("KYpFP",&YpFP);
 //tree->SetBranchAddress("KXt",&Xt);
 //tree->SetBranchAddress("KYt",&Yt);
 //tree->SetBranchAddress("KXpt",&Xpt);
 //tree->SetBranchAddress("KYpt",&Ypt);
 //tree->SetBranchAddress("KMom",&Mom);
 //tree->SetBranchAddress("LEDCFP",&LEDC);
 //tree->SetBranchAddress("ESPLTrig",&SPLTrig);
// for (int k=0;k<nLayerX;k++){
//   char str[100];
//   sprintf(str,"XEDC%d",k+1);
//   tree->SetBranchAddress(str,&XEDC[k]);
// }
// for (int k=0;k<nLayerY;k++){
//   char str[100];
//   if (k%2==0){
//     sprintf(str,"YEDC%d",2*k+3);
//   }
//   else {
//     sprintf(str,"YEDC%d",2*k+2);
//   }
//   tree->SetBranchAddress(str,&YEDC[k]);
// }
 bool goodflag = false;
 for (int i=0; i<nentries; i++){
   tree->GetEntry(i);
   goodflag = false;
   if(  //SPLTrig&& 
	jj < nMax 
	//&& Mom>0.0 
	//&& fabs(YFP)<40.0 // centselect (TG, 2017/8/24)
	//&& fabs(XFP)<4.0  // centselect (TG, 2017/8/24)
	/*&& XEDC[0]>-1000. 
	&& XEDC[1]>-1000. 
	&& XEDC[2]>-1000. 
	&& XEDC[3]>-1000. 
	&& XEDC[4]>-1000. 
	&& XEDC[5]>-1000. 
	&& XEDC[6]>-1000. 
	&& XEDC[7]>-1000. 
	&& XEDC[8]>-1000. 
	&& XEDC[9]>-1000. 
	&& XEDC[10]>-1000. 
	&& XEDC[11]>-1000. 
	&& XEDC[12]>-1000. 
	&& XEDC[13]>-1000. 
	&& XEDC[14]>-1000. 
	&& XEDC[15]>-1000.*/ 
      ){
//     Double_t z1[5];
//     Double_t z2[5];
//     Double_t z3[3];
//     Double_t z4[3];
//     for (int n=0; n<5 ; n++){
//       z1[n] =  -15. + 7.5*n -0.375 ;
//       z2[n] =  -15. + 7.5*n +0.375 ;
//     }
//     for (int n=0; n<3 ; n++){
//       z3[n] =  -1.905 + 0.635*n + 31.;
//       z4[n] =  0.3175 + 0.635*n + 31.;
//     }
//     Double_t zxdc[nLayerX]={z1[0],z2[0],z1[1],z2[1],z1[2],z2[2],
//			     z1[3],z2[3],z1[4],z2[4],
//			     z3[0],z3[1],z3[2],z4[0],z4[1],z4[2]};
//     Double_t zydc[nLayerY]={z1[1],z2[1],z1[3],z2[3],
//			     z3[0],z3[1],z4[1],z4[2]};
//     for (int k=0;k<nLayerX;k++){
//       Double_t s = sigma;
//       if (k==2 || k==3 || k==6 || k==7){
//	 s = sigma/cos(30.*3.14/180.);
//       }
//       else if (k==10 || k==11 || k==14 || k==15){
//	 s = sigma/cos(60.*3.14/180.);
//       }
//       XEDC[k]=gRandom->Gaus(XEDC[k],s);//[cm]
//     }	
//     for (int k=0;k<nLayerY;k++){
//       Double_t s = sigma;
//       if (k<=4){
//	 s = sigma/sin(30.*3.14/180.);
//       }
//       else {
//	 s = sigma/sin(60.*3.14/180.);
//       }
//       YEDC[k]=gRandom->Gaus(YEDC[k],s);//[cm]
//     }
//     
//     TVectorD x;x.Use(nLayerX,XEDC);
//     TVectorD y;y.Use(nLayerY,YEDC);
//     TVectorD zx;zx.Use(nLayerX,zxdc);
//     TVectorD zy;zy.Use(nLayerY,zydc);
//     TMatrixD Ax(nLayerX,nFit);
//     TMatrixDColumn(Ax,0) = 1.0;
//     TMatrixDColumn(Ax,1) = zx;
//     TDecompSVD svdx(Ax);
//     TMatrixD Ay(nLayerY,nFit);
//     TMatrixDColumn(Ay,0) = 1.0;
//     TMatrixDColumn(Ay,1) = zy;
//     TDecompSVD svdy(Ay);
//     Bool_t ok;
//     
//     const TVectorD x_svd = svdx.Solve(x,ok);
//     const TVectorD y_svd = svdy.Solve(y,ok);
//     /*XFP = x_svd(0)+FPpos*x_svd(1);
//       YFP = y_svd(0)+FPpos*y_svd(1);
//       XpFP = atan(x_svd(1));
//       YpFP = atan(y_svd(1));*/
//     
     
     
     // =================================== //
     // ======= Detector resolutions ====== //
     // =================================== //
     //XFP  = rnd->Gaus(XFP,0.010);   // sigma = 100 um
     //XpFP = rnd->Gaus(XpFP,0.5e-3); // sigma = 0.5 mrad
     //YFP  = rnd->Gaus(YFP,0.010);   // sigma = 100 um
     //YpFP = rnd->Gaus(YpFP,0.5e-3); // sigma = 0.5 mrad
     //XFP  = rnd->Gaus(XFP,0.010);   // sigma = 100 um
     //XpFP = rnd->Gaus(XpFP,0.05e-3); // sigma = 0.05 mrad
     //YFP  = rnd->Gaus(YFP,0.010);   // sigma = 100 um
     //YpFP = rnd->Gaus(YpFP,0.05e-3); // sigma = 0.05 mrad
     XFP   =  xdet[0];
     XpFP  =  xpdet[0];
     YFP   =  ydet[0];
     YpFP  =  ypdet[0];
     Mom   =  momentum[0];
     Xpt   =  xptar[0];
     if(LR_flag==2) Ypt   =  yptar[0] + 13.2/180. * 3.14159;
     else Ypt   =  yptar[0] - 13.2/180. * 3.14159;
     Zt    =  (ztarR[0]+ztarL[0])/2.0;
     //Zt    =  ztarR[0];
     //cout << "zzz" << Zt << endl;
     PL    =  pathl[0];
     
     if(Mom>1.3 && Mom<3.0 
	&& fabs(Zt)<0.25 
	&& fabs(Xpt)<0.1
	//&& s2nhit == 1
	&& fabs(ztarL[0]-ztarR[0])<0.05
	){
       goodflag = true;
     }
    
     else goodflag = false;
     
     if(goodflag==true){
       //cout << "zzz" << Zt << endl;
       // ---- y correction -------------------------------- //
       // ---- This corrections are used for p, x', y' ----- //
       //if(LR_flag == 2){
       //	 YFP  = YFP  - 0.11 * Zt;
       //	 YpFP = YpFP + 0.10 * Zt;
       //}
       //else{
       //	 YFP  = YFP  + 0.12 * Zt;
       //	 YpFP = YpFP - 0.09 * Zt;
       //}
       
       
       fVal[jj][0]=XFP;//[cm]
       fVal[jj][1]=XpFP;//[rad]
       fVal[jj][2]=YFP;//[cm]
       fVal[jj][3]=YpFP;//[rad]
       
       tVal[jj][0]=Xpt;  //[rad]
       tVal[jj][1]=Ypt;  //[rad]
       tVal[jj][2]=Mom;  //[GeV/c]
       tVal[jj][3]=PL; //[cm]
       tVal[jj][4]=Xt;   //[cm]
       tVal[jj][5]=Yt;   //[cm]
       tVal[jj][6]=Zt; //[cm]
       
       //sVal[jj][0]=XSS;//[cm]
       //sVal[jj][1]=XpSS;//[rad]
       //sVal[jj][2]=YSS;//[cm]
       //sVal[jj][3]=YpSS;//[rad]
       jj++;
     }
   }
 }
 nTrig=jj;
 
 return nTrig;
 
}

