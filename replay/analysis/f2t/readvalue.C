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
#include <TRandom.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <math.h>
using namespace std;

////////////////////////////////////////////////////////////
int readvalue(TFile *f, int nMax, double fVal[][4], double tVal[][8])
////////////////////////////////////////////////////////////
{
 int nentries=0;
 int nTrig=0;

 TTree *tree = (TTree*)f->Get("tree");
 nentries =(int)tree->GetEntries();
 TRandom rnd;
  
 int jj=0;
 double XFP=0.,XpFP=0.,YFP=0.,YpFP=0.,Mom=0.;
 double XSPLP=0.,YSPLP=0.;
 double Xpt=0.,Ypt=0.,LFP=0.,LEDC=0.;
 double XQ1i=0.,YQ1i=0.;
 double XQ2i=0.,YQ2i=0.;
 double XDi=0.,YDi=0.;
 double XDe=0.,YDe=0.;
 double XFP1=0.,YFP1=0.;
 double XpFP1=0.,YpFP1=0.;
 char SPLTrig=0,KTrig=0;
 double Zt = 0.0; // z at target 
 double PL = 0.0; // Path length
 const double  XFPm=-0.7, XpFPm=-0.15; 
 const double  YFPm=-0.05, YpFPm=-0.18; 
 const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74;
 const double  XFPr=1.3, XpFPr=0.27; 
 const double  YFPr=0.1, YpFPr=0.10; 
 const double  Xptr=0.15, Yptr=0.08, Momr=0.18;
 const double  PLm = 25.4, PLr=0.7;
 const double  Ztm = -0.15, Ztr=0.35;
 TCut trig="EDCTrig==1";
 
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
 double ztar[max], momentum[max], pathl[max];
 double s2nhit;
 
 // ======= Right HRS =====================
 //tree->SetBranchAddress("R.tr.x" ,&xdet); // dispersive plane
 //tree->SetBranchAddress("R.tr.th",&xpdet);
 //tree->SetBranchAddress("R.tr.y", &ydet); // non-dispersive plane
 //tree->SetBranchAddress("R.tr.ph",&ypdet);
 //tree->SetBranchAddress("R.tr.tg_th",&xptar);
 //tree->SetBranchAddress("R.tr.tg_ph",&yptar);
 //tree->SetBranchAddress("R.tr.p",&momentum);
 //tree->SetBranchAddress("R.tr.pathl",&pathl);
 //tree->SetBranchAddress("R.tr.vz",&ztar);
 //tree->SetBranchAddress("R.s2.nthit",&s2nhit);
 
 // ======= Left HRS =====================
 tree->SetBranchAddress("L.tr.x" ,&xdet); // dispersive plane
 tree->SetBranchAddress("L.tr.th",&xpdet);
 tree->SetBranchAddress("L.tr.y", &ydet); // non-dispersive plane
 tree->SetBranchAddress("L.tr.ph",&ypdet);
 tree->SetBranchAddress("L.tr.tg_th",&xptar);
 tree->SetBranchAddress("L.tr.tg_ph",&yptar);
 tree->SetBranchAddress("L.tr.p",&momentum);
 tree->SetBranchAddress("L.tr.pathl",&pathl);
 tree->SetBranchAddress("L.tr.vz",&ztar);
 tree->SetBranchAddress("L.s2.nthit",&s2nhit);
 
 

 bool goodflag = false;
 if (nMax>0){ 
   for (int i=0; i<nentries; i++){
     goodflag = false;
     tree->GetEntry(i);
     if( jj<nMax 
	 //&& SPLTrig == 1  
	 //&& fabs(XFP1)<53./2.
	 //&& fabs(YFP)<6.
	 //&& fabs(XFP)<15.
	 // && fabs(XFP)<4.0   // centselect (TG, 2017/8/24)
	 //&& fabs(YFP)<50.0  // centselect (TG, 2017/8/24)
	 ){
       //XFP=gRandom->Gaus(XFP,0.001);//[cm]
       //YFP=gRandom->Gaus(YFP,0.002);//[cm]
       //XpFP=gRandom->Gaus(XpFP,fabs(XpFP)*1e-4);//[rad]
       //YpFP=gRandom->Gaus(YpFP,fabs(YpFP)*1e-4);//[rad]
       XFP   =  xdet[0];
       XpFP  =  xpdet[0];
       YFP   =  ydet[0];
       YpFP  =  ypdet[0];
       Mom   =  momentum[0];
       Xpt   =  xptar[0];
       //Ypt   =  yptar[0];
       Ypt   =  yptar[0] + 13.2/180. * 3.14159;
       Zt    =  ztar[0];
       PL    =  pathl[0];

       if(Mom>1.7 && Mom<3.0 
	  && fabs(Zt)<0.25 
	  && fabs(Xpt)<0.1
	  && s2nhit == 1
	  ){
	 goodflag = true;
       }
       else goodflag = false;
       
       if(goodflag==true){
	 
	 fVal[jj][0]=(XFP-XFPm)/XFPr;//[m]
	 fVal[jj][1]=(XpFP-XpFPm)/XpFPr;//[rad]
	 fVal[jj][2]=(YFP-YFPm)/YFPr;//[m]
	 fVal[jj][3]=(YpFP-YpFPm)/YpFPr;//[rad]
	 
	 tVal[jj][0]=(Xpt-Xptm)/Xptr;//[rad]
	 tVal[jj][1]=(Ypt-Yptm)/Yptr;//[rad]
	 tVal[jj][2]=(Mom-Momm)/Momr;//[GeV/c]
	 //tVal[jj][3]=LFP;//[cm]
	 tVal[jj][3]=(PL-PLm)/PLr; // path length [m]
	 //tVal[jj][4]=1e-4;//expt // 
	 tVal[jj][4]=(Zt-Ztm)/Ztr; // z pos [m]
	 tVal[jj][5]=1e-4;//eypt
	 tVal[jj][6]=1e-5;//emom
	 tVal[jj][7]=1e-5;//eL
	 jj++;
       }
     }
   }
   nTrig=jj;
 }
 else{
   for (int i=0; i<nentries; i++){
     tree->GetEntry(i);
     if(  SPLTrig==1 
	  && XEDC[0]>-1000. 
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
	  && XEDC[15]>-1000. 
	  ){
       Double_t z1[5];
       Double_t z2[5];
       Double_t z3[3];
       Double_t z4[3];
       for (int n=0; n<5 ; n++){
	 z1[n] =  -15. + 7.5*i -0.375 ;
	 z2[n] =  -15. + 7.5*i +0.375 ;
       }
       for (int n=0; n<3 ; n++){
	 z3[n] =  -1.905 + 0.635*i + 26.2;
	 z4[n] =  0.3175 + 0.635*i + 26.2;
       }
       Double_t zxdc[nLayerX]={z1[0],z2[0],z1[1],z2[1],z1[2],z2[2],
			       z1[3],z2[3],z1[4],z2[4],
			       z3[0],z3[1],z3[2],z4[0],z4[1],z4[2]};
       Double_t zydc[nLayerY]={z1[1],z2[1],z1[3],z2[3],
			       z3[0],z3[1],z4[1],z4[2]};
       for (int k=0;k<nLayerX;k++){
	 Double_t s = sigma;
	 if (k==2 || k==3 || k==6 || k==7){
	   s = sigma/cos(30.*3.14/180.);
	 }
	 else if (k==10 || k==11 || k==14 || k==15){
	   s = sigma/cos(60.*3.14/180.);
	 }
	 XEDC[k]=gRandom->Gaus(XEDC[k],s);//[cm]
       }	
       for (int k=0;k<nLayerY;k++){
	 Double_t s = sigma;
	 if (k<=4){
	   s = sigma/sin(30.*3.14/180.);
	 }
	 else {
	   s = sigma/sin(60.*3.14/180.);
	 }
	 YEDC[k]=gRandom->Gaus(YEDC[k],s);//[cm]
       }
       
       TVectorD x;x.Use(nLayerX,XEDC);
       TVectorD y;y.Use(nLayerY,YEDC);
       TVectorD zx;zx.Use(nLayerX,zxdc);
       TVectorD zy;zy.Use(nLayerY,zydc);
       TMatrixD Ax(nLayerX,nFit);
       TMatrixDColumn(Ax,0) = 1.0;
       TMatrixDColumn(Ax,1) = zx;
       TDecompSVD svdx(Ax);
       TMatrixD Ay(nLayerY,nFit);
       TMatrixDColumn(Ay,0) = 1.0;
       TMatrixDColumn(Ay,1) = zy;
       TDecompSVD svdy(Ay);
       Bool_t ok;
       
       const TVectorD x_svd = svdx.Solve(x,ok);
       const TVectorD y_svd = svdy.Solve(y,ok);
       XFP = x_svd(0)+FPpos*x_svd(1);
       YFP = y_svd(0)+FPpos*y_svd(1);
       XpFP = atan(x_svd(1));
       YpFP = atan(y_svd(1));
       
       //X and Y in HES have - sign
       fVal[jj][0]=(XFP-XFPm)/XFPr;//[cm]
       fVal[jj][1]=(XpFP-XpFPm)/XpFPr;//[rad]
       fVal[jj][2]=(YFP-YFPm)/YFPr;//[cm]
       fVal[jj][3]=(YpFP-YpFPm)/YpFPr;//[rad]
       tVal[jj][0]=(Xpt-Xptm)/Xptr;//[rad]
       tVal[jj][1]=(Ypt-Yptm)/Yptr;//[rad]
       tVal[jj][2]=(Mom-Momm)/Momr;//[GeV/c]
       tVal[jj][3]=LEDC;//[cm]
       jj++;
     }
   }
   nTrig=jj;
 }
 
 return nTrig;
 
}

