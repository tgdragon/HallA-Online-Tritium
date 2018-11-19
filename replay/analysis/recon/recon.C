#include "recon.h"
#include <fstream>
#include <iostream>
#include <TMatrixD.h>
#include <TFile.h>
#include <TTree.h>
using namespace std;


///////////////////////////
int main()
///////////////////////////
{
 
  double Xt=0.0,Yt=0.0, Zt=0.0;
  double XFP=0.,XpFP=0.,YFP=0.,YpFP=0.;
  double Mom=0.,Xpt=0.,Ypt=0.,Len=0.;
  double MomR=0.,Xpt_R=0.,Ypt_R=0.,LenR=0.,Zt_R=0.0;
  double XSS=0.,XpSS=0.,YSS=0.,YpSS=0.;
  double XSSR=0.,XpSSR=0.,YSSR=0.,YpSSR=0.;
  double p=0.,xp=0.,yp=0.;
  double PL=0.0;
  double momres=0.;
  double Pmom[nParamT_4];
  double Pxpt[nParamT_4];
  double Pypt[nParamT_4];
  double Plen[nParamT];
  double Pzt[nParamT];
  double Pxss[nParamS];
  double Pxpss[nParamS];
  double Pyss[nParamS];
  double Pypss[nParamS];
  int XHoleID=100;
  int YHoleID=100;
  int nTrig=0;
  char name_Mmom[100];
  char name_Mxpt[100];
  char name_Mypt[100];
  char name_Mlen[100];
  char name_Mzt[100];
  char name_Mxss[100];
  char name_Mxpss[100];
  char name_Myss[100];
  char name_Mypss[100];
  char name_ifile[100];
  char name_ofile[100];
  //const double  XFPm=-10.0, XpFPm=-0.093; 
  //const double  YFPm=-135.0, YpFPm=-0.93; 
  //const double  Xptm=-0.245, Yptm=-0.082, Momm=1.67;
  //const double  XFPr=20.0, XpFPr=0.185; 
  //const double  YFPr=245.0, YpFPr=0.31; 
  //const double  Xptr=0.054, Yptr=0.163, Momr=0.25;
  const double  XFPm=-0.7, XpFPm=-0.15; 
  const double  YFPm=-0.05, YpFPm=-0.18; 
  const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74;
  const double  XFPr=1.3, XpFPr=0.27; 
  const double  YFPr=0.1, YpFPr=0.10; 
  const double  Xptr=0.15, Yptr=0.08, Momr=0.18;
  const double  PLm = 25.4, PLr=0.7;
  const double  Ztm = -0.15, Ztr=0.35;
  int LR_flag = 2; // 2: Right, Others: Left
  
  ifstream input("parRecon.dat");
  input >> name_Mmom;
  input >> name_Mxpt;
  input >> name_Mypt;
  input >> name_Mlen;
  input >> name_Mzt;
  input >> name_Mxss;
  input >> name_Mxpss;
  input >> name_Myss;
  input >> name_Mypss;
  input >> name_ifile;
  input >> name_ofile;
  input.close();
  
  ifstream Mmom(name_Mmom);
  for (int i=0;i<nParamT_4;i++){
    double par=0.;
    int p=0;
    Mmom >> par >> p >> p >> p >> p; 
    Pmom[i]=par;
  }
  Mmom.close();
  
  ifstream Mxpt(name_Mxpt);
  for (int i=0;i<nParamT_4;i++){
    double par=0.;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p; 
    Pxpt[i]=par;
  }
  Mxpt.close();
  
  ifstream Mypt(name_Mypt);
  for (int i=0;i<nParamT_4;i++){
    double par=0.;
    int p=0;
    Mypt >> par >> p >> p >> p >> p; 
    Pypt[i]=par;
  }
  Mypt.close();
  
  ifstream Mlen(name_Mlen);
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mlen >> par >> p >> p >> p >> p; 
    Plen[i]=par;
  }
  Mlen.close();
  
  ifstream Mzt(name_Mzt);
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p; 
    Pzt[i]=par;
  }
  Mlen.close();
  
  ifstream Mxss(name_Mxss);
  for (int i=0;i<nParamS;i++){
    double par=0.;
    int p=0;
    Mxss >> par >> p >> p >> p ; 
    Pxss[i]=par;
  }
  Mxss.close();
  
  ifstream Mxpss(name_Mxpss);
  for (int i=0;i<nParamS;i++){
    double par=0.;
    int p=0;
    Mxpss >> par >> p >> p >> p; 
    Pxpss[i]=par;
  }
  Mxpss.close();
  
  ifstream Myss(name_Myss);
   for (int i=0;i<nParamS;i++){
     double par=0.;
     int p=0;
     Myss >> par >> p >> p >> p ; 
     Pyss[i]=par;
   }
   Myss.close();
   
   ifstream Mypss(name_Mypss);
   for (int i=0;i<nParamS;i++){
     double par=0.;
     int p=0;
     Mypss >> par >> p >> p >> p; 
     Pypss[i]=par;
   }
   Mypss.close();
   
   
   TFile *ifile = new TFile(name_ifile);
   TTree *itree = (TTree*)ifile->Get("tree");
   TFile *ofile = new TFile(name_ofile,"recreate");
   TTree *otree = new TTree("tree","Tree of reconstructed event");
   int nnn = (int)itree->GetEntries();
   //const int nentries =(int)itree->GetEntries();
   nnn = 20000;
   const int nentries = nnn;
   double fVal[nentries][4];
   double tVal[nentries][7];
   double sVal[nentries][4];
   //cout << "open root file" << endl;
   
   //otree->Branch("Eloss",&Eloss,"Eloss/D");
   otree->Branch("KXt",&Xt,"KXt/D");
   otree->Branch("KYt",&Yt,"KYt/D");
   otree->Branch("KZt",&Zt,"KZt/D");
   otree->Branch("KXFP",&XFP,"KXFP/D");
   otree->Branch("KXpFP",&XpFP,"KXpFP/D");
   otree->Branch("KYFP",&YFP,"KYFP/D");
   otree->Branch("KYpFP",&YpFP,"KYpFP/D");
   otree->Branch("KXpt",&Xpt,"KXpt/D");
   otree->Branch("KYpt",&Ypt,"KYpt/D");
   otree->Branch("KLen",&Len,"KLen/D");
   otree->Branch("KMom",&Mom,"KMom/D");
   otree->Branch("KXptR",&Xpt_R,"KptR/D");
   otree->Branch("KYptR",&Ypt_R,"KYptR/D");
   otree->Branch("KLenR",&LenR,"KLenR/D");
   otree->Branch("KMomR",&MomR,"KMomR/D");
   otree->Branch("KZtR",&Zt_R,"KZtR/D");
   //otree->Branch("KXSS",&XSS,"KXSS/D");
   //otree->Branch("KYSS",&YSS,"KYSS/D");
   //otree->Branch("KXpSS",&XpSS,"KXpSS/D");
   //otree->Branch("KYpSS",&YpSS,"KYpSS/D");
   //otree->Branch("KXSSR",&XSSR,"KXSSR/D");
   //otree->Branch("KYSSR",&YSSR,"KYSSR/D");
   //otree->Branch("KXpSSR",&XpSSR,"KXpSSR/D");
   //otree->Branch("KYpSSR",&YpSSR,"KYpSSR/D");
   //otree->Branch("KXHoleID",&XHoleID,"KXHoleID/I");
   //otree->Branch("KYHoleID",&YHoleID,"KYHoleID/I");
   
   cout << "reading file..."<<endl;
   //nTrig=readvalue(ifile,0,fVal,tVal,sVal);
   nTrig=readvalue(ifile,nentries,fVal,tVal,sVal);
   cout <<"read "<< nTrig <<" events"<< endl;
   double YFP_cor, YpFP_cor;
   for (int i=0;i<nTrig;i++){
     int zzz=(int)((double)i/(double)nTrig*100.);
     if (i%1000==0){
       cout << zzz << "% finished"<<endl;
     }

     Zt  = tVal[i][6];
     //cout << Zt << endl;
     YFP = fVal[i][2];
     YpFP= fVal[i][3];
     if(LR_flag == 2){
       YFP_cor  = YFP  -( -0.11 * Zt);
       YpFP_cor = YpFP -(  0.10 * Zt);
     }
     else{
       YFP_cor  = YFP  -( 0.12 * Zt);
       YpFP_cor = YpFP -(-0.09 * Zt);
     }
     XFP=(fVal[i][0]-XFPm)/XFPr;
     XpFP=(fVal[i][1]-XpFPm)/XpFPr;
     YFP=(fVal[i][2]-YFPm)/YFPr;
     YpFP=(fVal[i][3]-YpFPm)/YpFPr;
     YFP_cor=(YFP_cor-YFPm)/YFPr;
     YpFP_cor=(YpFP_cor-YpFPm)/YpFPr;

     
     Xpt=(tVal[i][0]-Xptm)/Xptr;
     Ypt=(tVal[i][1]-Yptm)/Yptr;
     Mom=(tVal[i][2]-Momm)/Momr;
     Len=(tVal[i][3]-PLm)/PLr;
     Xt=tVal[i][4];
     Yt=tVal[i][5];
     //Zt=(tVal[i][6]-Ztm)/Ztr;
     
     MomR  = calcf2t_4th(Pmom,XFP,XpFP,YFP_cor,YpFP_cor);
     Xpt_R = calcf2t_4th(Pxpt,XFP,XpFP,YFP_cor,YpFP_cor);
     Ypt_R = calcf2t_4th(Pypt,XFP,XpFP,YFP_cor,YpFP_cor);
     LenR  = calcf2t_3rd(Plen,XFP,XpFP,YFP,YpFP);
     Zt_R  = calcf2t_3rd(Pzt,XFP,XpFP,YFP,YpFP);
     XSSR  = calct2s(Pxss,Xpt,Ypt,Mom);
     XpSSR = calct2s(Pxpss,Xpt,Ypt,Mom);
     YSSR  = calct2s(Pyss,Xpt,Ypt,Mom);
     YpSSR = calct2s(Pypss,Xpt,Ypt,Mom);
     //XSS=-sVal[i][0];
     //XpSS=-sVal[i][1];
     //YSS=-sVal[i][2];
     //YpSS=-sVal[i][3];
     
     
     //XFP = XFP*XFPr+XFPm;
     //XpFP = XpFP*XpFPr+XpFPm;
     //YFP = YFP*YFPr+YFPm;
     //YpFP = YpFP*YpFPr+YpFPm;
     XFP = XFP*XFPr+XFPm;
     XpFP = XpFP*XpFPr+XpFPm;
     YFP = YFP*YFPr+YFPm;
     YpFP = YpFP*YpFPr+YpFPm;
     Mom = Mom*Momr+Momm;
     Xpt = Xpt*Xptr+Xptm;
     Ypt = Ypt*Yptr+Yptm;
     //Zt = Zt*Ztr+Ztm;
     Len  = Len*PLr+PLm;
     
     MomR  = MomR*Momr+Momm;
     Xpt_R = Xpt_R*Xptr+Xptm;
     Ypt_R = Ypt_R*Yptr+Yptm;
     Zt_R  = Zt_R*Ztr+Ztm;
     LenR  = LenR*PLr+PLm;
     
     //cout << Zt << " " << Zt_R << endl;
     
     
     otree->Fill();
   }
   otree->Write();
   ofile->Close();
   
   return 0;
}

