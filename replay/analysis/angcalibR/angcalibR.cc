/*
  angcalibR.cc (RHRS)
  This is a macro to optimize angle parameters.
  
  Toshiyuki Gogami, Feb 20, 2019
*/

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TApplication.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TMarker.h"
#include"TMath.h"

#include <fstream>
#include <iostream>
using namespace std;


extern double calcf2t_zt(double* P, 
			 double xf, double xpf,
			 double yf, double ypf);
extern double calcf2t_4th_2(double*,
			    double, double, 
			    double, double,
			    double);
extern double tune(double* pa, int j, int angflag);
extern void fcn1(int &nPar, double* /*grad*/, 
		 double &fval, double* param, int /*iflag*/);
extern void fcn2(int &nPar, double* /*grad*/, 
		 double &fval, double* param, int /*iflag*/);


// =================================================== //
// ==== Offset and scaling factors for matrices ====== //
// =================================================== //
const double  XFPm=-0.7,  XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; 
const double  XFPr=1.3,   XpFPr=0.27; 
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; 
const double  PLm = 25.4, PLr=0.7; 
const double  Ztm = -0.15,Ztr=0.35; 

const int nfoil = 10;
double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125}; 
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};
//double selection_width = 0.0125; 
double selection_width = 0.005; // event selection width for z


const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SS pattern
const int ncol = 7;  // the number of column in SS pattern
const int nsshole = nrow*ncol; // the number of holes to consider 
double ssy_cent_real[nrow];
double ssx_cent_real[ncol];
double refx[nsshole];
double refy[nsshole];
//double selec_widthx = 0.60; // selection width in x (dispersive plane)
//double selec_widthy = 0.45; // selection width in y
//double selec_widthx = 0.65; // selection width in x (dispersive plane)
//double selec_widthy = 0.5; // selection width in y 
//double selec_widthx = 0.45;
//double selec_widthy = 0.35;
double selec_widthx = 0.50;
double selec_widthy = 0.35;


const int nParamT = 126;  // Number of parameters
const int nmax = 100000; // Number of events used for tuning
//const int nmax = 50000; // Number of events used for tuning
double x[nmax], y[nmax];
double xp[nmax], yp[nmax];
double z_recon[nmax]; // reconstructed z position
int foil_flag[nmax];  
int ntune_event = 0;
int holegroup[nmax];
double y_at_ss[nmax];

double l[nfoil];
//double l2[nfoil];
double dth[nfoil];
//double projectf[nfoil];
double OptPar1[nParamT];
double OptPar2[nParamT];


//const int nParamT2 = 4; 
//double parRaster[nParamT2];
//double Opt_Par[nParamT2];
//double Ras_curx[nmax];
//double Ras_cury[nmax];

const double hrs_ang = 13.2 * 3.14159 / 180.;

//void angcalib(){
int main(int argc, char** argv){
  TApplication* app = new TApplication("app", &argc, argv);
  
  int nite = 0; // The number of tuning iteration
  int flag = 1; // (tuning flag; 1=xp, others=yp)
  if(argc==1){
    cout << " nite=0: no tuning" << endl;
  }
  else if(argc==2){
    nite = atoi(argv[1]);
  }
  else if(argc>2){
    nite = atoi(argv[1]);
    flag = atoi(argv[2]);
  }
  else{
    cout << "Type ./angcalibR  #iteration  #flag(1=xp,others=yp)" << endl;
    return 100;
  }
  
  // =================================== //
  // ======== General conditions ======= //
  // =================================== //
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);


  // ================================== //
  // ======= Hole positions =========== //
  // ================================== //
  TMarker* mark[nsshole];
  int nhole = 0;
  for(int i=0 ; i<ncol ; i++){
    for(int j=0 ; j<nrow ; j++){
      ssy_cent_real[i] = -3.0*step + step*i;
      if(j%2==0) ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
      ssx_cent_real[j] = 5.0*step - step*j;
      //cout << " " << i << " " << j << ": "<<ssy_cent_real[i]  << " " << ssx_cent_real[j] << endl;

      refx[nhole] = ssx_cent_real[j];
      refy[nhole] = ssy_cent_real[i];
      
      //mark[nhole] = new TMarker(ssy_cent_real[i],ssx_cent_real[j],28);
      mark[nhole] = new TMarker(refy[nhole],refx[nhole],28);
      
      nhole++;
    }
  }

  
  // ======================================== //
  // ======= Opening a ROOT file ============ //
  // ======================================== //
  //TFile* f1 = new TFile("LHRS_SS_111721.root");
  TFile* f1 = new TFile("test_RHRS.root");
  TTree* t1 = (TTree*)f1->Get("tree");
  Double_t trig1;
  double ent = t1->GetEntries();
  //ent = 50000; // for test
  //if(nite>0) ent = 400000;
  
  cout << endl;
  cout << " The number of events read from a ROOT file: " << ent << endl;
  int evshift = 30000;
  evshift = 0;
  
  const int max = 100;
  double Xpt[max], Ypt[max];
  double Xt, Yt;
  double Zt[max];
  double XFP[max],  YFP[max];
  double XpFP[max], YpFP[max];
  t1->SetBranchAddress("Xt",   &Xt  );
  t1->SetBranchAddress("Yt",   &Yt  );
  t1->SetBranchAddress("Zt",   &Zt  ); // zt after raster correction
  t1->SetBranchAddress("Xpt",  &Xpt );
  t1->SetBranchAddress("Ypt",  &Ypt );
  t1->SetBranchAddress("XFP",  &XFP );
  t1->SetBranchAddress("YFP",  &YFP );
  t1->SetBranchAddress("XpFP", &XpFP);
  t1->SetBranchAddress("YpFP", &YpFP);
  
  TFile* fnew = new TFile("check_angcalib.root","recreate");
  //TTree* tnew = new TTree("tree","check for LHRS SS");
  //tnew->Branch("Xt",   &Xt,   "Xt/D");
  //tnew->Branch("Yt",   &Yt,   "Yt/D");
  //tnew->Branch("Zt",   &Zt,   "Zt[100]/D");
  //tnew->Branch("Xpt",  &Xpt,  "Xpt[100]/D");
  //tnew->Branch("Ypt",  &Ypt,  "Ypt[100]/D");
  //tnew->Branch("XFP",  &XFP,  "XFP[100]/D");
  //tnew->Branch("YFP",  &YFP,  "YFP[100]/D");
  //tnew->Branch("XpFP", &XpFP, "XpFP[100]/D");
  //tnew->Branch("YpFP", &YpFP, "YpFP[100]/D");
  

  char name_Mxt_R[500], name_Myt_R[500];
  //sprintf(name_Mxt_R,"../matrices/xpt_LHRS_4.dat");
  //sprintf(name_Myt_R,"../matrices/ypt_LHRS_4.dat");
  //sprintf(name_Mxt_R,"./newpar/newpar_xpt_1.dat"); 
  //sprintf(name_Myt_R,"./newpar/newpar_ypt_1.dat"); 
  sprintf(name_Mxt_R,"./newpar_xpt_2.dat");
  sprintf(name_Myt_R,"./newpar_ypt_2.dat");
  ifstream Mxt_R(name_Mxt_R);
  ifstream Myt_R(name_Myt_R);
  double Pxt_R[nParamT], Pyt_R[nParamT];
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxt_R >> par >> p >> p >> p >> p >> p;
    Pxt_R[i]  = par;
    OptPar1[i] = par;
    //cout << par << endl;
    Myt_R >> par >> p >> p >> p >> p >> p;
    Pyt_R[i]   = par;
    OptPar2[i] = par;
  }
  Mxt_R.close();
  Myt_R.close();

  istream* ifs = new ifstream("scale_offset_20190907.dat");
  double offs_xp[nfoil], offs_yp[nfoil];
  double scal_xp[nfoil], scal_yp[nfoil];
  bool offset_flag[nfoil];
  double temp;
  bool holethrough;
  for (int i=0 ; i<nfoil ; i++){
    *ifs >> temp
	 >> offs_xp[i] >> scal_xp[i]
	 >> offs_yp[i] >> scal_yp[i];
    //cout << offs_xp[i] << " " << offs_yp[i] << endl;
    if(scal_yp[i]>0.5) offset_flag[i] = true;
    else offset_flag[i] = false;
  }
  //ifs->std::close();
  
  ntune_event = 0;
  
  //TH2F* h1 = new TH2F("h1","",200,-5.0,5.0,200,-8.0,8.0); 
  TH2F* h1 = new TH2F("h1","",100,-5.0,5.0,100,-8.0,8.0); 
  h1->GetXaxis()->SetTitle("y (rad)");
  h1->GetXaxis()->CenterTitle();
  //h1->GetXaxis()->SetRangeUser(-0.2,0.2);
  h1->GetYaxis()->SetTitle("x");
  
  TH2F* h2[nfoil];
  TH2F* h2_new[nfoil];
  TH2F* h2_;
  TH2F* h3[nfoil];
  TH2F* h3_;
  char tempc[500];
  char tempc2[500];
  const double l0 = 1.03;
  //double dth[nfoil];
  int holeg_temp;
  int foilg_temp;
  bool filled=false;
   
  for(int i=0 ; i<nfoil ; i++){
    sprintf(tempc,"h2_%d",i);
    h2[i] = new TH2F(tempc, tempc,
		     h1->GetXaxis()->GetNbins(),
		     h1->GetXaxis()->GetXmin(),
		     h1->GetXaxis()->GetXmax(),
		     h1->GetYaxis()->GetNbins(),
		     h1->GetYaxis()->GetXmin(),
		     h1->GetYaxis()->GetXmax());
    //h2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
    h2[i]->GetXaxis()->SetTitle("y_SS (cm)");
    h2[i]->GetYaxis()->SetTitle("x_SS (cm)");
    l[i] = 0;
    //l[i] = sqrt(pow(l0,2.0) + pow(fcent_real[i]*100.,2.0) -2.0*l0*fcent_real[i]*100.*cos(hrs_ang));
    l[i] = (l0 - fcent_real[i]/cos(hrs_ang))*100.0;
    dth[i] = atan(l0*sin(hrs_ang)/(l0*cos(hrs_ang)-fcent_real[i])) - hrs_ang;
		  //projectf[i] = cos( dth[i] );
    //cout << i << " " << l[i] << " " << dth[i] << " " << projectf[i] << endl;
    sprintf(tempc,"h2_new_%d",i);
    h2_new[i] = (TH2F*)h2[i]->Clone(tempc);
    sprintf(tempc,"h3_%d",i);
    h3[i] = (TH2F*)h2[i]->Clone(tempc);
    h3[i]->SetMarkerColor(i+1);
  }
  h2_ = (TH2F*)h2[0]->Clone("h2_");
  h2_->SetTitle("SS pattern before event selection");
  h3_ = (TH2F*)h3[0]->Clone("h3_");
  h3_->SetTitle("SS pattern after event selection");
  
  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0; // x at FP
    y[i]    = -2222.0; // y at FP
    xp[i]   = -2222.0; // x' at FP
    yp[i]   = -2222.0; // y' at FP
    foil_flag[i] = -1; // group number of foil
    z_recon[ntune_event] = -2222.0;
    holegroup[ntune_event] = -1;
  }
  
  for (int i=0 ; i< ent ; i++){
    // ----- Initialization ------- //
    for(int j=0 ; j<max ; j++){
      XFP[j]  = -2222.0;
      XpFP[j] = -2222.0;
      YFP[j]  = -2222.0;
      YpFP[j] = -2222.0;
      Xpt[j]  = -2222.0;
      Ypt[j]  = -2222.0;
      Zt[j]   = -2222.0;
    }
    Xt = -2222.0;
    Yt = -2222.0;
    holeg_temp = -1;
    foilg_temp = -1;
    holethrough = false;
    filled = false;
    t1->GetEntry(i);
    //if(i+evshift<ent) t1->GetEntry(i+evshift); 
    //else t1->GetEntry(i-ent+evshift);

    // ------------------------------------------- //
    // ------- Event selection for tuning -------- //
    // ------------------------------------------- //
    XFP[0]  = (XFP[0]-XFPm)/XFPr;
    XpFP[0] = (XpFP[0]-XpFPm)/XpFPr;
    YFP[0]  = (YFP[0]-YFPm)/YFPr;
    YpFP[0] = (YpFP[0]-YpFPm)/YpFPr;
    Zt[0]   = (Zt[0]-Ztm)/Ztr;
    
    Xpt[0]  = calcf2t_4th_2(Pxt_R,
			    XFP[0], XpFP[0],
			    YFP[0], YpFP[0],
			    Zt[0]);
    Ypt[0] = calcf2t_4th_2(Pyt_R,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);
    
    Ypt[0]  = Ypt[0]*Yptr +Yptm;
    Xpt[0]  = Xpt[0]*Xptr +Xptm;
    Zt[0]   = Zt[0]*Ztr +Ztm;
    
    if(fabs(Xpt[0]) < 0.2 
       //&& fabs(Ypt[0]) < 0.04 ){
       && fabs(Ypt[0]) < 0.1 ){
      
      double ssx, ssy;
      for(int j=0 ; j<nfoil ; j++){
	
	if(fcent[j]-selection_width<Zt[0]
	   && Zt[0]<fcent[j]+selection_width){
	  
	  
	  ssy = l[j]*sin(atan(-Ypt[0]))/cos(dth[j]+atan(-Ypt[0])); // RHRS
	  //ssy = l[j]*sin(atan(-Ypt[0]))/cos(dth[j]-atan(-Ypt[0])); // LHRS
	  //ssx = -Xpt[0]*l[j];
	  double l2=0.0;
	  if(ssy>0){
	    l2 = sqrt(pow(l[j],2.0)+pow(ssy,2.0)+2.0*l[j]*sin(dth[j]));
	  }
	  else l2 = sqrt(pow(l[j],2.0)+pow(ssy,2.0)-2.0*l[j]*sin(dth[j]));
	  ssx = -Xpt[0]*l2;
	  
	  h2[j]->Fill(ssy,ssx);
	  
	  //if(offset_flag[j]==true){ 
	  if(offset_flag[j]==true || offset_flag[j]==false){ // (in case you don't need scale+offset for event selection)
	    
	    //ssx = (ssx + offs_xp[j])*scal_xp[j]*1.05;
	    //if(ssy>0){
	    //  ssy = (ssy + offs_yp[j])*scal_yp[j]*1.05;
	    //}
	    //else  ssy =(ssy + offs_yp[j])*scal_yp[j]*1.05;
	    
	    
	    //if (ssy>0) ssx = ssx * 1.08;

	    //if(j==8) ssy = ssy * 1.08;  // for second parameters
	    //if(j==9) ssy = ssy * 1.126; // for second parameters
	    
	    h2_new[j]->Fill(ssy,ssx);
	    h2_      ->Fill(ssy,ssx);
	    nhole=0;
	    for(int col=0 ; col<ncol ; col++){
	      for(int row=0 ; row<nrow ; row++){
		if(pow(ssx-refx[nhole],2.0)/pow(selec_widthx,2.0)
		   + pow(ssy-refy[nhole],2.0)/pow(selec_widthy,2.0)<1.0){
		  holeg_temp = nhole;
		  foilg_temp = j;
		  holethrough = true;
		}
		nhole++;
	      }
	    }
	  }
	  
	  if(ntune_event<nmax && holethrough==true
	     && filled==false ){
	    //foil_flag[ntune_event] = j;
	    foil_flag[ntune_event] = foilg_temp;
	    holegroup[ntune_event] = holeg_temp;
	    y_at_ss[ntune_event] = ssy;
	    x[ntune_event]  = XFP[0];  // scaled 
	    y[ntune_event]  = YFP[0];  // scaled
	    xp[ntune_event] = XpFP[0]; // scaled
	    yp[ntune_event] = YpFP[0]; // scaled
	    z_recon[ntune_event] = Zt[0]; // not scaled
	    ntune_event++;

	    filled=true;
	    h3[j]->Fill(ssy,ssx);
	    h3_  ->Fill(ssy,ssx);
	  }
	}
      }
    }
  }

  cout << " The number of events selected to be used for tuning: "
       << ntune_event << endl;

  
  
  double chi_sq1[nite];//, chi_sq2[nite]; // A chi square for each tuning
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;
  
  for(int i=0 ; i<nite ; i++){
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    x[i] = i+1;

    chi_sq1[i] = tune(OptPar1,i,flag);
    //chi_sq2[i] = tune(OptPar2,i,22);
    if(flag==1){
      sprintf(tempc,  "./newpar/newpar_xpt_%d.dat",i); 
      cout << " Tuning(xpt)# = " << i << ": chisq = "
	   << chi_sq1[i] << endl;
      cout << endl;
    }
    else {
      sprintf(tempc,  "./newpar/newpar_ypt_%d.dat",i); 
      cout << " Tuning(ypt)# = " << i << ": chisq = "
	   << chi_sq1[i] << endl;
      cout << endl;
    }
      
    //sprintf(tempc2, "./newpar/newpar_ypt_%d.dat",i); 

    ofstream * ofs = new ofstream(tempc);
    //ofstream * ofs2 = new ofstream(tempc2);
    int nppp = 0;
    const int nn = 4; // 4th order matrix using xf, xpf, y, ypf, and zt
    for(int i=0 ; i<nn+1 ; i++){
      for(int e=0 ; e<nn+1 ; e++){
	for(int d=0 ; d<nn+1 ; d++){
	  for(int c=0 ; c<nn+1 ; c++){
	    for(int b=0 ; b<nn+1 ; b++){
	      for(int a=0 ; a<nn+1 ; a++){  
		if(a+b+c+d+e==i){
		  if(flag==1){
		    *ofs << OptPar1[nppp] 
			 << " " << a 
			 << " " << b
			 << " " << c
			 << " " << d
			 << " " << e << endl;
		  }
		  else {
		    *ofs << OptPar2[nppp] 
			 << " " << a 
			 << " " << b
			 << " " << c
			 << " " << d
			 << " " << e << endl;
		  }
		  nppp++;
		  //	      cout << Plen_opt[nppp] 
		  //		   << " " << a 
		  //		   << " " << b
		  //		   << " " << c
		  //		   << " " << d << endl;
		}
	      }
	    }
	  }
	}
      }
    }
    ofs->close();
    ofs->clear();
    //ofs2->close();
    //ofs2->clear();

  }
  

  // =================================== //
  // ======== Draw histograms ========== //
  // =================================== //
  TCanvas* c1 = new TCanvas("c1","c1");
  TCanvas* c2 = new TCanvas("c2","c2");
  //h2[0]->Draw();
  c1->Divide(2,5);
  c2->Divide(2,5);
  //h1->Draw();
  for(int i=0 ; i<nfoil ; i++){
    c1->cd(i+1);h2[i]->Draw("col");
    c2->cd(i+1);h2_new[i]->Draw("col");
  }

  TCanvas* c3 = new TCanvas("c3","c3");
  h2_->Draw();
  for(int i=0 ; i<nsshole ; i++){
    mark[i]->SetMarkerColor(2);
    mark[i]->Draw("same");
  }

  TCanvas* c4 = new TCanvas("c4","c4");
  c4->cd();
  h3_->Draw();
  
  for(int i=0 ; i<nfoil ; i++){
    h3[i]->Draw("same");
  }
  for(int i=0 ; i<nsshole ; i++){
    //mark[i]->SetMarkerColor(3);
    mark[i]->Draw("same");
  }

  //fnew->Close(); 

  TObjArray h(1);
  for(int i=0 ; i<nfoil ; i++){
    h.Add(h2[i]);
    h.Add(h2_new[i]);
    h.Add(h3[i]);
  }
  h.Add(h2_);
  h.Add(h3_);
  h.Add(c1);
  h.Add(c2);
  
  if(nite>0){
    TGraph * gr1 = new TGraph(nite,x,chi_sq1);
    gr1->SetName("gr1");
    //TGraph * gr2 = new TGraph(nite,x,chi_sq2);
    //gr2->SetName("gr2");
    TCanvas * c3 = new TCanvas("c3","c3",600,600);
    //c3->Divide(1,2);
    c3->cd(1);gr1->Draw("*la");
    //c3->cd(2);gr2->Draw("*la");
    h.Add(gr1);
    //h.Add(gr2);
    h.Add(c3);
  }
  h.Write();

  app->Run();
  
  return 0;
  
} 

//////////////////////////////////////////////////
double calcf2t_zt(double* P, double xf, double xpf, 
		  double yf, double ypf)
//////////////////////////////////////////////////
{
  // -----3rd order ----- 
  // This is the third order claculation byb  using 35 parameter
  const int nMatT=3;  // These are for the RHRS use same for the LHRS
  const int nXf=3;
  const int nXpf=3;
  const int nYf=3;
  const int nYpf=3;
  
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

//////////////////////////////////////////////////
double calcf2t_4th_2(double* P, double xf, double xpf, 
		     double yf, double ypf, double zt)
//////////////////////////////////////////////////
{
  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
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
  }
  
  return Y; 
  
}

// #############################################################
double tune(double* pa, int j, int angflag) 
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = nParamT;
  //cout << allparam << endl;
  TMinuit* minuit = new TMinuit(allparam);
  if(angflag==1){
    minuit->SetFCN(fcn1);
  }
  else minuit->SetFCN(fcn2);


  double start[allparam];
  double step[allparam];
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=2; // The number of order is reduced for test (4-->2)
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  start[npar] = pa[npar];
		  //step[npar] = 1.0e-3;
		  step[npar] = pa[npar] * 5.0e-2;  
		}
		else{
		  start[npar] = 0.0;
		  step[npar] = 0.0;
		}
		npar++;
	      }
	    }
	  }
	}
      }    
    }
  }
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  
  double LLim[allparam];// Lower limit for all of the parameter
  double ULim[allparam];// Upper limit for all of the parameter
  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    //start[i] = pa[i]; 
    //step[i] = 1.0e-3;
    
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
    LLim[i] = pa[i] - 5.0; // temp
    ULim[i] = pa[i] + 5.0; // temp
    
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }
  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0; // original
  //arglist[0] = 1.0; // test
  //arglist[0] = 0.0;   // test
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    //minuit -> GetParameter(i,par[i],e);
    if(angflag==1){
      minuit -> GetParameter(i,OptPar1[i],er);
    }
    else minuit -> GetParameter(i,OptPar2[i],er);
  }
  
  return chi2;
}

// #############################################################
void fcn1(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  
  const double sigma = 1.0;
  double ztR      = 0.0;
  double refpos   = 0.0;
  double residual = 0.0;
  double ang      = 0.0;
  double sspos    = 0.0;
  double total_chi2 = 0.0;
  
  double nev[nfoil][nsshole];
  double chi2[nfoil][nsshole];
  double w[nfoil][nsshole];
  double ssy;

  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      nev[i][j]  = 0.0;
      chi2[i][j] = 0.0;
      w[i][j]    = 1.0;
    }
  }
  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    ang    = 0.0;
    sspos  = 0.0;
    refpos = 0.0;  refpos = refx[holegroup[i]];
    ztR    = 0.0;  ztR    = z_recon[i];
    ssy = y_at_ss[i];
    
    //if(foil_flag[i]==i) nev[i]++;

    ztR= (ztR-Ztm)/Ztr; // only zt was not scaled, so apply scaling here
    ang = calcf2t_4th_2(param,
			x[i], xp[i],
			y[i], yp[i],
			ztR);
    ang = -1.0 * (ang * Xptr +Xptm);
    //sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter
    //sspos = ang * l[foil_flag[i]]; // in centimeter
    double l2=0.0;
    if(ssy>0){
      l2 = sqrt(pow(l[foil_flag[i]],2.0)+pow(ssy,2.0)
		+2.0*l[foil_flag[i]]*sin(dth[foil_flag[i]]));
    }
    else {
      l2 = sqrt(pow(l[foil_flag[i]],2.0)+pow(ssy,2.0)
		-2.0*l[foil_flag[i]]*sin(dth[foil_flag[i]]));
    }
    sspos = ang * l2; // in centimeter
    
    // ------------------- //
    // --- Residual ------ //
    // ------------------- //
    residual = sspos - refpos;
    //cout << i << ": " << sspos << "-(" << refpos << ")=" << residual << endl;
    
    chi2[foil_flag[i]][holegroup[i]]
      = chi2[foil_flag[i]][holegroup[i]] + pow(residual,2.0);
    
    nev[foil_flag[i]][holegroup[i]]++;
  }

  double totevent = 0;
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      
      //if(nev[i][j]>0){
      if(nev[i][j]>5){ 
      //if(nev[i][j]>50){ // using only holes with more than 50 events
	//chi2[i][j] = chi2[i][j]/(nev[i][j]-126.0)/pow(sigma,2.0);
	chi2[i][j] = chi2[i][j]/pow(sigma,2.0);
      }
      else chi2[i][j] = 0.0;
      
      total_chi2 = total_chi2 + chi2[i][j]*w[i][j];
      totevent++;
    }
  }
  
  //fval = total_chi2/(double)nfoil/(double)nsshole;
  //if(totevent>125.){
  //  fval = total_chi2/(totevent-126.0);
  //}
  //else fval = total_chi2;

  fval = total_chi2;
}


// #############################################################
void fcn2(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  
  const double sigma = 1.0;
  double ztR      = 0.0;
  double refpos   = 0.0;
  double residual = 0.0;
  double ang      = 0.0;
  double sspos    = 0.0;
  double total_chi2 = 0.0;
  
  double nev[nfoil][nsshole];
  double chi2[nfoil][nsshole];
  double w[nfoil][nsshole];

  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      nev[i][j]  = 0.0;
      chi2[i][j] = 0.0;
      w[i][j]    = 1.0;
    }
  }
  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    ang    = 0.0;
    sspos  = 0.0;
    refpos = 0.0;  refpos = refy[holegroup[i]];
    ztR    = 0.0;  ztR    = z_recon[i];
    
    //if(foil_flag[i]==i) nev[i]++;

    ztR= (ztR-Ztm)/Ztr; // only zt was not scaled, so apply scaling here
    ang = calcf2t_4th_2(param,
			x[i], xp[i],
			y[i], yp[i],
			ztR);
    ang = -1.0 * (ang * Yptr +Yptm);
    //sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter
    sspos = l[foil_flag[i]] *sin(atan(ang))/cos(dth[foil_flag[i]]+atan(ang)); // RHRS
    //sspos = l[foil_flag[i]] *sin(atan(ang))/cos(dth[foil_flag[i]]-atan(ang)); // LHRS
    
    // ------------------- //
    // --- Residual ------ //
    // ------------------- //
    residual = sspos - refpos;
    //cout << i << ": " << sspos << "-(" << refpos << ")=" << residual << endl;
    
    chi2[foil_flag[i]][holegroup[i]]
      = chi2[foil_flag[i]][holegroup[i]] + pow(residual,2.0);
    
    nev[foil_flag[i]][holegroup[i]]++;
  }

  double totevent = 0;
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      
      //if(nev[i][j]>0){
      //if(nev[i][j]>10){ 
      if(nev[i][j]>5){ // using only holes with more than 50 events
	//chi2[i][j] = chi2[i][j]/(nev[i][j]-126.0)/pow(sigma,2.0);
	chi2[i][j] = chi2[i][j]/pow(sigma,2.0);
      }
      else chi2[i][j] = 0.0;
      
      total_chi2 = total_chi2 + chi2[i][j]*w[i][j];
      totevent++;
    }
  }
  
  //fval = total_chi2/(double)nfoil/(double)nsshole;
  //if(totevent>125.){
  //  fval = total_chi2/(totevent-126.0);
  //}
  //else fval = total_chi2;
  
  fval = total_chi2;
}



