/*
  pathl_opt.cc
  
  "pion time at the target + path length correction
  --> Path length optimization"
  
  Toshiyuki Gogami, Nov 7, 2018
*/

#include "pathl_opt.hh"
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <fstream>
#include <TMatrixD.h>
#include <iostream>
#include <math.h>
#include <TMinuit.h>
#include <TObjArray.h>
using namespace std;

//void cointime(){
int main(int argc, char** argv){
  // ===== General conditions ===== //
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  
  // ===== TApplication ======= //
  TApplication app("app",&argc,argv); 
  
  
  // ===== Opening a ROOT file ======= //
  //TFile* f1 = new TFile("coin_111170_111190.root"); // For the first H2 data (Oct31-Nov2, 2018)
  TFile* f1 = new TFile("coin_111480_111482.root"); // For the second H2 data (Nov10-, 2018)
  TTree* t1 = (TTree*)f1->Get("tree");
  
  // ===== 
  const double mpi = 0.13957;  // GeV/c2
  const double me  = 0.000511; // GeV/c2
  const int n   = 16;
  const int nf1 = 64;
  double ra[n],   la[n];
  double ra_R[n], la_R[n];
  double f1tdc[nf1], f1tdc_R[nf1];
  double ent = t1->GetEntries();
  double rf, rf_R;
  double timeL, timeR;
  double timeL_R, timeR_R;
  double nhit, nhit_R;
  const int max = 100;
  double hitseg[max];
  double hitseg_R[max];
  double pathl[max];
  double pathl_R[max];
  double pathl_s2[max];
  double pathl_s2_R[max];
  double mom[max];
  double mom_R[max];
  double xdet[max], ydet[max];
  double xpdet[max], ypdet[max];
  double xdet_L[max], ydet_L[max];
  double xpdet_L[max], ypdet_L[max];
  int seg, seg_R;
  double meantime;
  double meantime_R;
  double time_target;
  double time_rf, time_rf_R;
  double beta;
  double beta_R;
  double tref;
  double tref_R;
  double a1, a2;
  double msq;
  //const double toffset_R = -364.6-150.;
  const double toffset_R = -364.6-150.-15-0.788; // for second H2
  //double ra, la;
  //const double ch2time = 56.0e-12 * 1.0e+9; // 56 (ps/ch)
  const double ch2time = 58.0e-12 * 1.0e+9; // 58 (ps/ch) for second H2 data
  bool generalcut = false;
  bool adccut = false;
  bool pi_selection = false;
  //ent = 5000000;
  //ent = 1000000;
  //ent = 600000;
  ent = 300000;
  //ent = 150000;
  //ent = 50000;
  
  t1->SetBranchAddress("L.s2.ra_c",&ra);
  t1->SetBranchAddress("R.s2.ra_c",&ra_R);
  t1->SetBranchAddress("L.s2.la_c",&la);
  t1->SetBranchAddress("R.s2.la_c",&la_R);
  t1->SetBranchAddress("LTDC.F1FirstHit",&f1tdc);
  t1->SetBranchAddress("RTDC.F1FirstHit",&f1tdc_R);
  t1->SetBranchAddress("L.s2.nthit",&nhit);
  t1->SetBranchAddress("R.s2.nthit",&nhit_R);
  t1->SetBranchAddress("L.s2.t_pads",&hitseg);
  t1->SetBranchAddress("R.s2.t_pads",&hitseg_R);
  t1->SetBranchAddress("L.tr.p",&mom);
  t1->SetBranchAddress("R.tr.p",&mom_R);
  t1->SetBranchAddress("L.tr.pathl",&pathl);
  t1->SetBranchAddress("R.tr.pathl",&pathl_R);
  t1->SetBranchAddress("L.s2.trpath",&pathl_s2);
  t1->SetBranchAddress("R.s2.trpath",&pathl_s2_R);
  t1->SetBranchAddress("R.a1.asum_c",&a1);
  t1->SetBranchAddress("R.a2.asum_c",&a2);
  
  double XFP, XpFP;
  double YFP, YpFP;
  double XFP_L, XpFP_L;
  double YFP_L, YpFP_L;
  t1->SetBranchAddress("R.tr.x" ,&xdet); // dispersive plane
  t1->SetBranchAddress("R.tr.th",&xpdet);
  t1->SetBranchAddress("R.tr.y", &ydet); // non-dispersive plane
  t1->SetBranchAddress("R.tr.ph",&ypdet);
  t1->SetBranchAddress("L.tr.x" ,&xdet_L); // dispersive plane
  t1->SetBranchAddress("L.tr.th",&xpdet_L);
  t1->SetBranchAddress("L.tr.y", &ydet_L); // non-dispersive plane
  t1->SetBranchAddress("L.tr.ph",&ypdet_L);
  //const double  XFPm=-0.7, XpFPm=-0.15; 
  //const double  YFPm=-0.05, YpFPm=-0.18; 
  //const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74;
  //const double  XFPr=1.3, XpFPr=0.27; 
  //const double  YFPr=0.1, YpFPr=0.10; 
  //const double  Xptr=0.15, Yptr=0.08, Momr=0.18;
  //const double  PLm = 25.4, PLr=0.7;
  //const double  Ztm = -0.15, Ztr=0.35;
  double LenR, LenL;
  
  // ------- Matrix for path length reconstruction --------- //
  //double tval0[tunemax]; // x 
  //double tval1[tunemax]; // xp 
  //double tval2[tunemax]; // y 
  //double tval3[tunemax]; // yp
  //double betabeta[tunemax]; // beta
  //int ntune;
  char name_Mlen[100];
  sprintf(name_Mlen,"matrices/len_RHRS_1.dat"); // original
  //sprintf(name_Mlen,"outmat.dat");
  //sprintf(name_Mlen,"outmat2.dat");
  //sprintf(name_Mlen,"outmat3.dat");
  //sprintf(name_Mlen,"outmat4.dat");
  //sprintf(name_Mlen,"outmat5.dat");
  //sprintf(name_Mlen,"len_newmat_opt.dat");
  ifstream Mlen(name_Mlen);
  double Plen[nParamT];
  //double Plen_opt[nParamT];
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mlen >> par >> p >> p >> p >> p; 
    Plen[i]=par;
  }
  Mlen.close();
  
  //ifstream* s2_R_data = new ifstream("data/s2_t0_R.dat"); // H2_1
  ifstream* s2_R_data = new ifstream("data/s2_t0_R_2.dat"); // H2_2
  double s2_tzero[n];
  for(int i=0 ; i<n ; i++){
    *s2_R_data >> s2_tzero[i];
  }
  s2_R_data->close();

  //ifstream* s2_L_data = new ifstream("data/s2_t0_L.dat"); // H2_1
  ifstream* s2_L_data = new ifstream("data/s2_t0_L.dat"); // H2_2
  double s2_tzero_L[n];
  for(int i=0 ; i<n ; i++){
    *s2_L_data >> s2_tzero_L[i];
  }
  s2_L_data->close();
  

  ifstream* rtime_ycor_data = new ifstream("data/rtime_ycor.dat");
  for(int i=0 ; i<npar_rtime_ycor ; i++){
    *rtime_ycor_data >> par_rtime_ycor[i];
    //cout << par_rtime_ycor[i] << endl;
  }
  rtime_ycor_data->close();
  
  //char fname[500];
  //ifstream* rtime_ycor_data_L = new ifstream("data/rtime_ycor_L.dat"); // H2_1
  ifstream* rtime_ycor_data_L = new ifstream("data/rtime_ycor_L_2.dat"); // H2_2
  for(int i=0 ; i<npar_rtime_ycor_L ; i++){
    *rtime_ycor_data_L >> par_rtime_ycor_L[i];
    //cout << par_rtime_ycor_L[i] << endl;
  }
  rtime_ycor_data_L->close();
  
  ifstream* pathl_cor_data_L = new ifstream("data/pathl_L.dat");
  for(int i=0 ; i<npar_pathl_L_cor ; i++){
    *pathl_cor_data_L >> par_pathl_L_cor[i];
    //cout << par_rtime_ycor_L[i] << endl;
  }
  pathl_cor_data_L->close();

  
  
  // ======================================== //
  // ======== Histogram definitions ========= //
  // ======================================== //
  //TH1F* h1 = new TH1F("h1","",100,-10,10);
  //TH1F* h1 = new TH1F("h1","",100,0,4000);
  //TH1F* h1 = new TH1F("h1","",500,-50,50);
  //TH1F* h1 = new TH1F("h1","",200,370,420);
  TH1F* h1 = new TH1F("h1","",500,-500,500.);
  TH1F* h2 = new TH1F("h2","",300,-300,300.);
  //TH1F* h3 = new TH1F("h3","",100,-20.0,20.);
  TH1F* h3 = new TH1F("h3","",200,-20.0,20.);
  //TH1F* h3 = new TH1F("h3","",200,-200.0,200.);
  //TH2F* h4 = new TH2F("h4","",100,22.0,22.9,100,-10,10.);
  TH2F* h4 = new TH2F("h4","",100,21.9,22.8,100,-10,10.);
  h4->GetXaxis()->SetTitle("Path length to S2 (m)");
  h4->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  TH1F* h5 = new TH1F("h5","",200,-10.0,10.0);
  TH1F* h6 = (TH1F*)h3->Clone("h6");
  TH2F* h7 = new TH2F("h7","",100,-1.0,1.0,100,-10,10.);
  h7->GetXaxis()->SetTitle("x (m)");
  h7->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  TH2F* h8 = new TH2F("h8","",100,-0.15,0.15,100,-10,10.);
  h8->GetXaxis()->SetTitle("x' (rad)");
  h8->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  TH2F* h9 = new TH2F("h9","",100,-0.1,0.1,100,-10,10.);
  h9->GetXaxis()->SetTitle("y (m)");
  h9->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  TH2F* h10 = new TH2F("h10","",100,-0.06,0.05,100,-10,10.);
  h10->GetXaxis()->SetTitle("y' (rad)");
  h10->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  TH1F* h11 = new TH1F("h11","",300,-3000.,3000.);
  h11->GetXaxis()->SetTitle("RF (rad)");
  //TH1F* h12 = new TH1F("h12","",200,20.0,100.);
  TH1F* h12 = new TH1F("h12","",400,0.-20,20.);
  //TH2F* h13 = new TH2F("h13","",100,21.9,23.8,100,-10,10.);
  TH2F* h13 = new TH2F("h13","",100,21.9,23.5,100,-10,10.);
  h13->GetXaxis()->SetTitle("Path length to LHRS-S2 (m)");
  h13->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  
  TH2F* h14 = new TH2F("h14","",100,-1.0,1.0,100,-10,10.);
  h14->GetXaxis()->SetTitle("x (m)");
  h14->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  TH2F* h15 = new TH2F("h15","",100,-0.15,0.15,100,-10,10.);
  h15->GetXaxis()->SetTitle("x' (rad)");
  h15->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  TH2F* h16 = new TH2F("h16","",100,-0.1,0.1,100,-10,10.);
  h16->GetXaxis()->SetTitle("y (m)");
  h16->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  TH2F* h17 = new TH2F("h17","",100,-0.06,0.05,100,-10,10.);
  h17->GetXaxis()->SetTitle("y' (rad)");
  h17->GetYaxis()->SetTitle("#pi^{+} time at target (ns)");
  
  
  TH1F* hs2[n];
  char s2histname[500];
  char s2histtitle[500];
  for(int i=0 ; i<n ; i++){
    sprintf(s2histname,"hs2_%d",i);
    hs2[i] = new TH1F(s2histname,s2histname,100,-10,10);
    sprintf(s2histtitle,"time at target (RHRS-S2-%d) (ns) ",i);
    hs2[i]->GetXaxis()->SetTitle(s2histtitle);
  }
  
  TH1F* hs2_L[n];
  for(int i=0 ; i<n ; i++){
    sprintf(s2histname,"hs2_%d_L",i);
    hs2_L[i] = new TH1F(s2histname,s2histname,100,-10,10);
    sprintf(s2histtitle,"time at target (LHRS-S2-%d) (ns) ",i);
    hs2_L[i]->GetXaxis()->SetTitle(s2histtitle);
  }
  
  ntune = 0;
  double valval[max];
  
  // ========================================== //
  // ======= Main loop for event filling ====== //
  // ========================================== //
  for(int i=0 ; i<ent ; i++){
    // ----- Initialization ------ //
    for(int j=0 ; j<n ; j++){
      ra[j] = -2222.0;
      la[j] = -2222.0;
      ra_R[j] = -2222.0;
      la_R[j] = -2222.0;
    }
    for(int j=0; j<nf1 ; j++){
      f1tdc[j] = -2222.0;
      f1tdc_R[j] = -2222.0;
    }
    for(int j=0 ; j<max ; j++){
      hitseg[j]=-1;
      mom[j]   = -2222.0;
      pathl[j] = -2222.0;
      pathl_s2[j] = -2222.0;
      valval[j] = -2222.0;
      
      hitseg_R[j]=-1;
      mom_R[j]   = -2222.0;
      pathl_R[j] = -2222.0;
      pathl_s2_R[j] = -2222.0;
    }
    nhit     = 0.0;
    nhit_R   = 0.0;
    //hitseg   = 0;
    timeL    = -2222.0;
    timeR    = -2222.0;
    timeL_R  = -2222.0;
    timeR_R  = -2222.0;
    meantime = -2222.0;
    meantime_R = -2222.0;
    time_rf  = -2222.0;
    time_target = -2222.0;
    tref     = -2222.0;
    tref_R   = -2222.0;
    a1 = -2222.0;
    a2 = -2222.0;
    msq = -2222.0;
    
    generalcut = false;
    adccut = false;
    pi_selection = false;
    
    
    
    // ------ GetEntry ------- //
    t1->GetEntry(i);
    // ----------------------- //
    
    if(nhit == 1  && nhit_R == 1
       && (hitseg[0]>0 && hitseg[0]<17) 
       && (hitseg_R[0]>0 && hitseg_R[0]<17)
       && mom_R[0]>1.3 && mom[0]>1.3
       ) {
      generalcut = true;
    }
    else generalcut = false;
    
    if(generalcut==true){
      //seg = hitseg[0] - 1;
      seg   = hitseg[0];
      seg_R = hitseg_R[0];
      
      if(ra[seg]>2.0 && la[seg]>2.0 &&
	 ra_R[seg_R]>2.0 && la_R[seg_R]>2.0 ){
	adccut = true;
      }
      else adccut = false;
      
      if(a1>0.5 && a2>15.0) pi_selection = true;
      else pi_selection = false;
      
      //if(adccut==true || adccut==false){
      //if(adccut==true){
      //if(adccut==true && pi_selection==true){
      if(adccut==true){
      //if(adccut==true && seg==11 &&
	  //){
	//cout << ra[seg] << " " << la[seg] << endl;
	tref    = f1tdc[40] * ch2time;       // L1A  (LHRS)
	rf      = f1tdc[47] * ch2time;       // RF   (LHRS)
	timeL   = f1tdc[seg] * ch2time;      // Left PMT  (LHRS)
	timeR   = f1tdc[seg+48] * ch2time;   // Right PMT (LHRS)
	rf_R    = f1tdc[15] * ch2time;       // RF   (RHRS)
	tref_R  = f1tdc_R[9] * ch2time;      // L1A  (RHRS)
	timeL_R = f1tdc_R[seg_R+16] * ch2time; // Left PMT  (RHRS)
	timeR_R = f1tdc_R[seg_R+48] * ch2time; // Right PMT (RHRS)
	
	//cout << time_rf << endl;
	//h1->Fill(time_rf);
	if(timeL>0 && timeR>0 &&
	   timeL_R>0 && timeR_R>0){
	  
	  // ------- Path length reconstruction (RHRS) ------ //
	  XFP_L   = xdet_L[0];
	  YFP_L   = ydet_L[0];
	  XpFP_L  = xpdet_L[0];
	  YpFP_L  = ypdet_L[0];
	  XFP_L   =(XFP_L -XFPm)/XFPr;
	  XpFP_L  =(XpFP_L-XpFPm)/XpFPr;
	  YFP_L   =(YFP_L -YFPm)/YFPr;
	  YpFP_L  =(YpFP_L-YpFPm)/YpFPr;
	  LenL  = calcf2t(Plen,XFP_L,XpFP_L,YFP_L,YpFP_L);
	  LenL  = LenL*PLr+PLm;
	  //LenL  = pathl[0];
	  //cout << LenL << endl;
	  XFP_L   = XFP_L*XFPr + XFPm;
	  XpFP_L  = XpFP_L*XpFPr + XpFPm;
	  YFP_L   = YFP_L*YFPr + YFPm;
	  YpFP_L  = YpFP_L*YpFPr + YpFPm;
	  
	  // ------- Path length reconstruction (RHRS) ------ //
	  XFP = xdet[0];
	  YFP = ydet[0];
	  XpFP = xpdet[0];
	  YpFP = ypdet[0];
	  XFP   =(XFP -XFPm)/XFPr;
	  XpFP  =(XpFP-XpFPm)/XpFPr;
	  YFP   =(YFP -YFPm)/YFPr;
	  YpFP  =(YpFP-YpFPm)/YpFPr;
	  LenR  = calcf2t(Plen,XFP,XpFP,YFP,YpFP);
	  LenR  = LenR*PLr+PLm;
	  XFP   = XFP*XFPr + XFPm;
	  XpFP  = XpFP*XpFPr + XpFPm;
	  YFP   = YFP*YFPr + YFPm;
	  YpFP  = YpFP*YpFPr + YpFPm;
	  
	  //if(timeL>-50000 && timeR>-5000){	
	  beta = mom[0]/sqrt(pow(mom[0],2.0)+pow(me,2.0));
	  //cout << beta << endl;
	  //double cor   = (pathl[0]-pathl_s2[0])/3.0e+8/beta * 1.0e+9; // ns
	  //double cor   = (pathl[0]-pathl_s2[0])/3.0e+8/beta * 1.0e+9; // ns
	  //double cor   = (LenL-pathl_s2[0])/3.0e+8/beta * 1.0e+9; // ns
	  double cor   = (LenL-3.18)/3.0e+8/beta * 1.0e+9; // ns (test)
	  //cout << LenL << " " << pathl[0] << endl;
	  //cout << XFP_L << " " << XpFP_L << " " << YFP_L << " " << YpFP_L 
	  //cout     <<   " " << mom[0] << " " << LenL << " " << pathl_s2[0] << " " << beta << endl;
	  //cout << LenL-pathl_s2[0] << endl;
	  //cout << cor << endl;
	  beta_R = mom_R[0]/sqrt(pow(mom_R[0],2.0)+pow(mpi,2.0));
	  //double cor_R = (pathl_R[0]-pathl_s2_R[0])/3.0e+8/beta * 1.0e+9; // ns
	  double cor_R = (LenR-pathl_s2_R[0])/3.0e+8/beta_R * 1.0e+9; // ns
//	  cout << XFP << " " << XpFP << " " << YFP << " " << YpFP 
//	       <<   " " << mom_R[0] << " " << LenR << " " << pathl_s2_R[0] << " " << beta_R << endl;
	  //cout << cor << " " << cor_R << endl;
	  //cout << LenL << " " << LenR << endl;
	  meantime   = tref   - (timeL+timeR)/2.0 + 1762.+cor;
	  meantime_R = tref_R - (timeL_R+timeR_R)/2.0 + toffset_R + cor_R; // pathlength correction
	  meantime_R = meantime_R - s2_tzero[seg_R] - s2_tzero_L[seg]; // T0 correction 
	  double yfp_cor   = YFP * par_rtime_ycor[0] + YpFP * par_rtime_ycor[1];
	  //double yfp_cor_L = YFP_L * par_rtime_ycor_L[0] + YpFP_L * par_rtime_ycor_L[1];
	  valval[0] = XFP_L;
	  valval[1] = XpFP_L;
	  valval[2] = YFP_L;
	  valval[3] = YpFP_L;
	  double yfp_cor_L = Calc_FPcor(valval,par_rtime_ycor_L);
	  //cout << yfp_cor_L << endl;
	  meantime_R = meantime_R + yfp_cor + yfp_cor_L;
	  //double lenl_cor = -1.67*LenL + 43.08;
	  //double lenl_cor = -1.*(LenL-pathl_s2[0]) + 37;
	  //double lenl_cor = 1.80*LenL - 31.08;
	  //double lenl_cor = par_pathl_L_cor[0]*LenL + par_pathl_L_cor[1];
	  //double lenl_cor = par_pathl_L_cor[0]*(LenL-pathl_s2[0]) + par_pathl_L_cor[1];
	  //cout << lenl_cor << endl;
	  //lenl_cor = 0.0;
	  //meantime_R = meantime_R +lenl_cor;
	  //meantime_R = meantime_R;
	  meantime_R = meantime_R-cor+75.4;
	  //cout << cor << endl;
	  //meantime_R = meantime_R;
	  
	  time_rf   = rf-meantime;
	  time_rf_R = rf-meantime_R-tref_R;
	  //msq = (1./pow(beta,2.0))
	  //
	  
	  double ctime;
	  ctime = meantime - meantime_R;
	  //ctime = (meantime-cor) - (meantime_R-cor_R);
	  //ctime = meantime_R-cor_R;
	  //time_rf = rf-timeL;
	  //h1->Fill(meantime);
	  //h1->Fill(timeL);
	  //h1->Fill(time_rf-cor);
	  //h1->Fill(meantime-cor-150);
	  //h1->Fill(ctime+2430.);
	  //if(abs(meantime_R)<30.0 && 
	  if(abs(meantime_R)<3000.0 && 
	     //abs(meantime_R)<100.0 && 
	     abs(meantime)<300  ){
	    h1->Fill(ctime);
	    h2->Fill(meantime);
	    h3->Fill(meantime_R);
	    h4->Fill(LenR-pathl_s2_R[0],meantime_R);
	    h5->Fill(pathl_R[0]-LenR);
	    h7->Fill(xdet[0],meantime_R);
	    h8->Fill(xpdet[0],meantime_R);
	    h9->Fill(ydet[0],meantime_R);
	    h10->Fill(ypdet[0],meantime_R);
	    hs2[seg_R]->Fill(meantime_R);
	    hs2_L[seg]->Fill(meantime_R);
	    h11->Fill(time_rf_R);
	    h12->Fill(meantime_R);
	    //h13->Fill(LenL-pathl_s2[0],meantime_R);
	    h13->Fill(LenL-3.18,meantime_R); // (test)
	    h14->Fill(xdet_L[0],meantime_R);
	    h15->Fill(xpdet_L[0],meantime_R);
	    h16->Fill(ydet_L[0],meantime_R);
	    h17->Fill(ypdet_L[0],meantime_R);

	    
	    if(ntune<tunemax){
	      tval0[ntune] = XFP;
	      tval1[ntune] = XpFP;
	      tval2[ntune] = YFP;
	      tval3[ntune] = YpFP;
	      tval4[ntune] = LenR;
	      betabeta[ntune] = beta_R;
	      
	      tval0_L[ntune] = XFP_L;
	      tval1_L[ntune] = XpFP_L;
	      tval2_L[ntune] = YFP_L;
	      tval3_L[ntune] = YpFP_L;
	      tval4_L[ntune] = LenL;
	      betabeta_L[ntune] = beta;
	      //timet[ntune] = tref_R - (timeL_R+timeR_R)/2.0 + toffset_R;
	      //timet[ntune] = meantime_R - yfp_cor_L - yfp_cor;
	      
	      //timet[ntune] = meantime_R - yfp_cor_L;
	      timet[ntune] = meantime_R;
	      //timet[ntune] = meantime_R - cor_R;
	      s2plen[ntune]= pathl_s2_R[0];
	      s2plen_L[ntune]= pathl_s2[0];
	      ntune++;
	    }

	  }
	  //h1->Fill(meantime);
	  //h1->Fill(meantime-tref);
	}
      }
    }
  }

  double chi2 = 0.0;
  int nite = 0; // no tuning 
  //int nite = 1; //
  //int nite = 5; //
  for(int i=0 ; i<nite; i++){
    //chi2 = tune(Plen,i); // path length matrix
    //chi2 = tune2(par_rtime_ycor,i); // correction with RHRS-y and y'
    //chi2 = tune3(par_rtime_ycor_L,i); // correction with LHRS-y and y'
    chi2 = tune4(par_pathl_L_cor,i); // correction with LHRS-path length
    cout << endl;
    cout << "i: " << chi2 << endl;
  }
  
  /*
  for(int i=0 ; i<ent ; i++){
    t1->GetEntry(i);
    XFP = xdet[0];
    YFP = ydet[0];
    XpFP = xpdet[0];
    YpFP = xpdet[0];
    XFP   =(XFP -XFPm)/XFPr;
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    LenR  = calcf2t(Plen_opt,XFP,XpFP,YFP,YpFP);
    LenR  = LenR*PLr+PLm;
    beta_R = mom_R[0]/sqrt(pow(mom_R[0],2.0)+pow(mpi,2.0));
    //double cor_R = (pathl_R[0]-pathl_s2_R[0])/3.0e+8/beta * 1.0e+9; // ns
    double cor_R = (LenR-pathl_s2_R[0])/3.0e+8/beta_R * 1.0e+9; // ns
    //meantime   = tref   - (timeL+timeR)/2.0 + 1962.-cor;
    meantime_R = tref_R - (timeL_R+timeR_R)/2.0 + toffset_R - cor_R;
    h6->Fill(meantime_R);
  }
  */


  // ---- For path length matrix ----- 
  /*
  ofstream* ofs = new ofstream("len_newmat.dat");
  int nppp = 0;
  const int nn = 3; // 3rd order matrix
  for(int i=0 ; i<nn+1 ; i++){
    for(int d=0 ; d<nn+1 ; d++){
      for(int c=0 ; c<nn+1 ; c++){
	for(int b=0 ; b<nn+1 ; b++){
	  for(int a=0 ; a<nn+1 ; a++){  
	    if(a+b+c+d==i){
	      *ofs << Plen_opt[nppp] 
		   << " " << a 
		   << " " << b
		   << " " << c
		   << " " << d << endl;
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
  ofs->close();
  */
  
  /*
  ofstream* ofs = new ofstream("out_ycor.dat");
  for(int i=0 ; i<npar_rtime_ycor ; i++){
    *ofs << par_rtime_ycor[i] << endl;
    //cout << par_rtime_ycor[i] << endl;
  }
  ofs->close();
  */
  
  /*
  ofstream* ofs = new ofstream("out_ycor_L.dat");
  for(int i=0 ; i<npar_rtime_ycor_L ; i++){
    *ofs << par_rtime_ycor_L[i] << endl;
    //cout << par_rtime_ycor[i] << endl;
  }
  ofs->close();
  */
  
  ofstream* ofs = new ofstream("out_pathl_L.dat");
  for(int i=0 ; i<npar_pathl_L_cor ; i++){
    *ofs << par_pathl_L_cor[i] << endl;
    //cout << par_rtime_ycor[i] << endl;
  }
  ofs->close();
  
  
  // ================================= //
  // ====== Draw the histograms ====== //
  // ================================= //
  TCanvas* c1 = new TCanvas("c1","c1");
  h1->Draw();
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->Divide(1,2);
  c2->cd(1);h2->Draw();
  c2->cd(2);h3->Draw();
  TCanvas* c3 = new TCanvas("c3","c3");
  h4->Draw("col");
  TCanvas* c4 = new TCanvas("c4","c4");
  h5->Draw();
  TCanvas* c5 = new TCanvas("c5","c5");
  h6->Draw();
  TCanvas* c6 = new TCanvas("c6","c6");
  c6->Divide(2,2);
  c6->cd(1);h7->Draw();
  c6->cd(2);h8->Draw();
  c6->cd(3);h9->Draw();
  c6->cd(4);h10->Draw();
  TCanvas* c7 = new TCanvas("c7","c7");
  c7->Divide(4,4);
  for(int i=0 ; i<n ; i++){
    c7->cd(i+1);hs2[i]->Draw();
  }
  
  TFile* fnew = new TFile("out.root","recreate");
  //TTree* tnew = new TTree("tree","histfile");
  TObjArray h(0);
  for(int i=0; i<n ; i++){
    h.Add(hs2[i]);
    h.Add(hs2_L[i]);
  }
  h.Write();
  fnew->Close();
  
  TCanvas* c8 = new TCanvas("c8","c8");
  h11->Draw();
  TCanvas* c9 = new TCanvas("c9","c9");
  h12->Draw();
  TCanvas* c10 = new TCanvas("c10","c10");
  h13->Draw();
  TCanvas* c11 = new TCanvas("c11","c11");
  c11->Divide(2,2);
  c11->cd(1);h14->Draw();
  c11->cd(2);h15->Draw();
  c11->cd(3);h16->Draw();
  c11->cd(4);h17->Draw();
  TCanvas* c12 = new TCanvas("c12","c12");
  c12->Divide(4,4);
  for(int i=0 ; i<n ; i++){
    c12->cd(i+1);hs2_L[i]->Draw();
  }
  
  
  
  app.Run();// TApplication for display
  
}


// #############################################################
double tune(double* pa, int j)
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10];
  int ierflg = 0;
  int allparam = nParamT;
  
  //cout << allparam << endl;
  
  TMinuit* minuit = new TMinuit(allparam);
  minuit->SetFCN(fcn);
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];
  double ULim[allparam];
  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];
    step[i] = 1.0e-3;
    //if(i==4 || i==5) {
    //  step[i] = 0.0;
    //  start[i] = 0.0;
    //  pa[i]=0.0;
    //}
    //if(i==0 || i==1 || i==2) {
    //  step[i] = 0.0;
    //  start[i] = 0.0;
    //  pa[i]=0.0;
    //}
    //LLim[i] = pa[i] - pa[i]*0.05;
    //ULim[i] = pa[i] + pa[i]*0.05;
    LLim[i] = pa[i] - pa[i]*0.8;
    ULim[i] = pa[i] + pa[i]*0.8;
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }

  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  // ~~~~ Migrad + Simplex  ~~~~
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg);
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double e;
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);
  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    //minuit -> GetParameter(i,par[i],e);
    minuit -> GetParameter(i,Plen_opt[i],e);
  }
  
  return chi2;
}

// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  double chi2 = 0.0;
  //double chi2_tof = 0.0;
  //double chi2_1 = 0.0, chi2_2 = 0.0, chi2_3 = 0.0, chi2_4 = 0.0;
  //double newtime1, newtime2, newtime3, newtime4;
  //double newtof;
  //double res_tof; //residual
  //double res1, res2, res3, res4; //residual
  //const double reftof  = 0.0; // reference [ns]
  //const double reftime = 0.0; // reference [ns]
  //double n=0;
  //const double sigma_tof   = 0.120; // Sigma of TOF [ns]
  //const double sigma_time  = 0.070; // Sigma for normal counter [ns]
  //const double sigma_time0 = 0.040; // Sigma for timing counter [ns]
  double XFP, XpFP;
  double YFP, YpFP;
  double LenR_new;
  double time_tar;
  const double timeref = 0.0;
  double res_sq = 0.0;
  const double sigma = 0.3;
  
  for(int i=0 ; i<ntune ; i++){
    XFP   = tval0[i];
    XpFP  = tval1[i];
    YFP   = tval2[i];
    YpFP  = tval3[i];
    
    //cout << XFP << " " << YFP << " " << XpFP << " " << YpFP << endl;
    XFP   =(XFP -XFPm)/XFPr;
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    //LenR  = calcf2t(Plen,XFP,XpFP,YFP,YpFP);
    LenR_new  = calcf2t(param,XFP,XpFP,YFP,YpFP);
    //LenR_new  = calcf2t(Plen,XFP,XpFP,YFP,YpFP);
    LenR_new  = LenR_new*PLr+PLm;
    //cout << LenR_new << endl;
    
    double cor_R = (LenR_new - s2plen[i])/3.0e+8/betabeta[i] * 1.0e+9; // ns
    time_tar = timet[i] + cor_R;
    res_sq = res_sq + pow(time_tar-timeref,2.0);
  }

  chi2 = res_sq/(double)ntune/sigma;
  
  //chi2_1 = chi2_1/(double)n;
  //chi2_2 = chi2_2/(double)n;
  //chi2_3 = chi2_3/(double)n;
  //chi2_4 = chi2_4/(double)n;
  //chi2_tof = chi2_tof/(double)n;
  //
  //double w1, w2, w3, w4, wtof; // weight
  // double factor = 1./6.86; //
  //w1 = 0.0;
  //w2 = 1.0;
  //w3 = 1.0;
  //w4 = 1.0;
  //wtof = 3.0;
  
  //chi2_1 = chi2_1*w1/pow(sigma_time0 ,2.0);
  //chi2_2 = chi2_2*w2/pow(sigma_time ,2.0);
  //chi2_3 = chi2_3*w3/pow(sigma_time ,2.0);
  //chi2_4 = chi2_4*w4/pow(sigma_time ,2.0);
  //chi2_tof = chi2_tof*wtof/pow(sigma_tof ,2.0);
  //
  //chi2 = chi2_1 + chi2_2 + chi2_3 + chi2_4 + chi2_tof;
  //chi2 = chi2 / (w1+w2+w3+w4+wtof);

  fval = chi2;
}


// #############################################################
double tune2(double* pa, int j)
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10];
  int ierflg = 0;
  int allparam = npar_rtime_ycor;
  
  //cout << allparam << endl;
  
  TMinuit* minuit = new TMinuit(allparam);
  minuit->SetFCN(fcn2);
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];
  double ULim[allparam];
  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];
    step[i] = 1.0e-3;
    //if(i==4 || i==5) {
    //  step[i] = 0.0;
    //  start[i] = 0.0;
    //  pa[i]=0.0;
    //}
    //if(i==0 || i==1 || i==2) {
    //  step[i] = 0.0;
    //  start[i] = 0.0;
    //  pa[i]=0.0;
    //}
    //LLim[i] = pa[i] - pa[i]*0.05;
    //ULim[i] = pa[i] + pa[i]*0.05;
    LLim[i] = pa[i] - pa[i]*0.8;
    ULim[i] = pa[i] + pa[i]*0.8;
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }

  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  // ~~~~ Migrad + Simplex  ~~~~
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg);
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double e;
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);
  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    //minuit -> GetParameter(i,par[i],e);
    minuit -> GetParameter(i,par_rtime_ycor[i],e);
  }
  
  return chi2;
}

// #############################################################
void fcn2(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  double chi2 = 0.0;
  double XFP, XpFP;
  double YFP, YpFP;
  double Len;
  double time_tar;
  const double timeref = 0.0;
  double res_sq = 0.0;
  const double sigma = 0.3;
  
  for(int i=0 ; i<ntune ; i++){
    XFP  = tval0[i];
    XpFP = tval1[i];
    YFP  = tval2[i];
    YpFP = tval3[i];
    //Len  = tval4[i];
    
    //XFP   =(XFP -XFPm)/XFPr;
    //XpFP  =(XpFP-XpFPm)/XpFPr;
    //YFP   =(YFP -YFPm)/YFPr;
    //YpFP  =(YpFP-YpFPm)/YpFPr;
    
    //double cor_R = (Len - s2plen[i])/3.0e+8/betabeta[i] * 1.0e+9; // ns
    //time_tar = timet[i] + cor_R;
    time_tar = timet[i] + param[0]*YFP + param[1]*YpFP;
    res_sq = res_sq + pow(time_tar-timeref,2.0);
  }

  chi2 = res_sq/(double)ntune/sigma;
  
  //chi2_1 = chi2_1/(double)n;
  //chi2_2 = chi2_2/(double)n;
  //chi2_3 = chi2_3/(double)n;
  //chi2_4 = chi2_4/(double)n;
  //chi2_tof = chi2_tof/(double)n;
  //
  //double w1, w2, w3, w4, wtof; // weight
  // double factor = 1./6.86; //
  //w1 = 0.0;
  //w2 = 1.0;
  //w3 = 1.0;
  //w4 = 1.0;
  //wtof = 3.0;
  
  //chi2_1 = chi2_1*w1/pow(sigma_time0 ,2.0);
  //chi2_2 = chi2_2*w2/pow(sigma_time ,2.0);
  //chi2_3 = chi2_3*w3/pow(sigma_time ,2.0);
  //chi2_4 = chi2_4*w4/pow(sigma_time ,2.0);
  //chi2_tof = chi2_tof*wtof/pow(sigma_tof ,2.0);
  //
  //chi2 = chi2_1 + chi2_2 + chi2_3 + chi2_4 + chi2_tof;
  //chi2 = chi2 / (w1+w2+w3+w4+wtof);

  fval = chi2;
}

double tune3(double* pa, int j)
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10];
  int ierflg = 0;
  int allparam = npar_rtime_ycor_L;
  
  //cout << allparam << endl;
  
  TMinuit* minuit = new TMinuit(allparam);
  minuit->SetFCN(fcn3);
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];
  double ULim[allparam];
  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];
    //step[i] = 1.0e-3;
    //step[i] = 1.0e-2;
    step[i] = 1.0;

    //if(i==4 || i==5) {
    //  step[i] = 0.0;
    //  start[i] = 0.0;
    //  pa[i]=0.0;
    //}
    //if(i==0 || i==1 || i==2) {
    //  step[i] = 0.0;
    //  start[i] = 0.0;
    //  pa[i]=0.0;
    //}
    //LLim[i] = pa[i] - pa[i]*0.05;
    //ULim[i] = pa[i] + pa[i]*0.05;
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
    LLim[i] = pa[i] -100.0;
    ULim[i] = pa[i] +100.0;
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }

  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  // ~~~~ Migrad + Simplex  ~~~~
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg);
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double e;
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);
  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    //minuit -> GetParameter(i,par[i],e);
    minuit -> GetParameter(i,par_rtime_ycor_L[i],e);
  }
  
  return chi2;
}

// #############################################################
void fcn3(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  double chi2 = 0.0;
  double XFP, XpFP;
  double YFP, YpFP;
  double Len;
  double time_tar;
  const double timeref = 0.0;
  double res_sq = 0.0;
  const double sigma = 0.3;
  double yfp_cor_L;
  
  for(int i=0 ; i<ntune ; i++){
    XFP  = tval0_L[i];
    XpFP = tval1_L[i];
    YFP  = tval2_L[i];
    YpFP = tval3_L[i];
    //Len  = tval4[i];
    
    //cout 
    //XFP   =(XFP -XFPm)/XFPr;
    //XpFP  =(XpFP-XpFPm)/XpFPr;
    //YFP   =(YFP -YFPm)/YFPr;
    //YpFP  =(YpFP-YpFPm)/YpFPr;
    yfp_cor_L=0.0;
    //double cor_R = (Len - s2plen[i])/3.0e+8/betabeta[i] * 1.0e+9; // ns
    //time_tar = timet[i] + cor_R;
    yfp_cor_L = Calc_FPcor(tval0_L, param);
    //cout << tval0_L[0] << " " << tval1_L[1] << " " << param[0] << endl;
    //cout << param[0] << " " <<param[1] << " " << param[5] << endl;
    //cout << yfp_cor_L << endl;
    time_tar = timet[i] + yfp_cor_L;
    //time_tar = timet[i] + param[0]*YFP + param[1]*YpFP;
    res_sq = res_sq + pow(time_tar-timeref,2.0);
  }

  chi2 = res_sq/(double)ntune/sigma;
  
  //chi2_1 = chi2_1/(double)n;
  //chi2_2 = chi2_2/(double)n;
  //chi2_3 = chi2_3/(double)n;
  //chi2_4 = chi2_4/(double)n;
  //chi2_tof = chi2_tof/(double)n;
  //
  //double w1, w2, w3, w4, wtof; // weight
  // double factor = 1./6.86; //
  //w1 = 0.0;
  //w2 = 1.0;
  //w3 = 1.0;
  //w4 = 1.0;
  //wtof = 3.0;
  
  //chi2_1 = chi2_1*w1/pow(sigma_time0 ,2.0);
  //chi2_2 = chi2_2*w2/pow(sigma_time ,2.0);
  //chi2_3 = chi2_3*w3/pow(sigma_time ,2.0);
  //chi2_4 = chi2_4*w4/pow(sigma_time ,2.0);
  //chi2_tof = chi2_tof*wtof/pow(sigma_tof ,2.0);
  //
  //chi2 = chi2_1 + chi2_2 + chi2_3 + chi2_4 + chi2_tof;
  //chi2 = chi2 / (w1+w2+w3+w4+wtof);

  fval = chi2;
}


double Calc_FPcor(double* val, double* par){
  double x  = val[0];
  double xp = val[1];
  double y  = val[2];
  double yp = val[3];
  //const int nn = npar_rtime_ycor_L;
  //double par[nn];
  //for(int i=0 ; i<npar_rtime_ycor_L ; i++){
  //  par[i] = par_rtime_ycor_L[i];
  //}
  double cor = 0.0; 
  double cor1=0.0, cor2=0.0, cor3=0.0;
  //cout << yp << " " << par[3] << endl;
  
  cor1 = par[0]*y + par[1]*yp;
  //cor 1 
  //cor1 = par[0]*x + par[1]*xp + par[2]*y + par[3]*yp;
  //cor2 = par[4]*x*xp + par[5]*x*y + par[6]*x*yp + par[7]*xp*y + par[8]*xp*yp + par[9]*y*yp;
  //cor3 = par[10]*x*x + par[11]*xp*xp + par[12]*y*y + par[13]*yp*yp;
  //cout << cor << endl;
  
  cor = cor1+cor2+cor3;
  return cor;
  
}


double tune4(double* pa, int j)
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10];
  int ierflg = 0;
  int allparam = npar_pathl_L_cor;
  
  //cout << allparam << endl;
  
  TMinuit* minuit = new TMinuit(allparam);
  minuit->SetFCN(fcn4);
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];
  double ULim[allparam];
  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];
    step[i] = 1.0e-3;
    //step[i] = 1.0e-2;
    //step[i] = 0.1;

    //if(i==4 || i==5) {
    //  step[i] = 0.0;
    //  start[i] = 0.0;
    //  pa[i]=0.0;
    //}
    //if(i==0 || i==1 || i==2) {
    //  step[i] = 0.0;
    //  start[i] = 0.0;
    //  pa[i]=0.0;
    //}
    //LLim[i] = pa[i] - pa[i]*0.05;
    //ULim[i] = pa[i] + pa[i]*0.05;
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
    LLim[i] = pa[i] -100.0;
    ULim[i] = pa[i] +100.0;
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }

  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  // ~~~~ Migrad + Simplex  ~~~~
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg);
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double e;
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);
  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    //minuit -> GetParameter(i,par[i],e);
    minuit -> GetParameter(i,par_pathl_L_cor[i],e);
  }
  
  return chi2;
}

// #############################################################
void fcn4(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  double chi2 = 0.0;
  double XFP, XpFP;
  double YFP, YpFP;
  double Len;
  double time_tar;
  const double timeref = 0.0;
  double res_sq = 0.0;
  const double sigma = 0.5;
  double pathl_L_cor = 0.0;
  double s2len = 0.0;
  
  for(int i=0 ; i<ntune ; i++){
    XFP  = tval0_L[i];
    XpFP = tval1_L[i];
    YFP  = tval2_L[i];
    YpFP = tval3_L[i];
    Len  = tval4_L[i];
    s2len = s2plen_L[i];
    //cout << Len-s2len  << endl;
    cout << param[0] << endl;
    //cout 
    //XFP   =(XFP -XFPm)/XFPr;
    //XpFP  =(XpFP-XpFPm)/XpFPr;
    //YFP   =(YFP -YFPm)/YFPr;
    //YpFP  =(YpFP-YpFPm)/YpFPr;
    pathl_L_cor = 0.0;
    pathl_L_cor = (Len-s2len)*param[0] + param[1];
    time_tar = timet[i] + pathl_L_cor;
    res_sq = res_sq + pow(time_tar-timeref,2.0);
  }

  chi2 = sqrt(res_sq)/(double)ntune/sigma;
  
  //chi2_1 = chi2_1/(double)n;
  //chi2_2 = chi2_2/(double)n;
  //chi2_3 = chi2_3/(double)n;
  //chi2_4 = chi2_4/(double)n;
  //chi2_tof = chi2_tof/(double)n;
  //
  //double w1, w2, w3, w4, wtof; // weight
  // double factor = 1./6.86; //
  //w1 = 0.0;
  //w2 = 1.0;
  //w3 = 1.0;
  //w4 = 1.0;
  //wtof = 3.0;
  
  //chi2_1 = chi2_1*w1/pow(sigma_time0 ,2.0);
  //chi2_2 = chi2_2*w2/pow(sigma_time ,2.0);
  //chi2_3 = chi2_3*w3/pow(sigma_time ,2.0);
  //chi2_4 = chi2_4*w4/pow(sigma_time ,2.0);
  //chi2_tof = chi2_tof*wtof/pow(sigma_tof ,2.0);
  //
  //chi2 = chi2_1 + chi2_2 + chi2_3 + chi2_4 + chi2_tof;
  //chi2 = chi2 / (w1+w2+w3+w4+wtof);

  fval = chi2;
}
