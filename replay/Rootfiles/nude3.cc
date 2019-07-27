/*
  nude3.cc
  "Stripping unnecessary data from original ROOT file"
  
  Toshiyuki Gogami, Nov 9, 2018
*/

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
using namespace std;

const double  XFPm=-0.7, XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18; 
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74;
const double  XFPr=1.3, XpFPr=0.27; 
const double  YFPr=0.1, YpFPr=0.10; 
const double  Xptr=0.15, Yptr=0.08, Momr=0.18;
const double  PLm = 25.4, PLr=0.7;
const double  Ztm = -0.15, Ztr=0.35;
extern double calcf2t_plen(double* P, 
			   double xf, double xpf,
			   double yf, double ypf);
extern double calcf2t_4th_2(double* P, double xf, double xpf, 
			    double yf, double ypf, double zt);
extern double Calc_FPcor(double* val, double* par);

const int nParamT=35;     // For path length matrix (3rd order)
//double Plen_opt[nParamT]; // For path length matrix (3rd order)
const int nParamT_4=126;  // For xpt, ypt with z (4th order)

const int npar_rtime_ycor = 2;
double par_rtime_ycor[npar_rtime_ycor];
const int npar_rtime_ycor_L = 2;
double par_rtime_ycor_L[npar_rtime_ycor_L];
const int npar_pathl_L_cor = 2;
double par_pathl_L_cor[npar_pathl_L_cor];
//const double toffset_R = -364.6-150.;
double toffset_R = -364.6-150.; // for H2_1
const double toffset_L = 1762.0;
//const double ch2time = 56.0e-12 * 1.0e+9; // 56 (ps/ch); F1TDC 
double ch2time = 56.0e-12 * 1.0e+9; // 56 (ps/ch); F1TDC 
const double mpi = 0.13957;  // GeV/c2
const double me  = 0.000511; // GeV/c2

int nude3(int run=111179, int nf=0, int tflag=5){
  int dataflag = 1;
  if (run>111400){
    dataflag = 2;
    ch2time = 58.0e-12 * 1.0e+9;
    toffset_R = -364.6-150.-15-0.788-0.113; // for H2_2
  }
  //int main(int argc, char** argv){
  //int run=111179, nf=0;
  //if(argc==2){
  //  run = atoi(argv[1]);
  //}
  //else if(argc==3){
  //  run = atoi(argv[1]);
  //  nf  = atoi(argv[2]);
  //}
  //nf=0;
  
  //if(tflag==5) tflag=5;
  //if(tflag==4 ) tflag=5;
  //tflag: 4 (RHRS), 5 (coin)
  
  char inputfname[500], newfname[500];;
  if(nf<1){
    sprintf(inputfname,"./nnL/tritium_%d.root",run);
    sprintf(newfname,"./nnL/nude_dir2/nude_%d.root",run);
    if(tflag==5){
      //sprintf(newfname,"./nnL/coin_dragon2/tri_coin_%d.root",run);
      sprintf(newfname,"./nnL_/coin/tri_coin_%d.root",run);
    }
    else if(tflag==4){
      sprintf(newfname,"./nnL/RHRS_single_dragon/tri_Rsignle_%d.root",run);
    }
    else if(tflag==1){
      sprintf(newfname,"./nnL/LHRS_single_dragon/tri_Lsingle_%d.root",run);
    }
    else sprintf(newfname,"nude_dir2/nude_%d_%d.root",run,nf);
    //cout << "aaaaaaaaaaa" << endl;
  }
  else{
    sprintf(inputfname,"./nnL/tritium_%d_%d.root",run,nf);
    //sprintf(newfname,"nude_dir2/nude_%d_%d.root",run,nf);
    //sprintf(newfname,"nude_dir2/nude_%d_%d.root",run,nf);
    
    if(tflag==5){
      //sprintf(newfname,"./nnL/coin_dragon2/tri_coin_%d_%d.root",run,nf);
      sprintf(newfname,"./nnL_/coin/tri_coin_%d_%d.root",run,nf);
    }
    else if(tflag==4){
      sprintf(newfname,"./nnL/RHRS_single_dragon/tri_Rsingle_%d_%d.root",run,nf);
    }
    else if(tflag==1){
      sprintf(newfname,"./nnL/LHRS_single_dragon/tri_Lsingle_%d_%d.root",run,nf);
    }
    else sprintf(newfname,"./nnL/nude_dir2/nude_%d_%d.root",run,nf);
  }
  
  //sprintf(inputfname,"tritium_%d_%d.root",run,nf);
  //3565139
  TFile* f1 = new TFile(inputfname);
  if( !f1->IsOpen() ) {
    cout << " No root file: " << inputfname << endl;
    return 1;
  }
  
  TTree* t1 = (TTree*)f1->Get("T");
  const int max = 100;
  Double_t trig5[max];
  Double_t trig4[max];
  Double_t trig1[max];
  double ent = t1->GetEntries();
  //ent = 20000; // for test
  cout << endl;
  cout << " Stripping " << inputfname 
       << " (ev=" << ent 
       << ") --> " << newfname 
       << endl;
  
  double rtime_s0[max], ltime_s0[max];
  double rtime_s2[max], ltime_s2[max];
  double rtime[max], ltime[max];
  double rpathl[max], lpathl[max];
  double rpathl_s2[max], lpathl_s2[max];
  double a1, a2;
  double mom1[max], mom2[max];
  const int f1n = 64;
  double rf1tdc[f1n];
  double lf1tdc[f1n];
  double rvz[max], lvz[max];
  double vz_mean[max];
  double th1[max], ph1[max];
  double th2[max], ph2[max];
  Int_t runnum;
  double hallap;
  double r_s2_t_pads[max];
  double l_s2_t_pads[max];
  double r_s2_nthit;
  double l_s2_nthit;
  double r_th_fp[max];
  double l_th_fp[max];
  double r_ph_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double r_x_fp[max];
  double l_y_fp[max];
  double r_y_fp[max];
  const int n = 16;
  double r_s2_la_c[n];
  double r_s2_ra_c[n];
  double l_s2_la_c[n];
  double l_s2_ra_c[n];
  double rbeta[max];
  double lbeta[max];
  double nhit, nhit_R;
  double ps_asum;
  double a1_tdc[24];
  double a2_tdc[26];
  double rasterx, rastery;
  double dpp;
  UInt_t evid;
  
  //t1->SetBranchAddress("fEvtHdr.fRun", &runnum    );
  //t1->SetBranchAddress("fEvtHdr.fEvtNum", &evid    );
  t1->SetBranchAddress("HALLA_p", &hallap );
  t1->SetBranchAddress("HALLA_dpp", &dpp );
  t1->SetBranchAddress("DR.T1", &trig1    );
  t1->SetBranchAddress("DR.T4", &trig4    );
  t1->SetBranchAddress("DR.T5", &trig5    );
  t1->SetBranchAddress("R.tr.time", &rtime);
  t1->SetBranchAddress("L.tr.time", &ltime);
  t1->SetBranchAddress("R.tr.pathl", &rpathl);
  t1->SetBranchAddress("L.tr.pathl", &lpathl);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  t1->SetBranchAddress("R.tr.p", &mom1);
  t1->SetBranchAddress("L.tr.p", &mom2);
  t1->SetBranchAddress("RTDC.F1FirstHit", &rf1tdc);
  t1->SetBranchAddress("LTDC.F1FirstHit", &lf1tdc);
  t1->SetBranchAddress("R.tr.vz", &rvz);
  t1->SetBranchAddress("L.tr.vz", &lvz);
  t1->SetBranchAddress("R.tr.tg_th", &th1);
  t1->SetBranchAddress("R.tr.tg_ph", &ph1);

  t1->SetBranchAddress("L.tr.tg_th", &th2);
  t1->SetBranchAddress("L.tr.tg_ph", &ph2);
  t1->SetBranchAddress("R.s0.time", &rtime_s0);
  t1->SetBranchAddress("L.s0.time", &ltime_s0);
  t1->SetBranchAddress("R.s2.time", &rtime_s2);
  t1->SetBranchAddress("L.s2.time", &ltime_s2);
  t1->SetBranchAddress("R.s2.t_pads", &r_s2_t_pads);
  t1->SetBranchAddress("L.s2.t_pads", &l_s2_t_pads);
  t1->SetBranchAddress("R.s2.nthit",   &r_s2_nthit);
  t1->SetBranchAddress("L.s2.nthit",   &l_s2_nthit);
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
  t1->SetBranchAddress("R.s2.la_c",  &r_s2_la_c);
  t1->SetBranchAddress("R.s2.ra_c",  &r_s2_ra_c);
  t1->SetBranchAddress("L.s2.la_c",  &l_s2_la_c);
  t1->SetBranchAddress("L.s2.ra_c",  &l_s2_ra_c);
  t1->SetBranchAddress("R.tr.beta",  &rbeta);
  t1->SetBranchAddress("L.tr.beta",  &lbeta);
  t1->SetBranchAddress("R.s2.trpath",  &rpathl_s2);
  t1->SetBranchAddress("L.s2.trpath",  &lpathl_s2);
  t1->SetBranchAddress("L.s2.nthit",&nhit);
  t1->SetBranchAddress("R.s2.nthit",&nhit_R);
  t1->SetBranchAddress("R.ps.asum_c", &ps_asum);
  t1->SetBranchAddress("R.a1.t_fadc", &a1_tdc);
  t1->SetBranchAddress("R.a2.t_fadc", &a2_tdc);
  t1->SetBranchAddress("FbusRrb.Raster2.target.x", &rasterx);
  t1->SetBranchAddress("FbusRrb.Raster2.target.y", &rastery);
  double rast_curx, rast_cury;
  double rast_x, rast_y;
  double rast_x2; // raster x with new parameters
  double lcer_asum_c;
  t1->SetBranchAddress("Lrb.Raster2.rawcur.x", &rast_curx); // raster current
  t1->SetBranchAddress("Lrb.Raster2.rawcur.y", &rast_cury); // raster current
  t1->SetBranchAddress("Lrb.x", &rast_x);
  t1->SetBranchAddress("Lrb.y", &rast_y);
  t1->SetBranchAddress("L.cer.asum_c",  &lcer_asum_c);
  
  bool t5flag = false;
  bool t4flag = false;
  bool t1flag = false;
  bool genflag = false;
  bool acflag = false;
  bool trig_fire = false;
  double ctime[max];
  
  TFile* fnew = new TFile(newfname,"recreate");
  TTree* tnew = new TTree("tree","3H(e,e'K+)nnL experiment (2018)");
  
  tnew->Branch("DR.T1", &trig1, "DR.T1/D"   );
  tnew->Branch("DR.T4", &trig4, "DR.T4/D"  );
  tnew->Branch("DR.T5", &trig5, "DR.T5/D"   );
  
  //tnew->Branch("fEvtHdr.fRun", &runnum,   "fEvtHdr.fRun/I");
  tnew->Branch("runid", &run,   "runid/I");
  tnew->Branch("evid",  &evid,  "evid/I" );
  tnew->Branch("HALLA_p",  &hallap,    "HALLA_p/D");
  tnew->Branch("HALLA_dpp", &dpp ,     "HALLA_dpp/D");
  tnew->Branch("R.a1.asum_c", &a1,     "R.a1.asum_c/D");
  tnew->Branch("R.a2.asum_c", &a2,     "R.a2.asum_c/D");
  tnew->Branch("RTDC.F1FirstHit", &rf1tdc, "RTDC.F1FirstHit[64]/D");
  tnew->Branch("LTDC.F1FirstHit", &lf1tdc, "LTDC.F1FirstHit[64]/D");
  tnew->Branch("R.tr.tg_th", &th1,      "R.tr.tg_th[100]/D");
  tnew->Branch("R.tr.tg_ph", &ph1,      "R.tr.tg_ph[100]/D");
  tnew->Branch("L.tr.tg_th", &th2,      "L.tr.tg_th[100]/D");
  tnew->Branch("L.tr.tg_ph", &ph2,      "L.tr.tg_ph[100]/D");
  tnew->Branch("R.s0.time", &rtime_s0,  "R.s0.time[100]/D");
  tnew->Branch("L.s0.time", &ltime_s0,  "L.s0.time[100]/D");
  tnew->Branch("R.s2.time", &rtime_s2,  "R.s2.time[100]/D");
  tnew->Branch("L.s2.time", &ltime_s2,  "L.s2.time[100]/D");
  tnew->Branch("R.tr.time", &rtime,     "R.tr.time[100]/D");
  tnew->Branch("L.tr.time", &ltime,     "L.tr.time[100]/D");
  tnew->Branch("R.tr.vz", &rvz,         "R.tr.vz[100]/D");
  tnew->Branch("L.tr.vz", &lvz,         "L.tr.vz[100]/D");
  tnew->Branch("vz_mean", &vz_mean,     "vz_mean[100]/D");
  tnew->Branch("R.tr.p", &mom1,        "R.tr.p[100]/D");
  tnew->Branch("L.tr.p", &mom2,        "L.tr.p[100]/D");
  tnew->Branch("R.tr.pathl", &rpathl,  "R.tr.pathl[100]/D" );
  tnew->Branch("L.tr.pathl", &lpathl,  "L.tr.pathl[100]/D" );
  //tnew->Branch("R.s2.trpath", &rpathl_s2,  "R.s2.trpath[100]/D" );
  //tnew->Branch("L.s2.trpath", &lpathl_s2,  "L.s2.trpath[100]/D" );
  //tnew->Branch("R.s2.t_pads", &r_s2_t_pads, "R.s2.t_pads[100]/D" );
  //tnew->Branch("L.s2.t_pads", &l_s2_t_pads, "L.s2.t_pads[100]/D");
  //tnew->Branch("R.s2.nthit",   &r_s2_nthit,   "R.s2.nthit/D");
  //tnew->Branch("L.s2.nthit",   &l_s2_nthit,   "L.s2.nthit/D");
  tnew->Branch("R.tr.x",   &r_x_fp,  "R.tr.x[100]/D");
  tnew->Branch("L.tr.x",   &l_x_fp,  "L.tr.x[100]/D");
  tnew->Branch("R.tr.y",   &r_y_fp,  "R.tr.y[100]/D");
  tnew->Branch("L.tr.y",   &l_y_fp,  "L.tr.y[100]/D");
  tnew->Branch("R.tr.th",  &r_th_fp, "R.tr.th[100]/D");
  tnew->Branch("L.tr.th",  &l_th_fp, "L.tr.th[100]/D");
  tnew->Branch("R.tr.ph",  &r_ph_fp, "R.tr.ph[100]/D");
  tnew->Branch("L.tr.ph",  &l_ph_fp, "L.tr.ph[100]/D");
  //tnew->Branch("R.s2.la_c",  &r_s2_la_c, "R.s2.la_c[16]/D");
  //tnew->Branch("R.s2.ra_c",  &r_s2_ra_c, "R.s2.ra_c[16]/D");
  //tnew->Branch("L.s2.la_c",  &l_s2_la_c, "L.s2.la_c[16]/D");
  //tnew->Branch("L.s2.ra_c",  &l_s2_ra_c, "L.s2.ra_c[16]/D");
  tnew->Branch("R.tr.beta",  &rbeta, "R.tr.beta[100]/D");
  tnew->Branch("L.tr.beta",  &lbeta, "L.tr.beta[100]/D");
  tnew->Branch("ctime",      &ctime, "ctime[100]/D");
  tnew->Branch("R.ps.asum_c", &ps_asum, "R.ps.asum_c/D");
  //tnew->Branch("R.a1.t_fadc", &a1_tdc, "R.a1.t_fadc[24]/D");
  //tnew->Branch("R.a2.t_fadc", &a2_tdc, "R.a2.t_fadc[26]/D");
  tnew->Branch("FbusRrb.Raster2.target.x", &rasterx, "FbusRrb.Raster2.target.x/D");
  tnew->Branch("FbusRrb.Raster2.target.y", &rastery, "FbusRrb.Raster2.target.y/D");
  
  //tnew->Branch("Lrb.Raster2.rawcur.x", &rast_curx, "Lrb.Raster2.rawcur.x/D");
  //tnew->Branch("Lrb.Raster2.rawcur.y", &rast_cury, "Lrb.Raster2.rawcur.y/D");
  tnew->Branch("Lrb.x",  &rast_x, "Lrb.x/D");
  tnew->Branch("Lrb.x2", &rast_x2, "Lrb.x2/D");
  tnew->Branch("Lrb.y", &rast_y, "Lrb.y/D");
  tnew->Branch("L.cer.asum_c", &lcer_asum_c, "L.cer.asum_c/D");
  
  
  char name_Mlen[100];
  sprintf(name_Mlen,"matrices/len_RHRS_1.dat"); // original
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
  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,"matrices/zt_LHRS_opt.dat"); // optimized
  ifstream Mzt_L(name_Mzt_L);
  double Pzt_L[nParamT];
  //ent = 10000;
  //double Plenopt[nParamT];
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p; 
    Pzt_L[i]=par;
  }
  Mzt_L.close();
  
  char name_Mzt_R[500];
  sprintf(name_Mzt_R,"matrices/zt_RHRS_opt.dat");
  ifstream Mzt_R(name_Mzt_R);
  double Pzt_R[nParamT];
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mzt_R >> par >> p >> p >> p >> p; 
    Pzt_R[i]=par;
  }
  Mzt_R.close();
  
  char name_Mxp_L[500];
  sprintf(name_Mxp_L,"matrices/xp_LHRS_opt.dat");
  ifstream Mxp_L(name_Mxp_L);
  double Pxp_L[nParamT_4];
  for (int i=0;i<nParamT_4;i++){
    double par=0.;
    int p=0;
    Mxp_L >> par >> p >> p >> p >> p; 
    Pxp_L[i]=par;
  }
  Mxp_L.close();
  
  char name_Myp_L[500];
  sprintf(name_Myp_L,"matrices/yp_LHRS_opt.dat");
  ifstream Myp_L(name_Myp_L);
  double Pyp_L[nParamT_4];
  for (int i=0;i<nParamT_4;i++){
    double par=0.;
    int p=0;
    Myp_L >> par >> p >> p >> p >> p; 
    Pyp_L[i]=par;
  }
  Myp_L.close();
  
  ifstream* s2_R_data;
  ifstream* s2_L_data;
  ifstream* rtime_ycor_data_L;
  ifstream* rastx_data;
  if(dataflag == 1){
    s2_R_data = new ifstream("data/s2_t0_R.dat");
    s2_L_data = new ifstream("data/s2_t0_L.dat");
    rtime_ycor_data_L = new ifstream("data/rtime_ycor_L.dat"); 
  }
  else if(dataflag == 2){
    s2_R_data = new ifstream("data/s2_t0_R_2.dat");
    s2_L_data = new ifstream("data/s2_t0_L_2.dat");
    rtime_ycor_data_L = new ifstream("data/rtime_ycor_L_2.dat"); 
  }
  else{
    s2_R_data = new ifstream("data/s2_t0_R.dat");
    s2_L_data = new ifstream("data/s2_t0_L.dat");
    rtime_ycor_data_L = new ifstream("data/rtime_ycor_L.dat"); 
  }
  rastx_data = new ifstream("data/rasterx.dat");
  
  double s2_tzero_R[n];
  for(int i=0 ; i<n ; i++){
    *s2_R_data >> s2_tzero_R[i];
  }
  s2_R_data->close();

  
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
  
  
  double rastx_param[2];
  for(int i=0 ; i<2 ; i++){
    *rastx_data >> rastx_param[i];
  }
  rastx_data->close();
  
  
  int seg_L, seg_R;
  double XFP_R, XpFP_R;
  double YFP_R, YpFP_R;
  double XFP_L, XpFP_L;
  double YFP_L, YpFP_L;
  double LenL, LenR;
  double tref_L, tref_R;
  double timeL_L, timeR_L;
  double timeL_R, timeR_R;
  double valval[max];
  
  for(int i=0 ; i<ent ; i++){
    
    for(int j=0 ; j<max ;j++){
      trig1[j] = 0.0;
      trig4[j] = 0.0;
      trig5[j] = 0.0;
      rtime_s0[j] = -2222.0;
      ltime_s0[j] = -2222.0;
      rpathl[j]   = -2222.0;
      rtime_s2[j] = -2222.0;
      ltime_s2[j] = -2222.0;
      rtime[j]    = -2222.0;
      ltime[j]    = -2222.0;
      mom1[j]   = -2222.0;
      mom2[j]   = -2222.0;
      th1[j]   = -2222.0;
      ph1[j]   = -2222.0;
      th2[j]   = -2222.0;
      ph2[j]   = -2222.0;
      r_s2_t_pads[j] = -2222.0;
      l_s2_t_pads[j] = -2222.0;
      valval[j] = -2222.0;
      ctime[j]  = -2222.0;
      rvz[j]   = -2222.0;
      lvz[j]   = -2222.0;
      vz_mean[j] = -2222.0;
    }
    a1 = -2222.0;
    a2 = -2222.0;
    ps_asum = -2222.0;
    rast_x  = -2222.0;
    rast_x2 = -2222.0;
    rast_y  = -2222.0;
    rast_curx = -2222.0;
    rast_cury = -2222.0;
    lcer_asum_c = -2222.0;
    
    
    t5flag = false;
    t4flag = false;
    t1flag = false;
    genflag= false;
    acflag = false;
    trig_fire = false;
    

    // ------------------------- //
    // -------- Get Entry ------ // 
    // ------------------------- //
    t1->GetEntry(i);
    evid = i;
    
    
    // ------------------------------------ //
    // ------- Trigger conditions  -------- //
    // ------------------------------------ //
    if(trig5[0] > 0 && tflag==5) t5flag = true;
    else t5flag = false;
    
    if(trig4[0] > 0 && tflag==4) t4flag = true;
    else t4flag = false;
    
    if(trig1[0] > 0 && tflag==1) t1flag = true;
    else t1flag = false;
    
    if(t1flag==true || t4flag==true || t5flag==true) trig_fire = true;
    else trig_fire = false;
    

    // ------------------------------------------------ //
    // ------- General event selection ---------------- //
    // ------------------------------------------------ //
    if( nhit == 1 
	&& nhit_R ==1  // Single hit
        && mom1[0]>1.5 && mom1[0]<2.0
	&& mom2[0]>1.5 && mom2[0]<3.0
	) {
      
      seg_L  = l_s2_t_pads[0];
      seg_R  = r_s2_t_pads[0];
      
      if(r_s2_la_c[seg_R]>2.0 && r_s2_ra_c[seg_R]>2.0
	 && l_s2_la_c[seg_L]>2.0 && l_s2_ra_c[seg_L]>2.0
	 ){
	genflag = true;
	XFP_R  = r_x_fp[0];
	YFP_R  = r_y_fp[0];
	XpFP_R = r_th_fp[0];
	YpFP_R = r_ph_fp[0];
	XFP_L  = l_x_fp[0];
	YFP_L  = l_y_fp[0];
	XpFP_L = l_th_fp[0];
	YpFP_L = l_ph_fp[0];
      }
      else genflag = false;
    }
    
    
    // ------------------------------------------- //
    // ------ Aerogel Cherenkov selection -------- //
    // ------------------------------------------- //
    if(a1 > -10.0 && a1 < 30.0 && a2 > -10.0 && a2 < 60.0){
      acflag = true;
    }
    else acflag = false;
    
    if (trig_fire==true	
	&& genflag==true 
	&& acflag==true 
	){
      tref_L  = lf1tdc[40]       * ch2time;
      timeL_L = lf1tdc[seg_L]    * ch2time;
      timeR_L = lf1tdc[seg_L+48] * ch2time;
      
      tref_R  = rf1tdc[9]        * ch2time;
      timeL_R = rf1tdc[seg_R+16] * ch2time;
      timeR_R = rf1tdc[seg_R+48] * ch2time;
      
      if(timeL_L>0.0 && timeR_L>0.0
	 && timeL_R>0.0 && timeR_R>0.0){
	
	// ------- Path length reconstruction (LHRS) ------ //
	XFP_L   = (XFP_L -XFPm)/XFPr;
	XpFP_L  = (XpFP_L-XpFPm)/XpFPr;
	YFP_L   = (YFP_L -YFPm)/YFPr;
	YpFP_L  = (YpFP_L-YpFPm)/YpFPr;
	
	LenL    = calcf2t_plen(Plen,XFP_L,XpFP_L,YFP_L,YpFP_L);
	lvz[0]  = calcf2t_plen(Pzt_L,XFP_L,XpFP_L,YFP_L,YpFP_L);
	
	LenL    = LenL*PLr+PLm;
	XFP_L   = XFP_L*XFPr + XFPm;
	XpFP_L  = XpFP_L*XpFPr + XpFPm;
	YFP_L   = YFP_L*YFPr + YFPm;
	YpFP_L  = YpFP_L*YpFPr + YpFPm;
	lvz[0]  = lvz[0]*Ztr + Ztm;
	
	
	// ------- Path length reconstruction (RHRS) ------ //
	XFP_R   = (XFP_R -XFPm)/XFPr;
	XpFP_R  = (XpFP_R-XpFPm)/XpFPr;
	YFP_R   = (YFP_R -YFPm)/YFPr;
	YpFP_R  = (YpFP_R-YpFPm)/YpFPr;
	
	LenR    = calcf2t_plen(Plen,XFP_R,XpFP_R,YFP_R,YpFP_R);
	rvz[0]  = calcf2t_plen(Pzt_R,XFP_R,XpFP_R,YFP_R,YpFP_R);
	
	LenR    = LenR*PLr+PLm;
	XFP_R   = XFP_R*XFPr + XFPm;
	XpFP_R  = XpFP_R*XpFPr + XpFPm;
	YFP_R   = YFP_R*YFPr + YFPm;
	YpFP_R  = YpFP_R*YpFPr + YpFPm;
	rvz[0]  = rvz[0]*Ztr + Ztm;
	
	
	//cout << LenL << " " << LenR << endl;
	
	rast_x2 = rast_curx * rastx_param[1] + rastx_param[0];
	double cor_rast = rast_x2/tan(13.2/180.*3.14159);
	double rvz_cor = rvz[0] - cor_rast;
	double lvz_cor = lvz[0] + cor_rast;
	vz_mean[0] = (rvz_cor + lvz_cor)/2.0;
	
	// ----- Left arm angle ----- //
	XFP_L   = (XFP_L -XFPm)/XFPr;
	XpFP_L  = (XpFP_L-XpFPm)/XpFPr;
	YFP_L   = (YFP_L -YFPm)/YFPr;
	YpFP_L  = (YpFP_L-YpFPm)/YpFPr;
	double vzt = (vz_mean[0] - Ztm)/Ztr;
	th2[0]  = calcf2t_4th_2(Pxp_L, XFP_R,XpFP_R,YFP_R,YpFP_R,vzt);
	ph2[0]  = calcf2t_4th_2(Pyp_L, XFP_R,XpFP_R,YFP_R,YpFP_R,vzt);
	XFP_L   = XFP_L*XFPr + XFPm;
	XpFP_L  = XpFP_L*XpFPr + XpFPm;
	YFP_L   = YFP_L*YFPr + YFPm;
	YpFP_L  = YpFP_L*YpFPr + YpFPm;
	
	
	double beta_L = mom2[0]/sqrt(pow(mom2[0],2.0)+pow(me,2.0));
	double cor_L   = (LenL-3.18)/3.0e+8/beta_L * 1.0e+9; // (3.18 m; test)
	double beta_R = mom1[0]/sqrt(pow(mom1[0],2.0)+pow(mpi,2.0));
	double cor_R = (LenR-rpathl_s2[0])/3.0e+8/beta_R * 1.0e+9;
	
	double meantime_L = tref_L - (timeL_L+timeR_L)/2.0 + toffset_L + cor_L;
	double meantime_R = tref_R - (timeL_R+timeR_R)/2.0 + toffset_R + cor_R;
	
	// --- T0 correction --- 
	meantime_R = meantime_R - s2_tzero_R[seg_R] - s2_tzero_L[seg_L]; 
	
	double yfp_cor_R   = YFP_R * par_rtime_ycor[0] + YpFP_R * par_rtime_ycor[1];
	valval[0] = XFP_L;
	valval[1] = XpFP_L;
	valval[2] = YFP_L;
	valval[3] = YpFP_L;
	double yfp_cor_L = Calc_FPcor(valval,par_rtime_ycor_L);
	meantime_R = meantime_R + yfp_cor_R + yfp_cor_L;
	meantime_R = meantime_R-cor_L+75.4;
	ctime[0] = -meantime_R;
	
	if(tflag==5){
	  if(fabs(ctime[0])<20.0){
	    // ---- Filling data ------ //
	    tnew->Fill(); // ---------- //
	    //------------------------- //
	  }
	}
	else{
	  tnew->Fill();
	}
      }
      
    }
    
  }
  tnew->Write();
  fnew->Close();
  
  

  return 0;
}


//////////////////////////////////////////////////
double calcf2t_plen(double* P, double xf, double xpf, 
                 double yf, double ypf)
//////////////////////////////////////////////////
{
  // -----3rd order -----
  const int nMatT=3; 
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



double Calc_FPcor(double* val, double* par){
  double x  = val[0];
  double xp = val[1];
  double y  = val[2];
  double yp = val[3];
  double cor = 0.0; 
  double cor1=0.0, cor2=0.0, cor3=0.0;
  
  cor1 = par[0]*y + par[1]*yp;
  //cor1 = par[0]*x + par[1]*xp + par[2]*y + par[3]*yp;
  //cor2 = par[4]*x*xp + par[5]*x*y + par[6]*x*yp + par[7]*xp*y + par[8]*xp*yp + par[9]*y*yp;
  //cor3 = par[10]*x*x + par[11]*xp*xp + par[12]*y*y + par[13]*yp*yp;
  
  cor = cor1+cor2+cor3;
  return cor;
  
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


