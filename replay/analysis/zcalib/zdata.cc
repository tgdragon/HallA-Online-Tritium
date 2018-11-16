/*
  zdata.cc
  
  Toshiyuki Gogami, Nov 13, 2018
*/

extern double calcf2t_zt(double* P, 
			 double xf, double xpf,
			 double yf, double ypf);

const double  XFPm=-0.7, XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18; 
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74;
const double  XFPr=1.3, XpFPr=0.27; 
const double  YFPr=0.1, YpFPr=0.10; 
const double  Xptr=0.15, Yptr=0.08, Momr=0.18;
const double  PLm = 25.4, PLr=0.7;
const double  Ztm = -0.15, Ztr=0.35;
extern void fcn(int &nPar, double* /*grad*/, 
		double &fval, double* param, int /*iflag*/);
extern double tune(double* pa, int j);

const int nfoil = 10;
//double fcent[nfoil] = {-0.120, -0.096, -0.070, -0.043, -0.019,
//		       0.005, 0.031, 0.056, 0.109, 0.135};
double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125};
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};
//double selection_width = 0.01;
double selection_width = 0.012;
const int nParamT = 35;
const int nmax = 5000;
double x[nmax], y[nmax];
double xp[nmax], yp[nmax];
double z_recon[nmax];
int foil_flag[nmax];
int ntune_event = 0;
double Pzt_opt[nmax];


void zdata(){
  
  // ======================================== //
  // ======= Opening a ROOT file ============ //
  // ======================================== //
  TFile* f1 = new TFile("../../Rootfiles/nnL/tritium_111132.root");
  TTree* t1 = (TTree*)f1->Get("T");
  
  //TFile* f1 = new TFile("coin_111170_111190.root");
  //TFile* f1 = new TFile("coin_111480_111507.root");
  //TFile* f1 = new TFile("coin_H2_upto-111507.root");
  //TFile* f1 = new TFile("coin_H2_upto-111524.root");
  //TFile* f1 = new TFile("coin_H2_upto-111542.root");
  //TFile* f1 = new TFile("coin_H2_1.root");
  //TTree* t1 = (TTree*)f1->Get("tree");
  Double_t trig5;
  Double_t trig4;
  Double_t trig1;
  double ent = t1->GetEntries();
  ent = 200000; // for test
  int evshift = 30000;
    
  const int max = 100;
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
  double ctime[max];
  const double hrs_ang = 13.2 * 3.14159 / 180.;
  
  t1->SetBranchAddress("fEvtHdr.fRun", &runnum    );
  t1->SetBranchAddress("HALLA_p", &hallap );
  t1->SetBranchAddress("DR.T1", &trig1    );
  t1->SetBranchAddress("DR.T4", &trig4    );
  t1->SetBranchAddress("DR.T5", &trig5    );
  //t1->SetBranchAddress("R.tr.time", &rtime);
  //t1->SetBranchAddress("L.tr.time", &ltime);
  //t1->SetBranchAddress("R.tr.pathl", &rpathl);
  //t1->SetBranchAddress("L.tr.pathl", &lpathl);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  t1->SetBranchAddress("R.tr.p", &mom1);
  t1->SetBranchAddress("L.tr.p", &mom2);
  //t1->SetBranchAddress("RTDC.F1FirstHit", &rf1tdc);
  //t1->SetBranchAddress("LTDC.F1FirstHit", &lf1tdc);
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
  //t1->SetBranchAddress("R.s2.t_pads", &r_s2_t_pads);
  //t1->SetBranchAddress("L.s2.t_pads", &l_s2_t_pads);
  //t1->SetBranchAddress("R.s2.nthit",   &r_s2_nthit);
  //t1->SetBranchAddress("L.s2.nthit",   &l_s2_nthit);
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
  t1->SetBranchAddress("R.tr.beta",  &rbeta);
  t1->SetBranchAddress("L.tr.beta",  &lbeta);
  //t1->SetBranchAddress("R.s2.trpath",  &rpathl_s2);
  //t1->SetBranchAddress("L.s2.trpath",  &lpathl_s2);
  //t1->SetBranchAddress("L.s2.nthit",&nhit);
  //t1->SetBranchAddress("R.s2.nthit",&nhit_R);
  t1->SetBranchAddress("R.ps.asum_c", &ps_asum);
  //t1->SetBranchAddress("ctime", &ctime);
  t1->SetBranchAddress("R.a1.t_fadc", &a1_tdc);
  t1->SetBranchAddress("R.a2.t_fadc", &a2_tdc);
  
  TFile* fnew = new TFile("zt_rhrs.root","recreate");
  TTree* tnew = new TTree("tree","For z calibration (RHRS)");
  double ztR[max];
  tnew->Branch("fEvtHdr.fRun", &runnum,"fEvtHdr.fRun/D"    );
  tnew->Branch("HALLA_p", &hallap,"HALLA_p/D" );
  tnew->Branch("DR.T1", &trig1, "DR.T1/D"    );
  tnew->Branch("DR.T4", &trig4, "DR.T4/D"   );
  tnew->Branch("DR.T5", &trig5,  "DR.T5/D"  );
  tnew->Branch("R.tr.vz", &rvz, "R.tr.vz[100]/D");
  tnew->Branch("L.tr.vz", &lvz, "L.tr.vz[100]/D");
  tnew->Branch("R.tr.x",   &r_x_fp, "R.tr.x[100]/D"  );
  tnew->Branch("L.tr.x",   &l_x_fp, "L.tr.x[100]/D"  );
  tnew->Branch("R.tr.y",   &r_y_fp, "R.tr.y[100]/D"  );
  tnew->Branch("L.tr.y",   &l_y_fp, "L.tr.y[100]/D"  );
  tnew->Branch("R.tr.th",  &r_th_fp,"R.tr.th[100]/D" );
  tnew->Branch("L.tr.th",  &l_th_fp,"L.tr.th[100]/D" );
  tnew->Branch("R.tr.ph",  &r_ph_fp,"R.tr.ph[100]/D" );
  tnew->Branch("L.tr.ph",  &l_ph_fp,"L.tr.ph[100]/D" );
  tnew->Branch("R.tr.vz_TG",  &ztR, "R.tr.vz_TG[100]/D" );
  tnew->Branch("R.a1.asum_c", &a1, "R.a1.asum_c/D");
  tnew->Branch("R.a2.asum_c", &a2, "R.a2.asum_c/D");
  //tnew->Branch("L.tr.vz_TG",  &ztR,   "L.tr.vz_TG[100]/D" );
  

  double XFP, XpFP;
  double YFP, YpFP;
  
  ntune_event = 0;
  for(int i=0 ; i<nmax ; i++){
    Pzt_opt[i] = -2222.0;
  }

  char name_Mzt[500];
  //sprintf(name_Mzt,"matrices/zt_RHRS_2.dat"); // original
  sprintf(name_Mzt,"matrices/zt_RHRS_opt.dat"); // optimized parameter
  //sprintf(name_Mzt,"tune1/newpar_9.dat");
  //sprintf(name_Mzt,"tune2/newpar_7.dat");
  //sprintf(name_Mzt,"tune3/newpar_9.dat");
  //sprintf(name_Mzt,"tune4/newpar_9.dat");
  //sprintf(name_Mzt,"tune5/newpar_9.dat");
  //sprintf(name_Mzt,"tune6/newpar_19.dat");
  //sprintf(name_Mzt,"tune7/newpar_19.dat");
  //sprintf(name_Mzt,"tune8/newpar_19.dat");
  //sprintf(name_Mzt,"tune9/newpar_9.dat");
  //sprintf(name_Mzt,"tune10/newpar_9.dat");
  ifstream Mzt(name_Mzt);
  double Pzt[nParamT];
  ent = 100000;
  //double Plenopt[nParamT];
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p; 
    Pzt[i]=par;
    Pzt_opt[i]=Pzt[i];
    //cout << Pzt[i] << endl;
  }
  Mzt.close();
  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,"matrices/zt_LHRS_opt.dat"); // optimized
  ifstream Mzt_L(name_Mzt_L);
  double Pzt_L[nParamT];
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p; 
    Pzt_L[i]=par;
  }
  Mzt_L.close();
  
  
  TH1F* h1 = new TH1F("h1","",400,-0.5,0.5);
  h1->GetXaxis()->SetTitle("z-RHRS (m)");
  h1->GetXaxis()->SetRangeUser(-0.2,0.2);
  
  TH1F* h2[nfoil];
  char tempc[500];
  for(int i=0 ; i<nfoil ; i++){
    sprintf(tempc,"h2_%d",i);
    h2[i] = new TH1F(tempc,tempc,
		     h1->GetXaxis()->GetNbins(),
		     h1->GetXaxis()->GetXmin(),
		     h1->GetXaxis()->GetXmax());
    h2[i]->GetXaxis()->SetTitle("z-RHRS (m)");
    h2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
  }
  
  //double ztR;
  bool rtrig = false;
  bool ltrig = false;
  
  
  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0;
    y[i]    = -2222.0;
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    foil_flag[i] = -1;
  }
  
  
  for (int i=0 ; i< ent ; i++){
    for(int j=0 ; j<max ; j++){
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0;
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
    }
    trig1 = 0.0;
    trig4 = 0.0;
    trig5 = 0.0;
    rtrig = false;
    ltrig = false;
    
    
    //t1->GetEntry(i);
    if(i+evshift<ent) t1->GetEntry(i+evshift);
    else t1->GetEntry(i-evshift);
    
    if(trig4>1.0) rtrig = true;
    else rtrig = false;
    if(trig1>1.0) ltrig = true;
    else ltrig = false;
   
    XFP   = r_x_fp[0];
    XpFP  = r_th_fp[0];
    YFP   = r_y_fp[0];
    YpFP  = r_ph_fp[0];
    
    if(rtrig==true
       && fabs(XFP)  <2.0
       && fabs(XpFP) <0.1
       && fabs(YFP)  <0.5
       && fabs(YpFP) <0.1
       ){
      
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      ztR[0] = calcf2t_zt(Pzt, XFP, XpFP, YFP, YpFP);
      
      ztR[0] = ztR[0] * Ztr + Ztm;
      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;
       
      h1->Fill(ztR[0]);
      
      tnew->Fill();
      
      for(int j=0 ; j<nfoil ; j++){
	if(fcent[j]-selection_width<ztR[0] 
	   && ztR[0]< fcent[j]+selection_width){
	  h2[j]->Fill(ztR[0]);
	  if(j+2!=10) h2[j]->SetLineColor(j+2);
	  else h2[j]->SetLineColor(2);
	  h2[j]->SetLineStyle(9);
	  if(ntune_event<nmax){
	    foil_flag[ntune_event] = j;
	    x[ntune_event]  = XFP;
	    y[ntune_event]  = YFP;
	    xp[ntune_event] = XpFP;
	    yp[ntune_event] = YpFP;
	    z_recon[ntune_event] = ztR[0];
	    ntune_event++;
	  }
	}
      }
      //if(ntune_event<nmax){
      //	x[ntune_event]  = XFP;
      //	y[ntune_event]  = YFP;
      //	xp[ntune_event] = XpFP;
      //	yp[ntune_event] = YpFP;
      //	z_recon[ntune_event] = ztR;
      //	ntune_event++;
      //}
    }
  }
  tnew->Write();
  //fnew->Close();
  
  //for(int i=0 ; i<ntune_event ; i++){
  //  cout << " " << foil_flag[i] << " " << z_recon[i] << endl;
  //}
  

  // =================================== //
  // ======== Draw histograms ========== //
  // =================================== //
  TCanvas* c1 = new TCanvas("c1","c1");
  h1->Draw();
  for(int i=0 ; i<nfoil ; i++){
    h2[i]->Draw("same");
  }
  
  
  
  //const int nite = 10;
  const int nite = 0; // no tuning
  
  if (nite>0) cout << " Tuning started: " << endl;
  for(int i=0 ; i<nite ; i++){
    //sprintf(tempc, "./tune1/newpar_%d.dat",i);
    //sprintf(tempc, "./tune2/newpar_%d.dat",i);
    //sprintf(tempc, "./tune3/newpar_%d.dat",i);
    //sprintf(tempc, "./tune4/newpar_%d.dat",i);
    //sprintf(tempc, "./tune5/newpar_%d.dat",i);
    //sprintf(tempc, "./tune6/newpar_%d.dat",i);
    //sprintf(tempc, "./tune7/newpar_%d.dat",i);
    //sprintf(tempc, "./tune8/newpar_%d.dat",i);
    //sprintf(tempc, "./tune9/newpar_%d.dat",i);
    //sprintf(tempc, "./tune10/newpar_%d.dat",i);
    sprintf(tempc, "./testmat/newpar_%d.dat",i);
    ofstream * ofs = new ofstream(tempc);
    int nppp = 0;
    const int nn = 3; // 3rd order matrix
    for(int i=0 ; i<nn+1 ; i++){
      for(int d=0 ; d<nn+1 ; d++){
	for(int c=0 ; c<nn+1 ; c++){
	  for(int b=0 ; b<nn+1 ; b++){
	    for(int a=0 ; a<nn+1 ; a++){  
	      if(a+b+c+d==i){
		*ofs << Pzt_opt[nppp] 
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
    ofs->clear();
    
    
    // ---- Parameter tuning ----- //
    tune(Pzt_opt,i);
    //tune(Pzt,i);
  }


  
//  ofstream * ofs = new ofstream("newpar.dat");
//  int nppp = 0;
//  const int nn = 3; // 3rd order matrix
//  for(int i=0 ; i<nn+1 ; i++){
//    for(int d=0 ; d<nn+1 ; d++){
//      for(int c=0 ; c<nn+1 ; c++){
//	for(int b=0 ; b<nn+1 ; b++){
//	  for(int a=0 ; a<nn+1 ; a++){  
//	    if(a+b+c+d==i){
//	      *ofs << Pzt_opt[nppp] 
//		   << " " << a 
//		   << " " << b
//		   << " " << c
//		   << " " << d << endl;
//	      nppp++;
////	      cout << Plen_opt[nppp] 
////		   << " " << a 
////		   << " " << b
////		   << " " << c
////		   << " " << d << endl;
//	    }
//	  }
//	}
//      }
//    }
//  }
//  ofs->close();


  
  
  
  
  
}


//////////////////////////////////////////////////
double calcf2t_zt(double* P, double xf, double xpf, 
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
    minuit -> GetParameter(i,Pzt_opt[i],e);
  }
  
  return chi2;
}

// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  double chi2 = 0.0;
  
  double XFP, XpFP;
  double YFP, YpFP;
  const double sigma = 0.0045;
  
  double ztR;
  double refz = 0.0;
  double residual = 0.0;
  
  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    refz = 0.0;
    ztR = 0.0;
    
    XFP   = x[i];
    XpFP  = xp[i];
    YFP   = y[i];
    YpFP  = yp[i];
    refz  = fcent_real[foil_flag[i]];
    
    //cout << XFP << " " << YFP << " " << XpFP << " " << YpFP << endl;
    XFP   =(XFP -XFPm)/XFPr;
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    ztR = calcf2t_zt(param, XFP, XpFP, YFP, YpFP);
    
    ztR = ztR * Ztr + Ztm;

    residual = ztR-refz;
    //cout << residual*100 << " cm" << endl;
    chi2 = chi2 + pow(residual,2.0);
    
    
  }

  chi2 = sqrt(chi2)/(double)ntune_event/sigma;

  //chi2_2 = chi2_2*w2/pow(sigma_time ,2.0);
  //chi2_3 = chi2_3*w3/pow(sigma_time ,2.0);
  //chi2_4 = chi2_4*w4/pow(sigma_time ,2.0);
  //chi2_tof = chi2_tof*wtof/pow(sigma_tof ,2.0);
  //
  //chi2 = chi2_1 + chi2_2 + chi2_3 + chi2_4 + chi2_tof;
  //chi2 = chi2 / (w1+w2+w3+w4+wtof);

  fval = chi2;
}
