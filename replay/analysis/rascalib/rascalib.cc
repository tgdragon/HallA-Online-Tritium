/*
  rascalib.cc
  This is a macro to optimize raster parameters.
  
  Toshiyuki Gogami, Jan 25, 2019
*/

extern double calcf2t_zt(double* P, 
			 double xf, double xpf,
			 double yf, double ypf);
extern double tune(double* pa, int j);
extern void fcn(int &nPar, double* /*grad*/, 
		double &fval, double* param, int /*iflag*/);
extern double calcRasterCor(double,double,double);

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
double selection_width = 0.008; 

const int nParamT = 35; // Number of parameters for the 3rd order matrix (z recon)
const int nmax = 100000; // Number of events used for tuning
double x[nmax], y[nmax];
double xp[nmax], yp[nmax];
double z_recon[nmax]; // reconstructed z position
int foil_flag[nmax];  
int ntune_event = 0;  
const int nParamT2 = 4; 
double parRaster[nParamT2]; 
double Opt_Par[nParamT2];
double Ras_curx[nmax];
double Ras_cury[nmax];

const double hrs_ang = 13.2 * 3.14159 / 180.; 

void rascalib(){
  
  const int nite = 0;  // The number of tuning iteration
  
  // ======== General conditions ======= //
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  // ======================================== //
  // ======= Opening a ROOT file ============ //
  // ======================================== //
  TFile* f1 = new TFile("./tritium_111721.root");
  TTree* t1 = (TTree*)f1->Get("T");
  Double_t trig1;
  double ent = t1->GetEntries();
  if(nite>0) ent = 200000;
  //ent = 50000; // for test
  //ent = 100000; // for test
  cout << endl;
  cout << " The number of events read from a ROOT file: " << ent << endl;
  int evshift = 30000;
  evshift = 30000; 
  
  const int max = 100;
  double ltime_s0[max]; // lhrs s0 time
  double ltime_s2[max];
  // double ltime[max];// not in use
  // double lpathl[max];// not in use
  // double lpathl_s2[max];// probably not in use
  double mom2[max]; // mom2 = L.tr.p
  const int f1n = 64;
  
  double lvz[max]; // R.tr.vz ...vertex z in (m)
  double th2[max], ph2[max]; //L.tr.tg_th
  Int_t runnum; 
  double hallap;
  
  double l_th_fp[max];// L.tr.th  which is the tangent track theta angle
  double l_ph_fp[max];// R.tr.ph .. which is the tangent track phi angle
  double l_x_fp[max];// L.tr.x track x coordinate in (m)
  double l_y_fp[max];
  const int n = 16;
  
  double lbeta[max];
  double ctime[max]; // not in use
  double Lrb; // For raster correction.. x position of the raster
  double rpr; // z position of the raster
  
  
  double Ras_x; //  for raster ====================
  double Ras_y; // ================================
  double cer_asum; // cernkkov sum to remove the pions 
  
  t1->SetBranchAddress("fEvtHdr.fRun", &runnum   );
  t1->SetBranchAddress("HALLA_p", &hallap   ); // beam momentum 
  t1->SetBranchAddress("DR.T1", &trig1  );     // LHRS Trigger
  t1->SetBranchAddress("L.tr.p", &mom2);
  t1->SetBranchAddress("Lrb.x", &Lrb);         // raster x
  t1->SetBranchAddress("rpr.z", &rpr);         
  
  t1->SetBranchAddress("Lrb.Raster2.rawcur.x", &Ras_x); // raster current
  t1->SetBranchAddress("Lrb.Raster2.rawcur.y", &Ras_y); // raster current
  
  t1->SetBranchAddress("L.tr.vz", &lvz);      // z obtained by reply
  t1->SetBranchAddress("L.tr.tg_th", &th2);   // x' at target
  t1->SetBranchAddress("L.tr.tg_ph", &ph2);   // y' at target
  t1->SetBranchAddress("L.s0.time", &ltime_s0);
  t1->SetBranchAddress("L.s2.time", &ltime_s2);
  
  
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
  t1->SetBranchAddress("L.tr.beta",  &lbeta);
  t1->SetBranchAddress("L.cer.asum_c",  &cer_asum); // Gas cherenkov
  
  TFile* fnew = new TFile("./Rootfiles/check.root","recreate");
  TTree* tnew = new TTree("tree","check for raster parameter optimnization (LHRS)");
  double ztR[max];      // Reconstructed z position using a 3rd order matrix
  double ztR_wRC[max];  // z position with a raster correction
  
  tnew->Branch("fEvtHdr.fRun", &runnum,"fEvtHdr.fRun/D"    );
  tnew->Branch("HALLA_p", &hallap,"HALLA_p/D" );
  tnew->Branch("DR.T1", &trig1, "DR.T1/D"    );
  tnew->Branch("L.tr.vz", &lvz, "L.tr.vz[100]/D");
  tnew->Branch("L.tr.x",   &l_x_fp, "L.tr.x[100]/D"  );
  tnew->Branch("L.tr.y",   &l_y_fp, "L.tr.y[100]/D"  );
  tnew->Branch("L.tr.th",  &l_th_fp,"L.tr.th[100]/D" );
  tnew->Branch("L.tr.ph",  &l_ph_fp,"L.tr.ph[100]/D" );
  tnew->Branch("L.tr.vz_TG",    &ztR,    "L.tr.vz_TG[100]/D" );
  tnew->Branch("L.tr.vz_TG2",  &ztR_wRC, "L.tr.vz_TG2[100]/D" );
  tnew->Branch("Lrb.x", &Lrb, "Lrb.x/D"    ); // raster x
  tnew->Branch("rpr.z", &rpr, "rpr.z/D"    ); 
  tnew->Branch("L.cer.asum_c", &cer_asum, "L.cer.asum_c/D"); 
  
  tnew->Branch("Lrb.Raster2.rawcur.x", &Ras_x, "Lrb.Raster2.rawcur.x/D"    ); // raster current
  tnew->Branch("Lrb.Raster2.rawcur.y", &Ras_y, "Lrb.Raster2.rawcur.y/D"    ); // raster current
  
  double XFP, XpFP; // focal plane information
  double YFP, YpFP;
  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,"../matrices/zt_LHRS_opt.dat"); // An optimized z matrix for LHRS
  ifstream Mzt_L(name_Mzt_L);
  double Pzt_L[nParamT];
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    Pzt_L[i]=par;
  }
  Mzt_L.close();
  
  ntune_event = 0;
  
  // ------ Initialization of raster parameters------- //
  for(int j=0 ; j<nParamT2 ; j++){ //
    parRaster[j] = -2222.0;   // Input parameters
    Opt_Par[j]   = -2222.0;   // Parameters after tuning
  }
  
  char name_Raster_L[500];
  //sprintf(name_Raster_L,"./data/raster_ini.dat"); // initial raster parameters
  sprintf(name_Raster_L,"./data/newpar.dat"); // initial raster parameters
  ifstream Raster_L(name_Raster_L);
  for(int j=0;j<nParamT2;j++){
    Raster_L >> parRaster[j];
    
    // ---- raster x is tuned by this macro ---- //
    if(j==1 || j==3) parRaster[j] = 0.0;
  }
  Raster_L.close();
  
  // --- Check input parameters by cout --- //
  cout << " Input raster parameters (only x is used here): "
       << parRaster[0] << " " << parRaster[1] << " " 
       << parRaster[2] << " " << parRaster[3] << endl; 
  
  //TH1F* h1 = new TH1F("h1","Multi-C-foil target(After first tune,Raster off )",400,-0.5,0.5); // -0.5, 0.5
  TH1F* h1 = new TH1F("h1","",400,-0.5,0.5); 
  h1->GetXaxis()->SetTitle("z-LHRS (m)");
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetRangeUser(-0.2,0.2);
  h1->GetYaxis()->SetTitle("number of counts/2.5mm");
  h1->GetYaxis()->CenterTitle();
  TH1F* h2[nfoil];
  char tempc[500];
  for(int i=0 ; i<nfoil ; i++){
    sprintf(tempc,"h2_%d",i);
    h2[i] = new TH1F(tempc, tempc,
		     h1->GetXaxis()->GetNbins(),
		     h1->GetXaxis()->GetXmin(),
		     h1->GetXaxis()->GetXmax());
    h2[i]->GetXaxis()->SetTitle("z-LHRS (m)");
    h2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
  }
  
  bool ltrig = false;
  
  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0; // x at FP
    y[i]    = -2222.0; // y at FP
    xp[i]   = -2222.0; // x' at FP
    yp[i]   = -2222.0; // y' at FP
    foil_flag[i] = -1; // group number of foil
  }
  
  for (int i=0 ; i< ent ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0;
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      ztR[j]     = -2222.0;
      ztR_wRC[j] = -2222.0;
      z_recon[j] = -2222.0;
    }
    
    trig1 = 0.0;
    ltrig = false;
    
    //t1->GetEntry(i);
    if(i+evshift<ent) t1->GetEntry(i+evshift); 
    else t1->GetEntry(i-ent+evshift);
    if(trig1>1.0) ltrig = true;
    else ltrig = false;
    
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    
    if(ltrig==true
       && fabs(XFP)  <2.0 
       && fabs(XpFP) <0.1
       && fabs(YFP)  <0.5
       && fabs(YpFP) <0.1
       && cer_asum>1800   // pi^{-} rejection with gas Cherenkov 
       ){
      XFP  = (XFP-XFPm)/XFPr; 
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      ztR[0] = calcf2t_zt(Pzt_L, XFP, XpFP, YFP, YpFP); // calculated z position
      
      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;
      ztR[0] = ztR[0] * Ztr + Ztm; 
      
      // ----------------------------------------------------------- //
      // ------ Converting raster current to raster position ------- //
      // ----------------------------------------------------------- //
      double RasterCor = calcRasterCor(Ras_x, parRaster[2], parRaster[0]);
      RasterCor = RasterCor/tan(hrs_ang);
      ztR_wRC[0] = ztR[0]+RasterCor;
      h1->Fill(ztR[0]+RasterCor);
      tnew->Fill();
      
      // ------------------------------------------- //
      // ------- Event selection for tuning -------- //
      // ------------------------------------------- //
      if(nite>0){
	for(int j=0 ; j<nfoil ; j++){
	  if(fcent[j]-selection_width<ztR[0]+RasterCor 
	     && ztR[0]+RasterCor<fcent[j]+selection_width){
	    h2[j]->Fill(ztR[0]+RasterCor);
	    if(j+2!=10) h2[j]->SetLineColor(j+2);
	    else h2[j]->SetLineColor(2);
	    h2[j]->SetLineStyle(9);
	    if(ntune_event<nmax){
	      foil_flag[ntune_event] = j;
	      x[ntune_event]  = XFP;
	      y[ntune_event]  = YFP;
	      xp[ntune_event] = XpFP;
	      yp[ntune_event] = YpFP;
	      z_recon[ntune_event] = ztR[0]; // Reconstructed z position
	      //ztR_wRC[ntune_event] = ztR[0] + RasterCor; // z position with raster correction
	      Ras_curx[ntune_event] = Ras_x;
	      Ras_cury[ntune_event] = Ras_y;
	      ntune_event++;
	    }
	  }
	}
      } // nite>0
    }
  }
  
  tnew->Write();
  //fnew->Close(); 
  
  // =================================== //
  // ======== Draw histograms ========== //
  // =================================== //
  TCanvas* c1 = new TCanvas("c1","c1");
  h1->Draw();
  for(int i=0 ; i<nfoil ; i++){
    h2[i]->Draw("same");
  }
  
  
  double chi_sq[nite]; // A chi square for each tuning
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;
  
  for(int i=0 ; i<nite ; i++){
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    x[i] = i+1;
    if(i==0){ 
      chi_sq[i] = tune(parRaster,i); 
    }
    else { 
      chi_sq[i] = tune(Opt_Par,i); 
    }
    cout << " Tuning# = " << i << ": chisq = " << chi_sq[i] << endl; // 
    cout << endl;
    
    sprintf(tempc, "./data/newpar_%d.dat",i); //new matrix will be stored here
    ofstream * ofs = new ofstream(tempc);
    *ofs << Opt_Par[0] << " " << Opt_Par[1] << " " << Opt_Par[2] << " "<< Opt_Par[3] << endl;
    cout << Opt_Par[0] << " " << Opt_Par[1] << " " << Opt_Par[2] << " "<< Opt_Par[3] << endl;
    
    ofs->close();
    ofs->clear();
  }
  
  if(nite>0){
    TGraph * gr = new TGraph(nite,x,chi_sq);
    gr->SetName("gr");
    TCanvas * c2 = new TCanvas("c2","c2",600,600);
    gr->Draw("*la");
  }
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

// #############################################################
double tune(double* pa, int j) // tune fun defn
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = nParamT2;
  //cout << allparam << endl;
  TMinuit* minuit = new TMinuit(allparam);
  minuit->SetFCN(fcn);
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];// Lower limit for all of the parameter
  double ULim[allparam];// Upper limit for all of the parameter
  char pname[500];
  
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i]; // initial parameters
    //step[i] = 1.0e-5;  
    
    if(i==0){       // offset parameter for raster x
      step[i] = 1.0e-5;
    }
    else if(i==2){  // gradient parameters for raster x
      step[i] = 1.0e-10;
    }
    else step[i] = 0.0; // raster y is not tuned here
    
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
    LLim[i] = pa[i] -1.0; // temporary 
    ULim[i] = pa[i] +1.0; // temporary 
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }
  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double e;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    //minuit -> GetParameter(i,par[i],e);
    minuit -> GetParameter(i,Opt_Par[i],e);
  }
  
  return chi2;
}

// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  
  const double sigma = 0.00035;
  double ztR  = 0.0;
  double refz = 0.0;
  double nev[nfoil];
  double residual;
  double chi2[nfoil];
  double w[nfoil];
  double total_chi2 = 0.0;
  
  for(int i=0 ; i<nfoil ; i++){
    nev[i]  = 0.0;
    chi2[i] = 0.0;
    w[i]    = 1.0;
  }
  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    refz = 0.0;  refz = fcent_real[foil_flag[i]];
    ztR  = 0.0;  ztR  = z_recon[i];

    if(foil_flag[i]==i) nev[i]++;
    
    // ----------------------------- //
    // ----- Raster correction ----- //
    // ----------------------------- //
    double RasterCor = calcRasterCor(Ras_curx[i], param[2], param[0]); 
    RasterCor = RasterCor/tan(hrs_ang);
    double ztR_wRC; 
    ztR_wRC = ztR + RasterCor; // z with raster correction
    
    // ------------------- //
    // --- Residual ------ //
    // ------------------- //
    residual = ztR_wRC - refz; 
    chi2[foil_flag[i]] = chi2[foil_flag[i]] + pow(residual,2.0);
  }
  
  for(int i=0 ; i<nfoil ; i++){
    if(nev[i]!=0){
      chi2[i] = sqrt(chi2[i])/nev[i]/sigma;
    }
    total_chi2 = total_chi2 + chi2[i]*w[i];
  }
  
  fval = total_chi2/(double)nfoil;
}


double calcRasterCor(double a, double b, double c){
  return a * b + c;
}
