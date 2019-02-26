/*
  check2.cc
  This is a macro to optimize raster parameters.
  
  Toshiyuki Gogami, Jan 25, 2019
*/

extern double calcf2t_zt(double* P, 
			 double xf, double xpf,
			 double yf, double ypf);
extern double calcf2t_4th_2(double*,
			    double, double, 
			    double, double,
			    double);
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


const double step = 0.492 * 2.54;
const int nrow = 11;
const int ncol = 7;
const int nsshole = nrow*ncol;
//double ssy_cent[nrow];
//double ssx_cent[ncol];
double ssy_cent_real[nrow];
double ssx_cent_real[ncol];
double refx[nsshole];
double refy[nsshole];
double selec_widthx = 0.60;
double selec_widthy = 0.45;


const int nParamT = 126;  // Number of parameters
const int nmax = 100000; // Number of events used for tuning
double x[nmax], y[nmax];
double xp[nmax], yp[nmax];
double z_recon[nmax]; // reconstructed z position
int foil_flag[nmax];  
int ntune_event = 0;
int holegroup[nmax];
int holeg_temp;

const int nParamT2 = 4; 
double parRaster[nParamT2];
double Opt_Par[nParamT2];
double Ras_curx[nmax];
double Ras_cury[nmax];

const double hrs_ang = 13.2 * 3.14159 / 180.;

//#include<stdlibs>
//using namespace std;



void check2(){

  const int nite = 0;  // The number of tuning iteration
  
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
  TFile* f1 = new TFile("LHRS_SS_111721.root");
  TTree* t1 = (TTree*)f1->Get("tree");
  Double_t trig1;
  double ent = t1->GetEntries();
  if(nite>0) ent = 200000;
  ent = 50000; // for test
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
  
  TFile* fnew = new TFile("check.root","recreate");
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
  

  char name_Mxt_L[500], name_Myt_L[500];
  //sprintf(name_Mxt_L,"../matrices/xpt_LHRS_4.dat");
  //sprintf(name_Myt_L,"../matrices/ypt_LHRS_4.dat");
  sprintf(name_Mxt_L,"./sample_matrix/newpar_xpt_2.dat");
  sprintf(name_Myt_L,"./sample_matrix/newpar_ypt_2.dat");
  ifstream Mxt_L(name_Mxt_L);
  ifstream Myt_L(name_Myt_L);
  double Pxt_L[nParamT], Pyt_L[nParamT];
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxt_L >> par >> p >> p >> p >> p >> p;
    Pxt_L[i]=par;
    cout << par << endl;
    Myt_L >> par >> p >> p >> p >> p >> p;
    Pyt_L[i]=par;
  }
  Mxt_L.close();
  Myt_L.close();

  istream* ifs = new ifstream("scale_offset_20190210.dat");
  double offs_xp[nfoil], offs_yp[nfoil];
  double scal_xp[nfoil], scal_yp[nfoil];
  bool offset_flag[nfoil];
  double temp;
  bool holethrough;
  for (int i=0 ; i<nfoil ; i++){
    *ifs >> temp
	 >> offs_xp[i] >> scal_xp[i]
	 >> offs_yp[i] >> scal_yp[i];
    if(scal_yp[i]>0.5) offset_flag[i] = true;
    else offset_flag[i] = false;
  }
  //ifs->std::close();
  
  ntune_event = 0;
  
  TH2F* h1 = new TH2F("h1","",200,-5.0,5.0,200,-8.0,8.0); 
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
  const double l0 = 100.3;
  double dth[nfoil];
  double l[nfoil];
  double projectf[nfoil];
   
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
    l[i] = sqrt(pow(l0,2.0) + pow(fcent_real[i]*100.,2.0) -2.0*l0*fcent_real[i]*100.*cos(hrs_ang));
    dth[i] = asin(l0/l[i]*sin(hrs_ang)) - hrs_ang;
    projectf[i] = cos( dth[i] );
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
    holethrough = false;
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

    Xpt[0]  = calcf2t_4th_2(Pxt_L,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);
    Ypt[0] = calcf2t_4th_2(Pyt_L,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);

    Ypt[0]  = Ypt[0]*Yptr +Yptm;
    Xpt[0]  = Xpt[0]*Xptr +Xptm;
    Zt[0]   = Zt[0]*Ztr +Ztm;
    
    if(fabs(Xpt[0]) < 0.08 
       //&& fabs(Ypt[0]) < 0.04 ){
       && fabs(Ypt[0]) < 0.06 ){
      
      double ssx, ssy;
      for(int j=0 ; j<nfoil ; j++){
	
	if(fcent[j]-selection_width<Zt[0]
	   && Zt[0]<fcent[j]+selection_width){
	  
	  h2[j]->Fill(-Ypt[0]*l[j]*projectf[j],-Xpt[0]*l[j]*projectf[j]);
	  
	  //if(offset_flag[j]==true){
	  if(offset_flag[j]==true || offset_flag[j]==false){
	    
	    //ssx = (-Xpt[0]*l[j]*projectf[j] + offs_xp[j])*scal_xp[j];
	    //ssy = (-Ypt[0]*l[j]*projectf[j] + offs_yp[j])*scal_yp[j];
	    ssx = -Xpt[0]*l[j]*projectf[j];
	    ssy = -Ypt[0]*l[j]*projectf[j];
	    //if(j==8) ssy = ssy * 1.08;
	    //if(j==9) ssy = ssy * 1.126;
	    
	    h2_new[j]->Fill(ssy,ssx);
	    h2_      ->Fill(ssy,ssx);
	    nhole=0;
	    for(int col=0 ; col<ncol ; col++){
	      for(int row=0 ; row<nrow ; row++){
		if(pow(ssx-refx[nhole],2.0)/pow(selec_widthx,2.0)
		   + pow(ssy-refy[nhole],2.0)/pow(selec_widthy,2.0)<1.0){
		  holeg_temp = nhole;
		  holethrough = true;
		}
		nhole++;
	      }
	    }
	  }
	  
	  if(ntune_event<nmax && holethrough==true){
	    foil_flag[ntune_event] = j;
	    x[ntune_event]  = XFP[0];
	    y[ntune_event]  = YFP[0];
	    xp[ntune_event] = XpFP[0];
	    yp[ntune_event] = YpFP[0];
	    z_recon[ntune_event] = Zt[0];
	    holegroup[ntune_event] = holeg_temp;
	    ntune_event++;

	    h3[j]->Fill(ssy,ssx);
	    h3_  ->Fill(ssy,ssx);
	  }
	}
      }
    }
  }

  cout << " The number of events selected to be used for tuning: "
       << ntune_event << endl;

  TObjArray h(1);
  for(int i=0 ; i<nfoil ; i++){
    h.Add(h2[i]);
    h.Add(h2_new[i]);
    h.Add(h3[i]);
  }
  h.Add(h2_);
  h.Add(h3_);
  h.Write();
  
  
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
    c1->cd(i+1);h2[i]->Draw();
    c2->cd(i+1);h2_new[i]->Draw();
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
  
  
  //double chi_sq[nite]; // A chi square for each tuning
  //double x[nite];
  //if (nite>0) cout << " Tuning started: " << endl;
  
  //for(int i=0 ; i<nite ; i++){
  //  // --------------------------- //
  //  // ---- Parameter tuning ----- //
  //  // --------------------------- //
  //  x[i] = i+1;
  //  if(i==0){ 
  //    chi_sq[i] = tune(parRaster,i); 
  //  }
  //  else { 
  //    chi_sq[i] = tune(Opt_Par,i); 
  //  }
  //  cout << " Tuning# = " << i << ": chisq = " << chi_sq[i] << endl; // 
  //  cout << endl;
  //  
  //  sprintf(tempc, "./data/newpar_%d.dat",i); //new matrix will be stored here
  //  ofstream * ofs = new ofstream(tempc);
  //  *ofs << Opt_Par[0] << " " << Opt_Par[1] << " " << Opt_Par[2] << " "<< Opt_Par[3] << endl;
  //  cout << Opt_Par[0] << " " << Opt_Par[1] << " " << Opt_Par[2] << " "<< Opt_Par[3] << endl;
  //  
  //  ofs->close();
  //  ofs->clear();
  //}
  //
  //if(nite>0){
  //  TGraph * gr = new TGraph(nite,x,chi_sq);
  //  gr->SetName("gr");
  //  TCanvas * c2 = new TCanvas("c2","c2",600,600);
  //  gr->Draw("*la");
  //}
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
  // -----4rd order ----- 
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
