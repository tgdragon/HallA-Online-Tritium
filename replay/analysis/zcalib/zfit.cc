/*
  zfit.cc
  Fitting to z at the target 
  for the multifoil target data
  
  Toshiyuki Gogami, Oct 11, 2018
*/

double teng(double* x, double* par){
  double g = 0.0;
  double bg = 0.0;
  //double xstart = -12.5;
  for(int i=0 ; i<10 ; i++){
    g = g + par[i*3]*exp(-0.5*pow((x[0]-par[(i*3)+1])/par[(i*3)+2],2.0));
  }
  bg = par[30];
  return g+bg;
}

void zfit(int runnum = 93821){
  // ================================== //
  // ======== Opening ROOT file ======= //
  // ================================== //
  char fname[500];
  sprintf(fname,"./root/tritium_%d.root",runnum);
  TFile* f1 = new TFile(fname);
  TTree* t1 = (TTree*)f1->Get("T");

  const int zmax = 20;
  double vz[zmax];
  double t5;
  int nvz;
  t1->SetBranchAddress("R.tr.vz",&vz);
  t1->SetBranchAddress("Ndata.R.tr.vz",&nvz);
  t1->SetBranchAddress("DR.T5",&t5);
  double ent = t1->GetEntries();


  // ================================== //
  // ======== Empty histogram ========= //
  // ================================== //
  //TH1F* h1 = new TH1F("h1","",250,-40.,40.0);
  TH1F* h1 = new TH1F("h1","",400,-40.,40.0);
  
  for(int i=0 ; i<ent ; i++){
    for(int j=0 ; j<zmax ; j++){
      vz[j] = -2222.0;
    }
    t5 = -2222.0;
    
    t1->GetEntry(i);
    for(int j=0 ; j<nvz ; j++){
      vz[j] = vz[j] * 100.0; //m --> cm
    }

    if( nvz==1 // Only single hit event
        && t5 >1630.0 ){  // T5 (S0&&S2) && GC
      h1->Fill(vz[0]);
    }
  }

  // ===================================== //
  // ======= Drawing the histogram ======= //
  // ===================================== //
  TCanvas* c1 = new TCanvas("c1","c1");
  h1->GetXaxis()->SetRangeUser(-20.0,20.0);
  h1->Draw();

  
  // ============================ //
  // ======== Fitting =========== //
  // ============================ //
  const int nfoil = 10;  // The number of foils
  const int npar = 31;   // The number of paramters of fitting function
  double temp;
  int tempi;
  TF1* func1 = new TF1("func1",teng,-20.0,20.0,npar);
  double par[npar];
  double par_er[npar];
  func1->SetNpx(3000);
  func1->SetLineColor(9);
  for(int i=0 ; i<nfoil ; i++){
    temp = -12.5+2.5*i;
    if(i>7) temp = -12.5+2.5*i+2.5; // no 9th foil
    func1->SetParameter(0+(i*3), 30.0);
    func1->SetParameter(1+(i*3), temp);
    func1->SetParameter(2+(i*3), 0.5);
    func1->SetParLimits(2+(i*3),0.2,1.5);
    
  }
  func1->SetParameter(30, 4.0);
  //func1->FixParameter(30, 5.0);

  h1->Fit("func1","BN","",-20.0,20.0);
  func1->Draw("same");
  for(int i=0 ; i<npar ; i++){
    par[i]    = func1->GetParameter(i);
    par_er[i] = func1->GetParError(i);
  }
  //func1->Draw();

  // ============================= //
  // ====== Event selection ====== //
  // ============================= //
  double ssigma = 1.5;
  double width[nfoil];
  double width_er[nfoil];
  double mean[nfoil];
  double mean_er[nfoil];
  double meandiff[nfoil];
  double meandiff_er[nfoil];
  double x[nfoil];
  double x_er[nfoil];
  const double sigma2fwhm = 2.0*sqrt(2.0*log(2.0));
  cout << endl;
  cout << " " << "1.5 sigma event selection:" << endl;
  for(int i=0 ; i<nfoil ; i++){
    mean[i]      = par[(i*3)+1];
    mean_er[i]   = par_er[(i*3)+1];
    //width[i]     = par[(i*3)+2];
    //width_er[i]  = par_er[(i*3)+2];
    width[i]     = par[(i*3)+2] * sigma2fwhm;
    width_er[i]  = par_er[(i*3)+2] * sigma2fwhm;
    if(i>7) meandiff[i]  = mean[i] -(- 12.5 + i*2.5+2.5);
    else meandiff[i]  = mean[i] -(- 12.5 + i*2.5);
    meandiff_er[i]  = mean_er[i];
    x[i] = i;
    x_er[i]     = 0.0;
    //cout << i << " " << par[(i*3)+1] - par[(i*3)+2]*ssigna << " " << par[(i*3)+1] + par[(i*3)+2]*ssigma << endl;
    cout << " " << i << " " << mean[i] - width[i]*ssigma << " " << mean[i] + width[i]*ssigma << endl;
  }
  cout << endl;
  //cout << sigma2fwhm << endl; // for check

  // ==================~============================== //
  // ====== Creating graphs for resutls ============== //
  // ==================~============================== //
  TGraphErrors* gr1 = new TGraphErrors(nfoil,x,meandiff,x_er,meandiff_er);
  gr1->SetName("gr1");
  TGraphErrors* gr2 = new TGraphErrors(nfoil,x,width,x_er,width_er);
  gr2->SetName("gr2");

  TH2F* hemp1 = new TH2F("hemp1","",100,-1.0,11.0,100,-1.0,1.0);
  hemp1->GetXaxis()->SetTitle("Foil #");
  hemp1->GetYaxis()->SetTitle("z^{recon} - z^{0} (cm)");
  TH2F* hemp2 = new TH2F("hemp2","",100,-1.0,11.0,100,0.6,1.3);
  hemp2->GetXaxis()->SetTitle("Foil #");
  hemp2->GetYaxis()->SetTitle("FWHM (cm)");
  
  TCanvas* c2 = new TCanvas("c2","c2");
  hemp1->Draw();
  gr1->Draw("same*l");
  TCanvas* c3 = new TCanvas("c3","c3");
  hemp2->Draw();
  gr2->Draw("same*l");
  
  
  
}
