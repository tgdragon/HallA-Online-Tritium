/*
  nLambda.cc
  Counting the nubmer of events of Lambda/Sigma0
  
  Toshiyuki Gogami, Nov 10, 2018
*/


void nLambda(){
  TFile* f1 = new TFile("hist.root");
  TH1F* hh1     = (TH1F*)f1->Get("h1");
  TH1F* h_kaon  = (TH1F*)f1->Get("h1");
  TH1F* hh1_    = (TH1F*)f1->Get("h1_");
  TH1F* hh1_acc = (TH1F*)f1->Get("h1_acc");
  TH1F* hh1_shift = (TH1F*)f1->Get("h1_shift");
  TH1F* hh2       = (TH1F*)f1->Get("h2");
  TH1F* hh2_acc   = (TH1F*)f1->Get("h2_acc");
  TH1F* hh2_noacc = (TH1F*)f1->Get("h2");
  //hh2_noacc->Add(hh2_acc,-1);
  hh2->GetXaxis()->SetRangeUser(-80,150.0);
  
  TCanvas* c1 = new TCanvas("c1","c1");
  hh1_->Draw();
  hh1->Draw("same");
  hh1_acc->Draw("same");
  hh1_shift->Draw("same");
  
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->cd();
  h_kaon->Add(hh1_shift,-1);
  h_kaon->Draw();
  h_kaon->SetLineColor(1);
  double nkaon = h_kaon->GetEntries();
  
  TF1* func1 = new TF1("func1","[0]*exp(-0.5*pow((x-[1])/[2],2.0))",-5,5.3);
  func1->SetParameters(180,0.0,0.4);
  func1->SetNpx(1000);
  func1->SetLineColor(9);
  h_kaon->Fit("func1","NB","",-5.0,5.0);
  double sigma, sigma_er;
  sigma    = func1->GetParameter(2);
  sigma_er = func1->GetParError(2);
  func1->Draw("same");
  
  TCanvas* c3 = new TCanvas("c3","c3");
  int nbins = hh2->GetXaxis()->GetNbins();
  double min = hh2->GetXaxis()->GetXmin();
  double max = hh2->GetXaxis()->GetXmax();
  double dbin = (max-min)/(double)nbins;
  for(int i=0 ; i<nbins ; i++){
    double a    = hh2->GetBinContent(i+1);
    double a_er = sqrt(a);
    double b    = hh2_acc->GetBinContent(i+1);
    double b_er = hh2_acc->GetBinError(i+1);
    hh2_noacc->SetBinContent(i+1,a-b);
    hh2_noacc->SetBinError(i+1,sqrt(a_er*a_er+b_er*b_er));
  }
  c3->cd();
  hh2_noacc->Draw();
  hh2_noacc->SetLineColor(1);
  TF1* func2 = new TF1("func2","[0]*exp(-0.5*pow((x-[1])/[2],2.0))",-30,30.0);
  func2->SetParameters(100,0.0,2.0);
  func2->SetNpx(1000);
  //hh2_noacc->Fit("func2","NB","",-10.0,10.0);
  hh2_noacc->Fit("func2","NB","",-8.0,8.0);
  func2->SetLineColor(9);
  func2->Draw("same");
  
  TF1* func3 = new TF1("func3","[0]*exp(-0.5*pow((x-[1])/[2],2.0))",53,100.0);
  func3->SetParameters(20,77.0,3.0);
  func3->FixParameter(2,func2->GetParameter(2));
  func3->SetNpx(1000);
  hh2_noacc->Fit("func3","NB","",-70.0,85.0);
  func3->Draw("same");
  func2->SetLineColor(9);
  func3->SetLineColor(8);
  
  const double sigma2fwhm = 2.0*sqrt(2.0*log(2.0));
  double mmwidth    = func2->GetParameter(2);
  double mmwidth_er = func2->GetParError(2);  
  mmwidth = mmwidth * sigma2fwhm;
  mmwidth_er = mmwidth_er * sigma2fwhm;
  
  double areaL, areaS;
  double nL, nS;
  double convf = 1./(double)dbin;
  //cout << dbin << " " << convf << endl;
  
  areaL = func2->Integral(-1000,1000);
  nL = areaL * convf;
  areaS = func3->Integral(-1000,1000);
  nS = areaS * convf;
  
  cout << endl;
  cout << " Kaons: " << endl;
  cout << " sigma = " << sigma*1000. << " +/- " << sigma_er*1000. << " ps" << endl;
  
  cout << endl;
  cout << " Missing mass: " << endl;
  cout << " nK = " << nkaon << endl;
  cout << " nL = " << nL << endl;
  cout << " nS = " << nS 
       << " (nS/nL = " << nS/nL << ")"<< endl;
  cout << " FWHM = " << mmwidth << " +/- " << mmwidth_er << " MeV/c^{2}" << endl;
    

}

/* // ifarm1401: ~/workdir/HallA-Online-Tritium/replay/analysis/mm/figures_20181112/
   Processing nLambda.cc...
 FCN=11.235 FROM MIGRAD    STATUS=CONVERGED      72 CALLS          73 TOTAL
                     EDM=1.71588e-11    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.6 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           3.13792e+02   1.24943e+01   6.95479e-03  -1.62458e-07
   2  p1           1.67166e-01   1.27163e-02   6.12112e-06  -2.30502e-04
   3  p2           3.63377e-01   1.26549e-02   5.37446e-06   4.74549e-04
 FCN=29.3245 FROM MIGRAD    STATUS=CONVERGED      89 CALLS          90 TOTAL
                     EDM=1.26591e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.61830e+02   9.10653e+00   1.67862e-02   9.35505e-06
   2  p1          -1.72921e+00   1.58466e-01   4.22270e-04  -4.64615e-04
   3  p2           3.60004e+00   1.86404e-01   3.44784e-04   1.04856e-03
 FCN=999.95 FROM MIGRAD    STATUS=CONVERGED      39 CALLS          40 TOTAL
                     EDM=1.74615e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           4.34669e+01   3.53131e+00   5.45230e-02   2.28987e-05
   2  p1           7.57285e+01   4.05406e-01   6.27532e-03   1.45044e-03
   3  p2           3.60004e+00     fixed    

 Kaons: 
 sigma = 363.377 +/- 12.6549 ps

 Missing mass: 
 nK = 1059.16
 nL = 730.173
 nS = 196.122 (nS/nL = 0.268597)
 FWHM = 8.47744 +/- 0.438948 MeV/c^{2}

*/

