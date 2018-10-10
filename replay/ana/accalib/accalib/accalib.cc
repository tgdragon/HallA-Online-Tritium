/*
  "accalib.cc"
  Pedestal and gain parameters calibration for A1
  
  Toshiyuki Gogami, 6 Oct 2018
*/

const int n = 24; // Number of channels for A1
extern void CreateHistrootF();

double func1(double* x, double* par){
  double g1,g2,expo;
  g1 = par[0] * exp(-0.5*pow((x[0]-par[1])/par[2],2.0));
  g2 = par[3] * exp(-0.5*pow((x[0]-par[4])/par[5],2.0));
  expo = par[6]*exp(par[7]*(x[0]-par[1]));
  if(x[0]<0) expo = 0.0;
  //if(x[0]<0 || x[0]>500.0) expo = 0.0;
  //if(x[0]<0 || x[0]>400.0) expo = 0.0;
  return g1+g2+expo;
}
double func1_g(double* x, double* par){
  double g1;
  g1 = par[0] * exp(-0.5*pow((x[0]-par[1])/par[2],2.0));
  return g1;
}
double func1_e(double* x, double* par){
  double expo;
  expo = par[0]*exp(par[1]*(x[0]-par[2]));
  if(x[0]<0) expo = 0.0;
  return expo;
}



// ####################################### //
//   Main code to get pedestal and SPE     //
//   positions                             //
// ####################################### // 
int accalib(int ch = 10){
  
  if(ch<0 || ch>23) return 1;
  
  //bool iwanttoprint = true; // (true to print out)
  bool iwanttoprint = false; // (true to print out)
  
  // ==============~======================== //
  // ====== Opening a ROOT file ============ //
  // ======================================= //
  TFile* f1 = new TFile("hist.root");
  if(f1->IsOpen()==false) {
    cout << " hist.root was not found" << endl;
    cout << " --> Now, creating..." << endl;
    CreateHistrootF();
    f1 = new TFile("hist.root");
  }
  cout << endl;
  cout << " Channel = " << ch << " will be analyzed." << endl;
  cout << endl;
  char hname[500];
  sprintf(hname,"a1adc%d",ch);
  TH1F* h1 = (TH1F*)f1->Get(hname);
  h1->GetXaxis()->SetRangeUser(-300.0,1800.0);

  TCanvas* c1 = new TCanvas("c1","c1");
  gPad->SetLogy(1);
  h1->Draw();

  // ====================== //
  // ====== Fitting ======= //
  // ====================== //
  TF1* ffunc1 = new TF1("ffunc1",func1,-300.0,1500.0,8);
  ffunc1->SetNpx(1000);
  ffunc1->SetParameters(1.0e+4,0.0,30.0,
			1.0e+2,400.0,100.0,
			1.0e+2,-0.005);
  //ffunc1->SetParLimits(6,0.0,1.0+3);
  ffunc1->SetParLimits(7,-0.05,0);
  //ffunc1->FixParameter(6,0.0);
  h1->Fit("ffunc1","NB","",-200.0,800.0);
  ffunc1->Draw("same");

  TF1* ff1 = new TF1("ff1",func1_g,-300.,300.0,3);
  ff1->SetParameters(ffunc1->GetParameter(0),
		     ffunc1->GetParameter(1),
		     ffunc1->GetParameter(2));
  ff1->SetLineColor(1);
  ff1->SetLineStyle(7);
  ff1->Draw("same");
  
  TF1* ff2 = new TF1("ff2",func1_g,-300.0,1500.0,3);
  ff2->SetParameters(ffunc1->GetParameter(3),
		     ffunc1->GetParameter(4),
		     ffunc1->GetParameter(5));
  ff2->SetLineColor(1);
  ff2->SetLineStyle(7);
  ff2->Draw("same");
  
  TF1* ff3 = new TF1("ff3",func1_e,0.0,1000.0,3);
  ff3->SetParameters(ffunc1->GetParameter(6),
		     ffunc1->GetParameter(7),
		     ffunc1->GetParameter(1));
  ff3->SetLineColor(1);
  ff3->SetLineStyle(7);
  ff3->Draw("same");

  double ped, ped_er;
  double pedw, pedw_er;
  double spe, spe_er;
  ped     = ffunc1->GetParameter(1); // pedestal position
  ped_er  = ffunc1->GetParError(1);  // error on pedestal position 
  pedw    = ffunc1->GetParameter(2); // pedestal width
  pedw_er = ffunc1->GetParError(2);  // error on pedestal width 
  spe     = ffunc1->GetParameter(4); // spe peak position 
  spe_er  = ffunc1->GetParError(4);  // error on spe peak position
  cout << " ch  ped  ped_er  pedw  pedw_er  spe  spe_er ==> " << endl;
  cout << " "
       << ch << " "
       << ped << " " << ped_er << " "
       << pedw << " " << pedw_er << " "
       << spe-ped << " " << sqrt(pow(spe_er,2.0)+pow(ped_er,2.0))
       << endl;
  cout << endl;

  if (iwanttoprint==true){
    char pname[500];
    sprintf(pname,"./figures/fit_a1adc_%d.png",ch);
    c1->Print(pname);
    sprintf(pname,"./figures/fit_a1adc_%d.eps",ch);
    c1->Print(pname);
  }

  return 0;
}




// ####################################### //
//   A creation of a new ROOT file for     //
//   histograms to be fit                  //
// ####################################### // 

void CreateHistrootF(){ 
  // ==============~======================== //
  // ====== Opening a ROOT file ============ //
  // ======================================= //
  TFile* f1 = new TFile("tritium_93740.root");
  TTree* t1 = (TTree*)f1->Get("T");

  //const int n = 24;
  double adc[n];
  t1->SetBranchAddress("R.a1.a_c", &adc);
  double ent = t1->GetEntries();
  //ent = 1000;

  // ================~~=============== //
  // ======= Creating a hitogram ===== //
  // ================================= //
  TH1F* h0 = new TH1F("h0","",300,-300,2700);
  h0->GetXaxis()->SetTitle("ch");
  h0->GetYaxis()->SetTitle("Counts");
  TH1F* hist_adc[n];
  char hname[500];
  for(int i=0 ; i<n ; i++){
    sprintf(hname,"a1adc%d",i);
    hist_adc[i] = (TH1F*)h0->Clone(hname);
  }
  

  // ================================================= //
  // ======== Filling the histograms with data ======= //
  // ================================================= //
  for(int i=0 ; i<ent ; i++){

    for(int j=0 ; j<n ; j++){
      adc[j] = -2222.0;
    }
    
    t1->GetEntry(i);
    
    for(int j=0 ; j<n ; j++){
      hist_adc[j]->Fill(adc[j]);
    }

  }

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(5,5);
  for(int i=0 ; i<n ; i++){
    c1->cd(i+1);
    gPad->SetLogy(1);
    hist_adc[i]->Draw();
  }

  
  // ====================================================== //
  // ======== Creating a new ROOT file for histograms ===== //
  // ====================================================== //
  TFile* fnew = new TFile("hist.root","recreate");
  TObjArray h(1);
  h.Add(c1);
  h.Add(hist_adc[0]);
  h.Add(hist_adc[1]);
  h.Add(hist_adc[2]);
  h.Add(hist_adc[3]);
  h.Add(hist_adc[4]);
  h.Add(hist_adc[5]);
  h.Add(hist_adc[6]);
  h.Add(hist_adc[7]);
  h.Add(hist_adc[8]);
  h.Add(hist_adc[9]);
  h.Add(hist_adc[10]);
  h.Add(hist_adc[11]);
  h.Add(hist_adc[12]);
  h.Add(hist_adc[13]);
  h.Add(hist_adc[14]);
  h.Add(hist_adc[15]);
  h.Add(hist_adc[16]);
  h.Add(hist_adc[17]);
  h.Add(hist_adc[18]);
  h.Add(hist_adc[19]);
  h.Add(hist_adc[20]);
  h.Add(hist_adc[21]);
  h.Add(hist_adc[22]);
  h.Add(hist_adc[23]);
  h.Write();
  fnew->Close();
}
