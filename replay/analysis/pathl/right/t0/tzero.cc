/*
  tzero.cc
  
  Toshiyuki Gogami, Nov 7, 2018
*/

void tzero (){
  TFile* f1 = new TFile("out_20181111.root");
  const int n=16; // The number of segment in S2
  TH1F* hist[n];  
  char hname[500];
  char fname[500];
  char fff[500];
  double par[n];
  par[0] = 0.0;
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(4,4);
  TF1* func[n];
  
  //ofstream* ofs = new ofstream("s2_t0_L.dat");
  ofstream* ofs = new ofstream("s2_t0_R.dat");

  for(int i=0 ; i<n ; i++){
    sprintf(fname,"ffunc_%d",i);
    //func[i] = new TF1(fname,fff,-20.0,20.0,3);
    func[i] = new TF1(fname,"[0]*exp(-0.5*pow((x-[1])/[2],2.0))",-20.0,20.0);
    func[i]->SetParameters(100,0.0,1.3);
    sprintf(hname,"hs2_%d",i);
    hist[i] = (TH1F*)f1->Get(hname);
    c1->cd(i+1); 
    hist[i]->Draw();
    if(i!=0){
      hist[i]->Fit(fname,"NQB","",-3.0,3.0);
      //hist[i]->Fit(fname,"NQB","",-2.0,2.0);
      func[i]->Draw("same");
      par[i] = func[i]->GetParameter(1);
    }
    *ofs << par[i] << " " << endl;
  }
  ofs->close();
  
  
  
  
}
