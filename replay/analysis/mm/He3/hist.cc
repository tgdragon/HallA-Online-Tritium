/*
  h2.root
  
  T. Gogami, Aug 9, 2019
*/


void hist(){
  // ============================== //
  // ===== Open ROOT file ======== //
  // ============================== //
  TFile* f1 = new TFile("T2.root");
  TTree* t1 = (TTree*)f1->Get("tree");
  
  
  // ============================ //
  // ===== 
  // ============================ //
  const int max = 100;
  double a1, a2;
  double mm[max];
  double ctime[max];
  
  t1->SetBranchAddress("ctime", &ctime);
  t1->SetBranchAddress("mm", &mm);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);

  TH1F* hist1 = new TH1F("hist1","",100,-20.0,80.0);
  
  double ent = t1->GetEntries();
  
  for(int i=0 ; i<ent ; i++){
    for(int j=0 ; j<max ; j++){
      ctime[j] = -2222.0;
      mm[j]    = -2222.0;
    }
    a1 = -22222.0;
    a2 = -22222.0;
    
    
    t1->GetEntry(i);
    
    if(abs(ctime[0])<1.0
       && a1<3.0
       && a2>2.0
       && a2<15.0){
      hist1->Fill(mm[0]);
    }
    
  }
  
  TCanvas* c1 = new TCanvas("c1","c1");
  hist1->Draw();
  
}
  

  
