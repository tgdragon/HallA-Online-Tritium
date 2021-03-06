/*
  accheck.cc
  
  Toshiyuki Gogami, 7 Oct 2018
*/

const int n1 = 24; // the number of segments for A1
const int n2 = 26; // the number of segments for A2

int accheck(int runnum = 93740){
  // =================================== //
  // ======= General conditions ======== //
  // =================================== //
  gROOT->SetStyle("Plain");

  
  // ==================================== //
  // ======= Opening a ROOT file ======== //
  // ==================================== //
  char fname[500];
  sprintf(fname,"tritium_%d.root",runnum);
  TFile* f1 = new TFile(fname);
  TTree* t1 = (TTree*)f1->Get("T");

  if(!f1->IsOpen()){
    cout << "  NO  FILE:  " << fname << endl;
    return 1;
  }


  // ====================================  //
  // ======= Creating histograms ========  //
  // ====================================  //
  TH1F* h_npe1 = new TH1F("h_npe1","",400,-10.0,90.0);
  h_npe1->GetXaxis()->SetTitle("NPE");
  h_npe1->GetYaxis()->SetTitle("Counts / 0.25 NPE");
  t1->Project("h_npe1","R.a1.asum_c");
  h_npe1->GetXaxis()->SetRangeUser(-3.0,15.0);

  TH1F* h_npe2 = new TH1F("h_npe2","",800,-10.0,190.0);
  h_npe2->GetXaxis()->SetTitle("NPE");
  h_npe2->GetYaxis()->SetTitle("Counts / 0.25 NPE");
  t1->Project("h_npe2","R.a2.asum_c");
  h_npe2->GetXaxis()->SetRangeUser(-3.0,80.0);

  
  //---- A1 each PMT ----- //
  TH1F* h_a1[n1];
  double xmin1 = -10.0, xmax1 = 90.0;
  int xbin1 = 800;
  char hname[500], nname[500];
  double hist_max1[n1];
  double hmax1 = 0.0;
  double showmin1 = -1.8;
  double showmax1 =  7.0;
  for(int i=0 ; i<n1 ; i++){
    sprintf(hname,"a1_npe%d",i);
    h_a1[i] = new TH1F(hname,"",xbin1,xmin1,xmax1);
    sprintf(nname,"R.a1.a_c[%d]",i);
    t1->Project(hname,nname);
    h_a1[i]->GetXaxis()->SetRangeUser(showmin1,showmax1);
    h_a1[i]->SetLineWidth(2);
    hist_max1[i] = h_a1[i]->GetMaximum() * 1.2;
    //cout << i << ": "  << hist_max1[i] << endl;
    //h_a1[i]->SetMaximum(hist_max1[i]);
    if(hmax1<hist_max1[i]) hmax1 = hist_max1[i];
    //cout << i << " " << hmax1 << endl;
  }

  double temp = 0.0;
  //---- A2 each PMT ----- //
  TH1F* h_a2[n1];
  double xmin2 = -10.0, xmax2 = 190.0;
  int xbin2 = 1000;
  double hist_max2[n2];
  double hmax2 = 0.0;
  double showmin2 = showmin1;
  double showmax2 =  35.0;
  //cout << "aaa " << hmax1 << endl;
  for(int i=0 ; i<n2 ; i++){
    sprintf(hname,"a2_npe%d",i);
    h_a2[i] = new TH1F(hname,"",xbin2,xmin2,xmax2);
    sprintf(nname,"R.a2.a_c[%d]",i);
    t1->Project(hname,nname);
    h_a2[i]->GetXaxis()->SetRangeUser(showmin2,showmax2);
    h_a2[i]->SetLineWidth(2);
    hist_max2[i] = h_a2[i]->GetMaximum() * 1.2;
    //h_a1[i]->SetMaximum(hist_max1[i]);
    if(hmax2<hist_max2[i]) hmax2 = hist_max2[i];
  }

  //---- A1 vs. A2 ----- //
  TH2F* hcor = new TH2F("hcor","",
			h_npe1->GetXaxis()->GetNbins(),
			h_npe1->GetXaxis()->GetXmin(),
			h_npe1->GetXaxis()->GetXmax(),
			h_npe2->GetXaxis()->GetNbins(),
			h_npe2->GetXaxis()->GetXmin(),
			h_npe2->GetXaxis()->GetXmax());
  hcor->GetXaxis()->SetTitle("A1 NPE sum (n = 1.015)");
  hcor->GetYaxis()->SetTitle("A2 NPE sum (n = 1.055)");
  t1->Project("hcor","R.a2.asum_c:R.a1.asum_c");
  

  //cout << "bbb" << hmax1 << " "<< temp<< endl;
  //temp = hmax1;

  // ==================================== //
  // ======= Drawing histograms ========= //
  // ==================================== //
  TCanvas* c1 = new TCanvas("c1","NPE sum");
  c1->Divide(2,1);
  c1->cd(1); h_npe1->Draw(); // NPE sum A1
  c1->cd(2); h_npe2->Draw(); // NPE sum A2

  TCanvas* c2 = new TCanvas("c2","NPE A1");
  TLine* l1_a1 = new TLine(0.0,0.0,0.0,hmax1);
  l1_a1->SetLineColor(9);
  l1_a1->SetLineStyle(9);
  TLine* l2_a1 = new TLine(1.0,0.0,1.0,hmax1);
  l2_a1->SetLineColor(2);
  l2_a1->SetLineStyle(9);
  TLine* l3_a1 = new TLine(2.0,0.0,2.0,hmax1);
  l3_a1->SetLineColor(2);
  l3_a1->SetLineStyle(9);
  c2->Divide(6,4);
  for(int i=0 ; i<n1 ; i++){
    h_a1[i]->SetMaximum(hmax1);
    c2->cd(i+1);
    gPad->SetLogy(1);
    h_a1[i]->Draw();
    l1_a1->Draw("same");
    l2_a1->Draw("same");
    l3_a1->Draw("same");
  }
  c2->Update();

  
  TCanvas* c3 = new TCanvas("c3","NPE A2");
  c3->Divide(7,4);
  TLine* l1_a2 = new TLine(0.0,0.0,0.0,hmax2);
  l1_a2->SetLineColor(9);
  l1_a2->SetLineStyle(9);
  TLine* l2_a2 = new TLine(1.0,0.0,1.0,hmax2);
  l2_a2->SetLineColor(2);
  l2_a2->SetLineStyle(9);
  TLine* l3_a2 = new TLine(2.0,0.0,2.0,hmax2);
  l3_a2->SetLineColor(2);
  l3_a2->SetLineStyle(9);
  for(int i=0 ; i<n2 ; i++){
    h_a2[i]->SetMaximum(hmax2);
    c3->cd(i+1);
    gPad->SetLogy(1);
    h_a2[i]->Draw();
    l1_a2->Draw("same");
    l2_a2->Draw("same");
    //l3_a2->Draw("same");
  }
  c3->Update();

  TCanvas* c4 = new TCanvas("c4","A1 vs. A2");
  gPad->SetLogz(1);
  hcor->GetXaxis()->SetRangeUser(showmin1,showmax1);
  hcor->GetYaxis()->SetRangeUser(showmin2,showmax2);
  hcor->Draw("colz");


  char pname1[500];
  char pname2[500];
  char pname3[500];
  sprintf(pname1,"accheck_%d.pdf(",runnum);
  sprintf(pname2,"accheck_%d.pdf",runnum);
  sprintf(pname3,"accheck_%d.pdf)",runnum);
  c4->Print(pname1);
  c2->Print(pname2);
  c1->Print(pname2);
  c3->Print(pname3);
  

  return 0;

  
}
