/*
  mm.cc
  Missing mass reconstuction 
  
  Toshiyuki Gogami, Nov 9, 2018
*/

const double me = 0.000511;
const double mk = 0.493677;
//const double mp = 0.938272;
//const double mL = 1.115683;
//const double mp = 2.80839133; // He3 mass
//const double mL = 1.875612762 + 1.115683; // d+L
//
const double mp = 2.80892086;
const double mL = 0.939565379*2.0 + 1.115683; // nn+L


const double hrs_ang = 13.2 * 3.14159 / 180.;
extern double CalcMM(double ee, double* par_ep, double* par_k, double mt);

double nmix = 1;


void mixed_event(){
  TFile* f1 = new TFile("T2_20190920.root");
  TTree* t1 = (TTree*)f1->Get("tree");
  Double_t trig5;
  Double_t trig4;
  Double_t trig1;
  double ent = t1->GetEntries();
  //ent = 100000; // for test
    
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
  double vz_mean[max];
  double l_cer;
  
  
  t1->SetBranchAddress("fEvtHdr.fRun", &runnum    );
  t1->SetBranchAddress("HALLA_p", &hallap );
  //t1->SetBranchAddress("DR.T1", &trig1    );
  //t1->SetBranchAddress("DR.T4", &trig4    );
  //t1->SetBranchAddress("DR.T5", &trig5    );
  //t1->SetBranchAddress("R.tr.time", &rtime);
  //t1->SetBranchAddress("L.tr.time", &ltime);
  //t1->SetBranchAddress("R.tr.pathl", &rpathl);
  //t1->SetBranchAddress("L.tr.pathl", &lpathl);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  t1->SetBranchAddress("L.cer.asum_c", &l_cer);
  t1->SetBranchAddress("vz_mean", &vz_mean);
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
  //t1->SetBranchAddress("R.s0.time", &rtime_s0);
  //t1->SetBranchAddress("L.s0.time", &ltime_s0);
  //t1->SetBranchAddress("R.s2.time", &rtime_s2);
  //t1->SetBranchAddress("L.s2.time", &ltime_s2);
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
  //t1->SetBranchAddress("R.tr.beta",  &rbeta);
  //t1->SetBranchAddress("L.tr.beta",  &lbeta);
  //t1->SetBranchAddress("R.s2.trpath",  &rpathl_s2);
  //t1->SetBranchAddress("L.s2.trpath",  &lpathl_s2);
  //t1->SetBranchAddress("L.s2.nthit",&nhit);
  //t1->SetBranchAddress("R.s2.nthit",&nhit_R);
  t1->SetBranchAddress("R.ps.asum_c", &ps_asum);
  t1->SetBranchAddress("ctime", &ctime);
  //t1->SetBranchAddress("R.a1.t_fadc", &a1_tdc);
  //t1->SetBranchAddress("R.a2.t_fadc", &a2_tdc);
  
  
  //TCanvas* c1 = new TCanvas("c1","c1");
  TH1F* h1  = new TH1F("h1","",200,-20.,20.0); // e'K+ coincidence
  h1->GetXaxis()->SetTitle("coin time (ns)");
  h1->GetYaxis()->SetTitle("Counts / 200 ps");
  h1->GetXaxis()->SetRangeUser(-14.0,17.);
  TH1F* h1_ = (TH1F*)h1->Clone("h1_"); 
  TH1F* h1_acc = (TH1F*)h1->Clone("h1_acc"); // accidentals 
  TH1F* h1_acc_shift = (TH1F*)h1->Clone("h1_shift"); 
  
  h1 ->SetLineColor(2);
  h1->SetLineWidth(2);
  h1_->SetLineColor(1);
  h1_acc->SetLineColor(1);
  h1_acc->SetFillColor(1);
  h1_acc->SetFillStyle(3001);
  h1_acc_shift->SetLineColor(9);

  double xmin = -300.0, xmax = 300.0; int xbin = 300; // 2 MeV / bin
  //double xmin = -300.0, xmax = 300.0; int xbin = 600; // 1 MeV / bin
  //double xmin = -300.0, xmax = 300.0; int xbin = 200; // 3 MeV / bin
  TH1F* h2  = new TH1F("h2","",xbin,xmin,xmax);
  //TH1F* h2  = new TH1F("h2","",600,-300,300.);
  //h2->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (MeV/c^{2})");
  //h2->GetYaxis()->SetTitle("Counts / (2 MeV/c^{2})");
  h2->GetXaxis()->SetTitle("-B_{#Lambda} (MeV)");
  h2->GetYaxis()->SetTitle("Counts / 2 MeV");
  h2->GetXaxis()->SetRangeUser(-150.0,100.0);
  h2->SetLineColor(1);
  TH1F* h2_acc = (TH1F*)h2->Clone("h2_acc");
    
  TH2F* h3 = new TH2F("h3","",200,-3,30,200,-3,50.);
  h3->GetXaxis()->SetTitle("A1 NPE");
  h3->GetYaxis()->SetTitle("A2 NPE");
  h3->GetXaxis()->SetRangeUser(-1.0,12.0);
  h3->GetYaxis()->SetRangeUser(-1.0,60.0);
  TH2F* h3_k = (TH2F*)h3->Clone("h3_k");
  
  TH2F* h4 = new TH2F("h4","",200,-20,20,400,-3,82.);
  h4->GetXaxis()->SetTitle("Coin time (ns)");
  h4->GetYaxis()->SetTitle("A1 NPE");
  h4->GetXaxis()->SetRangeUser(-10.0,15.0);
  h4->GetYaxis()->SetRangeUser(-0.5,12.0);
  
  TH2F* h5 = new TH2F("h5","",200,-20,20,200,-3,82.);
  h5->GetXaxis()->SetTitle("Coin time (ns)");
  h5->GetYaxis()->SetTitle("A2 NPE");
  h5->GetXaxis()->SetRangeUser(-10.0,15.0);
  h5->GetYaxis()->SetRangeUser(-0.5,60.0);
  
  TH2F* h6 = new TH2F("h6","",300,-100,200,100,-1.5,1.5);
  h6->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h6->GetYaxis()->SetTitle("x_FP (e') (m)");
  h6->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h7 = new TH2F("h7","",300,-100,200,100,-0.5,0.5);
  h7->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h7->GetYaxis()->SetTitle("x'_FP (e') (rad)");
  h7->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h8 = new TH2F("h8","",300,-100,200,100,-0.1,0.1);
  h8->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h8->GetYaxis()->SetTitle("y_FP (e') (m)");
  h8->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h9 = new TH2F("h9","",300,-100,200,100,-0.1,0.1);
  h9->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h9->GetYaxis()->SetTitle("y'_FP (e') (rad)");
  h9->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h10 = new TH2F("h10","",300,-100,200,100,-1.5,1.5);
  h10->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h10->GetYaxis()->SetTitle("x_FP (K^{+}) (m)");
  h10->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h11 = new TH2F("h11","",300,-100,200,100,-0.5,0.5);
  h11->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h11->GetYaxis()->SetTitle("x'_FP (K^{+}) (rad)");
  h11->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h12 = new TH2F("h12","",300,-100,200,100,-0.1,0.1);
  h12->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h12->GetYaxis()->SetTitle("y_FP (K^{+}) (m)");
  h12->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h13 = new TH2F("h13","",300,-100,200,100,-0.1,0.1);
  h13->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h13->GetYaxis()->SetTitle("y'_FP (K^{+}) (rad)");
  h13->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h14 = new TH2F("h14","",300,-100,200,100,-0.3,0.3);
  h14->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h14->GetYaxis()->SetTitle("z (e') (m)");
  h14->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  TH2F* h15 = (TH2F*)h14->Clone("h15");
  h15->GetXaxis()->SetTitle("M_{x}-M_{#Lambda} (GeV/c^{2})");
  h15->GetYaxis()->SetTitle("z (K^{+}) (m)");
  //h15->GetXaxis()->SetRangeUser(-50.0,200.0); 
  
  
  bool ctime_selection = false;
  bool ac1_selection = false;
  bool ac2_selection = false;
  bool lcer_selection = false;
  bool zR_selection = false;
  bool zL_selection = false;
  bool vz_selection =false;
  bool acc1_selection = false;
  bool acc2_selection = false;
  bool acc3_selection = false;
  bool acc4_selection = false;
  bool acc5_selection = false;
  //const double kcenter = 3.1;
  //const double kcenter = 3.3; // h22.root
  const double kcenter = 0.0;
  const double rf_sep = 2.0;
  double par_ep[3];
  double par_k[3];
  double mm, mm_1st_cor;
  double shift;

  //ent = 50000; // for test
  // ================================== //
  // ======== Main loop =============== //
  // ================================== //
  for(int i=0 ; i<ent ; i++){
    
    for(int j=0 ; j<max ; j++){
      ctime[j] = -2222.0;
      rvz[j]   = -2222.0;
      lvz[j]   = -2222.0;
    }
    a1 = -2222.0;
    a2 = -2222.0;
    hallap = -2222.0;
    zR_selection    = false;
    zL_selection    = false;
    vz_selection    = false;
    ac1_selection   = false;
    ac2_selection   = false;
    lcer_selection  = false;
    
    ctime_selection = false;
    acc1_selection  = false;
    acc2_selection  = false;
    acc3_selection  = false;
    acc4_selection  = false;
    acc5_selection  = false;
    shift = 0.0;
    

    t1->GetEntry(i);
    ctime[0] = ctime[0] - kcenter;
    
    if(fabs(rvz[0])<0.1) zR_selection = true;
    else zR_selection = false;
    
    if(fabs(lvz[0])<0.1) zL_selection = true;
    else zL_selection = false;
    
    if (fabs(rvz[0]-lvz[0])<0.025
	&& fabs(vz_mean[0])<0.1){
      vz_selection =true;
    }
    else vz_selection = false;
    
    if(a1 < 3.0) ac1_selection = true;
    else ac1_selection = false;
    
    if(a2 > 3.0 && a2 <18.0 ) ac2_selection = true;
    else ac2_selection = false;
    
    if(l_cer>1500.0) lcer_selection = true;
    else lcer_selection = false;
    
    if(fabs(ctime[0])<1.0) ctime_selection = true;
    else  ctime_selection = false;
    
    if(abs(ctime[0]+5.0*rf_sep)<1.0) acc1_selection = true;
    else acc1_selection = false;
    if(abs(ctime[0]-1.0*rf_sep)<1.0) acc2_selection = true;
    else acc2_selection = false;
    if(abs(ctime[0]-2.0*rf_sep)<1.0) acc3_selection = true;
    else acc3_selection = false;
    if(abs(ctime[0]-7.0*rf_sep)<1.0) acc4_selection = true;
    else acc4_selection = false;
    if(abs(ctime[0]+4.0*rf_sep)<1.0) acc5_selection = true;
    else acc5_selection = false;
    
    
    //if(zR_selection==true && zL_selection==true){
    if(vz_selection==true){
      h3->Fill(a1,a2);
      h4->Fill(ctime[0],a1);
      h5->Fill(ctime[0],a2);
      if(ctime_selection){
	h3_k->Fill(a1,a2);
      }
    }
    
    
    if( //zR_selection     == true
	//&& zL_selection  == true
       //vz_selection == true
       //&& ac1_selection == true
       //&& ac2_selection == true 
       //&& lcer_selection==true
       1 == 1){

      
      par_ep[0] = mom2[0];
      par_ep[1] = th2[0];// right handed system
      par_ep[2] = ph2[0];// right handed system
      // ---- 400 um thick target -----
      //double dpe  = 184.3e-6; // GeV/c
      double dpep = 0.0; // GeV/c
      //double dpk  = 0.0; // GeV/c
      if(vz_mean[0]<8.0e-2){
	dpep = -1.35758 * sin(-4.59571*par_ep[2]) + 2.09;   // MeV/c
	//dpk  = -1.31749 * sin(-4.61513*par_k[2] ) + 2.0368; // MeV/c
      }
      else {
	dpep = 6.23e-3 * par_ep[2] + 0.403; // MeV/c
	//dpk  = 3.158e-2* par_k[2]  + 0.4058;// MeV/c
      }
      dpep = dpep / 1000.0; // MeV/c --> GeV/c
      //dpk  = dpk  / 1000.0; // MeV/c --> GeV/c
      //hallap = hallap - dpe;
      par_ep[0] = par_ep[0] + dpep;
      // par_k[0]  = par_k[0]  + dpk;
      
      if (vz_selection == true
	  && ac1_selection == true
	  && ac2_selection == true ){
	if(ctime_selection==true){
	  par_ep[0] = mom2[0];
	  par_ep[1] = th2[0];// right handed system
	  par_ep[2] = ph2[0];// right handed system
	  
	  par_k[0] = mom1[0];
	  par_k[1] = th1[0]; // right handed system
	  par_k[2] = ph1[0]; // right handed system
	  
	  // ---- 400 um thick target -----
	  double dpe  = 184.3e-6; // GeV/c
	  //double dpep = 0.0; // GeV/c
	  dpep = 0.0; // GeV/c
	  double dpk  = 0.0; // GeV/c
	  
	  if(vz_mean[0]<8.0e-2){
	    dpep = -1.35758 * sin(-4.59571*par_ep[2]) + 2.09;   // MeV/c
	    dpk  = -1.31749 * sin(-4.61513*par_k[2] ) + 2.0368; // MeV/c
	    
	  }
	  else {
	    dpep = 6.23e-3 * par_ep[2] + 0.403; // MeV/c
	    dpk  = 3.158e-2* par_k[2]  + 0.4058;// MeV/c
	  }
	  
	  dpep = dpep / 1000.0; // MeV/c --> GeV/c
	  dpk  = dpk  / 1000.0; // MeV/c --> GeV/c
	  
	  hallap = hallap - dpe;
	  par_ep[0] = par_ep[0] + dpep;
	  par_k[0]  = par_k[0]  + dpk;
	  
	  //      cout << hallap << " " 
	  //	   << par_ep[0] << " " << par_ep[1] << " " << par_ep[2] << " " << dpep << " "
	  //	   << par_k[0]  << " " << par_k[1]  << " " << par_k[2] << " " << dpk
	  //	   << endl;
	  
	  mm = CalcMM(hallap, par_ep, par_k, mp);
	  
	  
	  mm = (mm-mL)*1000.; // Gen --> MeV
	  
	  h1->Fill(ctime[0]);
	  
	  h2->Fill(mm);
	  h6->Fill(mm, l_x_fp[0]);
	  h7->Fill(mm, l_th_fp[0]);
	  h8->Fill(mm, l_y_fp[0]);
	  h9->Fill(mm, l_ph_fp[0]);
	  h14->Fill(mm, lvz[0]);
	  
	  h10->Fill(mm, r_x_fp[0]);
	  h11->Fill(mm, r_th_fp[0]);
	  h12->Fill(mm, r_y_fp[0]);
	  h13->Fill(mm, r_ph_fp[0]);
	  h15->Fill(mm, rvz[0]);
	}
	else h1_->Fill(ctime[0]);	
      }
      
      
      if(acc1_selection==true || 
	 acc2_selection==true ||
	 acc3_selection==true ||
	 acc4_selection==true || 
	 acc5_selection==true
	 //-8.0 < ctime[0] && ctime[0]<16.0
	 ){
	
	//if(acc1_selection ==true)      shift= 5.0*rf_sep;
	//else if(acc2_selection ==true) shift=-1.0*rf_sep;
	//else if(acc3_selection ==true) shift=-2.0*rf_sep;
	//else if(acc4_selection ==true) shift=-7.0*rf_sep;
	//else if(acc5_selection ==true) shift= 4.0*rf_sep;
	//else shift=-2222.0;
	
	//h1_acc       -> Fill(ctime[0]);
	//h1_acc_shift -> Fill(ctime[0]+shift);

	
	for(int j=0 ; j<(int)nmix ; j++){
	  
	  int mixed_ent = i+j;

	  if(mixed_ent<ent){
	    t1->GetEntry(mixed_ent);
	  }
	  else {
	    t1->GetEntry(mixed_ent-ent);
	  }

	  if(a1<3.0 && a2>3.0 && a2<18.0
	     && fabs(rvz[0]-lvz[0])<0.025
	     && fabs((rvz[0]+lvz[0])/2.0)<0.1 ){
	    par_k[0] = mom1[0];
	    par_k[1] = th1[0]; 
	    par_k[2] = ph1[0]; 
	    
	    // ---- 400 um thick target -----
	  double dpk  = 0.0; // GeV/c
	  
	  if(rvz[0] < 8.0e-2){
	    dpk  = -1.31749 * sin(-4.61513*par_k[2] ) + 2.0368; // MeV/c
	  }
	  else {
	    dpk  = 3.158e-2* par_k[2]  + 0.4058;// MeV/c
	  }
	  
	  dpk  = dpk  / 1000.0; // MeV/c --> GeV/c
	  
	  par_k[0]  = par_k[0]  + dpk;
	  
	  double mixed_mm = 0.0;
	  mixed_mm = CalcMM(hallap, par_ep, par_k, mp);
	  mixed_mm = (mixed_mm - mL) * 1000.; // Gen --> MeV
	  h2_acc   -> Fill(mixed_mm);
	  //cout << mixed_mm << endl;
	  }
	}
      } 
    }
  }
  
  h1_acc_shift->Scale(1./5.);
  //h2_acc->Scale(1./5.);
  //h2_acc->Scale(1./nmix/5.);
  //h2_acc->Scale(1./nmix/11.5);

  double dbin = (xmax-xmin)/(double)xbin;
  double minx=-100.0,maxx=-20.0;
  int fitmin = (minx-xmin)/dbin;
  int fitmax = (maxx-xmin)/dbin;
  double num1 = h2->Integral(fitmin,fitmax);
  double num2 = h2_acc->Integral(fitmin,fitmax);
  double mixscale = num1/num2;
  //cout << fitmin << "," << fitmax << ": "<< num1 << "/" << num2 << "= " << mixscale<< endl;
  cout << endl;
  cout << " " << minx << "-" << maxx << " MeV: mixscale = "<< num1 << "/" << num2 << "= " << mixscale<< endl;
  cout << " " << nmix << " x 5 bunches" << " =  " << nmix * 5.0 << " (" << "effective scale = "
       << num2 / num1 << ")" << endl;
    
  cout << endl;
  TH1F* h2_acc2 = (TH1F*)h2_acc->Clone();
  for(int i=0 ; i<xbin ; i++){
    h2_acc2->SetBinContent( i+1,
			    h2_acc->GetBinContent(i+1)*mixscale);
    h2_acc2->SetBinError( i+1,
			  h2_acc->GetBinError(i+1)*mixscale);
  }
  
  TCanvas* c1 = new TCanvas("c1","c1");
  h1_->Draw();
  h1 ->Draw("same");
  h1_acc->Draw("same");
  h1_acc_shift->Draw("same");
  
  TCanvas* c2 = new TCanvas("c2","c2");
  h2->Draw(); // Missing mass
  h2_acc2->SetLineColor(9);
  h2_acc2->Draw("same");
  
  //TCanvas* c3 = new TCanvas("c3","c3");
  //gPad->SetLogz(1);
  //h3->Draw("colz");
  
  //TCanvas* c3_k = new TCanvas("c3_k","c3_k");
  //gPad->SetLogz(1);
  //h3_k->Draw("colz");
  
  //TCanvas* c4= new TCanvas("c4","c4");
  //c4->Divide(1,2);
  //c4->cd(1);gPad->SetLogz(1);h4->Draw("colz");
  //c4->cd(2);gPad->SetLogz(1);h5->Draw("colz");
  
  //TCanvas* c5= new TCanvas("c5","c5");
  //c5->Divide(2,3);
  //c5->cd(1);gPad->SetLogz(0);h6->Draw("colz");
  //c5->cd(2);gPad->SetLogz(0);h7->Draw("colz");
  //c5->cd(3);gPad->SetLogz(0);h8->Draw("colz");
  //c5->cd(4);gPad->SetLogz(0);h9->Draw("colz");
  //c5->cd(5);gPad->SetLogz(0);h14->Draw("colz");
  
  //TCanvas* c6= new TCanvas("c6","c6");
  //c6->Divide(2,3);
  //c6->cd(1);gPad->SetLogz(0);h6->Draw("colz");
  //c6->cd(2);gPad->SetLogz(0);h7->Draw("colz");
  //c6->cd(3);gPad->SetLogz(0);h8->Draw("colz");
  //c6->cd(4);gPad->SetLogz(0);h9->Draw("colz");
  //c6->cd(5);gPad->SetLogz(0);h15->Draw("colz");
  
  TFile* fnew = new TFile("hist.root","recreate");
  TObjArray h(0);
  h.Add(h1);
  h.Add(h1_);
  h.Add(h1_acc);
  h.Add(h1_acc_shift);
  h.Add(h2);
  h.Add(h2_acc);
  h.Add(h3);
  h.Add(h3_k);
  h.Add(h4);
  h.Add(h5);
  h.Add(h6);
  h.Add(h7);
  h.Add(h8);
  h.Add(h9);
  h.Add(h10);
  h.Add(h11);
  h.Add(h12);
  h.Add(h13);
  h.Add(h14);
  h.Add(h15);
  h.Add(c1);
  h.Add(c2);
  //h.Add(c3);
  //h.Add(c4);
  //h.Add(c5);
  //h.Add(c6);
  h.Write();
  fnew->Close();
  
  //====== Print ===== //
  //c1->Print("cointime_H2_20181110.png","png");
  //c2->Print("MissigMass_H2_20181110.png","png");
  //c3->Print("a1a2_H2_20181110.png","png");
  //c3_k->Print("a1a2_kaon_H2_20181110.png","png");
  //c4->Print("cointime_a1a2_H2_20181110.png","png");
  //c5->Print("MM_vs_pare_H2_20181110.png","png");  
  //c6->Print("MM_vs_park_H2_20181110.png","png");  
  
}


double CalcMM(double ee, double* par_ep, double* par_k, double mt){
  
  double pe = ee;
  double Ee = sqrt(me*me + pe*pe);
  TVector3 vec_e (0.0, 0.0, pe);
  
  double pep  = par_ep[0];
  double xpep = par_ep[1];
  double ypep = par_ep[2];
  double px_ep, py_ep, pz_ep;
  pz_ep = pep / sqrt(1.0 + xpep*xpep + ypep*ypep);
  px_ep = xpep * pz_ep;
  py_ep = ypep * pz_ep;
  TVector3 vec_ep (px_ep, py_ep, pz_ep);
  vec_ep.RotateY(hrs_ang);
  //double Eep = sqrt(vec_ep * vec_ep);
  double Eep = sqrt(pep*pep + me*me);
  
  double pk  = par_k[0];
  double xpk = par_k[1];
  double ypk = par_k[2];
  double px_k, py_k, pz_k;
  pz_k = pk / sqrt(1.0 + xpk*xpk + ypk*ypk);
  px_k = xpk * pz_k;
  py_k = ypk * pz_k;
  TVector3 vec_k (px_k, py_k, pz_k);
  vec_k.RotateY(-hrs_ang);
  //double Ek = sqrt(vec_k * vec_k);
  double Ek = sqrt(pk*pk + mk*mk);
  
  double missingE2 = 0.0, missingP2 = 0.0, missingM2 = 0.0;
  missingE2 = pow(Ee + mt - Ek - Eep, 2.0);
  missingP2 = (vec_e - vec_ep - vec_k) * (vec_e - vec_ep - vec_k);
  missingM2 = missingE2 - missingP2;
  
  double MissingMass = 0.0;
  MissingMass = sqrt(missingM2);

  return MissingMass;
  
}
