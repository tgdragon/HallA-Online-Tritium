/*
  newpar.cc
  
  Toshiyuki Gogami, Nov 9, 2018
*/

void newpar(){
  //ifstream* ifs1 = new ifstream("s2_t0_L.dat");
  //ifstream* ifs2 = new ifstream("s2_t0_L_old.dat");
  //ofstream* ofs = new ofstream("s2_t0_L_2.dat");
  ifstream* ifs1 = new ifstream("s2_t0_R.dat");
  ifstream* ifs2 = new ifstream("s2_t0_R_old.dat");
  ofstream* ofs = new ofstream("s2_t0_R_2.dat");
  const int n=16;
  double par1[n];
  double par2[n];
  double newpar[n];
  for(int i=0 ; i<n ; i++){
    *ifs1 >> par1[i];
    *ifs2 >> par2[i];
    newpar[i] = par1[i] + par2[i];
    cout << " " << par1[i] << "+" << par2[i] << " = " << newpar[i] << endl;
    *ofs << newpar[i] << endl;
  }
  ofs->close();
  ifs1->close();
  ifs2->close();
  
  
  
}
