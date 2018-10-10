/*
  CreateParameters.cc
  
  This macro creates new parameters
  using original parameters and parameters
  obtained by the fitting (A1)

  Toshiyuki Gogami, Oct 6, 2018
*/

const int n = 24; // The number of channels in A1


void CreateParameters(){
  ifstream* org_param = new ifstream("org.dat");
  ifstream* new_param = new ifstream("a1_fitresults.dat");

  double opar1[n]; // pedestal
  double opar2[n]; // scale factor
  double npar1[n]; // pedestal (new)
  double npar2[n]; // scale factor (new)
  double temp;
  char tempc[500];
  cout << "READING ORGINAL PARAMETERS...." << endl;
  *org_param >> tempc;
  cout << tempc << endl;
  for(int i=0 ; i<n ; i++){
    *org_param >> opar1[i];
    cout << " " << i << ": "<<  opar1[i] << endl;
  }
  *org_param >> tempc;
  cout << tempc << endl;
  for(int i=0 ; i<n ; i++){
    *org_param >> opar2[i];
    cout << " " << i << ": "<<  opar2[i] << endl;
  }
  org_param->close();
  cout << endl;

  cout << "READING NEW PARAMETERS...." << endl;
  cout << "pedestal: " << endl;
  for(int i=0 ; i<n ; i++){
    *new_param >> temp >> npar1[i]  >> temp >> temp >> temp >> npar2[i] >> temp;
    cout << " " << i << ": "<<  npar1[i] << endl;
  }
  cout << "spe separation: " << endl;
  for(int i=0 ; i<n ; i++){
    cout << " " << i << ": "<<  npar2[i] << endl;
  }
  new_param->close();
  cout << endl;
  
  ofstream* ofs = new ofstream("db_R.a1.dat_newpedgain");
  *ofs << "R.a1.adc.pedestals = ";
  for(int i=0 ; i<n ; i++){
    *ofs << (int)opar1[i]-(int)npar1[i] << " ";
  }
  *ofs << endl;
  *ofs << endl;
  
  ofs->precision(4);
  *ofs << "R.a1.adc.gains = ";
  for(int i=0 ; i<n ; i++){
    //*ofs << npar2[i]*opar2[i] << " ";
    *ofs << opar2[i]/npar2[i] << " ";
  }
  *ofs << endl;
  
}
