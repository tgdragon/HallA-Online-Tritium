/*
  "parconv.cc"
  A sample code to create new parameter file.
  param[x,xp,y,yp,zt] --> param[x,xp,y,yp,zt,xt]
  
  Toshiyuki Gogami, March 31, 2019
*/


void parconv(){
  
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;
  const int nXt=0; // Note: this is specified to be zero 
  
  //double Y=0.;  // no use
  //double x=1.;  // no use
  int npar    =0; // The number of new parameters
  int npar_int=0; // The number of original parameters
  int temp;
  int a=0,b=0,c=0,d=0,e=0,f=0;
  double par = 0.0;
  
  char InputFileName[500];  // Input file 
  sprintf(InputFileName,"./sample_matrix/newpar_ypt_2.dat");
  ifstream* ifs = new ifstream(InputFileName);
  
  char OutputFileName[500]; // Output file
  sprintf(OutputFileName,"./sample_matrix/testtest.dat");
  ofstream* ofs = new ofstream(OutputFileName);

  for (int n=0;n<nMatT+1;n++){
    for(f=0;f<n+1;f++){
      for(e=0;e<n+1;e++){
	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){ 
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){
		
		if (a+b+c+d+e+f==n){
		  par_int = 0.0;
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt && f<=nXt){
		    *ifs >> par >> temp >> temp >> temp >> temp >> temp;
		    //cout << par << endl;
		    npar_int++;
		  }
		  else{
		    //x = 0.;
		  }
		  *ofs << par_int << " "
		       << a << " " << b << " " << c << " "
		       << d << " " << e << " " << f << endl;
		  npar++;
		}
		    
	      }
	    }
	  }
	}    
      }
    }
  }
  ifs->close();
  ofs->close();

  cout << endl;
  cout << " Input  file: " << InputFileName  << " " << npar_int<< " parameters" << endl;
  cout << " Output file: " << OutputFileName << " " << npar    << " parameters" << endl;
}
