/*
  "npar_matrix.cc"
  
  Checking the number of parameters
  in a matrix.
  
  Toshiyuki Gogami, Nov 8, 2018
*/


void npar_matrix(){
  
  int npar=0;
  //const int n=3; // order
  const int n=4; // order
  
  for(int i=0; i<n+1 ; i++){
    for(int d=0 ; d<n+1 ; d++){
      for(int c=0 ; c<n+1 ; c++){
	for(int b=0; b<n+1 ; b++){
	  for(int a=0; a<n+1 ; a++){
	    if(a+b+c+d == i) npar++;
	  }
	}
      }
    }
  }
  cout << n << ": " << npar << endl;
}
