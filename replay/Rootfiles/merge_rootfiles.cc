/*
  merge_rootfiles.cc
  Merging ROOT files that are made by nude_coin.cc(py)
  
  Toshiyuki Gogami, July 28,2019
*/

#include <TFile.h>
#include <TChain.h>

#include <fstream>
#include <iostream>
using namespace std;

int main(int argc, char** argv){
  char inputfile[300];
  char file_dir[300];
  //char out_dir[500];
  sprintf(file_dir,"./nnL/coin_dragon2");
  if (argc==2){
    sprintf(inputfile,"%s",argv[1]);
    
  }
  else if (argc==3){
    sprintf(inputfile,"%s",argv[1]);
    sprintf(file_dir, "%s",argv[2]);
  }
  else{
    cout << " Please specify data file (merge_rootfiles.cc)" << endl;
    cout << " --> quit." << endl;
    cout << endl;
    cout << "(a) ./merge_rootfiles (data_file)" << endl;
    cout << "(b) ./merge_rootfiles (data_file) (output_dir)" << endl;
    return 99;
  }
 
  char outfile[500];
  //if(inputfile=="h2.dat"){
  if(!strcmp(inputfile,"h2.dat")){
    //cout << " This is a list of H2 w/ H kinematics: " << endl;
    sprintf(outfile,"%s/h2.root",file_dir);
  }
  else if (!strcmp(inputfile,"h22.dat")){
    sprintf(outfile,"%s/h22.root",file_dir);
  }
  else if (!strcmp(inputfile,"T2.dat")){
    sprintf(outfile,"%s/T2.root",file_dir);
  }
  else if (!strcmp(inputfile,"He3.dat")){
    sprintf(outfile,"%s/He3.root",file_dir);
  }
  else if (!strcmp(inputfile,"Al_He3.dat")){
    sprintf(outfile,"%s/Al_He3.root",file_dir);
  }
  else if (!strcmp(inputfile,"Al_T2.dat")){
    sprintf(outfile,"%s/Al_T2.root",file_dir);
  }
  else if (!strcmp(inputfile,"Al_h2.dat")){
    sprintf(outfile,"%s/Al_h2.root",file_dir);
  }
  else if (!strcmp(inputfile,"Al_h22.dat")){
    sprintf(outfile,"%s/Al_h22.root",file_dir);
  }
  else {
    sprintf(outfile,"%s/temp.root",file_dir);
  }
  
  ifstream* dragon_data_file = new ifstream(inputfile);
  double run, num, flag;
  char tempc[900];
  int nfile = 0;
  if(dragon_data_file->is_open()==true){
    TChain* c1 = new TChain("tree");
    while (!dragon_data_file->eof()){
      *dragon_data_file >>run >> num >> flag;
      if(num==0){
	sprintf(tempc,"%s/tri_coin_%d.root",file_dir,(int)run);
      }
      else {
	sprintf(tempc,"%s/tri_coin_%d_%d.root",file_dir,(int)run,(int)num);
      }
      //cout << " " << run << " " << num << ": " << tempc << endl;
      cout << " Adding " << tempc << endl;
      c1->Add(tempc);
      nfile++;
    }
    cout << endl;
    cout << " Data list: " << inputfile << " (" 
	 << nfile << " files)" << endl;
    cout << " --> " << outfile << endl;
    cout << " Now merging....." << endl;
    c1->Merge(outfile);
    cout << endl;
    cout << " Done." << endl;
  }
  else return 88;
  
}
