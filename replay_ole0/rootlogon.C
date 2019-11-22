//////////////////////////////////////////////////////////////////////////
//                                                                          
// rootlogon.C                                                            
//                                                                          
// Load Lib, paths and possible scripts to the analyzer upon start
//                                                                          
//////////////////////////////////////////////////////////////////////////
//	
// Author : Jin Huang <mailto:jinhuang@jlab.org>    Mar 2008
//          Chao Gu         Update for g2p use      Mar 2012
//          Kalyan Allada   Update for GMp expt.    Dec 2013
//	    Florian Hauenstein   Update for Tritum expt.    Feb 2017
//          Tyler Hague     1.6 compatability update Sep 2017
//
//////////////////////////////////////////////////////////////////////////

{
  //Used for setting paths. %s is necessary so that it can be used to Form the paths.
  //This allows changing the directory in only one place for ease of portability.
  const char* thisdir = gSystem->Getenv("PWD");
  TString replay_dir_prefix(thisdir ? thisdir : "."); replay_dir_prefix.Append("/%s");
  TString Arch(gSystem->GetBuildArch());
  TString Arch32("linux");
  TString Arch64("linuxx8664gcc");
  TString ArchOSX("macosx64");
  if(Arch==Arch32){
    printf("\nrootlogon.C: Loading Linux 32-bit Replay Core Library...\n"); 
  //  gSystem->Load(Form(replay_dir_prefix,"old/ReplayCore_C.so"));
  }

  else if(Arch==Arch64){
     printf("\nrootlogon.C: Loading Linux 64-bit Replay Core Library...\n"); 
     gSystem->Load(Form(replay_dir_prefix,"ReplayCore64_C.so"));
     gSystem->Load(Form(replay_dir_prefix,"libTriton/libTriton.so"));
  }

  else if(Arch==ArchOSX){
     printf("\nrootlogon.C: Loading macOS 64-bit Replay Core Library...\n"); 
     gSystem->Load(Form(replay_dir_prefix,"ReplayCore64_C.so"));
     gSystem->Load(Form(replay_dir_prefix,"libTriton/libTriton.dylib"));
  }

  //Load more libs here, if necessary. 
  //Make sure it's in path of $LD_LIBRARY_PATH

  printf("rootlogon.C: Adding include directories...\n");

  gSystem->AddIncludePath("-I$ANALYZER");
  gInterpreter->AddIncludePath("$ANALYZER/include");

  // gSystem->AddIncludePath("-Ilibraries/include");
  // gInterpreter->AddIncludePath("libraries/include");

  printf("rootlogon.C: Done!\n\n");
}
