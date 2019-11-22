#ifndef Podd_ReplayCore
#define Podd_ReplayCore

void ReplayCore(
		Int_t runnumber=0,            //run #
		Int_t all=0,                  //-1=replay all;0=ask for a number
		Int_t DefReplayNum=-1,        //default replay event num
		const char* OutFileFormat="%s/tritium_%d.root", //output file format
		const char* OutDefineFile="HRS.odef",       //out define
		const char* CutDefineFile="HRS.cdef",       //cut define
		Bool_t EnableScalar=false,                    //Enable Scalar?
		Bool_t EnableHelicity=false,                  //Enable Helicity?
		Int_t FirstEventNum=0,         //First Event To Replay
		Bool_t QuietRun = kFALSE      //whether not ask question?
   );
#endif
