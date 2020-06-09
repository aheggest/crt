//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 18 16:50:09 2018 by ROOT version 6.12/06
// from TTree DetTree/MyCRTHit
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef TrueHit_h
#define TrueHit_h

#define DEBUG_truehit 0

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TNamed.h>
#include <TObject.h>

#include <vector>
#include <iostream>
#include <utility>
#include <limits>

class TrueHit : public TObject {

public :
   // Declaration of leaf types
   int            event;
   int            run;
   int            subRun;
   int            nHit;
   float          x;
   float          y;
   float          z;
   float          xErr;
   float          yErr;
   float          zErr;
   float          t0;
   float          t1;
   int            region;
   int            subSys;
   int            nFeb;
   vector<int>*   trackID;
   vector<int>*   pdg;
   vector<int>*   modID;
   vector<int>*   stripID;
   float          totPe;
   float          rmsPe;
   vector<float>* hitPe;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;
   TBranch        *b_subrun;
   TBranch        *b_nhit;   //!
   TBranch        *b_x; 
   TBranch        *b_y;
   TBranch        *b_z;
   TBranch        *b_xErr;
   TBranch        *b_yErr;
   TBranch        *b_zErr;
   TBranch        *b_t0;   //!
   TBranch        *b_t1;
   TBranch        *b_trackID;   //!
   TBranch        *b_region;   //!
   TBranch        *b_subSys;   //!
   TBranch        *b_nfeb;
   TBranch        *b_pdg;
   TBranch        *b_modID;
   TBranch        *b_stripID;
   TBranch        *b_totPe;
   TBranch        *b_rmsPe;
   TBranch        *b_hitPe;

   TrueHit();
   ~TrueHit() {};
   static TrueHit* giveThis();
   static TrueHit* giveThis(TTree* aTree, const std::string& option);
   static void releaseThis();
   void setupRead(const std::string& option = "");
   bool SetBranchAddresses();
   TrueHit(TTree* aTree, const std::string& option);
   TTree* GetInputTree();
   void GetHit(Long64_t entry);

   int   Event()        { return m_eve;        };
   int   Run()          { return m_run;        };
   int   SubRun()       { return m_subrun;     };
   int   NHit()         { return m_nhit;       };
   float X()            { return m_x;          };
   float Y()            { return m_y;          };
   float Z()            { return m_z;          };
   float XErr()         { return m_xerr;       };
   float YErr()         { return m_yerr;       };
   float ZErr()         { return m_zerr;       };
   float T0()           { return m_t0;         };
   float T1()           { return m_t1;         };
   int   TrackID(int i) { return m_trackid[i]; };
   int   Region()       { return m_region;     };
   int   SubSys()       { return m_subsys;     };
   int   NFeb()         { return m_nfeb;       };
   int   PDG(int i)     { return m_pdg[i];     };
   int   ModID(int i)   { return m_modid[i];   };
   int   StripID(int i) { return m_stripid[i]; };
   float TotPE()        { return m_totpe;      };
   float RmsPE()        { return m_rmspe;      };
   float HitPE(int i)   { return m_hitpe[i];   };
   int   NStrip()       { return m_stripid.size(); };

private:
   TTree*   m_treeIn;

   int           m_eve;
   int           m_run;
   int           m_subrun;
   int           m_nhit;
   float         m_x;
   float         m_y;
   float         m_z;
   float         m_xerr;
   float         m_yerr;
   float         m_zerr;
   int           m_t0;
   int           m_t1;
   int           m_region;
   int           m_subsys;
   int           m_nfeb;
   vector<int>   m_trackid;
   vector<int>   m_pdg; 
   vector<int>   m_modid;
   vector<int>   m_stripid;
   float         m_totpe;
   float         m_rmspe;
   vector<float> m_hitpe;

   static TrueHit* m_instance;
   ClassDef(TrueHit,0);
};

// initialisation of the Simulation poInt_ter
TrueHit* TrueHit::m_instance = 0;
TrueHit::TrueHit() {}
//to get a unique instance
inline TrueHit* TrueHit::giveThis()
{
     if(DEBUG_truehit) cout << "debug::m_instance TrueHit " << m_instance << endl;
     if (0 == m_instance){
        cout << "TrueHit::giveThis error not constructed properly " << endl;
     }
     return m_instance;
}

// Open a TTree for reading
inline TrueHit* TrueHit::giveThis(TTree* aTree, const std::string& option)
{
     if(DEBUG_truehit) cout << "debug::m_instance TrueHit " << m_instance << endl;
     if (0 == m_instance){
        m_instance = new TrueHit(aTree, option);
     }
     else{
        cout << "TrueHit::giveThis Warning " << aTree->GetTitle() << endl;
     }
     return m_instance;
}

// Delete unique instance
inline void TrueHit::releaseThis() {
     if ( m_instance != 0 ) {
        delete m_instance;
        m_instance = 0;
     }
}

// constructor for one TTree with read/write option
TrueHit::TrueHit(TTree* aTree, const std::string& option){
     if(option=="read"){
        m_treeIn = aTree;
        if(DEBUG_truehit) cout << " ....before setupread.... " << endl;
        setupRead();
     }
}

// setup the input tree for reading
void TrueHit::setupRead(const std::string& option){
     if(DEBUG_truehit) cout << " ....inside setupread.... " << endl;
     if(DEBUG_truehit) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     if(!m_treeIn){
        cout << "setupRead error: m_treeIn undefined " << endl;
        exit(1);
     }
     if(SetBranchAddresses()){}
     else{
        cerr << "TBranch error.."<< endl;
     }
}

bool TrueHit::SetBranchAddresses()
{
     if(DEBUG_truehit) cout << " ....inside setbranchaddresses.... " << endl;
     trackID = nullptr;
     pdg     = nullptr;
     hitPe   = nullptr;
     modID   = nullptr;
     stripID = nullptr;

     m_treeIn->SetMakeClass(1);
     m_treeIn->SetBranchAddress("event",   &event,   &b_event);
     m_treeIn->SetBranchAddress("run",     &run,     &b_run);
     m_treeIn->SetBranchAddress("subRun",  &subRun,  &b_subrun);
     m_treeIn->SetBranchAddress("nHit",    &nHit,    &b_nhit);
     m_treeIn->SetBranchAddress("x"   ,    &x,       &b_x);
     m_treeIn->SetBranchAddress("y"   ,    &y,       &b_y);
     m_treeIn->SetBranchAddress("z"   ,    &z,       &b_z);
     m_treeIn->SetBranchAddress("xErr",    &xErr,    &b_xErr);
     m_treeIn->SetBranchAddress("yErr",    &yErr,    &b_yErr);
     m_treeIn->SetBranchAddress("zErr",    &zErr,    &b_zErr);
     m_treeIn->SetBranchAddress("t0",      &t0,      &b_t0);
     m_treeIn->SetBranchAddress("t1",      &t1,      &b_t1);
     m_treeIn->SetBranchAddress("trackID", &trackID, &b_trackID);
     m_treeIn->SetBranchAddress("pdg",     &pdg,     &b_pdg);
     m_treeIn->SetBranchAddress("region",  &region,  &b_region);
     m_treeIn->SetBranchAddress("subSys",  &subSys,  &b_subSys);
     m_treeIn->SetBranchAddress("nFeb",    &nFeb,    &b_nfeb);
     m_treeIn->SetBranchAddress("modID",   &modID,   &b_modID);
     m_treeIn->SetBranchAddress("stripID", &stripID, &b_stripID);
     m_treeIn->SetBranchAddress("totPe",   &totPe,   &b_totPe);
     m_treeIn->SetBranchAddress("rmsPe",   &rmsPe,   &b_rmsPe);
     m_treeIn->SetBranchAddress("hitPe",   &hitPe,   &b_hitPe);

     if(DEBUG_truehit) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     return true;
}

// accessors for the input TTree
TTree* TrueHit::GetInputTree() {return m_treeIn;};

// get all branch contents of input TTree for entry i
void TrueHit::GetHit(Long64_t entry)
{
     if (!m_treeIn) {
        cout << "TrueHit::getEntry error" << endl;
        exit(1);
     }
     m_treeIn->GetEntry(entry);
     if(DEBUG_truehit) m_treeIn->Show(entry);

     m_eve       = event;
     m_run       = run;
     m_subrun    = subRun;
     m_nhit      = nHit;
     m_x         = x; 
     m_y         = y;
     m_z         = z;
     m_xerr      = xErr;
     m_yerr      = yErr;
     m_zerr      = zErr;
     m_region    = region;
     m_subsys    = subSys;
     m_nfeb      = nFeb;
     m_t0        = t0;
     m_t1        = t1;
     m_totpe     = totPe;
     m_rmspe     = rmsPe;

     m_trackid.clear();
     m_pdg.clear();
     m_hitpe.clear();
     m_modid.clear();
     m_stripid.clear();

     if(trackID->size()!=pdg->size())
	     cout << "trackID - pdg size mismatch!" << endl;

     for(size_t i=0; i<trackID->size(); i++){
         m_trackid.push_back(trackID->at(i));
	 m_pdg.push_back(pdg->at(i));
     }

     if(hitPe->size()!=modID->size())
	 cout << "WARNING: hitPe - modID size mismatch!" << endl;
     if(hitPe->size()!=stripID->size())
	 cout << "WARNING: hitPe - stripID size mismatch!" << endl;
     for(size_t i=0; i<hitPe->size(); i++){
         m_hitpe.push_back(hitPe->at(i));
	 m_modid.push_back(modID->at(i));
	 m_stripid.push_back(stripID->at(i));
     }

     if(DEBUG_truehit) m_treeIn->Print();
}
#endif
