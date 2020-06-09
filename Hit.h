//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 18 16:50:09 2018 by ROOT version 6.12/06
// from TTree DetTree/MyCRTHit
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef Hit_h
#define Hit_h

#define DEBUG_hit 0

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

class Hit : public TObject {

public :
   // Declaration of leaf types
   int            event;
   int            run;
   int            subRun;
   int            nHit;
   int            nFeb;
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
   float          totPe;
   float          rmsPe;
   vector<float>* hitPe;
   vector<int>*   modID;
   vector<int>*   stripID;
   vector<int>*   trackID;
   vector<int>*   pdg;

   // List of branches
   TBranch        *b_event;   
   TBranch        *b_run;
   TBranch        *b_subrun;
   TBranch        *b_nHit;   
   TBranch        *b_nFeb;
   TBranch        *b_x;
   TBranch        *b_y;
   TBranch        *b_z;
   TBranch        *b_xErr;
   TBranch        *b_yErr;
   TBranch        *b_zErr;
   TBranch        *b_t0;   
   TBranch        *b_t1;   
   TBranch        *b_trackID;   
   TBranch        *b_pdg;
   TBranch        *b_region;   
   TBranch        *b_subSys;   
   TBranch        *b_modID;
   TBranch        *b_stripID;
   TBranch        *b_totPe;
   TBranch        *b_rmsPe;
   TBranch        *b_hitPe;

   Hit();
   ~Hit() {};
   static Hit* giveThis();
   static Hit* giveThis(TTree* aTree, const std::string& option);
   static void releaseThis();
   void setupRead(const std::string& option = "");
   bool SetBranchAddresses();
   Hit(TTree* aTree, const std::string& option);
   TTree* GetInputTree();
   void GetHit(Long64_t entry);

   int         Event()            { return m_eve;            };
   int         Run()              { return m_run;            };
   int         SubRun()           { return m_subrun;         };
   int         NHit()             { return m_nhit;           };
   int         NFeb()             { return m_nfeb;           };
   float       X()                { return m_x;              };
   float       Y()                { return m_y;              };
   float       Z()                { return m_z;              };
   float       XErr()             { return m_xerr;           };
   float       YErr()             { return m_yerr;           };
   float       ZErr()             { return m_zerr;           };
   float       T0()               { return m_t0;             };
   float       T1()               { return m_t1;             };
   float       TotPe()            { return m_totpe;          };
   float       RMSPe()            { return m_rmspe;          };
   float       HitPe(int index)   { return m_hitpe[index];   };
   int         NTrack()           { return m_trackid.size(); };
   int         TrackID(int index) { return m_trackid[index]; };
   int         PDG(int index)     { return m_pdg[index];     };
   int         Region()           { return m_region;         };
   int         SubSys()           { return m_subsys;         };
   int         ModID(int index)   { return m_modid[index];   };
   int         StripID(int index) { return m_stripid[index]; };
   int         NMod()             { return m_modid.size();   };
   int         NStrip()           { return m_stripid.size(); };

private:
   TTree*        m_treeIn;
   int           m_eve;
   int           m_run;
   int           m_subrun;
   int           m_nhit;
   int           m_nfeb;
   float         m_x;
   float         m_y;
   float         m_z;
   float         m_xerr;
   float         m_yerr;
   float         m_zerr;
   float         m_t0;
   float         m_t1;
   float         m_totpe;
   float         m_rmspe;
   int           m_region;
   int           m_subsys;
   vector<float> m_hitpe;
   vector<int>   m_trackid;
   vector<int>   m_pdg;
   vector<int>   m_modid;
   vector<int>   m_stripid;

   static Hit* m_instance;
   ClassDef(Hit,0);
};

// initialisation of the Simulation poInt_ter
Hit* Hit::m_instance = 0;
Hit::Hit() {}
//to get a unique instance
inline Hit* Hit::giveThis()
{
     if(DEBUG_hit) cout << "debug::m_instance Hit " << m_instance << endl;
     if (0 == m_instance){
        cout << "Hit::giveThis error not constructed properly " << endl;
     }
     return m_instance;
}

// Open a TTree for reading
inline Hit* Hit::giveThis(TTree* aTree, const std::string& option)
{
     if(DEBUG_hit) cout << "debug::m_instance Hit " << m_instance << endl;
     if (0 == m_instance){
        m_instance = new Hit(aTree, option);
     }
     else{
        cout << "Hit::giveThis Warning " << aTree->GetTitle() << endl;
     }
     return m_instance;
}

// Delete unique instance
inline void Hit::releaseThis() {
     if ( m_instance != 0 ) {
        delete m_instance;
        m_instance = 0;
     }
}

// constructor for one TTree with read/write option
Hit::Hit(TTree* aTree, const std::string& option){
     if(option=="read"){
        m_treeIn = aTree;
        if(DEBUG_hit) cout << " ....before setupread.... " << endl;
        setupRead();
     }
}

// setup the input tree for reading
void Hit::setupRead(const std::string& option){
     if(DEBUG_hit) cout << " ....inside setupread.... " << endl;
     if(DEBUG_hit) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     if(!m_treeIn){
        cout << "setupRead error: m_treeIn undefined " << endl;
        exit(1);
     }
     if(SetBranchAddresses()){}
     else{
        cerr << "TBranch error.."<< endl;
     }
}

bool Hit::SetBranchAddresses()
{
     if(DEBUG_hit) cout << " ....inside setbranchaddresses.... " << endl;
     event=0;
     run=0;
     subRun=0;
     nHit = 0;
     nFeb = 0;
     x=0.;
     y=0.;
     z=0.;
     xErr=0.;
     yErr=0.;
     zErr=0.;
     t0=0.;
     t1=0.;
     region=0;
     subSys=0;
     totPe = 0.;
     rmsPe = 0.;
     hitPe=nullptr;
     modID=nullptr;
     stripID=nullptr;
     trackID=nullptr;
     pdg=nullptr;

     m_treeIn->SetMakeClass(1);
     m_treeIn->SetBranchAddress("event",   &event,   &b_event);
     m_treeIn->SetBranchAddress("run",     &run,     &b_run);
     m_treeIn->SetBranchAddress("subRun",  &subRun,  &b_subrun);
     m_treeIn->SetBranchAddress("nHit",    &nHit,    &b_nHit);
     m_treeIn->SetBranchAddress("nFeb",    &nFeb,   &b_nFeb);
     m_treeIn->SetBranchAddress("x"   ,    &x,       &b_x);
     m_treeIn->SetBranchAddress("y"   ,    &y,       &b_y);
     m_treeIn->SetBranchAddress("z"   ,    &z,       &b_z);
     m_treeIn->SetBranchAddress("xErr",    &xErr,    &b_xErr);
     m_treeIn->SetBranchAddress("yErr",    &yErr,    &b_yErr);
     m_treeIn->SetBranchAddress("zErr",    &zErr,    &b_zErr);
     m_treeIn->SetBranchAddress("t0",      &t0,      &b_t0);
     m_treeIn->SetBranchAddress("t1",      &t1,      &b_t1);
     m_treeIn->SetBranchAddress("totPe",   &totPe,   &b_totPe);
     m_treeIn->SetBranchAddress("rmsPe",   &rmsPe,   &b_rmsPe);
     m_treeIn->SetBranchAddress("hitPe",   &hitPe,   &b_hitPe);
     m_treeIn->SetBranchAddress("trackID", &trackID, &b_trackID);
     m_treeIn->SetBranchAddress("pdg",     &pdg,     &b_pdg);
     m_treeIn->SetBranchAddress("region",  &region,  &b_region);
     m_treeIn->SetBranchAddress("subSys",  &subSys,  &b_subSys);
     m_treeIn->SetBranchAddress("modID",   &modID,   &b_modID);
     m_treeIn->SetBranchAddress("stripID", &stripID, &b_stripID);

     if(DEBUG_hit) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     return true;
}

// accessors for the input TTree
TTree* Hit::GetInputTree() {return m_treeIn;};

// get all branch contents of input TTree for entry i
void Hit::GetHit(Long64_t entry)
{
     if (!m_treeIn) {
        cout << "Hit::getEntry error" << endl;
        exit(1);
     }
     m_treeIn->GetEntry(entry);
     if(DEBUG_hit) m_treeIn->Show(entry);

     m_hitpe.clear();
     m_trackid.clear();
     m_pdg.clear();
     m_modid.clear();
     m_stripid.clear();

     m_eve       = event;
     m_run       = run;
     m_subrun    = subRun;
     m_nhit      = nHit;
     m_nfeb      = nFeb;
     m_x         = x;
     m_y         = y;
     m_z         = z;
     m_xerr      = xErr;
     m_yerr      = yErr;
     m_zerr      = zErr;
     m_t0        = t0;
     m_t1        = t1;
     m_totpe     = totPe;
     m_rmspe     = rmsPe;
     m_region    = region;
     m_subsys    = subSys;

     if(trackID->size()!=pdg->size())
         cout << "WARNING: track and pdg arrays have different size!" << endl;
     for (size_t i=0; i<trackID->size(); i++) {
         m_trackid.push_back(trackID->at(i));
         m_pdg.push_back(pdg->at(i));
     }

     if(modID->size()!=stripID->size())
	 cout << "WARNING: modID and stripID are of different size!" << endl;
     if(modID->size()!=hitPe->size())
	 cout << "WARNING: modiD and hitPe are of different size!" << endl;
     for (size_t i=0; i<modID->size(); i++) {
         m_modid.push_back(modID->at(i));
         m_stripid.push_back(stripID->at(i));
	 m_hitpe.push_back(hitPe->at(i));
     }

     if(DEBUG_hit) m_treeIn->Print();
} 
#endif
