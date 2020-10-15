//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 18 16:50:09 2018 by ROOT version 6.12/06
// from TTree DetTree/MyCRTDetSim
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef DetSim_h
#define DetSim_h

#define DEBUG_det 0

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

class DetSim : public TObject {

public :
   // Declaration of leaf types
   Int_t           event;
   Int_t           run;
   Int_t           subRun;
   Int_t           entry;
   Int_t           nAbove;
   Int_t           t0;
   Int_t           t1;
   Int_t           adc[64];
   Int_t           maxAdc;
   Int_t           maxChan;
   vector<int>*    trackID;
   vector<int>*    pdg;
   Int_t           mac5;
   Int_t           region;
   Int_t           subSys;

   // List of branches
   TBranch        *b_event;
   TBranch        *b_run;
   TBranch        *b_subrun;
   TBranch        *b_entry;
   TBranch        *b_nAbove; 
   TBranch        *b_t0;  
   TBranch        *b_t1; 
   TBranch        *b_adc;
   TBranch        *b_maxChan;
   TBranch        *b_maxAdc;
   TBranch        *b_trackID;
   TBranch        *b_pdg;
   TBranch        *b_Mac5; 
   TBranch        *b_Region; 
   TBranch        *b_subSys;

   DetSim();
   ~DetSim() {};
   static DetSim* giveThis();
   static DetSim* giveThis(TChain* aTree, const std::string& option);
   static void releaseThis();
   void setupRead(const std::string& option = "");
   bool SetBranchAddresses();
   DetSim(TChain* aTree, const std::string& option);
   TChain* GetInputTree();
   void GetDet(Long64_t entry);

   Int_t Event()                { return m_eve;          };
   Int_t Run()                  { return m_run;          };
   Int_t SubRun()               { return m_subrun;       };
   Int_t NAbove()               { return m_nabove;       };
   Double_t T0()                { return m_t0;           };
   Double_t T1()                { return m_t1;           };
   Int_t ADC(Int_t index)       { return m_adc[index];   };
   Int_t MaxAdc()               { return m_maxadc;       };
   Int_t MaxChan()              { return m_maxchan;      };
   Int_t TrackID(Int_t index)   { return m_trkid[index]; };
   Int_t PDG(Int_t index)       { return m_pdg[index];   };
   Int_t NTrack()               { return m_trkid.size(); };
   Int_t Entry()                { return m_entry;        };
   Int_t Mac5()                 { return m_mac5;         };
   Int_t Region()               { return m_region;       };
   Int_t SubSys()               { return m_subsys;       };

private:
   TChain* m_treeIn;
   Int_t  m_eve;
   Int_t  m_run;
   Int_t  m_subrun;
   Int_t  m_nabove;
   Int_t  m_t0;
   Int_t  m_t1;
   Int_t  m_adc[64];
   Int_t  m_maxadc;
   Int_t  m_maxchan;
   vector<Int_t>  m_trkid;
   vector<Int_t>  m_pdg;
   Int_t  m_entry;
   Int_t  m_mac5;
   Int_t  m_region;
   Int_t  m_subsys;
   static DetSim* m_instance;
   ClassDef(DetSim,0);
};

// initialisation of the Simulation poInt_ter
DetSim* DetSim::m_instance = 0;
DetSim::DetSim() {}
//to get a unique instance
inline DetSim* DetSim::giveThis()
{
     if(DEBUG_det) cout << "debug::m_instance DetSim " << m_instance << endl;
     if (0 == m_instance){
        cout << "DetSim::giveThis error not constructed properly " << endl;
     }
     return m_instance;
}

// Open a TTree for reading
inline DetSim* DetSim::giveThis(TChain* aTree, const std::string& option)
{
     if(DEBUG_det) cout << "debug::m_instance DetSim " << m_instance << endl;
     if (0 == m_instance){
        m_instance = new DetSim(aTree, option);
     }
     else{
        cout << "DetSim::giveThis Warning " << aTree->GetTitle() << endl;
     }
     return m_instance;
}

// Delete unique instance
inline void DetSim::releaseThis() {
     if ( m_instance != 0 ) {
        delete m_instance;
        m_instance = 0;
     }
}

// constructor for one TTree with read/write option
DetSim::DetSim(TChain* aTree, const std::string& option){
     if(option=="read"){
        m_treeIn = aTree;
        if(DEBUG_det) cout << " ....before setupread.... " << endl;
        setupRead();
     }
}

// setup the input tree for reading
void DetSim::setupRead(const std::string& option){
     if(DEBUG_det) cout << " ....inside setupread.... " << endl;
     if(DEBUG_det) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     if(!m_treeIn){
        cout << "setupRead error: m_treeIn undefined " << endl;
        exit(1);
     }
     if(SetBranchAddresses()){}
     else{
        cerr << "TBranch error.."<< endl;
     }
}

bool DetSim::SetBranchAddresses()
{
     if(DEBUG_det) cout << " ....inside setbranchaddresses.... " << endl;
     event=0;
     run=0;
     subRun=0;
     entry=0;
     trackID=nullptr;
     pdg=nullptr;
     nAbove=0; 
     t0=0; 
     t1=0; 
     maxAdc=0;
     maxChan=0;
     subSys=0;
     region=0;

     m_treeIn->SetBranchAddress("event", &event, &b_event);
     m_treeIn->SetBranchAddress("run",   &run,   &b_run);
     m_treeIn->SetBranchAddress("subRun",&subRun,&b_subrun);
     m_treeIn->SetBranchAddress("nAbove", &nAbove, &b_nAbove);
     m_treeIn->SetBranchAddress("t0", &t0, &b_t0);
     m_treeIn->SetBranchAddress("t1", &t1, &b_t1);
     m_treeIn->SetBranchAddress("adc", adc, &b_adc);
     m_treeIn->SetBranchAddress("trackID", &trackID, &b_trackID);
     m_treeIn->SetBranchAddress("detPDG", &pdg, &b_pdg);
     m_treeIn->SetBranchAddress("entry", &entry, &b_entry);
     m_treeIn->SetBranchAddress("mac5", &mac5, &b_Mac5);
     m_treeIn->SetBranchAddress("region", &region, &b_Region);
     m_treeIn->SetBranchAddress("subSys", &subSys, &b_subSys);
     m_treeIn->SetBranchAddress("maxAdc", &maxAdc, &b_maxAdc);
     m_treeIn->SetBranchAddress("maxChan",&maxChan,&b_maxChan);
     if(DEBUG_det) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     return true;
}

// accessors for the input TTree
TChain* DetSim::GetInputTree() {return m_treeIn;};

// get all branch contents of input TTree for entry i
void DetSim::GetDet(Long64_t entry)
{
     if (!m_treeIn) {
        cout << "DetSim::getEntry error" << endl;
        exit(1);
     }
     m_treeIn->GetEntry(entry);
     if(DEBUG_det) m_treeIn->Show(entry);

     m_trkid.clear();
     m_pdg.clear();

     m_eve         = event;
     m_run         = run;
     m_subrun      = subRun;
     m_nabove      = nAbove;
     m_entry       = entry;
     m_mac5        = mac5;
     m_region      = region;
     m_subsys      = subSys;
     m_maxadc      = maxAdc;
     m_maxchan     = maxChan;
     m_t0          = t0;
     m_t1          = t1;
     for (size_t i=0; i<trackID->size(); i++) {
         m_trkid.push_back(trackID->at(i));
         m_pdg.push_back(pdg->at(i));
     }
     for(size_t ch=0; ch<64; ch++)
         m_adc[ch] = adc[ch];

     if(DEBUG_det) m_treeIn->Print();
}
#endif
