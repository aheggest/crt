//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug 30 17:51:48 2019 by ROOT version 6.08/06
// from TTree GenTree/truth information from the generator
// found on file: /home/chris/Data/icarus_simulation/genie/genie_detsim.root
//////////////////////////////////////////////////////////

#ifndef Gen_h
#define Gen_h

#define DEBUG_gen 0

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream> 

class Gen : public TObject {

public :

   // Declaration of leaf types
   Int_t           event;
   Int_t           run;
   Int_t           subRun;
   Int_t           nGen;
   vector<int>     *trackID;
   vector<int>     *pdg;
   vector<vector<double> > *startXYZT;
   vector<vector<double> > *endXYZT;
   vector<vector<double> > *startPE;
   vector<vector<double> > *endPE;
   //UChar_t         hasNu;
   //Int_t           nuPDG;
   //Int_t           nuInt;
   //Int_t           ccnc;
   //Double_t        nuTheta;
   //Int_t           nuNucleon;
   //Int_t           nuTarget;
   //Int_t           nuLeptonID;
   //Int_t           nuReg;
   //UChar_t         nuFid;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_subRun;   //!
   TBranch        *b_nGen;   //!
   TBranch        *b_trackID;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_startXYZT;   //!
   TBranch        *b_endXYZT;   //!
   TBranch        *b_startPE;   //!
   TBranch        *b_endPE;   //!
   //TBranch        *b_hasNu;   //!
   //TBranch        *b_nuPDG;   //!
   //TBranch        *b_nuInt;   //!
   //TBranch        *b_ccnc;   //!
   //TBranch        *b_nuTheta;   //!
   //TBranch        *b_nuNucleon;   //!
   //TBranch        *b_nuTarget;   //!
   //TBranch        *b_nuLeptonID;
   //TBranch        *b_nuReg;
   //TBranch        *b_nuFid;


   Gen();
   ~Gen() {};

   static Gen* giveThis();
   static Gen* giveThis(TTree* aTree, const std::string& option);
   static void releaseThis();
   void setupRead(const std::string& option = "");
   bool SetBranchAddresses();
   Gen(TTree* aTree, const std::string& option);
   TTree* GetInputTree();
   void GetGen(Long64_t entry);


   Int_t    Event()          { return m_event;           }
   Int_t    Run()            { return m_run;             }
   Int_t    SubRun()         { return subRun;            }
   Int_t    NGen()           { return m_nGen;            }
   Int_t    TrackID(Int_t i) { return m_trackID[i];      }
   Int_t    PDG(Int_t i)     { return m_pdg[i];          }
   Double_t StartX(Int_t i)  { return m_startXYZT[i][0]; }
   Double_t StartY(Int_t i)  { return m_startXYZT[i][1]; }
   Double_t StartZ(Int_t i)  { return m_startXYZT[i][2]; }
   Double_t StartT(Int_t i)  { return m_startXYZT[i][3]; }
   Double_t EndX(Int_t i)    { return m_endXYZT[i][0];   }
   Double_t EndY(Int_t i)    { return m_endXYZT[i][1];   }
   Double_t EndZ(Int_t i)    { return m_endXYZT[i][2];   }
   Double_t EndT(Int_t i)    { return m_endXYZT[i][3];   }
   Double_t StartPx(Int_t i) { return m_startPE[i][0];   }
   Double_t StartPy(Int_t i) { return m_startPE[i][1];   }
   Double_t StartPz(Int_t i) { return m_startPE[i][2];   }
   Double_t StartE(Int_t i)  { return m_startPE[i][3];   }
   Double_t EndPx(Int_t i)   { return m_endPE[i][0];     }
   Double_t EndPy(Int_t i)   { return m_endPE[i][1];     }
   Double_t EndPz(Int_t i)   { return m_endPE[i][2];     }
   Double_t EndE(Int_t i)    { return m_endPE[i][3];     }
   //Bool_t   HasNu()          { return m_hasNu;           }
   //Int_t    NuPDG()          { return m_nuPDG;           }
   //Int_t    NuInt()          { return m_nuInt;           }
   //Int_t    CCNC()           { return m_ccnc;            }
   //Double_t NuTheta()        { return m_nuTheta;         }
   //Int_t    NuNucleon()      { return m_nuNucleon;       }
   //Int_t    NuTarget()       { return m_nuTarget;        }
   //Int_t    NuLeptonID()     { return m_nuLeptonID;      }
   //Int_t    NuReg()          { return m_nuReg;           }
   //Bool_t   NuFid()          { return m_nuFid;           }

private:
   TTree*   m_treeIn;

   Int_t                   m_event;
   Int_t                   m_run;
   Int_t                   m_subRun;
   Int_t                   m_nGen;
   vector<int>             m_trackID;
   vector<int>             m_pdg;
   vector<vector<double> > m_startXYZT;
   vector<vector<double> > m_endXYZT;
   vector<vector<double> > m_startPE;
   vector<vector<double> > m_endPE;
   //Bool_t                  m_hasNu;
   //Int_t                   m_nuPDG;
   //Int_t                   m_nuInt;
   //Int_t                   m_ccnc;
   //Double_t                m_nuTheta;
   //Int_t                   m_nuNucleon;
   //Int_t                   m_nuTarget;
   //Int_t                   m_nuLeptonID;
   //Int_t                   m_nuReg;
   //Bool_t                  m_nuFid;

   static Gen* m_instance;
   ClassDef(Gen,0);
};

// initialisation of the Gen pointer
Gen* Gen::m_instance = 0;
Gen::Gen() {}

//to get a unique instance
inline Gen* Gen::giveThis()
{
     if(DEBUG_gen) cout << "debug::m_instance TrueHit " << m_instance << endl;
     if (0 == m_instance){
        cout << "Gen::giveThis error not constructed properly " << endl;
     }
     return m_instance;
}

// Open a TTree for reading
inline Gen* Gen::giveThis(TTree* aTree, const std::string& option)
{
     if(DEBUG_gen) cout << "debug::m_instance Gen " << m_instance << endl;
     if (0 == m_instance){
        m_instance = new Gen(aTree, option);
     }
     else{
        cout << "Gen::giveThis Warning " << aTree->GetTitle() << endl;
     }
     return m_instance;
}

// Delete unique instance
inline void Gen::releaseThis() {
     if ( m_instance != 0 ) {
        delete m_instance;
        m_instance = 0;
     }
}

// constructor for one TTree with read/write option
Gen::Gen(TTree* aTree, const std::string& option){
     if(option=="read"){
        m_treeIn = aTree;
        if(DEBUG_gen) cout << " ....before setupread.... " << endl;
        setupRead();
     }
}

// setup the input tree for reading
void Gen::setupRead(const std::string& option){
     if(DEBUG_gen) cout << " ....inside setupread.... " << endl;
     if(DEBUG_gen) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     if(!m_treeIn){
        cout << "setupRead error: m_treeIn undefined " << endl;
        exit(1);
     }
     if(SetBranchAddresses()){}
     else{
        cerr << "TBranch error.."<< endl;
     }
}

bool Gen::SetBranchAddresses()
{
   if(DEBUG_gen) cout << " ....inside setbranchaddresses.... " << endl;

   // Set object pointer
   trackID = 0;
   pdg = 0;
   startXYZT = 0;
   endXYZT = 0;
   startPE = 0;
   endPE = 0;

   // Set branch addresses and branch pointers
   m_treeIn->SetBranchAddress("event", &event, &b_event);
   m_treeIn->SetBranchAddress("run", &run, &b_run);
   m_treeIn->SetBranchAddress("subRun", &subRun, &b_subRun);
   m_treeIn->SetBranchAddress("nGen", &nGen, &b_nGen);
   m_treeIn->SetBranchAddress("trackID", &trackID, &b_trackID);
   m_treeIn->SetBranchAddress("pdg", &pdg, &b_pdg);
   m_treeIn->SetBranchAddress("startXYZT", &startXYZT, &b_startXYZT);
   m_treeIn->SetBranchAddress("endXYZT", &endXYZT, &b_endXYZT);
   m_treeIn->SetBranchAddress("startPE", &startPE, &b_startPE);
   m_treeIn->SetBranchAddress("endPE", &endPE, &b_endPE);
   //m_treeIn->SetBranchAddress("hasNu", &hasNu, &b_hasNu);
   //m_treeIn->SetBranchAddress("nuPDG", &nuPDG, &b_nuPDG);
   //m_treeIn->SetBranchAddress("nuInt", &nuInt, &b_nuInt);
   //m_treeIn->SetBranchAddress("ccnc", &ccnc, &b_ccnc);
   //m_treeIn->SetBranchAddress("nuTheta", &nuTheta, &b_nuTheta);
   //m_treeIn->SetBranchAddress("nuNucleon", &nuNucleon, &b_nuNucleon);
   //m_treeIn->SetBranchAddress("nuTarget", &nuTarget, &b_nuTarget);
   //m_treeIn->SetBranchAddress("nuLeptonID", &nuLeptonID, &b_nuLeptonID);
   //m_treeIn->SetBranchAddress("nuReg", &nuReg, &b_nuReg);
   //m_treeIn->SetBranchAddress("nuFid", &nuFid, &b_nuFid);

   if(DEBUG_gen) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
   return true;
}

// accessors for the input TTree
TTree* Gen::GetInputTree() {return m_treeIn;};

// get all branch contents of input TTree for entry i
void Gen::GetGen(Long64_t entry)
{
     if (!m_treeIn) {
        cout << "Gen::getEntry error" << endl;
        exit(1);
     }
     m_treeIn->GetEntry(entry);
     if(DEBUG_gen) m_treeIn->Show(entry);

     m_event      = event;         
     m_run        = run;           
     m_subRun     = subRun;
     m_nGen       = nGen;
     //m_hasNu      = hasNu;
     //m_nuPDG      = nuPDG;
     //m_nuInt      = nuInt;
     //m_ccnc       = ccnc;
     //m_nuTheta    = nuTheta;
     //m_nuNucleon  = nuNucleon;
     //m_nuTarget   = nuTarget;
     //m_nuLeptonID = nuLeptonID;
     //m_nuReg      = nuReg;
     //m_nuFid      = nuFid;

     m_trackID.clear(); 
     m_pdg.clear();
     m_startXYZT.clear();
     m_endXYZT.clear();
     m_startPE.clear();
     m_endPE.clear();

     for(Int_t i=0; i<nGen; i++) {
         m_trackID.push_back(trackID->at(i));
         m_pdg.push_back(pdg->at(i));
         m_startXYZT.push_back(startXYZT->at(i));
	 m_endXYZT.push_back(endXYZT->at(i));
         m_startPE.push_back(startPE->at(i));
         m_endPE.push_back(endPE->at(i));
     }

     if(DEBUG_gen) m_treeIn->Print();
}
#endif // #ifdef Gen_h
 
