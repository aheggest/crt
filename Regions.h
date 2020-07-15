//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 18 16:50:09 2018 by ROOT version 6.12/06
// from TTree RegionTree
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef Regions_h
#define Regions_h

#define DEBUG_region 0

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

class Regions : public TObject {

public :
   // Declaration of leaf types
   Int_t           event;
   Int_t           run;
   Int_t           subRun;
   Int_t           nReg;
   Int_t           fiducial;
   Int_t           active;
   Int_t           inactive;
   Int_t           crts;
   vector<int>*     regions;
   Int_t           pdg;
   Int_t           trackID;
   vector<double>*  eDep;
   vector<double>*  dL;
   vector<vector<double>>*  entryXYZT;
   vector<vector<double>>*  exitXYZT;
   vector<vector<double>>*  entrySlope;
   vector<vector<double>>*  exitSlope;
   vector<vector<double>>*  entryPE;
   vector<vector<double>>*  exitPE;
   vector<double>*          distToOpDet;//Enter;
   vector<vector<double>>*  opDetXYZT;//Enter;
   vector<int>*             opDetID;//Enter;
   //vector<double>*          distToOpDetExit;
   //vector<vector<double>>*  opDetXYZTExit;
   //vector<int>*             opDetIDExit;

   // List of branches
   TBranch        *b_event;  
   TBranch        *b_run;
   TBranch        *b_subrun;
   TBranch        *b_nReg;  
   TBranch        *b_fiducial;
   TBranch        *b_active;
   TBranch        *b_inactive;
   TBranch        *b_crts;
   TBranch        *b_regions;
   TBranch        *b_pdg;
   TBranch        *b_trackID;
   TBranch        *b_eDep;
   TBranch        *b_dL;
   TBranch        *b_entryXYZT; 
   TBranch        *b_exitXYZT; 
   TBranch        *b_entrySlope;
   TBranch        *b_exitSlope;
   TBranch        *b_entryPE; 
   TBranch        *b_exitPE; 
   TBranch        *b_distToOpDet;//Enter;
   TBranch        *b_opDetXYZT;//Enter;
   TBranch        *b_opDetID;//Enter;
   //TBranch        *b_distToOpDetExit;
   //TBranch        *b_opDetXYZTExit;
   //TBranch        *b_opDetIDExit;

   Regions();
   ~Regions() {};
   static Regions* giveThis();
   static Regions* giveThis(TChain* aTree, const std::string& option);
   static void releaseThis();
   void setupRead(const std::string& option = "");
   bool SetBranchAddresses();
   Regions(TChain* aTree, const std::string& option);
   TChain* GetInputTree();
   void GetRegion(Long64_t entry);

   //accessors
   Int_t    Event()                       { return m_eve;                };
   Int_t    Run()                         { return m_run;                };
   Int_t    SubRun()                      { return m_subrun;             };
   Int_t    NReg()                        { return m_nreg;               };
   Int_t    Fiducial()                    { return m_fiducial;           };
   Int_t    Active()                      { return m_active;             };
   Int_t    Inactive()                    { return m_inactive;           };
   Int_t    CRTs()                        { return m_crts;               };
   Double_t EDep(Int_t index)             { return m_edep[index];        };
   Double_t DL(Int_t index)               { return m_dl[index];          };
   Double_t EntryX(Int_t index)           { return m_entryx[index];      };
   Double_t EntryY(Int_t index)           { return m_entryy[index];      };
   Double_t EntryZ(Int_t index)           { return m_entryz[index];      };
   Double_t EntryT(Int_t index)           { return m_entryt[index];      };
   Double_t ExitX(Int_t index)            { return m_exitx[index];       };
   Double_t ExitY(Int_t index)            { return m_exity[index];       };
   Double_t ExitZ(Int_t index)            { return m_exitz[index];       };
   Double_t ExitT(Int_t index)            { return m_exitt[index];       };
   Double_t EntrySlopeX(Int_t index)      { return m_entryslopex[index]; };
   Double_t EntrySlopeY(Int_t index)      { return m_entryslopey[index]; };
   Double_t EntrySlopeZ(Int_t index)      { return m_entryslopez[index]; };
   Double_t ExitSlopeX(Int_t index)       { return m_exitslopex[index];  };
   Double_t ExitSlopeY(Int_t index)       { return m_exitslopey[index];  };
   Double_t ExitSlopeZ(Int_t index)       { return m_exitslopez[index];  };
   Double_t EntryPx(Int_t index)          { return m_entrypx[index];     };
   Double_t EntryPy(Int_t index)          { return m_entrypy[index];     };
   Double_t EntryPz(Int_t index)          { return m_entrypz[index];     };
   Double_t EntryE(Int_t index)           { return m_entrye[index];      };
   Double_t ExitPx(Int_t index)           { return m_exitpx[index];      };
   Double_t ExitPy(Int_t index)           { return m_exitpy[index];      };
   Double_t ExitPz(Int_t index)           { return m_exitpz[index];      };
   Double_t ExitE(Int_t index)            { return m_exite[index];       };
   Int_t    Region(Int_t index)           { return m_regions[index];     };
   Int_t    TrackID()                     { return m_trackid;            };
   Int_t    PDG()                         { return m_pdg;                };
   Double_t DistToOpDet(Int_t index) { return m_disttoopdet[index]; };
   Double_t OpDetX(Int_t index)      { return m_opdetx[index];      };
   Double_t OpDetY(Int_t index)      { return m_opdety[index];      };
   Double_t OpDetZ(Int_t index)      { return m_opdetz[index];      };
   Double_t OpDetT(Int_t index)      { return m_opdett[index];      };
   Int_t    OpDetID(Int_t index)     { return m_opdetid[index];     };
   //Double_t DistToEnterOpDet(Int_t index) { return m_disttoopdetenter[index]; };
   //Double_t OpDetEnterX(Int_t index)      { return m_opdetenterx[index];      };
   //Double_t OpDetEnterY(Int_t index)      { return m_opdetentery[index];      };
   //Double_t OpDetEnterZ(Int_t index)      { return m_opdetenterz[index];      };
   //Double_t OpDetEnterT(Int_t index)      { return m_opdetentert[index];      };
   //Int_t    OpDetEnterID(Int_t index)     { return m_opdetidenter[index];     };
   //Double_t DistToExitOpDet(Int_t index)  { return m_disttoopdetexit[index]; };
   //Double_t OpDetExitX(Int_t index)       { return m_opdetexitx[index];      };
   //Double_t OpDetExitY(Int_t index)       { return m_opdetexity[index];      };
   //Double_t OpDetExitZ(Int_t index)       { return m_opdetexitz[index];      };
   //Double_t OpDetExitT(Int_t index)       { return m_opdetexitt[index];      };
   //Int_t    OpDetExitID(Int_t index)      { return m_opdetidexit[index];     };

private:
   TChain*   m_treeIn;
   Int_t    m_eve;
   Int_t    m_run;
   Int_t    m_subrun;
   Int_t    m_nreg;
   Int_t    m_fiducial;
   Int_t    m_active;
   Int_t    m_inactive;
   Int_t    m_crts;
   vector<double> m_entryx; 
   vector<double> m_entryy; 
   vector<double> m_entryz; 
   vector<double> m_entryt; 
   vector<double> m_exitx; 
   vector<double> m_exity; 
   vector<double> m_exitz; 
   vector<double> m_exitt; 
   vector<double> m_entrypx; 
   vector<double> m_entrypy; 
   vector<double> m_entrypz; 
   vector<double> m_entrye; 
   vector<double> m_exitpx; 
   vector<double> m_exitpy; 
   vector<double> m_exitpz; 
   vector<double> m_exite; 
   vector<double> m_entryslopex; 
   vector<double> m_entryslopey; 
   vector<double> m_entryslopez; 
   vector<double> m_exitslopex; 
   vector<double> m_exitslopey; 
   vector<double> m_exitslopez; 
   vector<int> m_regions; 
   Int_t    m_trackid;
   Int_t    m_pdg;
   vector<double> m_edep; 
   vector<double> m_dl; 
   vector<double> m_disttoopdet;//enter; 
   vector<double> m_opdetx;//enterx; 
   vector<double> m_opdety;//entery; 
   vector<double> m_opdetz;//enterz; 
   vector<double> m_opdett;//entert; 
   vector<int>    m_opdetid;//enter;
   //vector<double> m_disttoopdetexit; 
   //vector<double> m_opdetexitx;
   //vector<double> m_opdetexity;
   //vector<double> m_opdetexitz;
   //vector<double> m_opdetexitt;
   //vector<int>    m_opdetidexit;
   
   static Regions* m_instance;
   ClassDef(Regions,0);
};

// initialisation of the Simulation poInt_ter
Regions* Regions::m_instance = 0;
Regions::Regions() {}
//to get a unique instance
inline Regions* Regions::giveThis()
{
     if(DEBUG_region) cout << "debug::m_instance Regions " << m_instance << endl;
     if (0 == m_instance){
        cout << "Region::giveThis error not constructed properly " << endl;
     }
     return m_instance;
}

// Open a TTree for reading
inline Regions* Regions::giveThis(TChain* aTree, const std::string& option)
{
     if(DEBUG_region) cout << "debug::m_instance Regions " << m_instance << endl;
     if (0 == m_instance){
        m_instance = new Regions(aTree, option);
     }
     else{
        cout << "Region::giveThis Warning " << aTree->GetTitle() << endl;
     }
     return m_instance;
}

// Delete unique instance
inline void Regions::releaseThis() {
     if ( m_instance != 0 ) {
        delete m_instance;
        m_instance = 0;
     }
}

// constructor for one TTree with read/write option
Regions::Regions(TChain* aTree, const std::string& option){
     if(option=="read"){
        m_treeIn = aTree;
        if(DEBUG_region) cout << " ....before setupread.... " << endl;
        setupRead();
     }
}

// setup the input tree for reading
void Regions::setupRead(const std::string& option){
     if(DEBUG_region) cout << " ....inside setupread.... " << endl;
     if(DEBUG_region) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     if(!m_treeIn){
        cout << "setupRead error: m_treeIn undefined " << endl;
        exit(1);
     }
     if(SetBranchAddresses()){}
     else{
        cerr << "TBranch error.."<< endl;
     }
}

bool Regions::SetBranchAddresses()
{
     if(DEBUG_region) cout << " ....inside setbranchaddresses.... " << endl;
     regions=0;
     eDep=0;
     dL=0;
     entryXYZT=0;
     exitXYZT=0;
     entrySlope=0;
     exitSlope=0;
     entryPE=0;
     exitPE=0;
     distToOpDet=nullptr;//Enter=0;
     opDetXYZT=nullptr;//Enter=0;
     opDetID=nullptr;//Enter=0;
     //distToOpDetExit=0;
     //opDetXYZTExit=0;
     //opDetIDExit=0;

     m_treeIn->SetBranchAddress("event",       &event,       &b_event);
     m_treeIn->SetBranchAddress("run",         &run,         &b_run);
     m_treeIn->SetBranchAddress("subRun",      &subRun,      &b_subrun);
     m_treeIn->SetBranchAddress("nReg",        &nReg,        &b_nReg);
     m_treeIn->SetBranchAddress("fiducial",    &fiducial,    &b_fiducial);
     m_treeIn->SetBranchAddress("active",      &active,      &b_active);
     m_treeIn->SetBranchAddress("inactive",    &inactive,    &b_inactive);
     m_treeIn->SetBranchAddress("crts",        &crts,        &b_crts);
     m_treeIn->SetBranchAddress("entryXYZT",   &entryXYZT,   &b_entryXYZT);
     m_treeIn->SetBranchAddress("exitXYZT",    &exitXYZT,    &b_exitXYZT);
     m_treeIn->SetBranchAddress("entryPE",     &entryPE,     &b_entryPE);
     m_treeIn->SetBranchAddress("exitPE",      &exitPE,      &b_exitPE);
     m_treeIn->SetBranchAddress("entrySlope",  &entrySlope,  &b_entrySlope);
     m_treeIn->SetBranchAddress("exitSlope",   &exitSlope,   &b_exitSlope);
     m_treeIn->SetBranchAddress("regions",     &regions,     &b_regions);
     m_treeIn->SetBranchAddress("pdg",         &pdg,         &b_pdg);
     m_treeIn->SetBranchAddress("trackID",     &trackID,     &b_trackID);
     m_treeIn->SetBranchAddress("eDep",        &eDep,        &b_eDep);
     m_treeIn->SetBranchAddress("dL",          &dL,          &b_dL);
     m_treeIn->SetBranchAddress("distToOpDet", &distToOpDet, &b_distToOpDet);
     m_treeIn->SetBranchAddress("opDetXYZT",   &opDetXYZT,   &b_opDetXYZT);
     m_treeIn->SetBranchAddress("opDetID",     &opDetID,     &b_opDetID);
     //m_treeIn->SetBranchAddress("distToOpDetEnter", &distToOpDetEnter, &b_distToOpDetEnter);
     //m_treeIn->SetBranchAddress("opDetXYZTEnter",   &opDetXYZTEnter,   &b_opDetXYZTEnter);
     //m_treeIn->SetBranchAddress("opDetIDEnter",     &opDetIDEnter,     &b_opDetIDEnter);
     //m_treeIn->SetBranchAddress("distToOpDetExit",  &distToOpDetExit, &b_distToOpDetExit);
     //m_treeIn->SetBranchAddress("opDetXYZTExit",    &opDetXYZTExit,   &b_opDetXYZTExit);
     //m_treeIn->SetBranchAddress("opDetIDExit",      &opDetIDExit,     &b_opDetIDExit);

     if(DEBUG_region) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     return true;
}

// accessors for the input TTree
TChain* Regions::GetInputTree() {return m_treeIn;};

// get all branch contents of input TTree for entry i
void Regions::GetRegion(Long64_t entry)
{
     if (!m_treeIn) {
        cout << "Regions::getEntry error" << endl;
        exit(1);
     }
     m_treeIn->GetEntry(entry);
     if(DEBUG_region) m_treeIn->Show(entry);

     m_eve       = event;
     m_run       = run;
     m_subrun    = subRun;
     m_nreg      = nReg;
     m_fiducial  = fiducial;
     m_active    = active;
     m_inactive  = inactive;
     m_crts      = crts;
     m_pdg       = pdg;
     m_trackid   = trackID;

     m_entryx.clear();
     m_entryy.clear();
     m_entryz.clear();
     m_entryt.clear();
     m_exitx.clear();
     m_exity.clear();
     m_exitz.clear();
     m_exitt.clear();
     m_entrypx.clear();
     m_entrypy.clear();
     m_entrypz.clear();
     m_entrye.clear();
     m_exitpx.clear();
     m_exitpy.clear();
     m_exitpz.clear();
     m_exite.clear();
     m_entryslopex.clear();
     m_entryslopey.clear();
     m_entryslopez.clear();
     m_exitslopex.clear();
     m_exitslopey.clear();
     m_exitslopez.clear();
     m_regions.clear();
     m_edep.clear();
     m_dl.clear();
     m_disttoopdet.clear();
     m_opdetx.clear();
     m_opdety.clear();
     m_opdetz.clear();
     m_opdett.clear();
     m_opdetid.clear();
     //m_disttoopdetenter.clear();
     //m_opdetenterx.clear();
     //m_opdetentery.clear();
     //m_opdetenterz.clear();
     //m_opdetentert.clear();
     //m_opdetidenter.clear();
     //m_disttoopdetexit.clear();
     //m_opdetexitx.clear();
     //m_opdetexity.clear();
     //m_opdetexitz.clear();
     //m_opdetexitt.clear();
     //m_opdetidexit.clear();

     if(DEBUG_region) {
         if((int)entryXYZT->size()!=nReg) cout << "entryXYZT size mismatch!" << endl;
         if((int)exitXYZT->size()!=nReg) cout << "exitXYZT size mismatch!" << endl;
         if((int)entryPE->size()!=nReg) cout << "entryPE size mismatch!" << endl;
         if((int)exitPE->size()!=nReg) cout << "exitPE size mismatch!" << endl;
         if((int)entrySlope->size()!=nReg) cout << "entrySlope size mismatch!" << endl;
         if((int)exitSlope->size()!=nReg) cout << "exitSlope size mismatch!" << endl;
         if((int)regions->size()!=nReg) cout << "regions size mismatch!" << endl;
         if((int)eDep->size()!=nReg) cout << "eDep size mismatch!" << endl;
         if((int)dL->size()!=nReg) cout << "dL size mismatch! #present - #expected = " << endl;
         //if((int)distToOpDetEnter->size()!=nReg) cout << "distToOpDetEnter size mismatch!" << endl;
         //if((int)opDetXYZTEnter->size()!=nReg) cout << "opDetXYZTEnter size mismatch!" << endl;
         //if((int)opDetIDEnter->size()!=nReg) cout << "opDetIDEnter size mismatch!" << endl;
         //if((int)distToOpDetExit->size()!=nReg) cout << "distToOpDetExit size mismatch!" << endl;
         //if((int)opDetXYZTExit->size()!=nReg) cout << "opDetXYZTExit size mismatch!" << endl;
         //if((int)opDetIDExit->size()!=nReg) cout << "opDetIDExit size mismatch!" << endl;
     }
   
     for (Int_t i=0; i<m_nreg; i++) {
         m_entryx.push_back((*entryXYZT)[i][0]); 
         m_entryy.push_back((*entryXYZT)[i][1]);
         m_entryz.push_back((*entryXYZT)[i][2]);
         m_entryt.push_back((*entryXYZT)[i][3]);
         m_exitx.push_back((*exitXYZT)[i][0]);
         m_exity.push_back((*exitXYZT)[i][1]);
         m_exitz.push_back((*exitXYZT)[i][2]);
         m_exitt.push_back((*exitXYZT)[i][3]);
         m_entrypx.push_back((*entryPE)[i][0]);
         m_entrypy.push_back((*entryPE)[i][1]);
         m_entrypz.push_back((*entryPE)[i][2]);
         m_entrye.push_back((*entryPE)[i][3]);
         m_exitpx.push_back((*exitPE)[i][0]);
         m_exitpy.push_back((*exitPE)[i][1]);
         m_exitpz.push_back((*exitPE)[i][2]);
         m_exite.push_back((*exitPE)[i][3]);
         m_entryslopex.push_back((*entrySlope)[i][0]);
         m_entryslopey.push_back((*entrySlope)[i][1]);
         m_entryslopez.push_back((*entrySlope)[i][2]);
         m_exitslopex.push_back((*exitSlope)[i][0]);
         m_exitslopey.push_back((*exitSlope)[i][1]);
         m_exitslopez.push_back((*exitSlope)[i][2]);
         m_regions.push_back((*regions)[i]);
         m_edep.push_back((*eDep)[i]);
         m_dl.push_back((*dL)[i]);
         m_disttoopdet.push_back((*distToOpDet)[i]);
         m_opdetx.push_back((*opDetXYZT)[i][0]);
         m_opdety.push_back((*opDetXYZT)[i][1]);
         m_opdetz.push_back((*opDetXYZT)[i][2]);
         m_opdett.push_back((*opDetXYZT)[i][3]);
         m_opdetid.push_back((*opDetID)[i]);
         //m_disttoopdetenter.push_back((*distToOpDetEnter)[i]);
         //m_opdetenterx.push_back((*opDetXYZTEnter)[i][0]);
         //m_opdetentery.push_back((*opDetXYZTEnter)[i][1]);
         //m_opdetenterz.push_back((*opDetXYZTEnter)[i][2]);
         //m_opdetentert.push_back((*opDetXYZTEnter)[i][3]);
         //m_opdetidenter.push_back((*opDetIDEnter)[i]);
         //m_disttoopdetexit.push_back((*distToOpDetExit)[i]);
         //m_opdetexitx.push_back((*opDetXYZTExit)[i][0]);
         //m_opdetexity.push_back((*opDetXYZTExit)[i][1]);
         //m_opdetexitz.push_back((*opDetXYZTExit)[i][2]);
         //m_opdetexitt.push_back((*opDetXYZTExit)[i][3]);
         //m_opdetidexit.push_back((*opDetIDExit)[i]);
     }

     if(DEBUG_region) m_treeIn->Print();
}
#endif
