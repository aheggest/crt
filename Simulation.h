//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 17 12:41:28 2018 by ROOT version 6.12/06
// from TTree SimTree/MyCRTSimulation
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef Simulation_h
#define Simulation_h

#define DEBUG_sim 0

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TNamed.h>
#include <TObject.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <utility>

class Simulation : public TObject {

public :
   // Declaration of leaf types
   Int_t                   event;
   Int_t                   subRun;
   Int_t                   run;
   Int_t                   trackID;
   Int_t                   pdg;
   Int_t                   nPoints;
   Float_t                 trackLength;
   Int_t                   process;
   Int_t                   endProcess;
   Int_t                   parentPDG;
   Float_t                 parentE;
   Int_t                   progenitor;
   vector<int>*            auxDetSensitiveID;
   vector<int>*            auxDetID;
   vector<double>*         auxDetEDep;
   vector<double>*         auxDetdEdx;
   vector<double>*         auxDetTrackLength;
   vector<vector<double>>* auxDetEnterXYZT;
   vector<vector<double>>* auxDetExitXYZT;
   vector<int>*            auxDetRegion;
   vector<int>*            mac5;
   vector<int>*            adType;
   double                  startXYZT[4];
   double                  endXYZT[4];
   double                  startPE[4];
   double                  endPE[4];
   int                     nChan;
   Int_t                   mother;
   Int_t                   nDaught;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_subRun;   //!
   TBranch        *b_run;   //!
   TBranch        *b_trackID;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_nPoints;   //!
   TBranch        *b_trackLenth;   //!
   TBranch        *b_process;   //!
   TBranch        *b_endProcess;   //!
   TBranch        *b_parentPDG;   //!
   TBranch        *b_parentE;   //!
   TBranch        *b_progenitor;   //!
   TBranch        *b_auxDetSensitiveID;   //!
   TBranch        *b_auxDetID;   //!
   TBranch        *b_auxDetEDep;   //!
   TBranch        *b_auxDetdEdx;   //!
   TBranch        *b_auxDetTrackLength;   //!
   TBranch        *b_auxDetEnterXYZT;   //!
   TBranch        *b_auxDetExitXYZT;   //!
   TBranch        *b_auxDetRegion;   //!
   TBranch        *b_mac5;   //!
   TBranch        *b_adType;
   TBranch        *b_startXYZT;   //!
   TBranch        *b_endXYZT;   //!
   TBranch        *b_startPE;   //!
   TBranch        *b_endPE;   //!
   TBranch        *b_nChan;   //!
   TBranch        *b_mother;   //!
   TBranch        *b_nDaught;   //!

   Simulation();
   ~Simulation() {};
   static Simulation* giveThis();
   static Simulation* giveThis(TChain* aTree, const std::string& option);
   static void releaseThis();
   void setupRead(const std::string& option = "");
   bool SetBranchAddresses();
   Simulation(TChain* aTree, const std::string& option);
   TChain* GetInputTree();
   void GetSim(Long64_t entry);

   Int_t Event()                  { return m_eve;     };
   Int_t Run()                    { return m_run;     };
   Int_t Subrun()                 { return m_subrun;  };
   Int_t NPoints()                { return m_npoints; };
   Int_t TrackID()                { return m_trkid;   };
   Int_t PDG()                    { return m_pdg;     };
   double TrackLength()           { return m_trklen;  };
   Int_t StartProcess()           { return m_process; };
   Int_t EndProcess()             { return m_endproc; };
   Int_t ParentPDG()              { return m_parpdg;  };
   Int_t Progenitor()             { return m_progen;  };
   Int_t Mother()                 { return m_mother;  };
   Int_t NDaughters()             { return m_ndaught; };
   double StartX()                { return m_startxt[0]; };
   double StartY()                { return m_startxt[1]; };
   double StartZ()                { return m_startxt[2]; };
   double StartT()                { return m_startxt[3]; };
   double EndX()                  { return m_endxt[0]; };
   double EndY()                  { return m_endxt[1]; };
   double EndZ()                  { return m_endxt[2]; };
   double EndT()                  { return m_endxt[3]; };
   double StartPx()               { return m_startpe[0]; };
   double StartPy()               { return m_startpe[1]; };
   double StartPz()               { return m_startpe[2]; };
   double StartE()                { return m_startpe[3]; };
   double EndPx()                 { return m_endpe[0]; };
   double EndPy()                 { return m_endpe[1]; };
   double EndPz()                 { return m_endpe[2]; };
   double EndE()                  { return m_endpe[3]; };
   int     NCRTHit()              { return m_nchan;    };
   Int_t   CRTModule(Int_t index) { return m_adid[index]; };
   Int_t   CRTStrip(Int_t index)  { return m_adsid[index]; };
   double CRTEnterX(Int_t index)  { return m_adentxt[index][0]; };
   double CRTEnterY(Int_t index)  { return m_adentxt[index][1]; };
   double CRTEnterZ(Int_t index)  { return m_adentxt[index][2]; };
   double CRTEnterT(Int_t index)  { return m_adentxt[index][3]; };
   double CRTExitX(Int_t index)   { return m_adextxt[index][0]; };
   double CRTExitY(Int_t index)   { return m_adextxt[index][1]; };
   double CRTExitZ(Int_t index)   { return m_adextxt[index][2]; };
   double CRTExitT(Int_t index)   { return m_adextxt[index][3]; };
   double CRTEDep(Int_t index)    { return m_adedep[index]; };
   double CRTdEdx(Int_t index)    { return m_addedx[index]; };
   double CRTTrackLength(Int_t index) { return m_adtrkl[index]; };
   Int_t   CRTRegion(Int_t index) { return m_adreg[index]; };
   Int_t   Mac5(Int_t index)      { return m_mac5[index]; };
   Int_t   CRTType(Int_t index)   { return m_adtype[index]; };

private:
   TChain*                 m_treeIn;     //pointer to input tree
   Int_t                  m_eve;        //event number
   Int_t                  m_subrun;     //subrun number
   Int_t                  m_run;        //run number
   Int_t                  m_trkid;      //track ID from G4
   Int_t                  m_pdg;        //PDG code
   Int_t                  m_npoints;    //no. of trajectory points
   double                 m_trklen;     //track length (|r_final - r_init|)
   Int_t                  m_process;    //process (code) that produced the particle
   Int_t                  m_endproc;    //process (code) that killed the particle
   Int_t                  m_parpdg;     //PDG code of parent of this particle
   Int_t                  m_pareng;     //energy of parent particle
   Int_t                  m_progen;     //progenitor of particle (oldest traceable ancestor)
   vector<int>            m_adsid;      //CRT strip ID no.
   vector<int>            m_adid;       //CRT module ID no.
   vector<double>         m_adedep;     //energy deposited in CRT strip
   vector<double>         m_addedx;     // <dE/dx> for this deposition
   vector<double>         m_adtrkl;     //path length of track in strip
   vector<vector<double>> m_adentxt;    //point and time of entry into CRT strip 
   vector<vector<double>> m_adextxt;    //point and time of entry from CRT strip
   vector<int>            m_adreg;      //region (code) of CRT hit
   vector<int>            m_mac5;       //mac5 ID of front-end board assigned to CRT module
   vector<int>            m_adtype;     //CRT module type code (0=CERN, 1=MINOS, 2=DC)
   double                 m_startxt[4]; //start position and time of particle
   double                 m_endxt[4];   //end position and time of particle
   double                 m_startpe[4]; //start momentum and energy of particle
   double                 m_endpe[4];   //end momentum and energy of particle
   int                    m_nchan;      //no. CRT strips hit channels by particle
   Int_t                  m_mother;     //track ID of mother particle
   Int_t                  m_ndaught;    //number of daughters produced by this particle
   static Simulation* m_instance;
   ClassDef(Simulation,0);

};//end class

// initialisation of the Simulation pointer
Simulation* Simulation::m_instance = 0;
Simulation::Simulation() {}
//to get a unique instance
inline Simulation* Simulation::giveThis()
{
     if(DEBUG_sim) cout << "debug::m_instance Simulation " << m_instance << endl;
     if (0 == m_instance){
        cout << "Simulation::giveThis error not constructed properly " << endl;
     }
     return m_instance;

}

// Open a TTree for reading
inline Simulation* Simulation::giveThis(TChain* aTree, const std::string& option)
{
     if(DEBUG_sim) cout << "debug::m_instance Simulation " << m_instance << endl;
     if (0 == m_instance){
        m_instance = new Simulation(aTree, option);
     } 
     else{
        cout << "Simulation::giveThis Warning " << aTree->GetTitle() << endl;
     }
     return m_instance;
}

// Delete unique instance
inline void Simulation::releaseThis() {
     if ( m_instance != 0 ) {
        delete m_instance;
        m_instance = 0;
     }
}

// constructor for one TTree with read/write option
Simulation::Simulation(TChain* aTree, const std::string& option){
     if(option=="read"){
        m_treeIn = aTree;
        if(DEBUG_sim) cout << " ....before setupread.... " << endl;
        setupRead();
     }
}

// setup the input tree for reading
void Simulation::setupRead(const std::string& option){
     if(DEBUG_sim) cout << " ....inside setupread.... " << endl;
     if(DEBUG_sim) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     if(!m_treeIn){
        cout << "setupRead error: m_treeIn undefined " << endl;
        exit(1);
     }
     if(SetBranchAddresses()){}
     else{
        cerr << "TBranch error.."<< endl;
     }
}

// Setting up Branch content of input TTree
bool Simulation::SetBranchAddresses() {
     if(DEBUG_sim) cout << " ....inside setbranchaddresses.... " << endl;

     auxDetSensitiveID=0;
     auxDetID=nullptr;
     auxDetEDep=nullptr;
     auxDetdEdx=nullptr;
     auxDetTrackLength=nullptr;
     auxDetEnterXYZT=nullptr;
     auxDetExitXYZT=nullptr;
     auxDetRegion=nullptr;
     mac5=nullptr;
     adType=nullptr;

     m_treeIn->SetMakeClass(1);
     m_treeIn->SetBranchAddress("event",             &event,             &b_event);
     m_treeIn->SetBranchAddress("subRun",            &subRun,            &b_subRun);
     m_treeIn->SetBranchAddress("run",               &run,               &b_run);
     m_treeIn->SetBranchAddress("trackID",           &trackID,           &b_trackID);
     m_treeIn->SetBranchAddress("pdg",               &pdg,               &b_pdg);
     m_treeIn->SetBranchAddress("nPoints",           &nPoints,           &b_nPoints);
     m_treeIn->SetBranchAddress("trackLength",       &trackLength,       &b_trackLenth);
     m_treeIn->SetBranchAddress("process",           &process,           &b_process);
     m_treeIn->SetBranchAddress("endProcess",        &endProcess,        &b_endProcess);
     m_treeIn->SetBranchAddress("parentPDG",         &parentPDG,         &b_parentPDG);
     m_treeIn->SetBranchAddress("parentE",           &parentE,           &b_parentE);
     m_treeIn->SetBranchAddress("progenitor",        &progenitor,        &b_progenitor);
     m_treeIn->SetBranchAddress("auxDetSensitiveID", &auxDetSensitiveID, &b_auxDetSensitiveID);
     m_treeIn->SetBranchAddress("auxDetID",          &auxDetID,          &b_auxDetID);
     m_treeIn->SetBranchAddress("auxDetEDep",        &auxDetEDep,        &b_auxDetEDep);
     m_treeIn->SetBranchAddress("auxDetdEdx",        &auxDetdEdx,        &b_auxDetdEdx);
     m_treeIn->SetBranchAddress("auxDetTrackLength", &auxDetTrackLength, &b_auxDetTrackLength);
     m_treeIn->SetBranchAddress("auxDetEnterXYZT",   &auxDetEnterXYZT,    &b_auxDetEnterXYZT);
     m_treeIn->SetBranchAddress("auxDetExitXYZT",    &auxDetExitXYZT,     &b_auxDetExitXYZT);
     m_treeIn->SetBranchAddress("auxDetRegion",      &auxDetRegion,      &b_auxDetRegion);
     m_treeIn->SetBranchAddress("mac5",              &mac5,              &b_mac5);
     m_treeIn->SetBranchAddress("adType",            &adType,            &b_adType);
     m_treeIn->SetBranchAddress("startXYZT",         startXYZT,          &b_startXYZT);
     m_treeIn->SetBranchAddress("endXYZT",           endXYZT,            &b_endXYZT);
     m_treeIn->SetBranchAddress("startPE",           startPE,            &b_startPE);
     m_treeIn->SetBranchAddress("endPE",             endPE,              &b_endPE);
     m_treeIn->SetBranchAddress("nChan",             &nChan,             &b_nChan);
     m_treeIn->SetBranchAddress("mother",            &mother,            &b_mother);
     m_treeIn->SetBranchAddress("nDaught",           &nDaught,           &b_nDaught);
     if(DEBUG_sim) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
     return true;
}
 
// accessors for the input TTree
TChain* Simulation::GetInputTree() {return m_treeIn;};

// get all branch contents of input TTree for entry i
void Simulation::GetSim(Long64_t entry)
{ 
     if (!m_treeIn) {
        cout << "Simulation::getEntry error" << endl;
        exit(1);
     }
     m_treeIn->GetEntry(entry);
     if(DEBUG_sim) m_treeIn->Show(entry);
     m_eve     = event;
     m_subrun  = subRun;
     m_run     = run;
     m_trkid   = trackID;
     m_pdg     = pdg;
     m_npoints = nPoints;
     m_trklen  = trackLength;
     m_process = process;
     m_endproc = endProcess;
     m_parpdg  = parentPDG;
     m_pareng  = parentE;
     m_progen  = progenitor;
     m_nchan   = nChan;
     m_mother  = mother;
     m_ndaught = nDaught;

     m_adsid.clear();
     m_adid.clear();
     m_adedep.clear();
     m_addedx.clear();
     m_adtrkl.clear();
     m_adreg.clear();
     m_mac5.clear();
     m_adtype.clear();
     m_adentxt.clear();
     m_adextxt.clear();

     for (int i=0; i<m_nchan; i++) {
        m_adsid.push_back(auxDetSensitiveID->at(i));
        m_adid.push_back(auxDetID->at(i));
        m_adedep.push_back(auxDetEDep->at(i));
        m_addedx.push_back(auxDetdEdx->at(i));
        m_adtrkl.push_back(auxDetTrackLength->at(i));
        m_adreg.push_back(auxDetRegion->at(i));
        m_mac5.push_back(mac5->at(i));
        m_adtype.push_back(adType->at(i));
        vector<double> tmpentxt, tmpextxt;
        for (int j=0; j<4; j++) {
           tmpentxt.push_back((*auxDetEnterXYZT)[i][j]);
           tmpextxt.push_back((*auxDetExitXYZT)[i][j]);
        }
        m_adentxt.push_back(tmpentxt);
        m_adextxt.push_back(tmpextxt);
     }

     for (int i=0; i<4; i++) {
        m_startpe[i] = startPE[i];
        m_endpe[i]   = endPE[i];
        m_startxt[i] = startXYZT[i];
        m_endxt[i]   = endXYZT[i];
     }

     if(DEBUG_sim) m_treeIn->Print();
}

#endif
