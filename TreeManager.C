#include "TreeManager.h"

#include <iostream>
#include <stdlib.h>
#include <list>

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif
#ifdef __MAKECINT__
#pragma link off all class;
#pragma link C++ class CosmicDisplay;
#pragma link C++ class Simulation;
#pragma link C++ class Regions;
#pragma link C++ class DetSim;
#pragma link C++ class Hit;
#pragma link C++ class TrueHit;
#pragma link C++ class Gen;

#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;
#endif

class CosmicDisplay;
class Simulation;
class Regions;
class DetSim;
class Hit;
class TrueHit;
class Gen;

TreeManager::TreeManager(std::string infilename){
        m_inFilename = infilename;
        cout << "Processing file: " << m_inFilename.c_str() << endl;

        //m_inFile = new TFile(m_inFilename.c_str(), "READ");
        //if( m_inFile->IsOpen() == kFALSE ) return;
        //cout << "File opened!" << endl;

        m_nameToPtr["SimTree"] = m_treeSim;
        m_nameToPtr["RegTree"] = m_treeReg;
        m_nameToPtr["DetTree"] = m_treeDet;
        m_nameToPtr["HitTree"] = m_treeHit;
        m_nameToPtr["TrueCRTHitTree"] = m_treeTrueHit;
        m_nameToPtr["DisplayTree"] = m_treeDis;
        m_nameToPtr["GenTree"] = m_treeGen;

        std::string suff = infilename.substr(infilename.length()-4,4);
        char type = 'e';
        if(suff=="root") type = 'f';
        if(suff=="list") type = 'l';

        init(type);
}

TreeManager::~TreeManager(){
        delete  m_treeDis;
        delete  m_treeSim;
        delete  m_treeReg;
        delete  m_treeDet;
        delete  m_treeHit;
        delete  m_treeTrueHit;
        delete  m_treeGen;
        delete  m_inFile;
}

void TreeManager::nullify(){
        m_treeDis = NULL;
        m_treeSim = NULL;
        m_treeReg = NULL;
        m_treeDet = NULL;
        m_treeHit = NULL;
        m_treeTrueHit = NULL;
        m_treeGen = NULL;
        CosmicDisplay::releaseThis();
        Simulation::releaseThis();
        Regions::releaseThis();
        DetSim::releaseThis();
        Hit::releaseThis();
        TrueHit::releaseThis();
        Gen::releaseThis();
}

void TreeManager::fillchain(const char type, const std::string tname) {

    if(tname=="SimTree") m_treeSim = new TChain(tname.c_str());
    if(tname=="RegTree") m_treeReg = new TChain(tname.c_str());
    if(tname=="DetTree") m_treeDet = new TChain(tname.c_str());
    if(tname=="HitTree") m_treeHit = new TChain(tname.c_str());
    if(tname=="TrueCRTHitTree") m_treeTrueHit = new TChain(tname.c_str());
    if(tname=="DisplayTree") m_treeDis = new TChain(tname.c_str());
    if(tname=="GenTree") m_treeGen = new TChain(tname.c_str());
    //m_nameToPtr[tname] = new TChain(tname.c_str());
    std::string tpath = "CRTSimAnalysis/"+tname;

    if(type=='l') {

        ifstream fin;
        fin.open(m_inFilename.c_str());
        std::string line;
        //m_nameToPtr[tname] = new TChain(tname.c_str());
        //std::string tpath = "CRTSimAnalysis/"+tname;

        while(getline(fin,line)){
           // m_nameToPtr[tname]->AddFile(line.c_str(),TTree::kMaxEntries,tpath.c_str());
            if(tname=="SimTree") m_treeSim->AddFile(line.c_str(),TTree::kMaxEntries,tpath.c_str());
            if(tname=="RegTree") m_treeReg->AddFile(line.c_str(),TTree::kMaxEntries,tpath.c_str());
            if(tname=="DetTree") m_treeDet->AddFile(line.c_str(),TTree::kMaxEntries,tpath.c_str());
            if(tname=="HitTree") m_treeHit->AddFile(line.c_str(),TTree::kMaxEntries,tpath.c_str());
            if(tname=="TrueCRTHitTree") m_treeTrueHit->AddFile(line.c_str(),TTree::kMaxEntries,tpath.c_str());
            if(tname=="DisplayTree") m_treeDis->AddFile(line.c_str(),TTree::kMaxEntries,tpath.c_str());
            if(tname=="GenTree") m_treeGen->AddFile(line.c_str(),TTree::kMaxEntries,tpath.c_str());
        }

    }

    else if(type=='f'){
        //m_nameToPtr[tname]->AddFile(m_inFilename.c_str(),TTree::kMaxEntries,tpath.c_str()); 
        if(tname=="SimTree") m_treeSim->AddFile(m_inFilename.c_str(),TTree::kMaxEntries,tpath.c_str());
        if(tname=="RegTree") m_treeReg->AddFile(m_inFilename.c_str(),TTree::kMaxEntries,tpath.c_str());
        if(tname=="DetTree") m_treeDet->AddFile(m_inFilename.c_str(),TTree::kMaxEntries,tpath.c_str());
        if(tname=="HitTree") m_treeHit->AddFile(m_inFilename.c_str(),TTree::kMaxEntries,tpath.c_str());
        if(tname=="TrueCRTHitTree") m_treeTrueHit->AddFile(m_inFilename.c_str(),TTree::kMaxEntries,tpath.c_str());
        if(tname=="DisplayTree") m_treeDis->AddFile(m_inFilename.c_str(),TTree::kMaxEntries,tpath.c_str());
        if(tname=="GenTree") m_treeGen->AddFile(m_inFilename.c_str(),TTree::kMaxEntries,tpath.c_str());
        //cout << "new chain with " << m_nameToPtr[tname]->GetEntries() << endl;
    }

    else {
        //throw cet::exception("TreeManger") << "uknown input file type" << std::endl;
        std::cout << "CRITICAL ERROR: UNKNOWN INPUT FILE TYPE! SUFFIX MUST BE '.list' OR '.root'" << std::endl;
    }

}

void TreeManager::init(const char type){
        nullify();
        int ntree=0;
        for(auto const& name : m_nameToPtr){
            fillchain(type,name.first);
            ntree++;
        }
        cout << "filled " << ntree << " TChains"  << endl;

        //m_treeSim = (TChain*) m_inFile->FindObjectAny("SimTree");
        if( m_treeSim != NULL ) {
                cout << "m_treeSim found!" << endl;
                m_tsim = Simulation::giveThis(m_treeSim,"read");
                if( m_tsim != NULL ) {
                        cerr << " [1] it is ok!!" << endl;
                }
        }
        else {std::cout << "null SimTree!" << std::endl;}

        //m_treeReg = (TChain*) m_inFile->FindObjectAny("RegTree");
        if( m_treeReg != NULL ) {
                cout << "m_treeReg found!" << endl;
                m_treg = Regions::giveThis(m_treeReg,"read");
                if( m_treg != NULL ) {
                        cerr << " [2] it is ok!!" << std::endl;
                }
        }
        else{ std::cout << "null RegTree!" << std::endl;}

        //m_treeDet = (TChain*) m_inFile->FindObjectAny("DetTree");
        if( m_treeDet != NULL ) {
                cout << "m_treeDet found!" << endl;
                m_tdet = DetSim::giveThis(m_treeDet,"read");
                if( m_tdet != NULL ) {
                        cerr << " [3] it is ok!!" << endl;
                }
        }
        else{ std::cout << "null DetTree!" << std::endl;}

        //m_treeHit = (TChain*) m_inFile->FindObjectAny("HitTree");
        if( m_treeHit != NULL ) {
                cout << "m_treeHit found!" << endl;
                m_thit = Hit::giveThis(m_treeHit,"read");
                if( m_thit != NULL ) {
                        cerr << " [4] it is ok!!" << endl;
                }
        }
        else {std::cout << "null HitTree!" << std::endl;}

        //m_treeTrueHit = (TChain*) m_inFile->FindObjectAny("TrueCRTHitTree");
        if( m_treeTrueHit != NULL ) {
                cout << "m_treeTrueHit found!" << endl;
                m_ttruehit = TrueHit::giveThis(m_treeTrueHit,"read");
                if( m_ttruehit != NULL ) {
                        cerr << " [5] it is ok!!" << endl;
                }
        }
        else {std::cout << "null TrueHitTree!" << std::endl;}

        //m_treeDis = (TChain*) m_inFile->FindObjectAny("DisplayTree");
        if( m_treeHit != NULL ) {
                cout << "m_treeDis found!" << endl;
                m_tdis = CosmicDisplay::giveThis(m_treeDis,"read");
                if( m_tdis != NULL ) {
                        cerr << " [6] it is ok!!" << endl;
                }
        }
        else {std::cout << "null DisplayTree!" << std::endl;}

        //m_treeGen = (TChain*) m_inFile->FindObjectAny("GenTree");
        if( m_treeGen != NULL ) {
                cout << "m_treeGen found!" << endl;
                m_tgen = Gen::giveThis(m_treeGen,"read");
                if( m_tgen != NULL ) {
                        cerr << " [7] it is ok!!" << endl;
                }
        }
        else {std::cout << "null GenTree!" << std::endl;}
}
