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
        init();
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

void TreeManager::init(){
        nullify();
        m_inFile = new TFile(m_inFilename.c_str(), "READ");
        if( m_inFile->IsOpen() == kFALSE ) return;
        cout << "File opened!" << endl;

        m_treeSim = (TTree*) m_inFile->FindObjectAny("SimTree");
        if( m_treeSim != NULL ) {
                cout << "m_treeSim found!" << endl;
                m_tsim = Simulation::giveThis(m_treeSim,"read");
                if( m_tsim != NULL ) {
                        cerr << " [1] it is ok!!" << endl;
                }
        }

        m_treeReg = (TTree*) m_inFile->FindObjectAny("RegTree");
        if( m_treeReg != NULL ) {
                cout << "m_treeReg found!" << endl;
                m_treg = Regions::giveThis(m_treeReg,"read");
                if( m_treg != NULL ) {
                        cerr << " [2] it is ok!!" << endl;
                }
        }

        m_treeDet = (TTree*) m_inFile->FindObjectAny("DetTree");
        if( m_treeDet != NULL ) {
                cout << "m_treeDet found!" << endl;
                m_tdet = DetSim::giveThis(m_treeDet,"read");
                if( m_tdet != NULL ) {
                        cerr << " [3] it is ok!!" << endl;
                }
        }

        m_treeHit = (TTree*) m_inFile->FindObjectAny("HitTree");
        if( m_treeHit != NULL ) {
                cout << "m_treeHit found!" << endl;
                m_thit = Hit::giveThis(m_treeHit,"read");
                if( m_thit != NULL ) {
                        cerr << " [4] it is ok!!" << endl;
                }
        }

        m_treeTrueHit = (TTree*) m_inFile->FindObjectAny("TrueCRTHitTree");
        if( m_treeTrueHit != NULL ) {
                cout << "m_treeTrueHit found!" << endl;
                m_ttruehit = TrueHit::giveThis(m_treeTrueHit,"read");
                if( m_ttruehit != NULL ) {
                        cerr << " [5] it is ok!!" << endl;
                }
        }

        m_treeDis = (TTree*) m_inFile->FindObjectAny("DisplayTree");
        if( m_treeHit != NULL ) {
                cout << "m_treeDis found!" << endl;
                m_tdis = CosmicDisplay::giveThis(m_treeDis,"read");
                if( m_tdis != NULL ) {
                        cerr << " [6] it is ok!!" << endl;
                }
        }

        m_treeGen = (TTree*) m_inFile->FindObjectAny("GenTree");
        if( m_treeGen != NULL ) {
                cout << "m_treeGen found!" << endl;
                m_tgen = Gen::giveThis(m_treeGen,"read");
                if( m_tgen != NULL ) {
                        cerr << " [7] it is ok!!" << endl;
                }
        }
}
