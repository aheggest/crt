#ifndef TREEMANAGER_h
#define TREEMANAGER_h

#include "CosmicDisplay.h"
#include "Simulation.h"
#include "Regions.h"
#include "DetSim.h"
#include "Hit.h"
#include "TrueHit.h"
#include "Gen.h"

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

//#include "cetlib/exception.h"
#include <string>

#ifdef __CINT__
#pragma link C++ class vector<vector<double> >;
#pragma link C++ class vector<vector<float> >;
#pragma link C++ class vector<vector<int> >;
#else
template class std::vector<std::vector<double> >;
template class std::vector<std::vector<float> >;
template class std::vector<std::vector<int> >;
#endif


class TreeManager {
        public :
                TreeManager(std::string infilename);
                virtual ~TreeManager();

                CosmicDisplay *tmCD()   { cout << "..tmCD.."  << endl; return m_tdis; };
                Simulation *tmSim()     { cout << "..tmSim.." << endl; return m_tsim; };
                Regions    *tmReg()     { cout << "..tmReg.." << endl; return m_treg; };
                DetSim     *tmDet()     { cout << "..tmDet.." << endl; return m_tdet; };
                Hit        *tmHit()     { cout << "..tmHit.." << endl; return m_thit; };
                TrueHit    *tmTrueHit() { cout << "..tmTrueHit.." << endl; return m_ttruehit; };
                Gen        *tmGen()     { cout << "..tmGen.." << endl; return m_tgen; };

        protected:
                void          init(const char type);
                void          fillchain(const char type, const std::string tname);
                void          nullify();

        private:
                std::string   m_inFilename;

                TFile         *m_inFile;

                TChain         *m_treeDis;
                TChain         *m_treeSim;
                TChain         *m_treeReg;
                TChain         *m_treeDet;
                TChain         *m_treeHit;
                TChain         *m_treeTrueHit;
                TChain         *m_treeGen;

                CosmicDisplay *m_tdis;
                Simulation    *m_tsim;
                Regions       *m_treg;
                DetSim        *m_tdet;
                Hit           *m_thit;
                TrueHit       *m_ttruehit;
                Gen           *m_tgen;

                std::map<std::string,TChain*> m_nameToPtr;
};
#endif //TreeManager_h

