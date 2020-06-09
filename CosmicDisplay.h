#ifndef COSMICDISPLAY_H
#define COSMICDISPLAY_H

#define DEBUG_cd 0

#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>
#include <TNamed.h>
#include <TObject.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <vector>


	class CosmicDisplay : public TObject {

	  public :

		// Declaration of leaf types
		Int_t        event;
                Int_t        run;
                Int_t        subRun;
                Int_t        trackID;
		Int_t        nSeg;
		//vector<int>        regions;
		Int_t        pdg;
                vector<vector<double>>*    slopes; //direction cosines
		vector<vector<double>>*    pe;   //4-momentum
		vector<vector<double>>*    xyzt;   //4-position

		// List of branches

		TBranch      *b_event;   //!
                TBranch      *b_run;
                TBranch      *b_subrun;
                TBranch      *b_trackID;
		TBranch      *b_nSeg;   //!
		TBranch      *b_regions;   //!
		TBranch      *b_pdg;   //!
                TBranch      *b_slopes;
		TBranch      *b_pe;   //!
		TBranch      *b_xyzt;   //!

		CosmicDisplay();
		~CosmicDisplay(){};

		static CosmicDisplay* giveThis();
		static CosmicDisplay* giveThis(TTree* aTree, const std::string& option);
		static void releaseThis();
		void setupRead(const std::string& option = "");
		bool SetBranchAddresses();
		CosmicDisplay(TTree* aTree, const std::string& option);
		TTree* GetInputTree();// {return m_treeIn;};

		void GetCosmic(Long64_t entry);

		Int_t       Event()                   { return m_event;         };
                Int_t       Run()                     { return m_run;           };
                Int_t       SubRun()                  { return m_subrun;        };
                Int_t       TrackID()                 { return m_trackid;       };
		Int_t       NSeg()                    { return m_nseg;          };
		//Int_t       Region(Int_t index)       { return m_regions[index];};
		Int_t       PDG()                     { return m_pdg;           };
		Double_t    SlopeX(Int_t index)       { return m_slopex[index]; };
		Double_t    SlopeY(Int_t index)       { return m_slopey[index]; };
		Double_t    SlopeZ(Int_t index)       { return m_slopez[index]; };
		Double_t    E(Int_t index)            { return m_e[index];      };
		Double_t    Px(Int_t index)           { return m_px[index];     };
                Double_t    Py(Int_t index)           { return m_py[index];     };
                Double_t    Pz(Int_t index)           { return m_pz[index];     };
                Double_t    P(Int_t index)            { return m_p[index];      }
		Double_t    X(Int_t index)            { return m_x[index];      };
		Double_t    Y(Int_t index)            { return m_y[index];      };
		Double_t    Z(Int_t index)            { return m_z[index];      };
		Double_t    T(Int_t index)            { return m_t[index];      };

	private:

		TTree* m_treeIn;

		Int_t    m_event;
                Int_t    m_run;
                Int_t    m_subrun;
                Int_t    m_trackid;
                Int_t    m_nseg; 
		//vector<int>    m_regions;
                Int_t    m_pdg;
		vector<double> m_slopex;
                vector<double> m_slopey;
                vector<double> m_slopez;
		vector<double> m_x; 
                vector<double> m_y;
                vector<double> m_z;
                vector<double> m_t;
		vector<double> m_e;
                vector<double> m_px;
                vector<double> m_py;
                vector<double> m_pz;
                vector<double> m_p;

		static CosmicDisplay* m_instance;
		ClassDef(CosmicDisplay,0)

	};

	// initialisation of the CosmicDisplay pointer
	CosmicDisplay* CosmicDisplay::m_instance = 0;
	CosmicDisplay::CosmicDisplay(){	}

	// to get a unique instance 
	inline CosmicDisplay* CosmicDisplay::giveThis()
	{
		if(DEBUG_cd) cout << "debug::m_instance cosmicDisplay " << m_instance << endl;
		if (0 == m_instance){
			cout << "CosmicDisplay::giveThis error not constructed properly " << endl;
		}
		return m_instance;
	}

	// Open a TTree for reading (writing option not now)
	inline CosmicDisplay* CosmicDisplay::giveThis(TTree* aTree, const std::string& option)
	{
		if(DEBUG_cd) cout << "debug::m_instance cosmicDisplay " << m_instance << endl;
		if (0 == m_instance){
			m_instance = new CosmicDisplay(aTree, option);
		} 
                else{
			cout << "CosmicDisplay::giveThis Warning " << aTree->GetTitle() << endl;
		}
		return m_instance;
	}

	// Delete unique instance
	inline void CosmicDisplay::releaseThis() {
		if ( m_instance != 0 ) {
			delete m_instance;
			m_instance = 0;
		}
	}

	// constructor for one TTree with read/write option
	CosmicDisplay::CosmicDisplay(TTree* aTree, const std::string& option){
		if(option=="read"){
			m_treeIn = aTree;
			if(DEBUG_cd) cout << " ....before setupread.... " << endl;
			setupRead();
		}
	}

	// setup the input tree for reading
	void CosmicDisplay::setupRead(const std::string& option){
		if(DEBUG_cd) cout << " ....inside setupread.... " << endl;
		if(DEBUG_cd) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
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
	bool CosmicDisplay::SetBranchAddresses(){
		if(DEBUG_cd) cout << " ....inside setbranchaddresses.... " << endl;

                slopes = 0;
                pe = 0;
                xyzt = 0;


		m_treeIn->SetMakeClass(1);
		m_treeIn->SetBranchAddress("event",   &event,  &b_event);
                m_treeIn->SetBranchAddress("run",     &run,    &b_run);
                m_treeIn->SetBranchAddress("subRun",  &subRun, &b_subrun);
                m_treeIn->SetBranchAddress("trackID", &trackID,&b_trackID);
		m_treeIn->SetBranchAddress("nSeg",    &nSeg,   &b_nSeg);
		//m_treeIn->SetBranchAddress("regions", &regions, &b_regions);
		m_treeIn->SetBranchAddress("pdg",     &pdg,    &b_pdg);
                m_treeIn->SetBranchAddress("slopes",  &slopes,  &b_slopes);
		m_treeIn->SetBranchAddress("pe",      &pe,      &b_pe);
		m_treeIn->SetBranchAddress("xyzt",    &xyzt,    &b_xyzt);

		if(DEBUG_cd) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
		return true;
	}

	// accessors for the input TTree
	TTree* CosmicDisplay::GetInputTree() {return m_treeIn;};

	// get all branch cpntents of input TTree for entry i
	void CosmicDisplay::GetCosmic(Long64_t entry)
	{
		if (!m_treeIn) {
			cout << "CosmicDisplay::getEntry error" << endl;
			exit(1);
		}
		m_treeIn->GetEntry(entry);
		if(DEBUG_cd) m_treeIn->Show(entry);
		m_event=event;
                m_run = run;
                m_subrun = subRun;
                m_trackid=trackID;
                m_nseg=nSeg; 
                m_pdg        = pdg;

                m_slopex.clear();
                m_slopey.clear();
                m_slopez.clear();
                m_e.clear();
                m_px.clear();
                m_py.clear();
                m_pz.clear();
                m_p.clear();
                m_x.clear();
                m_y.clear();
                m_z.clear();
                m_t.clear();

		for(int i=0; i<nSeg;i++){
			//m_regions.push_back(regions[i]);
			m_slopex.push_back((*slopes)[i][0]);
			m_slopey.push_back((*slopes)[i][1]);
			m_slopez.push_back((*slopes)[i][2]);
			m_e.push_back((*pe)[i][3]);
			m_px.push_back((*pe)[i][0]);
			m_py.push_back((*pe)[i][1]);
			m_pz.push_back((*pe)[i][2]);
                        m_p.push_back(TMath::Sqrt(m_px[i]*m_px[i]
                                         +m_py[i]*m_py[i] + m_pz[i]*m_pz[i] ));
                        m_x.push_back((*xyzt)[i][0]);
                        m_y.push_back((*xyzt)[i][1]);
                        m_z.push_back((*xyzt)[i][2]);
			m_t.push_back((*xyzt)[i][3]);
		}
		if(DEBUG_cd) m_treeIn->Print();
	}

	
#endif


