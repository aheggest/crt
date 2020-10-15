//local includes
#include "TreeManager.h"

//ROOT includes
#include <TGeoManager.h>
#include <TPolyMarker3D.h>
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveElement.h>
#include <TEveUtil.h>
#include <TEveGeoNode.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEveBoxSet.h>
#include <TBranch.h>
#include <TEveViewer.h>
#include <TGLViewer.h>
#include <TEveTrans.h>
#include <TEveScene.h>
#include <TEveBrowser.h>
#include <TH3F.h>
#include <TObjArray.h>
#include <THashTable.h>
#include <TPaveText.h>
#include <TLegend.h>

//C++ includes
#include <string>
#include <iostream>
#include <list>
#include <algorithm>
#include <set>
#include <cmath>
#include <climits>

using namespace std;

//file to be analyzed
TreeManager filereader("./files.list");

//////////////////////////////////////
// Simulation Tree specific functions
//////////////////////////////////////

//Plot PDG codes found in tree
void PlotPDGs() {

    Simulation *sim = filereader.tmSim(); //pointer to Simulation class
    Int_t entriesSim = sim->GetInputTree()->GetEntries(); //number of entries in tree
    TH1F* h = new TH1F("h","G4 PDG Codes",10000,-5000,5000); //histogram for PDS codes

    //loop over entries
    for ( Int_t ientry = 0; ientry < entriesSim; ientry++ )
    {
        sim->GetSim(ientry); //get info for this entry
        h->Fill(sim->PDG()); //fill histo with PDG code
    }

    h->Draw();//draw the histogram

}

//Plot the number of CRT strips hit per track
void PlotModuleEntries() {
    Simulation *sim = filereader.tmSim();
    Int_t entriesSim = sim->GetInputTree()->GetEntries();

    TH1F* h = new TH1F("h","True CRT Hit Strip Multiplicity",10,0,10);

    for ( Int_t ientry = 0; ientry < entriesSim; ientry++ ) {
        sim->GetSim(ientry);
        h->Fill(sim->NCRTHit());
    }

    h->Draw();
}

//Plot the CRT module IDs
void PlotHitModules() {
    Simulation *sim = filereader.tmSim();
    Int_t entriesSim = sim->GetInputTree()->GetEntries();

    TH1F* h = new TH1F("h","Hit Modules",300,0,300);

    for ( Int_t ientry = 0; ientry < entriesSim; ientry++ ) {
        sim->GetSim(ientry);
        for ( Int_t i=0; i< sim->NCRTHit(); i++) {
            h->Fill(sim->CRTModule(i));
        }
    }

    h->Draw();

}

//Plot the number of CRT modules hit per track
void PlotModuleMult() {

    Simulation *sim = filereader.tmSim();
    Int_t entriesSim = sim->GetInputTree()->GetEntries();

    TH1F* h = new TH1F("h","Hit Module Multiplicity (per Track)",10,0,10);
    set<uint32_t> mods;//use set so as to not double count modules

    //loop over entries (one track per entry)
    for ( Int_t ientry = 0; ientry < entriesSim; ientry++ ) {
        sim->GetSim(ientry);
	//loop over CRT strips hit by the track
        for ( Int_t i=0; i< sim->NCRTHit(); i++) {
             mods.insert(sim->CRTModule(i));
        }
        if (mods.size() > 10 ) cout << "too many modules found!" << endl;
        if (mods.size()>0) h->Fill(mods.size()); //only fill histo when module(s) is(are) hit
        mods.clear();//clear for the next track
    }

    h->SetLineWidth(2);
    h->GetXaxis()->SetTitle("multiplicity");

    h->Draw();
    
}

//Plot the true hit module multiplicity per track per CRT subsystem
void PlotSubSystemModMult() {

    Simulation *sim = filereader.tmSim();
    Int_t entriesSim = sim->GetInputTree()->GetEntries();

    //one histogram per subsystm
    TH1F* hc = new TH1F("hc","Hit Module Multiplicity",10,0,10);
    TH1F* hm = new TH1F("hm","Hit Module Multiplicity",10,0,10);
    TH1F* hd = new TH1F("hd","Hit Module Multiplicity",10,0,10);
    set<Int_t> mods, modsC, modsM, modsD;

    // loop over events/trackIDs
    for ( Int_t ientry = 0; ientry < entriesSim; ientry++ ) {

        sim->GetSim(ientry);

        // loop over sim data for track
        for ( Int_t i=0; i< sim->NCRTHit(); i++) {

             Int_t id = sim->CRTModule(i);
             Int_t type = sim->CRTType(i);
             mods.insert(id);

	     switch(type) {
		 case 0:
                     modsC.insert(id);
                     break;
	         case 1:
                     modsM.insert(id);
		     break;
		 case 2:
		     modsD.insert(id);
		     break;
	     }
        } //for track hits
        
        // if at least one module hit
        if (mods.size()>0) {
            if(modsC.size()>0) hc->Fill(modsC.size());
            if(modsM.size()>0) hm->Fill(modsM.size());
            if(modsD.size()>0) hd->Fill(modsD.size());
        }

        mods.clear();
        modsC.clear();
        modsM.clear();
        modsD.clear();
    } // for events/tracks

    hc->SetLineColor(kBlue);
    hm->SetLineColor(kRed);
    hd->SetLineColor(kMagenta);

    hc->SetLineWidth(2);
    hm->SetLineWidth(2);
    hd->SetLineWidth(2);

    hc->GetXaxis()->SetTitle("multiplicity");

    hc->Draw();
    hm->Draw("same");
    hd->Draw("same");

}

//Plot the number of strips hit per module per track
// if scale = true, histograms are area normalized
void PlotStripMult(Bool_t splitlat=kFALSE, Bool_t scale=kTRUE) {
    Simulation *sim = filereader.tmSim();
    Int_t entriesSim = sim->GetInputTree()->GetEntries();

    TH1F* htop = new TH1F("htop","Hit Strips Multiplicity (per Module)",20,0,20);
    TH1F* hslat = new TH1F("hslat","Hit Strips Multiplicity (per Module)",20,0,20);
    TH1F* hslong = new TH1F("hslong","Hit Strips Multiplicity (per Module)",20,0,20);
    TH1F* hlat = new TH1F("hlat","Hit Strips Multiplicity (per Module)",20,0,20);
    TH1F* hlong = new TH1F("hlong","Hit Strips Multiplicity (per Module)",20,0,20);
    TH1F* hd = new TH1F("hd","Hit Strips Multiplicity (per Module)",20,0,20);

    map<uint32_t,set<uint32_t>> hitsTop;
    map<uint32_t,set<uint32_t>> hitsSLong;
    map<uint32_t,set<uint32_t>> hitsSLat;
    map<uint32_t,set<uint32_t>> hitsLat;
    map<uint32_t,set<uint32_t>> hitsLong;
    map<uint32_t,set<uint32_t>> hitsd;

    // 1 sim entry / trackID
    for ( Int_t ientry = 0; ientry < entriesSim; ientry++ ) {
        sim->GetSim(ientry);
        for ( Int_t i=0; i< sim->NCRTHit(); i++) {
            if (sim->CRTRegion(i)==30)
               hitsTop[sim->CRTModule(i)].insert(sim->CRTStrip(i));
            if (sim->CRTRegion(i)==31 || sim->CRTRegion(i)==32) 
               hitsSLat[sim->CRTModule(i)].insert(sim->CRTStrip(i));
            if (sim->CRTRegion(i)==33 || sim->CRTRegion(i)==34)
               hitsSLong[sim->CRTModule(i)].insert(sim->CRTStrip(i));
            if (sim->CRTRegion(i)>=40 && sim->CRTRegion(i)<=45) 
               hitsLat[sim->CRTModule(i)].insert(sim->CRTStrip(i));
            if (sim->CRTRegion(i)==46 || sim->CRTRegion(i)==47)
               hitsLong[sim->CRTModule(i)].insert(sim->CRTStrip(i));
            if (sim->CRTRegion(i)==50) hitsd[sim->CRTModule(i)].insert(sim->CRTStrip(i));
        }
        for (auto hit : hitsTop) {
            htop->Fill(hit.second.size());
        }
        for (auto hit : hitsSLat) {
            hslat->Fill(hit.second.size());
        }
        for (auto hit : hitsSLong) {
            hslong->Fill(hit.second.size());
        }
        for (auto hit : hitsLat) {
            hlat->Fill(hit.second.size());
        }
        for (auto hit : hitsLong) {
            hlong->Fill(hit.second.size());
        }
        for (auto hit : hitsd) {
            hd->Fill(hit.second.size());
        }

        hitsTop.clear();
        hitsSLat.clear();
        hitsSLong.clear();
        hitsLat.clear();
        hitsLong.clear();
        hitsd.clear();
    }

    if (!splitlat) {
      hslat->Add(hslong);
      hlat->Add(hlong);
    }

    htop->SetLineColor(kRed);
    hslat->SetLineColor(kOrange);
    if(splitlat) hslong->SetLineColor(kYellow);
    hlat->SetLineColor(kBlue);
    if(splitlat) hlong->SetLineColor(kCyan);
    hd->SetLineColor(kMagenta);

    htop->SetLineWidth(2);
    hslat->SetLineWidth(2);
    if(splitlat) hslong->SetLineWidth(2);
    hlat->SetLineWidth(2);
    if(splitlat) hlong->SetLineWidth(2);
    hd->SetLineWidth(2);

    htop->GetXaxis()->SetTitle("strip multiplicity");

    if(scale) {
      htop->Scale(1.0/htop->Integral());
      hslat->Scale(1.0/hslat->Integral());
      if(splitlat) hslong->Scale(1.0/hslong->Integral());
      hlat->Scale(1.0/hlat->Integral());
      if(splitlat) hlong->Scale(1.0/hlong->Integral());
      hd->Scale(1.0/hd->Integral());
    }

    htop->Draw("hist");
    hslat->Draw("same,hist");
    if(splitlat) hslong->Draw("same,hist");
    hlat->Draw("same,hist");
    if(splitlat) hlong->Draw("same,hist");
    hd->Draw("same,hist");

    TLegend *leg = new TLegend(0.75,0.7,0.95,0.95);
    leg->AddEntry(htop,"top","l");
    if(splitlat) {
      leg->AddEntry(hslat,"slope lat","l");
      leg->AddEntry(hslong,"slope long","l");
    }
    else
      leg->AddEntry(hslat,"slope","l");
    if(splitlat) {
      leg->AddEntry(hlat,"side lat","l");
      leg->AddEntry(hlong,"side long","l");
    }
    else
      leg->AddEntry(hlat,"side","l");
    leg->AddEntry(hd,"bottom","l");
    leg->Draw();
}

//Plot hit strip IDs vs. module IDs 
void PlotModuleStrips() {
    Simulation *sim = filereader.tmSim();
    Int_t entriesSim = sim->GetInputTree()->GetEntries();

    TH2F* h = new TH2F("h","Hit Strips",300,0,300,65,0,65);

    for ( Int_t ientry = 0; ientry < entriesSim; ientry++ ) {
        sim->GetSim(ientry);

        for ( Int_t i=0; i< sim->NCRTHit(); i++) {
            h->Fill(sim->CRTModule(i),sim->CRTStrip(i));
        }    
    }

    h->GetXaxis()->SetTitle("module ID");
    h->GetYaxis()->SetTitle("strip ID");
    h->Draw("colz");
}

//Plot deposited energy per CRT strip per track by CRT subsystem 
void PlotSubSystemEDep() {

    Simulation *sim = filereader.tmSim();
    Int_t entriesSim = sim->GetInputTree()->GetEntries();

    TH1F* hc = new TH1F("hc","True Deposited Energy",100,0,10);
    TH1F* hm = new TH1F("hm","True Deposited Energy",100,0,10);
    TH1F* hd = new TH1F("hd","True Deposited Energy",100,0,10);

    for ( Int_t ientry = 0; ientry < entriesSim; ientry++ ) {

    	    sim->GetSim(ientry);

        for ( Int_t i=0; i< sim->NCRTHit(); i++) {

            switch (sim->CRTType(i)) {
                case 0 : //c
                    hc->Fill(sim->CRTEDep(i)*1000); //fill with energy converted to MeV
                    break;
                case 1 : //m
                    hm->Fill(sim->CRTEDep(i)*1000);
                    break;
                case 2 : //d
                    hd->Fill(sim->CRTEDep(i)*1000);
                    break;
            }//switch type
        }//for crt hits
    }//for entries

    hc->GetXaxis()->SetTitle("deposited energy per strip (MeV)");

    hc->SetLineColor(kBlue);
    hm->SetLineColor(kRed);
    hd->SetLineColor(kMagenta);

    hc->SetLineWidth(2);
    hm->SetLineWidth(2);
    hd->SetLineWidth(2);

    hc->Draw();
    hm->Draw("sames");
    hd->Draw("sames");

}

///////////////////////////////////
// Regions Tree specific functions
///////////////////////////////////

//calculate the fraction of muons entering through a particular CRT region
void RegRates() {

    const bool printcounts = false;

    struct eveid{
        int eve, run, subrun;

    };
    auto cmp_eveid = [](const eveid& lhs, const eveid& rhs) {
        if( lhs.run != rhs.run ) return lhs.run < rhs.run;
        if( lhs.subrun != rhs.subrun ) return lhs.subrun < rhs.subrun;
        return lhs.eve < rhs.eve;
    };

    Regions* reg = filereader.tmReg();
    const Int_t nentries = reg->GetInputTree()->GetEntries();
    cout << "reading from RegionsTree with " << nentries << " entries..." << endl;

    Int_t n_crt_nolar=0, n_crt_noav=0, n_crt_iv=0, n_crt_av=0, n_crt_fv=0;
    Int_t n = 0, n_crt=0, n_noav = 0, n_av = 0, n_fv = 0;
    set<eveid,decltype(cmp_eveid)> events(cmp_eveid);
    map<int,int> regcounts_all;
    map<int,int> regcounts_nolar;
    map<int,int> regcounts_noav;
    map<int,int> regcounts_av;
    map<int,int> regcounts_fv;

    map<int,double> regrates_all  , regfrac_all;
    map<int,double> regrates_nolar, regfrac_nolar;
    map<int,double> regrates_noav , regfrac_noav;
    map<int,double> regrates_av   , regfrac_av;
    map<int,double> regrates_fv   , regfrac_fv;

    int regs[14] = {30,31,32,33,34,40,41,42,43,44,45,46,47,50};
    for(int i=0; i<14; i++) {
        regcounts_all[regs[i]] = 0;
        regcounts_nolar[regs[i]] = 0;
        regcounts_noav[regs[i]] = 0;
        regcounts_av[regs[i]] = 0;
        regcounts_fv[regs[i]] = 0;
    }

    // loop over tracks
    for ( Int_t ientry = 0; ientry < nentries; ientry++ ) {

        struct eveid eve = {reg->Event(),reg->Run(),reg->SubRun()};
        events.insert(eve);

        reg->GetRegion(ientry);
        // if muon doesn't hit CRT or LAr, then skip
        if (reg->Active()==0&&reg->CRTs()==0&&reg->Inactive()==0) continue;

        // if not a muon, skip
        if ( abs(reg->PDG()) != 13 ) continue;

        //loop over regions, get all CRT regions hit by track
        vector<int> muregs;
        for ( Int_t hit=0; hit<reg->NReg(); hit++ ) {

            // skip non-CRT regions
            if (reg->Region(hit) < 30) continue;

            muregs.push_back(reg->Region(hit));

        } // end loop over regions

        n++;
        if(reg->CRTs()>0){
            n_crt++;
            regcounts_all[muregs[0]]++;;
            if(reg->Active()==0 && reg->Inactive()==0){
                n_crt_nolar++;
                regcounts_nolar[muregs[0]]++;
            }
            if(reg->Active()==0 && reg->Inactive()>0){
                n_crt_noav++;
                regcounts_noav[muregs[0]]++;
            }
            if(reg->Active()>0){ 
                n_crt_av++;
                regcounts_av[muregs[0]]++;
            }
            if(reg->Fiducial()>0){
                n_crt_fv++;
                regcounts_fv[muregs[0]]++;
            }
        }
        if(reg->Active()==0 && reg->Inactive()>0) n_noav++;
        if(reg->Active()>0)   n_av++;
        if(reg->Fiducial()>0) n_fv++;

        muregs.clear();

    }// loop over events/tracks

    double sampleTime = events.size()*3.3e-3;
    int ntop[5], nside[5];
    double ntop_rate[5], ntop_frac[5];
    double nside_rate[5], nside_frac[5];
    for(int i=0; i<5; i++) {ntop[i]=0; nside[i]=0;}
    for(int i=30; i<35; i++){
        ntop[0]+=regcounts_all[i];
        ntop[1]+=regcounts_nolar[i];
        ntop[2]+=regcounts_noav[i];
        ntop[3]+=regcounts_av[i];
        ntop[4]+=regcounts_fv[i];
    }
    for(int i=40; i<48; i++){
        nside[0]+=regcounts_all[i];
        nside[1]+=regcounts_nolar[i];
        nside[2]+=regcounts_noav[i];
        nside[3]+=regcounts_av[i];
        nside[4]+=regcounts_fv[i];
    }
    for(int i=0; i<5; i++){
        ntop_rate[i] = 1.e-3*ntop[i]/sampleTime;
        nside_rate[i] = 1.e-3*nside[i]/sampleTime;
    }
    ntop_frac[0]  = 1.*ntop[0]/n_crt;
    ntop_frac[1]  = 1.*ntop[1]/n_crt_nolar;
    ntop_frac[2]  = 1.*ntop[2]/n_crt_noav;
    ntop_frac[3]  = 1.*ntop[3]/n_crt_av;
    ntop_frac[4]  = 1.*ntop[4]/n_crt_fv;
    nside_frac[0] = 1.*nside[0]/n_crt;
    nside_frac[1] = 1.*nside[1]/n_crt_nolar;
    nside_frac[2] = 1.*nside[2]/n_crt_noav;
    nside_frac[3] = 1.*nside[3]/n_crt_av;
    nside_frac[4] = 1.*nside[4]/n_crt_fv;

    for(int i=0; i<14; i++) {
        regrates_all[regs[i]]   = 1.e-3*regcounts_all[regs[i]]  /sampleTime; 
        regrates_nolar[regs[i]] = 1.e-3*regcounts_nolar[regs[i]]/sampleTime;
        regrates_noav[regs[i]]  = 1.e-3*regcounts_noav[regs[i]] /sampleTime;  
        regrates_av[regs[i]]    = 1.e-3*regcounts_av[regs[i]]   /sampleTime; 
        regrates_fv[regs[i]]    = 1.e-3*regcounts_fv[regs[i]]   /sampleTime; 

        regfrac_all[regs[i]]   = 1.*regcounts_all[regs[i]]  /n_crt;
        regfrac_nolar[regs[i]] = 1.*regcounts_nolar[regs[i]]/n_crt_nolar;
        regfrac_noav[regs[i]]  = 1.*regcounts_noav[regs[i]] /n_crt_noav;
        regfrac_av[regs[i]]    = 1.*regcounts_av[regs[i]]   /n_crt_av;
        regfrac_fv[regs[i]]    = 1.*regcounts_fv[regs[i]]   /n_crt_fv;
    }


    cout << '\n' << " found " << events.size() << " events amounting to " << sampleTime << " seconds" << '\n' 
     << endl;

    if(printcounts)
    cout << "  ---   counts   --- " << '\n'
         << "   muons hitting CRT and/or LAr:          " << n << '\n'
         << "   muons hitting CRT:                     " << n_crt << '\n'
         << "   muons hitting CRT but miss LAr:        " << n_crt_nolar << '\n'
         << "   muons hitting CRT and IV but miss AV:  " << n_crt_noav << '\n'
         << "   muons hitting IV but not AV:           " << n_noav << '\n'
         << "   muons hitting AV:                      " << n_av << '\n'
         << "   muons hitting FV:                      " << n_fv << '\n'
     << endl;

    cout << " --- tagging efficiencies ---" << '\n'
         << "   muons hitting IV but not AV:    " << 1.*n_crt_noav/n_noav << '\n'
         << "   muons hitting AV:               " << 1.*n_crt_av/n_av << '\n'
         << "   muons hitting FV:               " << 1.*n_crt_fv/n_fv << '\n'
     << endl;

    cout << " --- fraction of tracks hitting the CRT that... ---" << '\n'
         << "   hit the CRT only:  " << 1.*n_crt_nolar/n_crt << '\n'
         << "   hit IV but not AV: " << 1.*n_crt_noav/n_crt << '\n'
         << "   hit AV:            " << 1.*n_crt_av/n_crt << '\n'
         << "   hit FV:            " << 1.*n_crt_fv/n_crt << '\n'
     << endl;

    if(printcounts)
    cout << "  --- first hit CRT region counts ---"  << '\n'
         <<"            | all | CRT only | CRT+IV only | CRT+AV | CRT+FV " << '\n'
         << "   top:    " << ntop[0] << " | " << ntop[1] << " | " << ntop[2] << " | " << ntop[3] << " | " << ntop[4]  << '\n'
         << "      roof:   " << regcounts_all[30] << " | " << regcounts_nolar[30] << " | " << regcounts_noav[30] << " | " << regcounts_av[30] << " | " << regcounts_fv[30] << '\n'
         << "      N rim:  " << regcounts_all[34] << " | " << regcounts_nolar[34] << " | " << regcounts_noav[34] << " | " << regcounts_av[34] << " | " << regcounts_fv[34] << '\n'
         << "      S rim:  " << regcounts_all[33] << " | " << regcounts_nolar[33] << " | " << regcounts_noav[33] << " | " << regcounts_av[33] << " | " << regcounts_fv[33] << '\n'
         << "      W rim:  " << regcounts_all[31] << " | " << regcounts_nolar[31] << " | " << regcounts_noav[31] << " | " << regcounts_av[31] << " | " << regcounts_fv[31] << '\n'
         << "      E rim:  " << regcounts_all[32] << " | " << regcounts_nolar[32] << " | " << regcounts_noav[32] << " | " << regcounts_av[32] << " | " << regcounts_fv[32] << '\n'
         << "   side:   " << nside[0] << " | " << nside[1] << " | " << nside[2] << " | " << nside[3] << " | " << nside[4] << '\n'
         << "      S wall: " << regcounts_all[46] << " | " << regcounts_nolar[46] << " | " << regcounts_noav[46] << " | " << regcounts_av[46] << " | " << regcounts_fv[46] << '\n'
         << "      N wall: " << regcounts_all[47] << " | " << regcounts_nolar[47] << " | " << regcounts_noav[47] << " | " << regcounts_av[47] << " | " << regcounts_fv[47] << '\n'
         << "      W-N:    " << regcounts_all[42] << " | " << regcounts_nolar[42] << " | " << regcounts_noav[42] << " | " << regcounts_av[42] << " | " << regcounts_fv[42] << '\n'
         << "      W-C:    " << regcounts_all[41] << " | " << regcounts_nolar[41] << " | " << regcounts_noav[41] << " | " << regcounts_av[41] << " | " << regcounts_fv[41] << '\n'
         << "      W-S:    " << regcounts_all[40] << " | " << regcounts_nolar[40] << " | " << regcounts_noav[40] << " | " << regcounts_av[40] << " | " << regcounts_fv[40] << '\n'
         << "      E-N:    " << regcounts_all[45] << " | " << regcounts_nolar[45] << " | " << regcounts_noav[45] << " | " << regcounts_av[45] << " | " << regcounts_fv[45] << '\n'
         << "      E-C:    " << regcounts_all[44] << " | " << regcounts_nolar[44] << " | " << regcounts_noav[44] << " | " << regcounts_av[44] << " | " << regcounts_fv[44] << '\n'
         << "      E-S:    " << regcounts_all[43] << " | " << regcounts_nolar[43] << " | " << regcounts_noav[43] << " | " << regcounts_av[43] << " | " << regcounts_fv[43] << '\n'
         << "   bottom: "    << regcounts_all[50] << " | " << regcounts_nolar[50] << " | " << regcounts_noav[50] << " | " << regcounts_av[50] << " | " << regcounts_fv[50] << '\n'
     << endl;

    cout << " --- rates [kHz] --- " << '\n'
         << "   muons hitting CRT and/or LAr:   " << 1.e-3*n          /sampleTime << '\n'
         << "   muons hitting CRT:              " << 1.e-3*n_crt      /sampleTime << '\n'
         << "   muons hitting CRT but miss LAr: " << 1.e-3*n_crt_nolar/sampleTime << '\n'
         << "   muons hitting CRT but miss AV:  " << 1.e-3*n_crt_noav /sampleTime << '\n'
         << "   muons hitting IV but not AV:    " << 1.e-3*n_noav     /sampleTime << '\n'
         << "   muons hitting AV:               " << 1.e-3*n_av       /sampleTime << '\n'
         << "   muons hitting FV:               " << 1.e-3*n_fv       /sampleTime << '\n'
         <<"           | all | CRT only | CRT+IV only | CRT+AV | CRT+FV " << '\n'
         << "   top:    " << ntop_rate[0] << " | " << ntop_rate[1] << " | " << ntop_rate[2] << " | " << ntop_rate[3] << " | " << ntop_rate[4]  << '\n'
         << "      roof:   " << regrates_all[30] << " | " << regrates_nolar[30] << " | " << regrates_noav[30] << " | " << regrates_av[30] << " | " << regrates_fv[30] << '\n'
         << "      N rim:  " << regrates_all[34] << " | " << regrates_nolar[34] << " | " << regrates_noav[34] << " | " << regrates_av[34] << " | " << regrates_fv[34] << '\n'
         << "      S rim:  " << regrates_all[33] << " | " << regrates_nolar[33] << " | " << regrates_noav[33] << " | " << regrates_av[33] << " | " << regrates_fv[33] << '\n'
         << "      W rim:  " << regrates_all[31] << " | " << regrates_nolar[31] << " | " << regrates_noav[31] << " | " << regrates_av[31] << " | " << regrates_fv[31] << '\n'
         << "      E rim:  " << regrates_all[32] << " | " << regrates_nolar[32] << " | " << regrates_noav[32] << " | " << regrates_av[32] << " | " << regrates_fv[32] << '\n'
         << "   side:   " << nside_rate[0] << " | " << nside_rate[1] << " | " << nside_rate[2] << " | " << nside_rate[3] << " | " << nside_rate[4] << '\n'
         << "      S wall: " << regrates_all[46] << " | " << regrates_nolar[46] << " | " << regrates_noav[46] << " | " << regrates_av[46] << " | " << regrates_fv[46] << '\n'
         << "      N wall: " << regrates_all[47] << " | " << regrates_nolar[47] << " | " << regrates_noav[47] << " | " << regrates_av[47] << " | " << regrates_fv[47] << '\n'
         << "      W-N:    " << regrates_all[42] << " | " << regrates_nolar[42] << " | " << regrates_noav[42] << " | " << regrates_av[42] << " | " << regrates_fv[42] << '\n'
         << "      W-C:    " << regrates_all[41] << " | " << regrates_nolar[41] << " | " << regrates_noav[41] << " | " << regrates_av[41] << " | " << regrates_fv[41] << '\n'
         << "      W-S:    " << regrates_all[40] << " | " << regrates_nolar[40] << " | " << regrates_noav[40] << " | " << regrates_av[40] << " | " << regrates_fv[40] << '\n'
         << "      E-N:    " << regrates_all[45] << " | " << regrates_nolar[45] << " | " << regrates_noav[45] << " | " << regrates_av[45] << " | " << regrates_fv[45] << '\n'
         << "      E-C:    " << regrates_all[44] << " | " << regrates_nolar[44] << " | " << regrates_noav[44] << " | " << regrates_av[44] << " | " << regrates_fv[44] << '\n'
         << "      E-S:    " << regrates_all[43] << " | " << regrates_nolar[43] << " | " << regrates_noav[43] << " | " << regrates_av[43] << " | " << regrates_fv[43] << '\n'
         << "   bottom: "    << regrates_all[50] << " | " << regrates_nolar[50] << " | " << regrates_noav[50] << " | " << regrates_av[50] << " | " << regrates_fv[50] << '\n'
     << endl;

    cout << " --- entry fractions ---"<< '\n'
         <<"          | all | CRT only | CRT+IV only | CRT+AV | CRT+FV " << '\n'
         << "   top:    " << ntop_frac[0] << " | " << ntop_frac[1] << " | " << ntop_frac[2] << " | " << ntop_frac[3] << " | " << ntop_frac[4]  << '\n'
         << "      roof:   " << regfrac_all[30] << " | " << regfrac_nolar[30] << " | " << regfrac_noav[30] << " | " << regfrac_av[30] << " | " << regfrac_fv[30] << '\n'
         << "      N rim:  " << regfrac_all[34] << " | " << regfrac_nolar[34] << " | " << regfrac_noav[34] << " | " << regfrac_av[34] << " | " << regfrac_fv[34] << '\n'
         << "      S rim:  " << regfrac_all[33] << " | " << regfrac_nolar[33] << " | " << regfrac_noav[33] << " | " << regfrac_av[33] << " | " << regfrac_fv[33] << '\n'
         << "      W rim:  " << regfrac_all[31] << " | " << regfrac_nolar[31] << " | " << regfrac_noav[31] << " | " << regfrac_av[31] << " | " << regfrac_fv[31] << '\n'
         << "      E rim:  " << regfrac_all[32] << " | " << regfrac_nolar[32] << " | " << regfrac_noav[32] << " | " << regfrac_av[32] << " | " << regfrac_fv[32] << '\n'
         << "   side:   " << nside_frac[0] << " | " << nside_frac[1] << " | " << nside_frac[2] << " | " << nside_frac[3] << " | " << nside_frac[4] << '\n'
         << "      S wall: " << regfrac_all[46] << " | " << regfrac_nolar[46] << " | " << regfrac_noav[46] << " | " << regfrac_av[46] << " | " << regfrac_fv[46] << '\n'
         << "      N wall: " << regfrac_all[47] << " | " << regfrac_nolar[47] << " | " << regfrac_noav[47] << " | " << regfrac_av[47] << " | " << regfrac_fv[47] << '\n'
         << "      W-N:    " << regfrac_all[42] << " | " << regfrac_nolar[42] << " | " << regfrac_noav[42] << " | " << regfrac_av[42] << " | " << regfrac_fv[42] << '\n'
         << "      W-C:    " << regfrac_all[41] << " | " << regfrac_nolar[41] << " | " << regfrac_noav[41] << " | " << regfrac_av[41] << " | " << regfrac_fv[41] << '\n'
         << "      W-S:    " << regfrac_all[40] << " | " << regfrac_nolar[40] << " | " << regfrac_noav[40] << " | " << regfrac_av[40] << " | " << regfrac_fv[40] << '\n'
         << "      E-N:    " << regfrac_all[45] << " | " << regfrac_nolar[45] << " | " << regfrac_noav[45] << " | " << regfrac_av[45] << " | " << regfrac_fv[45] << '\n'
         << "      E-C:    " << regfrac_all[44] << " | " << regfrac_nolar[44] << " | " << regfrac_noav[44] << " | " << regfrac_av[44] << " | " << regfrac_fv[44] << '\n'
         << "      E-S:    " << regfrac_all[43] << " | " << regfrac_nolar[43] << " | " << regfrac_noav[43] << " | " << regfrac_av[43] << " | " << regfrac_fv[43] << '\n'
         << "   bottom: "    << regfrac_all[50] << " | " << regfrac_nolar[50] << " | " << regfrac_noav[50] << " | " << regfrac_av[50] << " | " << regfrac_fv[50] << '\n'
     << endl;


}

//Dump region info for a given event
void DumpRegionEvent(Int_t event=-1) {

        Regions *reg = filereader.tmReg();
        Long64_t entriesReg = reg->GetInputTree()->GetEntries();
        cout << "Reading all regions " << entriesReg << endl;

        cout << "-------------------------------------------------------" << '\n'
             << "| region | time |" << '\n'
             << "-------------------------------------------------------" << endl;

        int ntracks = 0, nav = 0;
        bool isav = false;
        for( Long64_t ientry=0; ientry<entriesReg; ++ientry )
        {
                reg->GetRegion(ientry);
                if (event!=-1&&reg->Event() != event) continue;
                cout << "track " << reg->TrackID() << " found with " << reg->NReg() << " entries" << endl;
                ntracks++;
                isav = false;
                for (Int_t hit=0; hit<reg->NReg(); hit++) {
                  cout << "-------------------------------------------------------" << '\n'
                       << " | " << reg->Region(hit)
                       << " | " << reg->EntryT(hit) << " |" << endl;
                  if (reg->Region(hit)>0&&reg->Region(hit)<10) isav = true;
                }
                if (isav) nav++;
        }
        cout << ntracks << " tracks in event " << event << endl;
        cout << " of which " << nav << " pass through LArAV" << endl;
}


////////////////////////////////////
///  Det Tree specific functions ///
////////////////////////////////////

//Plot charge for every channel, event
void PlotADC() {

    DetSim *det = filereader.tmDet();
    Int_t entriesDet = det->GetInputTree()->GetEntries();
    TH1F* h = new TH1F("h","ADC&",4100/10,0,4100);

    for (Int_t ientry = 0; ientry < entriesDet; ientry++) {

    	det->GetDet(ientry);

	size_t maxchan=64; //max no. chans for bottom
	if(det->SubSys()!=2) maxchan=32;//top and sides on have 32 channels

	//loop over all channels
        for (size_t chan=0; chan<maxchan;  chan++) {
            if(det->ADC(chan)!=0) //zero suppress
		h->Fill(det->ADC(chan));
        }
    }

    h->Draw();
}

//Plot CRT ADC spectra by subsystem
void PlotSubSysADC() {

    DetSim *det = filereader.tmDet();
    Int_t entriesDet = det->GetInputTree()->GetEntries();

    TH1F* hc = new TH1F("hc","Charge Spectra",4100/20,0,4100); //top
    TH1F* hm = new TH1F("hm","ADC by System", 4100/20,0,4100);  //side
    TH1F* hd = new TH1F("hd","ADC by System", 4100/20,0,4100);  //bottom

    //loop over entries
    for (Int_t ientry = 0; ientry < entriesDet; ientry++) {

        det->GetDet(ientry);

	size_t maxchan=64; //max no. chans for bottom
        if(det->SubSys()!=2) maxchan=32;//top and sides on have 32 channels

        for (size_t chan=0; chan<maxchan; chan++) {
	    if(det->ADC(chan)==0) continue; //zero suppress
            switch (det->SubSys()) {
                case 0 :
                    hc->Fill(det->ADC(chan));
                    break;
                case 1 :
                    hm->Fill(det->ADC(chan));
                    break;
                case 2 :
                    hd->Fill(det->ADC(chan));
                    break;
            }
        }
    }
   
    hc->GetXaxis()->SetTitle("ADC");

    hc->SetLineColor(kBlue);
    hm->SetLineColor(kRed);
    hd->SetLineColor(kMagenta);

    hc->SetLineWidth(2);
    hm->SetLineWidth(2);
    hd->SetLineWidth(2);

    auto legend = new TLegend(0.53,0.7,0.75,0.9);
    legend->SetHeader("CRT Subsystem","C"); // option "C" allows to center the header
    legend->AddEntry(hc,"CERN","lc");
    legend->AddEntry(hm,"MINOS","lc");
    legend->AddEntry(hd,"DC","lc");

    auto const* c = new TCanvas();
    hc->Draw();
    hm->Draw("sames");
    hd->Draw("sames");
    legend->Draw();

}

/////////////////////////////////////////////
// Hit Tree functions
/////////////////////////////////////////////

void PlotHitXY() {
    Hit *hit = filereader.tmHit();
    Int_t entriesHit = hit->GetInputTree()->GetEntries();

    TH2F *h = new TH2F("h","Reconstructed CRT Hits",600,-600,600,550,-400,700);

    for (Int_t ientry=0; ientry<entriesHit; ientry++) {
        hit->GetHit(ientry);
        h->Fill(hit->X(),hit->Y()); 
    }

    h->GetXaxis()->SetTitle("Reco X (cm)");
    h->GetYaxis()->SetTitle("Reco Y (cm)");
    h->SetMarkerStyle(8);
    h->SetMarkerColor(kBlue);
    h->Draw();

}

void PlotHitYZ() {
    Hit *hit = filereader.tmHit();
    Int_t entriesHit = hit->GetInputTree()->GetEntries();

    TH2F *h = new TH2F("h","Reconstructed CRT Hits",1000,-1000,1000,550,-400,700);

    for (Int_t ientry=0; ientry<entriesHit; ientry++) {
        hit->GetHit(ientry);
        h->Fill(hit->Z(),hit->Y());
    }

    h->GetXaxis()->SetTitle("Reco Z (cm)");
    h->GetYaxis()->SetTitle("Reco Y (cm)");
    h->SetMarkerStyle(8);
    h->SetMarkerColor(kBlue);
    h->Draw();

}

void PlotHitXZ() {
    Hit *hit = filereader.tmHit();
    Int_t entriesHit = hit->GetInputTree()->GetEntries();

    TH2F *h = new TH2F("h","Reconstructed CRT Hits",1000,-1000,1000,600,-600,600);

    for (Int_t ientry=0; ientry<entriesHit; ientry++) {
        hit->GetHit(ientry);
        h->Fill(hit->Z(),hit->X()); 
    }

    h->GetXaxis()->SetTitle("Reco Z (cm)");
    h->GetYaxis()->SetTitle("Reco X (cm)");
    h->SetMarkerStyle(8);
    h->SetMarkerColor(kBlue);
    h->Draw();

}

void PlotRecoHitByRegion(Int_t reg=0) {

    Hit* hit = filereader.tmHit();
    Int_t entries = hit->GetInputTree()->GetEntries();

    TH2F *hzx = new TH2F("hzx","Reco CRT Hits",1000,-1200,1200,600,-600,600);
    TH2F *hzy = new TH2F("hzy","Reco CRT Hits",1000,-1020,1200,550,-400,700);
    TH2F *hxy = new TH2F("hxy","Reco CRT Hits",600,-650,650,550,-400,700);

    for (Int_t ientry=0; ientry<entries; ientry++) {

        hit->GetHit(ientry);

        if (hit->Region()!=reg) continue;

        hzx->Fill(hit->Z(),hit->X());
        hzy->Fill(hit->Z(),hit->Y());
        hxy->Fill(hit->X(),hit->Y());

    }

    hxy->GetXaxis()->SetTitle("X (cm)");
    hxy->GetYaxis()->SetTitle("Y (cm)");

    hzy->GetXaxis()->SetTitle("Z (cm)");
    hzy->GetYaxis()->SetTitle("Y (cm)");

    hzx->GetXaxis()->SetTitle("Z (cm)");
    hzx->GetYaxis()->SetTitle("X (cm)");

    hxy->SetMarkerStyle(8);
    hzx->SetMarkerStyle(8);
    hzy->SetMarkerStyle(8);

    hxy->SetMarkerColor(kBlue);
    hzy->SetMarkerColor(kBlue);
    hzx->SetMarkerColor(kBlue);

    new TCanvas();
    hxy->Draw();

    new TCanvas();
    hzy->Draw();

    new TCanvas();
    hzx->Draw();
}

void HitTrackIDs() {

    Hit* hit = filereader.tmHit();
    Int_t nentries = hit->GetInputTree()->GetEntries();

    TH1F *h = new TH1F("h","multiplicity of tracks associated with reco hit",20,0,20);

    for ( Int_t ientry=0; ientry<nentries; ientry++ ) {

        hit->GetHit(ientry);
        h->Fill(hit->NTrack());
    }

    h->Draw();
}

void NonMuCRTHits() {

    Hit* hit = filereader.tmHit();
    Int_t nentries = hit->GetInputTree()->GetEntries();

    TH1F* h = new TH1F("h","non-mu PDGs in CRTHits",6000,-3000,3000);
    Int_t nNoMu = 0;

    for ( Int_t ientry = 0; ientry<64; ientry++ ) {

        hit->GetHit(ientry);

        bool noMu=true;
        for (size_t i=0; i<hit->NTrack(); i++){
            if ( abs(hit->PDG(i))==13 ) {
                noMu = false;
                break;
            }
        }
        if ( !noMu ) continue;

        for (size_t i=0; i<hit->NTrack(); i++){
            h->Fill(hit->PDG(i));
        }
        nNoMu++;

    }

    cout << " non-muon hits: " << nNoMu << " (" << 100.0*nNoMu/nentries << "%)" << endl;

    h->Draw();

}


//////////////////////////////////////
// TrueCRTHit tree functions
/////////////////////////////////////

void PlotTrueCRTHit() {

    TrueHit* trueHit = filereader.tmTrueHit();
    Int_t entriesTrue = trueHit->GetInputTree()->GetEntries();

    TH1F *hresx_c = new TH1F("hresx_c","True CRT CERN Hit Position RMS",1000, 0,100);
    TH1F *hresy_c = new TH1F("hresy_c","True CRT CERN Hit Position RMS",1000, 0,100);
    TH1F *hresz_c = new TH1F("hresz_c","True CRT CERN Hit Position RMS",1000, 0,100);

    TH1F *hresx_m = new TH1F("hresx_m","True CRT MINOS Hit Position RMS",1000, 0,100);
    TH1F *hresy_m = new TH1F("hresy_m","True CRT MINOS Hit Position RMS",1000, 0,100);
    TH1F *hresz_m = new TH1F("hresz_m","True CRT MINOS Hit Position RMS",1000, 0,100);

    TH1F *hresx_d = new TH1F("hresx_d","True CRT DC Hit Position RMS",1000, 0,100);
    TH1F *hresy_d = new TH1F("hresy_d","True CRT DC Hit Position RMS",1000, 0,100);
    TH1F *hresz_d = new TH1F("hresz_d","True CRT DC Hit Position RMS",1000, 0,100);

    TH2F *hzx = new TH2F("hzx","True CRT Hits",1000,-1000,1000,600,-600,600);
    TH2F *hzy = new TH2F("hzy","True CRT Hits",1000,-1000,1000,550,-400,700);
    TH2F *hxy = new TH2F("hxy","True CRT Hits",600,-600,600,550,-400,700);

    for (Int_t ientry=0; ientry<entriesTrue; ientry++) {
        trueHit->GetHit(ientry);

        if (trueHit->SubSys()==0) {
            hresx_c->Fill(trueHit->XErr());
            hresy_c->Fill(trueHit->YErr());
            hresz_c->Fill(trueHit->ZErr());
        }
        if (trueHit->SubSys()==1) {
            hresx_m->Fill(trueHit->XErr());
            hresy_m->Fill(trueHit->YErr());
            hresz_m->Fill(trueHit->ZErr());
        }
        if (trueHit->SubSys()==2) {
            hresx_d->Fill(trueHit->XErr());
            hresy_d->Fill(trueHit->YErr());
            hresz_d->Fill(trueHit->ZErr());
        }

        hzx->Fill(trueHit->Z(),trueHit->X());
        hzy->Fill(trueHit->Z(),trueHit->Y());
        hxy->Fill(trueHit->X(),trueHit->Y());
    }

    hresx_c->GetXaxis()->SetTitle("RMS (cm)");
    hresx_m->GetXaxis()->SetTitle("RMS (cm)");
    hresx_d->GetXaxis()->SetTitle("RMS (cm)");

    hresx_c->SetLineColor(kBlue);
    hresy_c->SetLineColor(kRed);
    hresz_c->SetLineColor(kMagenta);
    hresx_m->SetLineColor(kBlue);
    hresy_m->SetLineColor(kRed);
    hresz_m->SetLineColor(kMagenta);
    hresx_d->SetLineColor(kBlue);
    hresy_d->SetLineColor(kRed);
    hresz_d->SetLineColor(kMagenta);

    hxy->GetXaxis()->SetTitle("X (cm)");
    hxy->GetYaxis()->SetTitle("Y (cm)");

    hzy->GetXaxis()->SetTitle("Z (cm)");
    hzy->GetYaxis()->SetTitle("Y (cm)");

    hzx->GetXaxis()->SetTitle("Z (cm)");
    hzx->GetYaxis()->SetTitle("X (cm)");

    hxy->SetMarkerStyle(8);
    hzx->SetMarkerStyle(8);
    hzy->SetMarkerStyle(8);

    hxy->SetMarkerColor(kBlue);
    hzy->SetMarkerColor(kBlue);
    hzx->SetMarkerColor(kBlue);

    new TCanvas();
    hxy->Draw();

    new TCanvas();
    hzy->Draw();

    new TCanvas();
    hzx->Draw();

    new TCanvas();
    hresx_c->Draw();
    hresy_c->Draw("same");
    hresz_c->Draw("same");

    new TCanvas();
    hresx_m->Draw();
    hresy_m->Draw("same");
    hresz_m->Draw("same");

    new TCanvas();
    hresx_d->Draw();
    hresy_d->Draw("same");
    hresz_d->Draw("same");

}

void PlotTrueCRTHitByRegion(Int_t reg=0) {

    TrueHit* trueHit = filereader.tmTrueHit();
    Int_t entriesTrue = trueHit->GetInputTree()->GetEntries();

    TH1F *hresx = new TH1F("hresx","True CRT Hit Position RMS",5000, 0,1000);
    TH1F *hresy = new TH1F("hresy","True CRT Hit Position RMS",5000, 0,1000);
    TH1F *hresz = new TH1F("hresz","True CRT Hit Position RMS",5000, 0,1000);

    TH2F *hzx = new TH2F("hzx","True CRT Hits",1000,-1200,1200,600,-600,600);
    TH2F *hzy = new TH2F("hzy","True CRT Hits",1000,-1020,1200,550,-400,700);
    TH2F *hxy = new TH2F("hxy","True CRT Hits",600,-650,650,550,-400,700);

    for (Int_t ientry=0; ientry<entriesTrue; ientry++) {

        trueHit->GetHit(ientry);

        if (trueHit->Region()!=reg) continue;

        hresx->Fill(trueHit->XErr());
        hresy->Fill(trueHit->YErr());
        hresz->Fill(trueHit->ZErr());

        hzx->Fill(trueHit->Z(),trueHit->X());
        hzy->Fill(trueHit->Z(),trueHit->Y());
        hxy->Fill(trueHit->X(),trueHit->Y());

    }

    hresx->GetXaxis()->SetTitle("RMS (cm)");

    hresx->SetLineColor(kBlue);
    hresy->SetLineColor(kRed);
    hresz->SetLineColor(kMagenta);

    hxy->GetXaxis()->SetTitle("X (cm)");
    hxy->GetYaxis()->SetTitle("Y (cm)");

    hzy->GetXaxis()->SetTitle("Z (cm)");
    hzy->GetYaxis()->SetTitle("Y (cm)");

    hzx->GetXaxis()->SetTitle("Z (cm)");
    hzx->GetYaxis()->SetTitle("X (cm)");

    hxy->SetMarkerStyle(8);
    hzx->SetMarkerStyle(8);
    hzy->SetMarkerStyle(8);

    hxy->SetMarkerColor(kBlue);
    hzy->SetMarkerColor(kBlue);
    hzx->SetMarkerColor(kBlue);

    new TCanvas();
    hxy->Draw();

    new TCanvas();
    hzy->Draw();

    new TCanvas();
    hzx->Draw();

    new TCanvas();
    hresx->Draw();
    hresy->Draw("same");
    hresz->Draw("same");

}

/////////////////////////////////////////////////////////
// CosmicDisplay tree functions, including event display
/////////////////////////////////////////////////////////

//dump all info from CosmicDisplay tree
void DumpCosmicTracks()
{
        CosmicDisplay *track = filereader.tmCD();
        Long64_t entriesTrk = track->GetInputTree()->GetEntries();
        cout << "Reading all tracks " << entriesTrk << endl;
        for( Long64_t ientry=0; ientry<entriesTrk; ++ientry )
        {
                track->GetCosmic(ientry);
                //if( ientry>20 ) break;
                if ( DEBUG_cd ) cout << "tracks:: " << track->Event() << " " << track->NSeg() << endl;
                for( Long64_t eve=0; eve<track->NSeg(); ++eve )
                {
                        if ( DEBUG_cd )
                        cout << " segments["<< eve <<"] "  << track->PDG()
                                 << " " << track->SlopeX(eve) << " " << track->SlopeY(eve) << " " << track->SlopeZ(eve)
                                 << " " << track->X(eve) << " " << track->Y(eve) << " " << track->Z(eve)
                                 << " " << track->E(eve) << " " << track->P(eve)
                                 << " " << track->T(eve)  << endl;
                } // end of for eve loop
        } // end of for ientry loop
}//end def

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Simple Event Display 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
void DisplayEvent(Long64_t event = -1, Long64_t trackID = -1, Bool_t fromFile=kFALSE, Bool_t register=kTRUE)
{
    const double posDetEncl = 154.995;

    gStyle->SetPalette(-1);
    TEveManager* evem;
    if (gEve)
    {
        evem = gEve;
        gEve->GetViewers()->DeleteAnnotations();
        gEve->GetCurrentEvent()->DestroyElements();
    }
    else
    {
        evem = TEveManager::Create();
    }


    TEveRGBAPalette* pal = new TEveRGBAPalette(0, 50);
    TEveRGBAPalette* pal2 = new TEveRGBAPalette(0, 50);
    TEveViewer *ev;
    TGLViewer *gv;

    ev = gEve->GetDefaultViewer();
    gv = ev->GetGLViewer();
    gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);

    // Skip MultiView Options
    //frames

    TEveWindowSlot *slot = 0;
    TEveWindowPack *pack = 0;
    TEveViewer     *T3DView = 0;
    TEveViewer     *TRPhiView = 0;
    TEveViewer     *TRhoZView = 0;
    TEveViewer     *MuonView = 0;

    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    pack = slot->MakePack(); // slot destroyed
    pack->SetElementName("Multi View");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kFALSE);

    pack->NewSlot()->MakeCurrent(); // new slot created from pack
    T3DView = gEve->SpawnNewViewer("Y-Z View", "");
    T3DView->AddScene(gEve->GetGlobalScene());
    T3DView->AddScene(gEve->GetEventScene());
    gv = T3DView->GetGLViewer();
    gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
    pack = pack->NewSlot()->MakePack();
    pack->SetShowTitleBar(kFALSE);
    pack->NewSlot()->MakeCurrent();

    TRPhiView = gEve->SpawnNewViewer("X-Y View", "");
    TRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    TRPhiView->AddScene(gEve->GetGlobalScene());
    TRPhiView->AddScene(gEve->GetEventScene());
    gv = TRPhiView->GetGLViewer();
    gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
    pack->NewSlot()->MakeCurrent();

    TRhoZView = gEve->SpawnNewViewer("X-Z View", "");
    TRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
    TRhoZView->AddScene(gEve->GetGlobalScene());
    TRhoZView->AddScene(gEve->GetEventScene());
    gv = TRhoZView->GetGLViewer();
    gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);

    //TString file("./icarus_complete_test_crt_03aug_nowires.gdml");
    TString file("./icarus_complete_nowires.gdml");
    TGeoManager::Import(file);
    TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("volDetEnclosure_0");
    TEveGeoTopNode* evnode = new TEveGeoTopNode(gGeoManager,node);
    evnode->SetMainTransparency(95);
    gEve->AddElement(evnode);
    

    // Building .....
    double XmaxF = 500;
    double XminF = -500;
    double YmaxF = 1000;
    double YminF = -1000;
    double ZmaxF = 10000;
    double ZminF = -10000;
    double dx = XmaxF-XminF;
    double dy = YmaxF-YminF;
    double dz = ZmaxF-ZminF;
    TEveFrameBox* frm = new TEveFrameBox();
    frm->SetAABox(XminF, YminF, ZminF, dx, dy, dz);
    frm->SetFrameColor(kBlack);
    frm->SetBackColorRGBA(120,120,120,20);
    frm->SetDrawBack(kTRUE);

    cout << "I am going to extract track segments" << endl;
    map<TString,TEveBoxSet*> boxSets;
    // extract segments information from ResultsTTree
    CosmicDisplay *trk = filereader.tmCD();
    Int_t nentriesCD = trk->GetInputTree()->GetEntries();
    int charge;
    bool cdfind = false;

    for(Int_t ientry=0; ientry<nentriesCD; ientry++) {

        trk->GetCosmic(ientry);
        //cout << "retreived data for entry " << ientry << endl;
        if(trk->Event() != event && event!= -1) continue;
        if(trk->TrackID() != trackID && event!=-1 && trackID!=-1) continue;
        cdfind = true;

        TString trkname = "Track_";
        trkname+=event; trkname+="_"; trkname+=trackID;
        boxSets[trkname] = new TEveBoxSet(trkname);
        TEveBoxSet* seg = boxSets[trkname];
        seg->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
        seg->SetPalette(pal);
        seg->UseSingleColor();
        seg->SetMainColor(kRed);
        seg->SetMainTransparency(50);

        // loop over track segments
        for(Int_t i=0; i<trk->NSeg(); i++)
        {
            //if( trk->Region(i) == 0 ) continue;
            double size = 10;
            double xminT = trk->X(i) - size/2;
            double yminT = trk->Y(i) - size/2 - posDetEncl;
            double zminT = trk->Z(i) - size/2;
            charge = 1;//GetCharged(trk->PDG());
            if( charge ){ // drawing charged tracks only
                seg->AddBox(xminT, yminT, zminT, size, size, size);
                seg->SetMainColor(kRed);
                seg->DigitValue(15);
            }
        } // segments

        seg->RefitPlex();
        seg->SetPickable(kTRUE);
        seg->RefMainTrans().SetPos(0,0,0);
        gEve->AddElement(seg);

    } // trackIDs

    if(!cdfind) cout << "warning: event and/or trackID not found in CD tree" << endl;

    // Sim hits will be drawn 
    Hit* simHit = filereader.tmHit();
    Int_t nentriesHit = simHit->GetInputTree()->GetEntries();

    // loop over individual hits
    for (Int_t ientry=0; ientry<nentriesHit; ientry++) {

        simHit->GetHit(ientry);

        if(simHit->Event() != event && event!= -1) continue;
        if(simHit->TrackID(0) != trackID && event!=-1 && trackID!=-1) continue;
        //if(simHit->SubSys() != 0 ) continue;

        TString hitname = "SimHit_";
        hitname+=event; hitname+="_"; hitname+=trackID;
        boxSets[hitname] = new TEveBoxSet(hitname);
        TEveBoxSet* seg = boxSets[hitname];

        seg->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
        seg->SetPalette(pal);
        seg->UseSingleColor();
        seg->SetMainColor(kBlue);
        //seg->SetMainTransparency(50);

        double xminHit = simHit->X()-10;//simHit->XErr();
        double yminHit = simHit->Y()-10 - posDetEncl;//simHit->YErr();
        double zminHit = simHit->Z()-10;//simHit->ZErr();
        double sizexHit = 20;//2*simHit->XErr();
        double sizeyHit = 20;//2*simHit->YErr();
        double sizezHit = 20;//2*simHit->ZErr();

        seg->AddBox(xminHit, yminHit, zminHit, sizexHit, sizeyHit, sizezHit);
        //seg->SetMainColor(kBlue);
        seg->DigitValue(15);

        seg->RefitPlex();
        seg->SetPickable(kTRUE);
        seg->RefMainTrans().SetPos(0,0,0);
        gEve->AddElement(seg);

    }// for simHits

    //True hits will be drawn 
    TrueHit* trueHit = filereader.tmTrueHit();
    Int_t nentriesTrueHit = trueHit->GetInputTree()->GetEntries();

    // loop over individual hits
    for (Int_t ientry=0; ientry<nentriesTrueHit; ientry++) {

        trueHit->GetHit(ientry);

        if(trueHit->Event() != event && event!= -1) continue;
        if(trueHit->TrackID(0) != trackID && event!=-1) continue;
        //if(trueHit->SubSys() != 0 ) continue;

        TString hitname = "TrueHit_";
        hitname+=event; hitname+="_"; hitname+=trackID;
        boxSets[hitname] = new TEveBoxSet(hitname);
        TEveBoxSet* seg = boxSets[hitname];

        seg->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
        seg->SetPalette(pal);
        seg->UseSingleColor();
        seg->SetMainColor(kMagenta);
        //seg->SetMainTransparency(50);

        double xminHit = trueHit->X()-10;//simHit->XErr();
        double yminHit = trueHit->Y()-10 - posDetEncl;//simHit->YErr();
        double zminHit = trueHit->Z()-10;//simHit->ZErr();
        double sizexHit = 20;//2*simHit->XErr();
        double sizeyHit = 20;//2*simHit->YErr();
        double sizezHit = 20;//2*simHit->ZErr();

        seg->AddBox(xminHit, yminHit, zminHit, sizexHit, sizeyHit, sizezHit);
        //seg->SetMainColor(kBlue);
        seg->DigitValue(15);

        seg->RefitPlex();
        seg->SetPickable(kTRUE);
        seg->RefMainTrans().SetPos(0,0,0);
        gEve->AddElement(seg);

    }// for trueHits


    //gEve->GetBrowser()->GetTabRight()->SetTab(0);
    gEve->Redraw3D(kTRUE);
 
}//end def


//////////////////////////////////////////////////////////
/// functions mixing between trees
//////////////////////////////////////////////////////////

//find time of flight between CRT and earliest PMT hit
void PlotCRTToPMTTime() {

        double refIndex = 1.38; //index of refraction in LAr

        Regions *reg = filereader.tmReg();
        Long64_t entriesReg = reg->GetInputTree()->GetEntries();

        int nmiss = 0;
        TH1D* h = new TH1D("h","time delay bettween CRT and PMT",50,0,50);
        h->GetXaxis()->SetTitle("time difference (ns)");
        h->SetLineWidth(2);

        double tCRT = 0., tPMT = 0.;
        bool tagCRT = false, tagPMT = false;

        for( Long64_t ientry=0; ientry<entriesReg; ++ientry )
        {
                reg->GetRegion(ientry);
                if (reg->Active()==0) continue;

                tagCRT = false; tagPMT = false;
                for (Int_t hit=0; hit<reg->NReg(); hit++) {

                    if (reg->Region(hit)>30) 
                    {
                        tCRT = reg->EntryT(hit);
                        tagCRT = true;
                    }
                    if (reg->Region(hit)>0&&reg->Region(hit)<10)
                    {
                        tPMT = reg->OpDetT(hit);
                        tagPMT = true;
                    }
                    if (tagCRT && tagPMT) break;
                }
                if (!tagCRT||!tagPMT) nmiss++;
                if (tPMT-tCRT == 0 ) cout << "0 delay found!" << endl;
                cout << "uncorr tPMT: " << tPMT << '\n'
                     << "  corr tPMT: " << tPMT*refIndex << endl;
                cout << "dT: " << tPMT-tCRT << endl;
                h->Fill(tPMT-tCRT);
        }

        cout << "found " << h->Integral() << " matches and missed " << nmiss << endl;

        h->Draw();
}

//trying to match true CRT hits to tracks entering/leaving TPC regions
void TrueCRTTPCMatch() {

    TrueHit *hit = filereader.tmTrueHit();
    Int_t nentriesHit = hit->GetInputTree()->GetEntries();

    Regions* reg = filereader.tmReg();
    Int_t nentriesReg = reg->GetInputTree()->GetEntries();

    Int_t nMissCRT=0, nAV=0;
    TH1D *hdist = new TH1D("hdist","CRTHit to AV Point of Entry",10000/5,0,10000);
    TH1D *hproj = new TH1D("hproj","CRT Hit distance from AV back projection",250,0,1000);

    for ( Int_t entryReg=0; entryReg<nentriesReg; entryReg++ ) {

        reg->GetRegion(entryReg);

        if ( reg->Active() == 0 ) continue;
        if ( reg->CRTs() == 0 ) {
            nMissCRT++;
            continue;
        }
        nAV++;

        for ( Int_t i=0; i<reg->NReg(); i++ ) {

            if ( reg->Region(i) < 5 || reg->Region(i) > 8 ) continue;
            if ( reg->DL(i) < 20 ) continue;

            std::vector<double> hitX, hitY, hitZ;

            for ( Int_t entryHit=0; entryHit<nentriesHit; entryHit++ ) {

                hit->GetHit(entryHit);

                if (reg->Event() != hit->Event()) continue;
                if (reg->TrackID() != hit->TrackID(0)) continue;

                hitX.push_back(hit->X());
                hitY.push_back(hit->Y());
                hitZ.push_back(hit->Z());         

            } // loop over CRTHits

            if (hitX.size()==0) continue;

            double distance = -1.0, distance_tmp=-1.0;
            int minIndex = 0;
            //double 
            for ( size_t j=0; j<hitX.size(); j++) {
                distance_tmp = sqrt(pow(reg->EntryX(i)-hitX[j],2)
                      +pow(reg->EntryY(i)-hitY[j],2)+pow(reg->EntryZ(i)-hitZ[j],2));
                if (distance_tmp<distance||j==0) {
                    distance = distance_tmp; 
                    minIndex = j;
                }
            }

            hdist->Fill(distance);

            double dx, dy, dz;
            dx = reg->EntryX(i) - hitX[minIndex];
            dy = reg->EntryY(i) - hitY[minIndex];
            dz = reg->EntryX(i) - hitX[minIndex];
            double a, bx, by, bz, distMin;
            a = reg->EntrySlopeX(i) * dx;
            a+= reg->EntrySlopeY(i) * dy;
            a+= reg->EntrySlopeZ(i) * dz;

            bx = reg->EntrySlopeX(i) * a - dx;
            by = reg->EntrySlopeY(i) * a - dy;
            bz = reg->EntrySlopeZ(i) * a - dz;

            distMin = sqrt( bx*bx + by*by + bz*bz );
            if (distMin>99999)
            cout << "cosx: " << reg->EntrySlopeX(i) << ", cosy: " << reg->EntrySlopeY(i) << ", cosz: " << reg->EntrySlopeZ(i) << '\n'
                 << "bx: " << bx << ", by: " << by << ", bz: " << bz << '\n'
                 << "distMin: " << distMin << endl;
            hproj->Fill(distMin);
        } // loop over region entries
    } // loop over region events

    hdist->Draw();
    new TCanvas();
    hproj->Draw();

}


///////////////////// end file /////////////////////////
