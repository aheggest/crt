//local includes
#include "TreeManager.h"

//ROOT includes
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH3F.h>
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

//file to be analyzed//


//1. merged55 file, a merged file of 55 muonsonly files.
//TreeManager filereader("~/Downloads/research/analysis/crtanalysis/ROOTfiles/merged55.root");

//2. mergedNOoverburden55 file, a merged file of 55 muonsonly files with NO OVERBURDEN.
TreeManager filereader("~/Downloads/research/analysis/crtanalysis/ROOTfiles/mergedNOoverburden55.root");

//3. 1 no overburden file, for quicker compiling
//TreeManager filereader("~/Downloads/research/analysis/crtanalysis/ROOTfiles/prodcosmics_icarus_nooverburden_cosmics_nooverburden_simana.root");

//TreeManager filereader("~/Downloads/research/analysis/crtanalysis/ROOTfiles/simana.root");

//TreeManager filereader("/icarus/data/users/chilgenb/prod_v08_55_00/simana/prodcosmics_icarus_nooverburden/filesana.list");

/////////////////////////////////////
// Generator Tree Specific Functions
/////////////////////////////////////


//Draw PDG functions, using genTree instead of the simTree as done in Analyze.C
void PlotGenPDGs() {

  Gen *gen = filereader.tmGen(); //pointer to Generator class
  Int_t nentries = gen->GetInputTree()->GetEntries(); //number of entries in tree
  cout << "reading from GenTree with " << nentries << " entries..." << endl; //read number of entries
  TH1F* h = new TH1F("h","Generator PDG Codes",10000,-5000,5000); //intitialize histogram for PDG codes

  //loop over entries
  for ( Int_t ientry = 0; ientry < nentries; ientry++ ){
        gen->GetGen(ientry); //get info for ith entry
        h->Fill(gen->PDG(ientry)); //fill the histo with PDG code
  }
  h->Draw(); //draw the histogram
}

void PlotFlux() {
    //TFile fout("/Downloads/research/analysis/crtanalysis/ROOTfiles/outfile.root","RECREATE");
    Gen *gen = filereader.tmGen(); //pointer to Generator class
    Int_t nentries = gen->GetInputTree()->GetEntries(); //number of entries in generator tree
    cout << "reading from GenTree with " << nentries << " entries..." << endl; //read number of entries


    const size_t size = 101;
    const double mult = 1.18; //constant multiplier, smaller M = finer bins
    double bins[size];      //intitialize array "bins"
    bins[0] = 0.01;
    double nfiles = 55; //how many input files
    Int_t nevents = nentries / nfiles; //how many events per file
    cout << "reading from " << nevents << " events from " << nfiles << " files " << endl;


    for (size_t i = 1; i < size; i++) {
      bins[i] = bins[i-1]*mult;
      }

    //one histogrem per particle
    TH1F* hmun = new TH1F("hmun","Particle fluxes",size-1,bins); //muon(-) histo
    TH1F* hmup = new TH1F("hmup","Particle fluxes",size-1,bins); //muon(+) histo
    TH1F* hph = new TH1F("hph","Particle fluxes",size-1,bins); //photon histo
    TH1F* hn = new TH1F("hn","Particle fluxes",size-1,bins); //neutron histo
    TH1F* hpr = new TH1F("hpr","Particle fluxes",size-1,bins); //proton histo
    TH1F* he = new TH1F("he","Particle fluxes",size-1,bins); //electron histo
    TH1F* hpo = new TH1F("hpo","Particle fluxes",size-1,bins); //positron histo
    TH1F* hpip = new TH1F("hpip","Particle fluxes",size-1,bins); //pion(+) histo
    TH1F* hpin = new TH1F("hpin","Particle fluxes",size-1,bins); //pion(-) histo
    TH1F* hk = new TH1F("hk","Particle fluxes",size-1,bins); //kaon0 histo



    //loop over entries
   for ( Int_t ientry = 0; ientry < nentries; ientry++ ){

      gen->GetGen(ientry); //get info for ith entry

      for (Int_t i = 0; i < gen->NGen(); i++) {


        Int_t weight = 1.0*gen->StartE(i);
        if (gen->PDG(i) == 13) hmun->Fill(gen->StartE(i),weight); //Fill histo with muon(-) energies
        if (gen->PDG(i) == -13) hmup->Fill(gen->StartE(i),weight); //Fill histo with muon(+) energies
        if (gen->PDG(i) == 22) hph->Fill(gen->StartE(i),weight); //Fill histo with photon energies
        if (gen->PDG(i) == 2112) hn->Fill(gen->StartE(i),weight); //Fill histo with Neutron energies
        if (gen->PDG(i) == 2212) hpr->Fill(gen->StartE(i),weight); //Fill histo with proton energies
        if (gen->PDG(i) == 11) he->Fill(gen->StartE(i),weight); //Fill histo with electron energies
        if (gen->PDG(i) == -11) hpo->Fill(gen->StartE(i),weight); //Fill histo with positron energies
        if (gen->PDG(i) == 211) hpip->Fill(gen->StartE(i),weight); //Fill histo with pion(+) energies
        if (gen->PDG(i) == -211) hpin->Fill(gen->StartE(i),weight); //Fill histo with pion(-) energies
        if (gen->PDG(i) == 311) hk->Fill(gen->StartE(i),weight); //Fill histo with kaon0 energies

      } // end for ngen
    } //end for entries


  double simTime = nevents*3.3e-3; //.033[s]
  double simArea = 1.2e+6; //[120 m^2 = 1.2*10^6 cm^2]
  cout << "simTime " << simTime << endl;
  hmun->Scale(1/(simTime*simArea));
  hmup->Scale(1/(simTime*simArea));
  hph->Scale(1/(simTime*simArea));
  hn->Scale(1/(simTime*simArea));
  hpr->Scale(1/(simTime*simArea));
  he->Scale(1/(simTime*simArea));
  hpo->Scale(1/(simTime*simArea));
  hpip->Scale(1/(simTime*simArea));
  hpin->Scale(1/(simTime*simArea));
  hk->Scale(1/(simTime*simArea));

  hmun->GetXaxis()->SetTitle("Energy(GeV)");
  hmun->GetYaxis()->SetTitle("[Part/m^2/s]");

  hmun->SetLineColor(kRed);
  hmup->SetLineColor(kPink);
  hph->SetLineColor(kGreen);
  hn->SetLineColor(kOrange);
  hpr->SetLineColor(kBlack);
  he->SetLineColor(9);
  hpo->SetLineColor(kCyan);
  hpip->SetLineColor(6);
  hpin->SetLineColor(30);
  hk->SetLineColor(kGray);

  auto legend = new TLegend(0.75,0.7,0.95,0.95);
  legend->AddEntry(hn,"neutrons","lc");
  legend->AddEntry(hpr,"protons","lc");
  legend->AddEntry(hpip,"pion(+)","lc");
  legend->AddEntry(hpin,"pion(-)","lc");
  legend->AddEntry(hk,"kaon0","lc");
  legend->AddEntry(hmup,"muon(+)","lc");
  legend->AddEntry(hmun,"muon(-)","lc");
  legend->AddEntry(hph,"photons","lc");
  legend->AddEntry(he,"electron","lc");
  legend->AddEntry(hpo,"positron","lc");

  hmun->GetYaxis()->SetRangeUser(1e-12,10);

  //set log scale
  auto c = new TCanvas("c1");
  c->SetLogx();
  c->SetLogy();

  hmun->Draw("hist");
  hmup->Draw("samehist");
  hph->Draw("samehist");
  hn->Draw("samehist");
  hpr->Draw("samehist");
  he->Draw("samehist");
  hpo->Draw("samehist");
  hpip->Draw("samehist");
  hpin->Draw("samehist");
  hk->Draw("samehist");
  legend->Draw();
  //c->Write();
  //c->SaveAs("/Downloads/research/analysis/crtanalysis/ROOTfiles/c1.png");
  //fout.Close();
}

//Dump generator info
void DumpGeneratorInfo() {

    Gen *gen = filereader.tmGen(); //pointer to Generator class
    Int_t nentries = gen->GetInputTree()->GetEntries(); //number of entries in tree
    cout << "Total number of entries in Generator tree: " << nentries << endl;

    //for loop to get&dump info for each entry
    for( Int_t ientry = 0; ientry < nentries; ientry++)
    {
      gen->GetGen(ientry); //get info for ith entry

      //print things
      cout << "TrackID for ith entry: "<< gen->TrackID(ientry) << endl;
      cout << "PDG for ith entry: " << gen->PDG(ientry) << endl;
      cout << "-----------------------------" << endl;

    }
}


//////////////////////////////////
// Regions Tree Specific Functions
///////////////////////////////////

//count regions hit with reg->Active(), Reg->CRTs(),..
void CountRegionsHit(){

        Regions* reg = filereader.tmReg();
        Int_t nentries = reg->GetInputTree()->GetEntries();
        cout << "reading from RegionsTree with " << nentries << " entries..." << endl;

        //setup counters

        Int_t nCRTonly = 0, nAVonly = 0, nIVonly = 0; //counters for how many pass through only CRT, only AV, only IV
        Int_t nCRTAV = 0, nCRTIV = 0, nAVIV = 0, nCRTtagAV = 0;
        Int_t nAVmissCRT = 0, nIVmissCRT =0; //miss CRT and hit av only.
        Int_t nAV = 0;
        Int_t nALL = 0; // counter for those that hit all
        Int_t ntracks = 0;
        bool isCRTonly = false, isAVonly = false, isIVonly = false;
        bool isCRTAV = false, isCRTIV = false, isAVIV = false;
        bool isAVmissCRT = false, isIVmissCRT=false;
        bool isALL = false, isCRTtagAV = false;
        bool isAV = false;
        set<int> events;
        map<int,vector<int>> muHits;

        reg->GetRegion(0);
        int currentEvent = reg->Event();
        events.insert(currentEvent);

        for (Int_t ientry = 0; ientry < nentries; ientry++) {

          reg->GetRegion(ientry);

          if (abs(reg->PDG()) != 13 ) continue; //if not a muon, skip
          ntracks++;
          //pass through CRT only
          isCRTonly = false;
          if (reg->CRTs()>0 && reg->Active()==0  && reg->Inactive()==0) isCRTonly = true;
          if (isCRTonly) nCRTonly++;
          //pass through AV only
          isAVonly = false;
          if (reg->CRTs()==0 && reg->Active()>0  && reg->Inactive()==0) isAVonly = true;
          if (isAVonly) nAVonly++;
          //pass through IV only
          isIVonly = false;
          if (reg->CRTs()==0 && reg->Active()==0  && reg->Inactive()>0) isIVonly = true;
          if (isIVonly) nIVonly++;
          //pass through CRT and enter AV only
          isCRTAV = false;
          if (reg->CRTs()>0 && reg->Active()>0  && reg->Inactive()==0) isCRTAV = true;
          if (isCRTAV) nCRTAV++;
          //pass through CRT and enter IV only
          isCRTIV = false;
          if (reg->CRTs()>0 && reg->Active()==0  && reg->Inactive()>0) isCRTIV = true;
          if (isCRTIV) nCRTIV++;
          //pass through all CRT, IV and AV
          isALL = false;
          if (reg->CRTs()>0 && reg->Active()>0  && reg->Inactive()>0) isALL = true;
          if (isALL) nALL++;
          //miss CRT but enter AV and IV
          isAVIV = false;
          if (reg->CRTs()==0 && reg->Active()>0  && reg->Inactive()>0) isAVIV = true;
          if (isAVIV) nAVIV++;
          //tagged by CRT and AV
          isCRTtagAV = false;
          if (reg->CRTs()>0 && reg->Active()>0) isCRTtagAV = true;
          if (isCRTtagAV) nCRTtagAV++;
          //hit AV but miss CRT
          isAVmissCRT = false;
          if (reg->CRTs()==0 && reg->Active()>0) isAVmissCRT = true;
          if (isAVmissCRT) nAVmissCRT++;
          //hit IV but miss CRT
          isIVmissCRT = false;
          if (reg->CRTs()==0 && reg->Inactive()>0) isIVmissCRT = true;
          if (isIVmissCRT) nIVmissCRT++;
          //total that hit AV
          isAV = false;
          if (reg->Active()>0) isAV = true;
          if (isAV) nAV++;

        }

        cout << ntracks << " tracks " <<  endl;
        cout << " of which " << nCRTonly << " pass through CRT only" << endl;
        cout << " of which " << nAVonly << " pass through AV only" << endl;
        cout << " of which " << nIVonly << " pass through IV only" << endl;
        cout << " of which " << nCRTAV << " pass through the CRT and only enter into AV and IV=0" << endl;
        cout << " of which " << nCRTIV << " pass through the CRT and only enter into IV" << endl;
        cout << " of which " << nAVIV << " miss CRT but enter both AV and IV" << endl;
        cout << " of which " << nALL << " pass through all CRT, AV, IV" << endl;
        cout << " of which " << nCRTtagAV << " pass through AV and CRT" << endl;
        cout << " of which " << nAVmissCRT << " miss CRT but enter AV" << endl;
        cout << " of which " << nIVmissCRT << " miss CRT but enter IV" << endl;
        cout << " of which " << nAV << " pass through the AV." << endl;

        //tagging efficiencies
        cout << "Efficiencies of... " << endl;
        cout << "    CRT tagging muons entering AV: " << 1.0*nCRTtagAV/nAV << endl;

}

// Rate of cosmic muons passing through each CRT region
void DumpRates(){
  Regions* reg = filereader.tmReg();
  Int_t nentries = reg->GetInputTree()->GetEntries();
  cout << "reading from RegionsTree with " << nentries << " entries..." << endl;

  //setup counters

  Int_t nCRTonly = 0, nAVonly = 0, nIVonly = 0; //counters for how many pass through only CRT, only AV, only IV
  Int_t nCRTAV = 0, nCRTIV = 0, nAVIV = 0, nCRTtagAV = 0;
  Int_t nAVmissCRT = 0, nIVmissCRT =0; //miss CRT and hit av only.
  Int_t nAV = 0;
  Int_t nALL = 0; // counter for those that hit all
  Int_t ntracks = 0;
  bool isCRTonly = false, isAVonly = false, isIVonly = false;
  bool isCRTAV = false, isCRTIV = false, isAVIV = false;
  bool isAVmissCRT = false, isIVmissCRT=false;
  bool isALL = false, isCRTtagAV = false;
  bool isAV = false;
  set<int> events;
  map<int,vector<int>> muHits;

  reg->GetRegion(0);
  int currentEvent = reg->Event();
  events.insert(currentEvent);

  for (Int_t ientry = 0; ientry < nentries; ientry++) {

    reg->GetRegion(ientry);

    if (abs(reg->PDG()) != 13 ) continue; //if not a muon, skip
    if ( reg->Event()==currentEvent || ientry==nentries-1) {
      ntracks++;

      //pass through CRT only
      isCRTonly = false;
      if (reg->CRTs()>0 && reg->Active()==0  && reg->Inactive()==0) isCRTonly = true;
      if (isCRTonly) nCRTonly++;

      //pass through AV only
      isAVonly = false;
      if (reg->CRTs()==0 && reg->Active()>0  && reg->Inactive()==0) isAVonly = true;
      if (isAVonly) nAVonly++;

      //pass through IV only
      isIVonly = false;
      if (reg->CRTs()==0 && reg->Active()==0  && reg->Inactive()>0) isIVonly = true;
      if (isIVonly) nIVonly++;

      //pass through CRT and enter AV only
      isCRTAV = false;
      if (reg->CRTs()>0 && reg->Active()>0  && reg->Inactive()==0) isCRTAV = true;
      if (isCRTAV) nCRTAV++;

      //pass through CRT and enter IV only
      isCRTIV = false;
      if (reg->CRTs()>0 && reg->Active()==0  && reg->Inactive()>0) isCRTIV = true;
      if (isCRTIV) nCRTIV++;

      //pass through all CRT, IV and AV
      isALL = false;
      if (reg->CRTs()>0 && reg->Active()>0  && reg->Inactive()>0) isALL = true;
      if (isALL) nALL++;

      //miss CRT but enter AV and IV
      isAVIV = false;
      if (reg->CRTs()==0 && reg->Active()>0  && reg->Inactive()>0) isAVIV = true;
      if (isAVIV) nAVIV++;

      isCRTtagAV = false;
      if (reg->CRTs()>0 && reg->Active()>0) isCRTtagAV = true;
      if (isCRTtagAV) nCRTtagAV++;

      isAVmissCRT = false;
      if (reg->CRTs()==0 && reg->Active()>0) isAVmissCRT = true;
      if (isAVmissCRT) nAVmissCRT++;

      isIVmissCRT = false;
      if (reg->CRTs()==0 && reg->Inactive()>0) isIVmissCRT = true;
      if (isIVmissCRT) nIVmissCRT++;

      isAV = false;
      if (reg->Active()>0) isAV = true;
      if (isAV) nAV++;
    }

    muHits.clear();
    currentEvent=reg->Event();
    events.insert(reg->Event());
  }

  double sampleTime = events.size()*3.3e-3;

  cout << " entry counts in " << events.size() << " events: " << '\n'
      << "   ntracks:  " << 1.0*ntracks << '\n'
      << "   nCRTonly:  " << 1.0*nCRTonly << '\n'
      << "   nAVonly:  " << 1.0*nAVonly << '\n'
      << "   nIVonly:  " << 1.0*nIVonly << '\n'
      << "   nCRTAV:  " << 1.0*nCRTAV << '\n'
      << "   nCRTIV:  " << 1.0*nCRTAV << '\n'
      << "   nALL:  " << 1.0*nALL << '\n'
      << "   nAVIV:  " << 1.0*nAVIV << '\n'
      << "   nCRTtagAV: " << 1.0*nCRTtagAV << '\n'
      << "   nAVmissCRT: " << 1.0*nAVmissCRT << '\n'
      << "   nIVmissCRT: " << 1.0*nIVmissCRT << '\n'
  << endl;

  cout << " entry rates in " << events.size() << " events: " << '\n'
       << "   ntracks:  " << 1.0*ntracks/sampleTime << '\n'
       << "   nCRTonly:  " << 1.0*nCRTonly/sampleTime << '\n'
       << "   nAVonly:  " << 1.0*nAVonly/sampleTime << '\n'
       << "   nAVonly:  " << 1.0*nIVonly/sampleTime << '\n'
       << "   nCRTAV:  " << 1.0*nCRTAV/sampleTime << '\n'
       << "   nCRTIV:  " << 1.0*nCRTIV/sampleTime << '\n'
       << "   nALL:   " << 1.0*nALL/sampleTime << '\n'
       << "   nAVIV:  " << 1.0*nAVIV/sampleTime << '\n'
       << "   nCRTtagAV: " << 1.0*nCRTtagAV/sampleTime << '\n'
       << "   nAVmissCRT: " << 1.0*nAVmissCRT/sampleTime << '\n'
       << "   nIVmissCRT: " << 1.0*nIVmissCRT/sampleTime << '\n'
  << endl;


}

/*void PlotRates(){

  Regions* reg = filereader.tmReg();
  Int_t nentries = reg->GetInputTree()->GetEntries();
  cout << "reading from RegionsTree with " << nentries << " entries..." << endl;

  Int_t nmu = 0;
  Int_t nrate = 0;

  Int_t nmu_tmp = 0;
  Int_t nrate_tmp = 0;

  set<int> events;
  map<int,vector<int>> muHits;

//  TH1F* hmu = new TH1F("hmu","#mu's per event",300,0,300);
//  TH1F* hrate = new TH1F("h_rate","#mu's rate per event",1000000,20000,1000000);

  reg->GetRegion(0
  int currentEvent = reg->Event();
  events.insert(currentEvent);

  double sampleTime = events.size()*3.3e-3;


  for ( Int_t ientry = 0; ientry < nentries; ientry++ ) {

      reg->GetRegion(ientry);
      if (reg->Active()==0) continue;

      // if not a muon, skip
      if ( abs(reg->PDG()) != 13 ) continue;

      //if same event
      if ( reg->Event()==currentEvent){
        nmu_tmp++;
        nrate_tmp++;
      }// if same event

      //if new event or last entry
      if ( reg->Event()!=currentEvent || ientry==nentries-1) {
  //      hmu->Fill(nmu_tmp);
        nmu+=nmu_tmp; nmu_tmp=0;

        nrate_tmp = 1.0*nmu/sampleTime;
  //      hrate->Fill(nrate_tmp);
        nrate+=nrate_tmp; nrate_tmp=0;

        muHits.clear();
        currentEvent=reg->Event();
        events.insert(reg->Event());
      //  cout << "new event: " << currentEvent << endl;

      } //if new event or last entry
  } //for loop over events/tracks


  new TCanvas(); hmu->Draw();
  new TCanvas(); hrate->Draw();


}*/
// rate of cosmogenic muons entering the active volume(AV), inactive volume(IV) only
///////TO be Filled!!!//

// Rate of cosmic muons passing through CRT and enter into AV, pass through CRT and only enter IV
///////TO be Filled!!!//

//  miss CRT but enter AV, miss CRT and enter IV only
///////TO be Filled!!!//

///////////////////////////////////
// DetSim Tree Specific Functions//
///////////////////////////////////

//trigger rate per FEB, for a given region
void TriggerRates(Int_t reg = 0) {

    struct eveid{
        int eve, run, subrun;

    };

    auto cmp_eveid = [](const eveid& lhs, const eveid& rhs) {
        if( lhs.run != rhs.run ) return lhs.run < rhs.run;
        if( lhs.subrun != rhs.subrun ) return lhs.subrun < rhs.subrun;
        return lhs.eve < rhs.eve;
    };

    DetSim *det = filereader.tmDet();
    Int_t entriesDet = det->GetInputTree()->GetEntries();
    cout << "reading from DetSim Tree with " << entriesDet << " entries...." << endl;
    //TH1F* h = new TH1F("h","trigger counts on FEBs",1200,0,1201);
    set<eveid,decltype(cmp_eveid)> events(cmp_eveid);
    map<int,int> febcounts;
    //map<int,TH1F*> triggerMap;

    TObjArray Hlist(0);     //create an array of histograms
    TH1F* h;                //create a pointer to a histogram
    TH1F* hsingle = new TH1F("hsingle","trigger rates on FEBs",1000,200.5,1200.5);
    int nfebcounts = 0;


    for (size_t ientry = 0; ientry < entriesDet; ientry++) {

      struct eveid eve = {det->Event(),det->Run(),det->SubRun()};
      events.insert(eve);

      det->GetDet(ientry);
      if(det->Mac5()<0) cout << "negative mac5: " << det->Mac5() << endl;
      if(det->Region()!=reg) continue;


      if (febcounts.find(det->Mac5())==febcounts.end()){
          febcounts[det->Mac5()]=0;
      }
      else {
        febcounts[det->Mac5()]++;
      }
      //febcounts[det->Mac5()]++; //fill febcounts map with Mac5 (febID!)
    }//loop over events/tracks

    double sampleTime = events.size()*3.3e-3;
    //int nhist = 1;
    int nfebID = 1;
    cout << "triggers per FEB in region " << reg << "..." << endl;
    cout << "febID -> counts -> trigger rate" << endl;

    auto legend = new TLegend(0.75,0.7,0.95,0.95);

    for (auto const& febcount : febcounts) {
      //if(febcount.second==0) continue;
      if (febcount.first < 0) continue;
      double rate = febcount.second/sampleTime;
      cout << febcount.first << " -> " << febcount.second << " -> " << rate << " (Hz)" << endl;

      //construct a histogram for every febID
      string hname = "FEB"+to_string(febcount.first);
      string htitle = "Trigger Rate on FEBs in region "+to_string(reg);
      h = new TH1F(hname.c_str(),htitle.c_str(),1000,200.5,1200.5);
      Hlist.Add(h);
      h->SetLineColor(nfebID);
      h->Fill(rate);
      hsingle->Fill(rate);
      //new TCanvas;
      h->GetXaxis()->SetTitle("rate (Hz)");
      h->Draw("histsame");
      nfebcounts+=febcount.second;

      legend->AddEntry(h,hname.c_str(),"lc");
      nfebID++;
    } //loop over all febcounts
    Hlist.Add(hsingle);
    legend->Draw("same");
    //new TCanvas;
    //h->GetXaxis()->SetTitle("rate (Hz)");
    //h->Draw("hist");
    hsingle->GetXaxis()->SetTitle("rate (Hz)");
    new TCanvas; hsingle->Draw();

    cout << "number of febIDs: " << febcounts.size() << endl;
    cout << "nfebID: " << nfebID << endl;
    cout << "total number of counts: " << nfebcounts << endl;
    cout << "number of events:" << events.size() << endl;
    cout << "sampleTime: " << sampleTime << endl;
    cout << "trigger rate across all febs: " << 1.0*nfebcounts / sampleTime << "(Hz)" <<endl;

    //TFile f("/icarus/data/users/aheggest/crtanalysis/rootplots/outfile.root","recreate");
    //TFile f("~/Downloads/research/analysis/crtanalysis/ROOTfiles/outfileTrig.root","RECREATE");

    //TFile f("outfile.root");
    //gStyle->SetOptStat(0);
    Hlist.Write();
    //f.Close();

}

//plot spectra for a given region
void PlotSpectra(/*Int_t reg = 0*/) {
  //TFile fout("~/Downloads/research/analysis/crtanalysis/ROOTfiles/outfile.root","RECREATE");

  struct eveid{
      int eve, run, subrun;

  };

  auto cmp_eveid = [](const eveid& lhs, const eveid& rhs) {
      if( lhs.run != rhs.run ) return lhs.run < rhs.run;
      if( lhs.subrun != rhs.subrun ) return lhs.subrun < rhs.subrun;
      return lhs.eve < rhs.eve;
  };

  DetSim *det = filereader.tmDet();
  Int_t entriesDet = det->GetInputTree()->GetEntries();
  cout << "reading from DetSim Tree with " << entriesDet << " entries...." << endl;
  double pedestal = 63.6;
  double gain = 70.0;
  int ncounts = 0;
   int ncounts_maxchan_dyn = 0;
   int ncounts_full = 0;
   int ncounts_dyn = 0;
   int ncounts_maxchan_full = 0;
   int ncounts_regioncut = 0;

  set<eveid,decltype(cmp_eveid)> events(cmp_eveid);
  map<int,vector<TH1F*>> spectraMap;

  //loop over tracks
  for(Int_t ientry = 0; ientry < entriesDet; ientry++){
    struct eveid eve = {det->Event(),det->Run(),det->SubRun()};
    events.insert(eve);

    det->GetDet(ientry);
    if(det->Mac5()<0) cout << "negative mac5: " << det->Mac5() << endl;
    //if(det->Region()!=reg) continue; //for a given region
    if(det->Region()!= 41 && det->Region()!=47) continue;
    ncounts_regioncut ++;


    if(spectraMap.find(det->Mac5())==spectraMap.end()){
      for(int i=0; i<32; i++){
         string hname = to_string(det->Mac5())+"_"+to_string(i);
         string htitle = "FEB "+to_string(det->Mac5())+" in region "/*+to_string(reg)*/;
         spectraMap[det->Mac5()].push_back(new TH1F(hname.c_str(),htitle.c_str(),121,-0.5,120.5));
         ncounts++;

      }
      //full range: (121,-0.5,120.5), dynamic range: (44, 6.5, 50.5)
    }
    for(int i=0; i<32; i++){
        //if (det->MaxChan()!=i) continue; //if not maxchan, skip
        float pe = ((det->ADC(i) -63.6)/70.0);
        if(det->ADC(i)< 63.6+3.5*70.0) continue; //set pedestal
        ncounts_full++;
        if (6.5 < pe && pe < 50.5) ncounts_dyn++;
        if (det->MaxChan()==i) {
          ncounts_maxchan_full++;
          if (6.5 < pe && pe < 50.5) ncounts_maxchan_dyn++;
        }
        if(det->ADC(i)==0) continue;
        spectraMap[det->Mac5()][i]->Fill((det->ADC(i)-63.6)/70.0); //fill histo for FEBid convert to p.e.u
    }
  }//loop over events/tracks
  //cout << "ncounts: " << ncounts << endl;
  TH1F* havg = new TH1F("havg","Average spectra across FEBs",121,-0.5,120.5);

  int nhist = 1;
  for(auto const& febhists : spectraMap) {
    //cout << "Febhists.first: "<< febhists.first << endl;

    febhists.second[0]->GetXaxis()->SetTitle("p.e.u.");
    febhists.second[0]->SetLineColor(1);
    febhists.second[0]->Sumw2();
    febhists.second[0]->Scale(1.0/febhists.second[0]->Integral());
    //auto c1 = new TCanvas();
    //c1->SetLogy();
    //new TCanvas();
    //febhists.second[0]->Draw("e0hist");
    havg->Add(febhists.second[0],1);
    for(int i=1; i<32; i++){
        febhists.second[i]->GetXaxis()->SetTitle("p.e.u");
        febhists.second[i]->SetLineColor(i);
        febhists.second[i]->Sumw2();
        febhists.second[i]->Scale(1.0/febhists.second[i]->Integral());
        //febhists.second[i]->Draw("e0histsame");
        if (febhists.second[i]->GetEntries()==0) continue;
        int weight = febhists.second[i]->GetEntries();
        //wcout << "weight: " << weight << endl;
        havg->Add(febhists.second[i],weight);
        //havg->Add(febhists.second[i],1);
        nhist++;
        if(febhists.second[i]==0) cout << "Channel is empty!" << endl;
     }

    //c1->Write();
    //c1->SaveAs("~/Downloads/research/analysis/crtanalysis/ROOTfiles/c1.png");
  }
  double sampleTime = events.size()*3.3e-3;
  cout << "ncounts_regioncut: " << ncounts_regioncut << endl;
  cout << "ncounts_maxchan_dyn: " << ncounts_maxchan_dyn << endl;
  cout << "ncounts_maxchan_full: " << ncounts_maxchan_full << endl;
  cout << "ncounts_dyn: " << ncounts_dyn << endl;
  cout << "ncounts_full: " << ncounts_full << endl;
  cout << "nhist: " << nhist << endl;
  cout << "spectraMap.size(): " << spectraMap.size() << endl;
  cout << "number of events: " << events.size() << endl;
  cout << "sampleTime: " << sampleTime << " (s)" << endl;
  //cout << "rate (full): " << 1.0*ncounts_full/sampleTime << " (Hz)" << endl;
  //cout << "rate (dynamic): " << 1.0*ncounts_dyn/sampleTime << " (Hz)" << endl;


  TFile fout("~/Downloads/research/analysis/crtanalysis/ROOTfiles/outfile.root","RECREATE");
  //TFile fout("/icarus/data/users/aheggest/crtanalysis/rootplots/outfile.root","recreate");
  auto c = new TCanvas();
  c->SetLogy();
  havg->GetXaxis()->SetTitle("p.e.u.");
  havg->Scale(1.0/nhist);
  havg->Scale(1.0/havg->Integral());
  cout << "havg->Integral(): " << havg->Integral() << endl;
  havg->Write();
  havg->Draw("e0hist"); //draw average histo across all febs
  //gStyle->SetOptStat(0);
  //c->Write();
  c->SaveAs("~/Downloads/research/analysis/crtanalysis/ROOTfiles/c.png");
  //c->SaveAs("/icarus/data/users/aheggest/crtanalysis/rootplots/c.png");
  fout.Close();
}

void DumpDetSimInfo( ) {

    DetSim *det = filereader.tmDet();
    Int_t entriesDet = det->GetInputTree()->GetEntries();


    for (Int_t ientry = 0; ientry < entriesDet; ientry++) {
      det->GetDet(ientry);
      /*cout << "Event: " << det->Event() << endl;
      cout << "Run: " << det->Run() << endl;
      cout << "SubRun: " << det->SubRun() << endl;
      cout << "NAbove: " << det->NAbove() << endl;
      cout << "T0: " << det->T0() << endl;
      cout << "T1: " << det->T1() << endl;*/
      cout << "ADC for ith entry: "<< det->ADC(ientry) << endl;
      cout << "MaxAdc: " << det->MaxAdc() << endl;
      cout << "MaxChan: " << det->MaxChan() << endl;
    /*  cout << "TrackID for ith entry: "<< det->TrackID(ientry) << endl;
      cout << "PDG for ith entry: "<< det->PDG(ientry) << endl;
      cout << "NTrack: " << det->NTrack() << endl;
      cout << "Entry: " << det->Entry() << endl;
      cout << "Mac5: " << det->Mac5() << endl;
      cout << "Region: " << det->Region() << endl;
      cout << "SubSys: " << det->SubSys() << endl;*/
      cout << "-.-..--..---.-..--.-.--.-.--" << endl;


    }

}
