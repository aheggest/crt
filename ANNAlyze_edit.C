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

/////////////////////////////////////
// Generator Tree Specific Functions
/////////////////////////////////////


void PlotFlux() {

    Gen *gen = filereader.tmGen(); //pointer to Generator class
    Int_t nentries = gen->GetInputTree()->GetEntries(); //number of entries in generator tree

    const size_t size = 26;
    const double mult = 1.7; //constant multiplier, smaller M = finer bins
    double bins[size];      //intitialize array "bins"
    bins[0] = 0.1;

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

  hmun->GetXaxis()->SetTitle("Energy(GeV)");

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

  hmun->GetYaxis()->SetRangeUser(0.001,100000);


  //set log scale
  auto *c1 = new TCanvas("c1");
  c1->SetLogx();
  c1->SetLogy();

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
}
