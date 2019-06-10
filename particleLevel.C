#define particleLevel_cxx
#include "particleLevel.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>


#define GeV 0.001

void particleLevel::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L particleLevel.C
  //      root> particleLevel t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;

  fChain->SetBranchStatus("*",1);  // enable all branches
 
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;

  int n=0;

  TH1F* VLQ_T_Lep_Mass = new TH1F("VLQ_T_Lep_Mass","VLQ_T_Lep_Mass",100,0,3000);
  TH1F* VLQ_T_Had_Mass = new TH1F("VLQ_T_Had_Mass","VLQ_T_Had_Mass",100,0,3500);

  //Loop through events.
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    vector <int> ei, mi;
    bool byes;
    TLorentzVector met4, lepton4, T_lep4, T_lep4_0, T_lep4_1, T_had4, T_had4_0, T_had4_1;
    vector <TLorentzVector> other_jets4, btagged_jets4;  
    float T_had_mass, T_lep_mass;

    //loop to test electrons
    for(uint i=0; i < (*el_pt).size(); ++i){
      if((*el_pt)[i]*GeV > 30 && fabs((*el_eta)[i]) < 2.0){
	ei.push_back(i);
      }
    }

    //loop to test muons
    for(uint i=0; i < (*mu_pt).size(); ++i){
      if((*mu_pt)[i]*GeV > 30 && fabs((*mu_eta)[i]) < 2.0){
	mi.push_back(i);
      }
    }

    //loop to test jets
    for(uint i=0; i < (*jet_pt).size(); ++i){
      if((*jet_pt)[i]*GeV > 30 && fabs((*jet_eta)[i]) < 2.5){
	// this is a good jet

	TLorentzVector jet4;
	jet4.SetPtEtaPhiE((*jet_pt)[i], (*jet_eta)[i], (*jet_phi)[i], (*jet_e)[i]);
	
	byes = false;

	// does this jet also match bquark b? if so, increment n_bjets and fill TLorentzVectors
	for(uint b=0; b < (*bquark_pt).size(); ++b){
	  
	  TLorentzVector jet_b4;
	  jet_b4.SetPtEtaPhiE((*bquark_pt)[b], (*bquark_eta)[b], (*bquark_phi)[b], (*bquark_e)[b]);

	  float R = jet_b4.DeltaR(jet4); //calculate R^2 value = eta and phi difference.
	  
	  if(R < 0.2){
	    byes = true;
	  } 
	}
	//Set up b tagged and regular vectors
	if(byes == false){
	  other_jets4.push_back(jet4);
	}
	else{
	  btagged_jets4.push_back(jet4);
	}
      }
    }

    //test to make sure we get events with 1 lepton and met > 30 GeV.
    if(ei.size()+mi.size() != 1) continue;
    if(met_met < 30) continue;
    if(btagged_jets4.size() + other_jets4.size() < 4) continue;
    if(btagged_jets4.size() < 2) continue;
    
    //fill the lepton 4vector
    if(ei.size() == 1){
      lepton4.SetPtEtaPhiE((*el_pt)[0], (*el_eta)[0], (*el_phi)[0], (*el_e)[0]);
    }
    else{
      lepton4.SetPtEtaPhiE((*mu_pt)[0], (*mu_eta)[0], (*mu_phi)[0], (*mu_e)[0]);
    }

    //met 4 vector
    met4.SetPtEtaPhiE(met_met, 0, met_phi, met_met);

    //add together 4 vectors with different combos
    T_lep4_0 = met4 + lepton4 + btagged_jets4[0];
    T_lep4_1 = met4 + lepton4 + btagged_jets4[1];
    T_had4_0 = other_jets4[0] + other_jets4[1] + btagged_jets4[0];
    T_had4_1 = other_jets4[0] + other_jets4[1] + btagged_jets4[1];

    //  get mass and test to see which combo is best
    float D1, D2;

    D1 = fabs(T_lep4_0.M() - T_had4_1.M());
    D2 = fabs(T_lep4_1.M() - T_had4_0.M());

    T_lep_mass = T_lep4_0.M();
    T_had_mass = T_had4_1.M();

    if(D2 < D1){
      T_lep_mass = T_lep4_1.M();
      T_had_mass = T_had4_0.M();
    }
    
    //fill the histogram
    VLQ_T_Lep_Mass->Fill(T_lep_mass*GeV);
    VLQ_T_Had_Mass->Fill(T_had_mass*GeV);

    ++n;
  }

  //make canvases to draw onto
  TCanvas* c0 = new TCanvas("c0","c0"); 
  VLQ_T_Lep_Mass->Draw();
  TCanvas* c1 = new TCanvas("c1","c1");
  VLQ_T_Had_Mass->Draw();

  c0->SaveAs("VLQ_T_Lep_Mass.png");
  c1->SaveAs("VLQ_T_Had_Mass.png");

  cout << n << endl;
}

