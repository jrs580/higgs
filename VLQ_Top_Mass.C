// Training exercise given to Joshua Stewart by Joe Haley
// June 4th, 2019
// Use Assignment 3 to cut down to 1 lepton, then make the following cuts:
// met > 30 Gev
// 4 jets, 2 b tagged

// Author: Joshua Stewart

#include "particleLevel.h"
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>

void VLQ_Top_Mass()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
 
  gROOT->ProcessLine(".L particleLevel.C+"); //Loads particleLevel class.

  particleLevel lp;

  lp.Loop();  //calls upon the function Loop inside the class.
}
