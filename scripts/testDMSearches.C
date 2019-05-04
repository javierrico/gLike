// Macro testDMSearches.C
// Author: J. Rico
// Date: Jan 2019
// For beginners, to understand the basic usage of the Iact1dUnbinned and Iact1dBinned classes
// IMPORTANT NOTE: for some still-to-be-understood "feature", macros
// containing Lkl-based objects MUST be run in compiled mode, i.e.
// run this macro with:
// .x testDMSearches.C+

#include <iostream>
#include <fstream>
#include "Iact1dUnbinnedLkl.h"
#include "Iact1dBinnedLkl.h"
#include "JointLkl.h"
#include "TCanvas.h"

void testDMSearches()
{
  // input data
  const Double_t logJ         = 19.;   // [GeV^2 cm^-5] log_10 of J-factor of the assumed DM source
  const Double_t DlogJ        = 0;     // [GeV^2 cm^-5] statistical error in log_10 of J-factor of the assumed DM source
  const Double_t mass         = 1000.; // [GeV] mass of the DM particle
  const TString  dNdEFileName = TString(Form("./DM/dNdE/Cirelli/dNdESignal_bb_%.1fmass.root",mass)); // dN/dE input file  
  const TString  inputFile1   = "./data/genericIact_dataIRF_01.root";  // input file with event list and their associated IRF
  const TString  inputFile2   = "./data/genericIact_dataIRF_02.root";  // input file with event list and their associated IRF
  const Double_t errorDef     = 4;

  // create and configure an Iact1dUnbinnedLkl object for 1D unbinned likelihood analysis
  Iact1dUnbinnedLkl* unbn = new Iact1dUnbinnedLkl(Form("logJ=%.2f DlogJ=0 inputfile=%s",logJ,inputFile1.Data()));
  unbn->ReaddNdESignal(dNdEFileName);
  unbn->SetDMAnnihilationUnitsForG(mass);    // set units for DM annihilation <sv>

  // create and configure an Iact1dBinnedLkl object for 1D unbinned likelihood analysis
  Iact1dBinnedLkl* bn = new Iact1dBinnedLkl(Form("logJ=%.2f DlogJ=0 inputfile=%s",logJ,inputFile2.Data()));
  bn->ReaddNdESignal(dNdEFileName);
  bn->SetDMAnnihilationUnitsForG(mass);    // set units for DM annihilation <sv>

  // create and fill a JointLkl object for the combined analysis of both datasets
  JointLkl* jointLkl = new JointLkl(Form("DlogJ=%.2f",DlogJ));
  jointLkl->SetErrorDef(errorDef);
  jointLkl->AddSample(unbn);
  jointLkl->AddSample(bn);

   // print how JointLkl looks like
  jointLkl->PrintData();

  std::ofstream data;
  Double_t svStep = 0.;
  //Double_t svMin = 1.e-28;
  //Double_t svMax = 1.e-22;
  Double_t svMin = 1.e-28;
  Double_t svMax = 1.e-18;
  Int_t svNPoints = 1000;
  Double_t isv = svMax;
  Double_t isv2 = svMin;
  Double_t svScan[svNPoints+1];
  double svScan2[svNPoints+1];
  Double_t glow = -1000.;
  Double_t gupp = +1000.;

  cout << "units of G = " << jointLkl->GetUnitsOfG() << "   glow = " << svMin/jointLkl->GetUnitsOfG() << "     gupp = " << svMax/jointLkl->GetUnitsOfG() << "   gmin = " << jointLkl->GetLklMin() << endl;
  jointLkl->MinimizeLkl();
  //Double_t svMin_2 = (glow+jointLkl->GetLklMin())*jointLkl->GetUnitsOfG();
  //Double_t svMax_2 = (gupp+jointLkl->GetLklMin())*jointLkl->GetUnitsOfG();
  Double_t svMin_2 = svMin;
  Double_t svMax_2 = svMax;
  cout << "units of G = " << jointLkl->GetUnitsOfG() << "   glow = " << svMin_2/jointLkl->GetUnitsOfG() << "     gupp = " << svMax_2/jointLkl->GetUnitsOfG() << "   gmin = " << jointLkl->GetLklMin() << "   g for lkl 0 = " << jointLkl->GetGForLkl(0.,kTRUE) << endl;
  //cout << "units of G = " << jointLkl->GetUnitsOfG() << "   glow = " << svMin_2/jointLkl->GetUnitsOfG()-jointLkl->GetLklMin() << "     gupp = " << svMax_2/jointLkl->GetUnitsOfG()-jointLkl->GetLklMin() << "   gmin = " << jointLkl->GetLklMin() << "   g for lkl 0 = " << jointLkl->GetGForLkl(0.,kTRUE) << endl;
  cout << "units of G = " << jointLkl->GetUnitsOfG() << "   glow = " << svMin_2/jointLkl->GetUnitsOfG() << "     gupp = " << svMax_2/jointLkl->GetUnitsOfG() << "   gmin = " << jointLkl->GetLklMin() << "   g for lkl 0 = " << jointLkl->GetGForLkl(0.,kTRUE) << endl;

  //gSystem->Exec(Form("mkdir -p %s/Data/txt",fPlotsDir.Data()));
  //TString dataFile = fPlotsDir+"Data/txt/"+label+".txt";
  TString dataFile = "GloryDuck_MAGIC_MC_test_file.txt";
  data.open(dataFile);

  cout << "logJ = " << logJ << "   ";
  data << logJ << "        " << mass << endl;
  Int_t counter = 0;
  svStep = TMath::Exp(TMath::Log(TMath::Abs(svMax/svMin))/svNPoints);
  for(Int_t i=0; i<=svNPoints; i++)
    {
      //data << " " << isv << " ";
      //cout << " " << isv << " " << isv/jointLkl->GetUnitsOfG()/*-jointLkl->GetLklMin()*/ << " " << endl; //counter << endl;
      svScan[counter] = isv;
      svScan2[counter] = isv2/jointLkl->GetUnitsOfG()-jointLkl->GetLklMin();
      //cout << " " << isv2 << " " << isv2/jointLkl->GetUnitsOfG()-jointLkl->GetLklMin() << " " << endl; //counter << endl;
      counter++;
      isv/=svStep;
      isv2*=svStep;
    }

  //cout << mass << "      " ;
  //data << endl;

  Double_t diameter = (svMax-svMin)/(2.*jointLkl->GetUnitsOfG());
  cout << "Diam = " << diameter << endl;

  // call minimization and print/plot results
  jointLkl->ComputeLklVsG();
  //jointLkl->ComputeLklVsG(kFALSE,10000,svMin_2/jointLkl->GetUnitsOfG()-diameter,svMax_2/jointLkl->GetUnitsOfG()-diameter,/*svScan2,*/kTRUE,kTRUE);
  //jointLkl->ComputeLklVsG(kFALSE,1000,svMin_2/jointLkl->GetUnitsOfG()-jointLkl->GetLklMin(),svMax_2/jointLkl->GetUnitsOfG()-jointLkl->GetLklMin(),svScan2,kTRUE,kTRUE);
  //jointLkl->ComputeLklVsG(kFALSE,svNPoints,svMin_2/jointLkl->GetUnitsOfG()-jointLkl->GetLklMin(),svMax_2/jointLkl->GetUnitsOfG()-jointLkl->GetLklMin(),kFALSE,kTRUE);
  jointLkl->PrintOverview();     // print the details from the fit
  TCanvas* c1 = new TCanvas("c1","",700,500);
  jointLkl->GetLklVsG()->Draw(); // plot the -2logL vs g curve

  //data << mass << "      " ;
  for (Int_t isv=0; isv<=svNPoints; isv++)
    {
      //cout << "sv = " << svScan[isv] << " offset sv = " << jointLkl->GetGForLkl(0.,kTRUE) << " gsv = " << svScan[isv]+jointLkl->GetGForLkl(0.,kTRUE) << " g = " << svScan[isv]/jointLkl->GetUnitsOfG() << " offsetg = " << jointLkl->GetGForLkl(0.,kFALSE) << " gsv2 = " << svScan[isv]/jointLkl->GetUnitsOfG()+jointLkl->GetGForLkl(0.,kFALSE) << " res = " << jointLkl->GetLklVsG(kTRUE)->Eval(svScan[isv]+jointLkl->GetGForLkl(0.,kTRUE)) << " res2 = " << jointLkl->GetLklVsG(kFALSE)->Eval(svScan[isv]/jointLkl->GetUnitsOfG()+jointLkl->GetGForLkl(0.,kFALSE)) << endl;
      data << svScan[isv] << "   " << jointLkl->GetLklVsG()->Eval(svScan[isv]+jointLkl->GetGForLkl(0.,kTRUE)) << endl;
    }
  // go to next line in the file
  //cout << endl;
  //data << endl;
  data.close();
}
