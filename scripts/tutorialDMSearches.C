// tutorial illustrating how to compute a flux UL given a DM flux model
#include <iostream>

using namespace std;

void tutorialDMSearches() {
    // input data
    // log_10 of J-factor of the assumed DM source
    const Double_t logJ = 19.; // [GeV^2 cm^-5]
    // statistical error in log_10 on the J-factor of the assumed DM source
    const Double_t DlogJ = 0; // [GeV^2 cm^-5]
    // mass of the DM particle
    const Double_t mass = 1000.; // [GeV]
    // dN/dE input file
    const TString dNdEFileName =
        TString(Form("$GLIKESYS/DM/dNdE/Cirelli/dNdESignal_bb_%.1fmass.root", mass));
    // input file with event list and their associated IRFs
    const TString inputFile1 = "$GLIKESYS/data/root/genericIact_dataIRF_01.root";
    // input file with event list and their associated IRFs
    const TString inputFile2 = "$GLIKESYS/data/root/genericIact_dataIRF_02.root";
    const Double_t errorDef = 4;

    // create and configure an Iact1dUnbinnedLkl object for 1D unbinned
    // likelihood analysis
    Iact1dUnbinnedLkl *unbn = new Iact1dUnbinnedLkl(
        Form("logJ=%.2f DlogJ=0 inputfile=%s", logJ, inputFile1.Data()));
    unbn->ReaddNdESignal(dNdEFileName);
    // set units for DM annihilation <sv>
    unbn->SetDMAnnihilationUnitsForG(mass);

    // create and configure an Iact1dBinnedLkl object for 1D unbinned likelihood
    // analysis
    Iact1dBinnedLkl *bn = new Iact1dBinnedLkl(
        Form("logJ=%.2f DlogJ=0 inputfile=%s", logJ, inputFile2.Data()));
    bn->ReaddNdESignal(dNdEFileName);
    // set units for DM annihilation <sv>
    bn->SetDMAnnihilationUnitsForG(mass);

    // create and fill a JointLkl object for the combined analysis of both
    // datasets
    JointLkl *jointLkl = new JointLkl(Form("DlogJ=%.2f", DlogJ));
    jointLkl->SetErrorDef(errorDef);
    jointLkl->AddSample(unbn);
    jointLkl->AddSample(bn);

    // print how JointLkl looks like
    jointLkl->PrintData();

    // call minimization and print/plot results
    jointLkl->ComputeLklVsG();
    jointLkl->PrintOverview(); // print the details from the fit
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    jointLkl->GetLklVsG()->Draw(); // plot the -2logL vs g curve
}
