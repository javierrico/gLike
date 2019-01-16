// Macro testJointLkl.C
// Author: J. Rico
// Date: Jan 2019
// For beginners, to understand the basic usage of the JointLkl class
// IMPORTANT NOTE: for some still-to-be-understood "feature", macros
// containing JointLkl objects MUST be run in compiled mode, i.e.
// run this macro with:
// .x testJointLkl.C+

#include <iostream>
#include "PoissonLkl.h"
#include "JointLkl.h"
#include "TCanvas.h"

void testJointLkl()
{
  // define number of PoissonLkl objects that we will use in the JointLkl and their configuration parameters
  const Int_t nsamples = 2;
  const Int_t Non[nsamples]  = {130,90};
  const Int_t Noff[nsamples] = {105,110};

  // create the JointLkl object that will collect the PoissonLkl objects
  const Double_t errorDef = 4;
  JointLkl* j = new JointLkl;
  j->SetErrorDef(errorDef);
    
  // configure samples and introduce them into the JointLkl
  PoissonLkl* p[nsamples];
  for(Int_t isample=0;isample<nsamples;isample++)
    {
      p[isample] = new PoissonLkl(Non[isample],Noff[isample]);
      p[isample]->SetName(Form("PoissonLkl_%1d",isample));
      j->AddSample(p[isample]);
    }

  // print how JointLkl looks like
  j->PrintData();

  // call minimization and print/plot results
  j->ComputeLklVsG();
  j->PrintOverview();     // print the details from the fit
  TCanvas* c1 = new TCanvas("c1","",700,500);
  j->GetLklVsG()->Draw(); // plot the -2logL vs g curve

  // compare it to each PoissonLkl object separately
  for(Int_t isample=0;isample<nsamples;isample++)
    {
      p[isample]->SetErrorDef(errorDef);      
      p[isample]->ComputeLklVsG();
      p[isample]->PrintOverview();     // print the details from the fit
      TGraph* lklvsg = p[isample]->GetLklVsG();
      lklvsg->SetLineColor(isample+2); // plot the -2logL vs g curve
      lklvsg->Draw("same");            // plot the -2logL vs g curve
    }

  // We assume g is the number of gamma rays measured by some Cherenkov
  // telescope from a given steady source in the On region, during a given
  // known observation time, with known IRF and threshold.
  // The spectrum of the source is a power law with index -2,
  // and unknown normalization k = (dPhi/dE)(E=100GeV)
  // We consider the following simplifications:
  // - constant effective area with energy
  // - perfect energy resolution

  // these are the values of observation time, IRF, etc...
  const Double_t Emin[nsamples]  = {100,200};           // [GeV] energy threshold
  const Double_t Emax[nsamples]  = {10000,10000};       // [GeV] maximum measurable energy
  const Double_t Tobs[nsamples]  = {2*60*60,1.8*60*60}; // [s] observation times
  const Double_t Aeff[nsamples]  = {1e9,0.98e9};        // [cm^2] effective area (considered constant)
  Double_t IntegraldFdE[nsamples];                      // [GeV] integral of the spectral shape function betwen Emin and Emax
  Double_t unitsOfG[nsamples];                          // [s-1 cm-2 GeV-1] 1/(Tobs*Aeff*IntegraldFdE)
  for(Int_t isample=0;isample<nsamples;isample++)
    {
      IntegraldFdE[isample] = 100.*100.*(1./Emin[isample]-1./Emax[isample]);
      unitsOfG[isample]     = 1./(Tobs[isample]*Aeff[isample]*IntegraldFdE[isample]);
      p[isample]->SetUnitsOfG(unitsOfG[isample]);
    }

  // print how JointLkl looks like
  cout << endl << "##################################################" << endl;
  j->PrintData();

  // call minimization and print/plot results
  j->ComputeLklVsG();
  j->PrintOverview();     // print the details from the fit
  TCanvas* c2 = new TCanvas("c2","",700,500);
  j->GetLklVsG()->Draw(); // plot the -2logL vs g curve

   // compare it to each PoissonLkl object separately
  for(Int_t isample=0;isample<nsamples;isample++)
    {
      p[isample]->SetErrorDef(errorDef);      
      p[isample]->ComputeLklVsG();
      p[isample]->PrintOverview();     // print the details from the fit      
      TGraph* lklvsg = p[isample]->GetLklVsG();
      lklvsg->SetLineColor(isample+2); // plot the -2logL vs g curve
      lklvsg->Draw("same");            // plot the -2logL vs g curve
    }

  
}
