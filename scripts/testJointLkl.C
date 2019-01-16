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
