// Macro testPoissonLkl.C
// Author: J. Rico
// Date: Jan 2019
// For beginners, to understand the basic works of a simple Lkl-based object
// and, in particular, a standalone PoissonLkl object, configured with
// different level of complexities (i.e. number of nuisance parameters)
// IMPORTANT NOTE: for some still-to-be-understood "feature", macros
// containing Lkl-based objects MUST be run in compiled mode, i.e.
// run this macro with:
// .x testPoissonLkl.C+

#include <iostream>
#include "PoissonLkl.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRolke.h"

using namespace std;

void testPoissonLkl()
{
  // create and configure the simplest possible PoissonLkl object
  const Int_t Non = 130;
  const Int_t Noff = 90;
  const Double_t tau = 1.;
  PoissonLkl* p = new PoissonLkl(Non,Noff,tau);
  
  // configure the minimization algorithm
  const Double_t errorDef = 4;
  const Double_t unitsOfG = 2;
  p->SetErrorDef(errorDef);
  p->SetUnitsOfG(unitsOfG);

  // call the minimization
  p->ComputeLklVsG();

  // print the fit results
  cout << endl << "PrintOverview:" << endl;
  p->PrintOverview(); // print the details from the fit

  // plot fit results
  TCanvas* c1 = new TCanvas("canvas1", "canvas1", 800, 600);
  TGraph* lklVsG = p->GetLklVsG(); // the -2logL vs g curve
  TLegend* legend = new TLegend(0.5, 0.8, 0.9, 0.9);
  legend->SetTextSize(0.03);
  lklVsG->SetLineWidth(2);
  lklVsG->GetXaxis()->SetTitle("g");
  lklVsG->GetXaxis()->CenterTitle();
  lklVsG->GetYaxis()->SetTitle("-2 log L");
  lklVsG->GetYaxis()->CenterTitle();
  lklVsG->Draw();  
  legend->AddEntry(lklVsG, "Poisson Likelihood");

  // access the gLike results
  const Double_t gmin = p->GetGLklMin();                // get the value of g that minimizes -2logL
  const Double_t gerr = p->GetGLklMinErr();             // get the value of gerr such that -2logL(gmin+/-gerr)=-2logLmin+fErrorDef
  const Double_t sig  = sqrt(p->GetLklVsG()->Eval(0)) ; // compute significance of detection

  // compute significance with LiMa formula for comparison
  const Double_t LiMasig = sqrt(2.*(Non*log((tau+1.)*(1.*Non/(Non+Noff)))+Noff*log((tau+1)/tau*(1.*Noff/(Non+Noff)))));

  // compare confidence interval with Rolke method for comparison
  TRolke* r = new TRolke(0.9544);
  r->SetPoissonBkgKnownEff(Non,Noff,tau,1);
  Double_t low,upp;
  r->GetLimits(low,upp);
  Double_t grolke    = unitsOfG*(upp+low)/2;
  Double_t grolkeerr = unitsOfG*(upp-low)/2;

  // print and compare results
  cout << endl;
  cout << "The max lkl  value of g = " << gmin << ",\t with "<< sqrt(errorDef) << "-sigma error bar = " << gerr << endl;
  cout << "The Rolke    value of g = " << grolke << ",\t with 2-sigma error bar = " << grolkeerr << endl;
  cout << "The gLike significance of the detection of signal is " << sig     << " sigma" << endl;
  cout << "The Li&Ma significande of the detection of signal is " << LiMasig << " sigma" << endl;

  // let us consider, in the PoissonLkl, 1% uncertainties on the gamma-ray 
  // detection efficiency (Deff) and on the exposure ratio between the 
  // On and Off samples (Dtau) 
  Double_t Deff  = 0.1;
  Double_t Dtau  = 0.1;
  p->SetDEff(Deff);
  p->SetDTau(Dtau);
  // repeat the Likelihood computation
  p->ComputeLklVsG();
  p->PrintOverview();    
  // compare Likelihood profiles w/ and w/o syst. uncertainties
  TGraph* lklVsGUnc = p->GetLklVsG(); 
  lklVsGUnc->SetLineColor(2);
  lklVsGUnc->SetLineWidth(2);
  lklVsGUnc->Draw("same"); // plot the -2logL vs g curve
  legend->AddEntry(lklVsGUnc, "Poisson Likelihood + syst. unc.");
  legend->Draw("same");
}
