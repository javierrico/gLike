/* ======================================================================== *\
!
!   Author(s): Javier Rico         01/2015 <mailto:jrico@ifae.es>
!
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// IMPORTANT NOTE: THE USE OF THIS CODE TO PRODUCE PAPERS OF THE MAGIC
// AND/OR CTA COLLABORATIONS IS ALLOWED FOLLOWING THEIR RESPECTIVE
// PUBLICATION POLICIES FOR FULL-COLLABORATION PAPERS. FOR
// PUBLICATIONS OUTSIDE THOSE FRAMEWORKS PLEASE CONTACT FIRST THE
// AUTHORS (Jelena Aleksic <jelena@ifae.es> AND Javier Rico
// <mailto:jrico@ifae.es>), WHO COULD CLAIM AUTHORSHIP OF THE
// RESULTING PAPER.
//
// WHEN USING Iact1dUnbinnedLkl CLASS, A REFERENCE SHOULD BE MADE TO THE 
// FOLLOWING PUBLICATION:
// Aleksic, Rico & Martinez JCAP 10 (2012) 032
//
//
// PoissonLkl
// 
// Defines a double-Poisson joint likelihood for observed number of events
// in On and Off regions. 
//
// The free parameter is the number of signal events in the On region (g)
// The number of background events in the On region (b) is a nuisance 
// parameter, constrained by the measurement of the Off regions.
// The ratio of exposures between Off and On samples (tau) is also a 
// nuisance parameter with Gaussian distribution characterized by fTau
// (mean) and fDTau (sigma), that can be fixed making fDTau<=0.
// You can specify a Gaussian uncertainty (fDEff) in the gamma-ray
// detection efficiency (following Rolke et al. 2005) using SetDEff(deff).
//
// Additionally, you can specify a fraction of G to be spilled to the
// Off region (e.g. for very extended sources), and the number of
// foreground events in the On region (e.g. when there is gamma-ray 
// source close to the one of our interest).
// 
// This is the likelihood considered by, e.g. Rolke et al. 
// (NIM A551 (2005) 493–503) for computing upper limits
// or Li & Ma (ApJ 272 (1983) 317) for computing statistical significance
// of observed excess events (fixed tau).
//
// To combine several bins into a joint likelihood use class 
// Iact1dBinnedLkl.
//
//////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "TMath.h"

#include "PoissonLkl.h"

ClassImp(PoissonLkl);

using namespace std;

// static constants
static const TString  gName            = "PoissonLkl";
static const TString  gTitle           = "Double-poisson (classical) Likelihood";
static const Int_t    gNPars           = 4;                           // Number of free+nuisance parameters
const Char_t*  gParName[gNPars] = {"g","b","tau","eff"};                    // Name of parameters
 
// -2logL function for minuit
void poissonLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;

////////////////////////////////////////////////////////////////
// 
// Default constructor
// At least you need to provide number of on and off events
// Also you can provide tau (default = 1) and DeltaTau
// (default = 0, i.e. tau is a fixed parameter)
//
PoissonLkl::PoissonLkl(UInt_t non,UInt_t noff,Float_t tau,Float_t dTau,Float_t dEff,TString name,TString title) : 
  Lkl(gNPars," ",name,title), fNon(non), fNoff(noff), fTau(tau), fDTau(dTau), fDEff(dEff),
  fGFractionInOff(0), fFrgNEvents(0), fKnownBackground(kFALSE)
{
  // if no name/title are provided use default ones
  if(strcmp(name,"")==0)  SetName(gName);
  if(strcmp(title,"")==0) SetTitle(gTitle);
}

////////////////////////////////////////////////////////////////
// 
// consString constructor
// OJO: to be filled
//
PoissonLkl::PoissonLkl(TString confString,TString name,TString title) : 
  Lkl(gNPars," ",name,title), fNon(0), fNoff(0), fTau(1), fDTau(0), fDEff(0),
  fGFractionInOff(0), fFrgNEvents(0), fKnownBackground(kFALSE)
{
  // if no name/title are provided use default ones
  // if(strcmp(name,"")==0)  SetName(gName);
  // if(strcmp(title,"")==0) SetTitle(gTitle);
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
PoissonLkl::~PoissonLkl()
{
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of minuit subtleties)
//
void PoissonLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN fulllkl
}		   
////////////////////////////////////////////////////////////////
//
// Check that everything is correct before calling -2logL minimization
//
// Return 0 in case of success, 1 otherwise
//
Int_t PoissonLkl::MakeChecks()
{
  if(IsChecked()) return 0;
  if(fTau<0)
    {
     cout << "PoissonLkl::MakeChecks Warning: tau cannot be a negative number (it's " 
	  << fTau << ")" << endl;
     return 1;     
    }
  
  SetChecked();
  return 0;
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance),
// fix tau if requested
// fix eff if requested
// fix b   if requested (for tests mostly)
//
void  PoissonLkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(poissonLkl);
  fMinuit->SetName(Form("%s_Minuit",GetName()));
  Double_t pStart[gNPars] = {ginit,fNoff/fTau,fTau,1};
  Double_t pDelta[gNPars] = {0.01,0.01*TMath::Sqrt(fNoff/fTau),0.01,0.01};       // Precision of parameters during minimization
  
  SetParameters(gParName,pStart,pDelta);

  // Fix tau if requested
  fMinuit->Release(PoissonLkl::gTauIndex);
  FixPar(PoissonLkl::gTauIndex,kFALSE);
  if(GetDTau()<=0)
    { 
      fMinuit->FixParameter(PoissonLkl::gTauIndex);
      FixPar(PoissonLkl::gTauIndex);
    }

  // Fix eff if requested
  fMinuit->Release(PoissonLkl::gEffIndex);
  FixPar(PoissonLkl::gEffIndex,kFALSE);
  if(GetDEff()<=0)
    { 
      fMinuit->FixParameter(PoissonLkl::gEffIndex);
      FixPar(PoissonLkl::gEffIndex);
    }

  // Fix b if requested
  fMinuit->Release(PoissonLkl::gBIndex);
  FixPar(PoissonLkl::gBIndex,kFALSE);
  if(fKnownBackground)
    { 
      fMinuit->FixParameter(PoissonLkl::gBIndex);
      FixPar(PoissonLkl::gBIndex);
    }
  
}  

//////////////////////////////////////////////////////////////////
//
// Print info about the experimental data
//
void PoissonLkl::PrintData(Int_t level) 
{
  Lkl::PrintData(level);
  
  Margin(level); cout << "                     Non  = " << fNon  << endl;
  Margin(level); cout << "                     Noff = " << fNoff << endl;
  Margin(level); cout << "            Measured tau  = " << fTau << " +/- " << fDTau << endl;
  Margin(level); cout << "            Measured Deff = " << "1 +/- " << fDEff << endl;
  Margin(level); cout << "     Fraction of G in Off = " <<  fGFractionInOff << endl;
  Margin(level); cout << "  Foreground events in On = " <<  fFrgNEvents << endl;
  Margin(level); cout << "    Background in On (b) is " <<  (fKnownBackground? "FIXED" : "NUISANCE") << endl;
}

////////////////////////////////////////////////////////////////////////
//
// Poisson likelihood function (-2logL) 
// To be minimized by TMinuit
// Free+nuisance parameters:
// par[0] = g (total estimated number of signal events in On region)
// par[1] = b (total estimated number of background events in On region)
// par[2] = estimated value of tau
// par[3] = estimated value of efficiency (e in Rolke et al. 2005)
//
void poissonLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;
  
  // get internal object, histos, values, etc
  PoissonLkl*     mylkl         = dynamic_cast<PoissonLkl*>(minuit->GetObjectFit());
  UInt_t           Non           = mylkl->GetNon();
  UInt_t           Noff          = mylkl->GetNoff();
  Float_t          tau           = mylkl->GetTau();
  Float_t          dTau          = mylkl->GetDTau();
  Float_t          dEff          = mylkl->GetDEff();
  Bool_t           tauIsNuisance = !(mylkl->IsParFixed(PoissonLkl::gTauIndex));
  Bool_t           effIsNuisance = !(mylkl->IsParFixed(PoissonLkl::gEffIndex));
 
  // Estimated number of background events in signal and background regions
  Double_t effest  = par[3];
  Double_t g       = par[0]*effest;
  Double_t b       = par[1];
  Double_t trytau  = par[2];
  Double_t boff    = b*trytau;
  Double_t goff    = g*mylkl->GetGFractionInOff()*trytau;
  Double_t frg     = mylkl->GetFrgNEvents();
  if(frg<0) frg=0;
  
  // compute likelihood
  if(b<0 || g+b<0 || trytau<0)
    f = 1e99;
  else
    {
      f = -2*TMath::Log(TMath::Poisson(Non,g+b+frg))-2*TMath::Log(TMath::Poisson(Noff,goff+boff));
      if(tauIsNuisance)	
	f += -2*TMath::Log(TMath::Gaus(trytau, tau, dTau, kTRUE));
      if(effIsNuisance)
	f += -2*TMath::Log(TMath::Gaus(effest, 1, dEff, kTRUE));
    }
}
