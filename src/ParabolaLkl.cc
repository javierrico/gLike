/* ======================================================================== *\
!
!   Author(s): Javier Rico         01/2015 <mailto:jrico@ifae.es>
!
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// ParabolaLkl
//
// This is the simplest form of a Lkl daughter class: it stores the
// value of -2logL vs g, which has been computed elsewhere, and 
// return its value when queried. Used e.g. by FermiTables2016Lkl.
//
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>

#include "TMath.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TLegend.h"

#include "ParabolaLkl.h"

ClassImp(ParabolaLkl);

using namespace std;

// static constants
static const TString  gName            = "ParabolaLkl";
static const TString  gTitle           = "Likelihood vs G";
static const Int_t    gNPars           = 1;         // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"g"};     // Name of parameters
static       Double_t gPDelta[gNPars]  = {0.01};    // Precision of parameters during minimization

// -2logL function for minuit
void parabolaLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;

////////////////////////////////////////////////////////////////
//
// Constructor, just creates the TGraph
//
ParabolaLkl::ParabolaLkl(Int_t n,Double_t* x,Double_t* y,TString name,TString title) :
  Lkl(gNPars," ",name,title)
{
  // if no name/title are provided use default ones
  if(strcmp(name,"")==0)  SetName(gName);
  if(strcmp(title,"")==0) SetTitle(gTitle);

  SetGLklVsG(n,x,y);
}

////////////////////////////////////////////////////////////////
//
// confString constructor
// OJO: to be filled!!!
//
ParabolaLkl::ParabolaLkl(TString fileName,TString name,TString title) :
  Lkl(gNPars," ",name,title)
{
  // if no name/title are provided use default ones
  // if(strcmp(name,"")==0)  SetName(gName);
  // if(strcmp(title,"")==0) SetTitle(gTitle);

  // SetGLklVsG(n,x,y);
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
ParabolaLkl::~ParabolaLkl()
{
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of minuit subtleties)
//
void ParabolaLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN fulllkl
}
				   
////////////////////////////////////////////////////////////////
//
// Check that the samples are correct
//
// Return 0 in case of success, 1 otherwise
//
Int_t ParabolaLkl::MakeChecks()
{
  if(IsChecked()) return 0;

  // check that the TGraph to compute the -2logL values exists
  if(!GetLklVsG()) return 1;

  SetChecked();

  return 0;
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance),
//
void ParabolaLkl::SetFunctionAndPars(Double_t ginit)
{
  Double_t pStart[gNPars]  = {ginit};       // Precision of parameters during minimization
  fMinuit->SetFCN(parabolaLkl);
  fMinuit->SetName(Form("%s_Minuit",GetName()));
  SetParameters(gParName,pStart,gPDelta);
}  

////////////////////////////////////////////////////////////////////////
// parabola likelihood function (-2logL) 
// Free parameters:
// par[0] = g (x-axis of the parabola)
//
void parabolaLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;
  f = 0;

  Double_t g = par[0];

  // get internal object
  minuit->GetObjectFit();  
  ParabolaLkl* mylkl = dynamic_cast<ParabolaLkl*>(minuit->GetObjectFit());
  TGraph* lklVsG = mylkl->GetLklVsG(kFALSE);
  Double_t gmin  = lklVsG->GetX()[0];
  Double_t gmax  = lklVsG->GetX()[lklVsG->GetN()-1];
  Double_t fgmin = lklVsG->GetY()[0];
  Double_t fgmax = lklVsG->GetY()[lklVsG->GetN()-1];
  Double_t fmax  = TMath::Max(fgmin,fgmax);
  
  // return high value if out of bounds
  if(g<gmin)      f = fgmin+(gmin-g)*fmax+(gmin-g);
  else if(g>gmax) f = fgmax+(g-gmax)*fmax+(g-gmax);
  else            f = lklVsG->Eval(g,0,"S");
} 
