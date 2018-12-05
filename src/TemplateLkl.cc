/* ======================================================================== *\
!
!
!   Author(s): Javier Rico         12/2018 <mailto:jrico@ifae.es>
! 
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// TemplateLkl
//
// This is a template to create new likelihood classes
// It contains the minimum such class should contain, and a few guidelines
// to be able to create you own likelihood class
//
// Modify this and have a look to class e.g. Iact1dUnbinned to see a practical
// implementation
//
//////////////////////////////////////////////////////////////////////////////

// include c++ needed classes

// include Root needed classes
#include "TMath.h"

// include gLike needed classes
#include "TemplateLkl.h"

ClassImp(TemplateLkl);

using namespace std;

// class name and title
static const TString  gName            = "TemplateLkl";
static const TString  gTitle           = "Template Likelihood";

// List of free parameters.
// Here we consider the minimal case with just the "g" free parameter
// As many nuisance parameters as wanted can be added to the list
// (check also TemplateLkl::SetFunctionAndPars)

static const Int_t    gNPars           = 1;                      // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"g"};                  // Name of parameters

// -2logL function for minuit
// you can change its name but the format is required by ROOT, you should not change it
// if you change the name, do it consistently in all this file
void templateLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
// also, do not change this
static TMinuit* minuit = NULL;


//////////////////////////////////////////////////////////////////////////////
//
// String constructor
// The string contains the elements for the constructor
// which you can read from e.g. an rc input file
//
TemplateLkl::TemplateLkl(TString inputString) :
  Lkl(gNPars,inputString,gName,gTitle)
{
  if(InterpretInputString(inputString))
    cout << "TemplateLkl::TemplateLkl Warning: there were problems interpreting the input string" << endl;      
}

//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
// The string has been passed to Lkl::InterpretInputString in constructor
// normally read from an input rc file
//
// Return 0 in case checks are Ok, 1 otherwise
//
Int_t TemplateLkl::InterpretInputString(TString inputString)
{
  // decode your inputString and configure the class accordingly

  return 0;
}

////////////////////////////////////////////////////////////////
// 
// Destructor
// make sure you do not leave any memory leaks
//
TemplateLkl::~TemplateLkl()
{
  
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of TMinuit subtleties)
// Leave this unchanged
//
void TemplateLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance)
// in case of need, some of the nuisance parameters can be fixed here
//
void TemplateLkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(templateLkl);  
  fMinuit->SetName(Form("%s_Minuit",GetName()));

  // set and pars initial and step values
  Double_t pStart[gNPars] = {ginit};
  Double_t pDelta[gNPars] = {1};    // Precision of parameters during minimization

  // initialize the free (and nuisance) parameters
  SetParameters(gParName,pStart,pDelta);
}		
		

////////////////////////////////////////////////////////////////
//
// Check that all elements needed for the fit are present, otherwise
// try to create/compute them
//
// Return 0 in case checks are Ok, 1 otherwise
//
Int_t TemplateLkl::MakeChecks()
{
  if(IsChecked()) return 0;

  // add your checks here and try to mend whatever needs to be mended
  
  SetChecked();
  return 0;
}	      
			   

////////////////////////////////////////////////////////////////////////
// likelihood function (-2logL) 
// To be minimized by TMinuit
// Free parameters:
// par[0] = g 
//
// this is a trivial funcion that just returs 0, so probably not a good idea
// to try to minimize it. You must replace this but your likelihood function
//
void templateLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;

  // get internal object, histos, values, etc
  Double_t g       = par[0];
  f = -2*TMath::Log(1);
}		
