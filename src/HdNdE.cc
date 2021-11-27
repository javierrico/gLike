/* ======================================================================== *\
!
!   Author(s): Javier Rico         11/2021 <mailto:jrico@ifae.es>
! 
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// HdNdE
//
// Extension of the TH1F class with read methods for dNdE histograms, so that
// all read methods can be shared by all Lkl-based classes containing dNdE-type
// of histograms, in such a way that the methods do not need to be implemented
// individually for each of such Lkl-based classes
//
//////////////////////////////////////////////////////////////////////////////

// include c++ needed classes

// include Root needed classes
#include "TH1F.h"
#include "TFile.h"

// include gLike needed classes
#include "HdNdE.h"

ClassImp(HdNdE);

using namespace std;

// class name and title
static const TString  gName   = "HdNdE";
static const TString  gTitle  = "TH1F extended with read methods for dNdE";



//////////////////////////////////////////////////////////////////////////////
//
// Constructor
//
HdNdE::HdNdE(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup)
  : TH1F(name,title,nbinsx,xlow,xup)
{
}

//////////////////////////////////////////////////////////////////////////////
//
// Constructor
//
HdNdE::HdNdE(const char *name, const char *title, Int_t nbinsx, const Double_t *xbins)
  : TH1F(name,title,nbinsx,xbins)
{
}

//////////////////////////////////////////////////////////////////////////////
//
// Constructor
//
HdNdE::HdNdE(const char *name, const char *title, Int_t nbinsx, const Float_t *xbins)
  : TH1F(name,title,nbinsx,xbins)
{
}

//////////////////////////////////////////////////////////////////////////////
//
// Constructor
//
HdNdE::HdNdE(const TH1F &h1f)
  : TH1F(h1f)
{
}

////////////////////////////////////////////////////////////////
// 
// Destructor
// make sure you do not leave any memory leaks
//
HdNdE::~HdNdE()
{
}

