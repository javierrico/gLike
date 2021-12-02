/* ======================================================================== *\
!
!   Author(s): Javier Rico         11/2021 <mailto:jrico@ifae.es>
! 
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// HdNdE
//
// Abstract class to be inherited by all the Lkl classes that use a dNdE-type
// histogram, in such a way that the methods do not need to be implemented
// individually for each of such Lkl-based classes
//
//////////////////////////////////////////////////////////////////////////////

// include c++ needed classes
#include <iostream>

// include Root needed classes
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"

// include gLike needed classes
#include "HdNdE.h"

ClassImp(HdNdE);

using namespace std;

// class name and title
static const TString  gName      = "fHdNdESignal";
static const TString  gTitle     = "HdNdE-based abstract likelihood";
static const TString  gXTitle    = "log_{10}(E [GeV])";
static const TString  gYTitle    = "dN/dE [GeV^{-1}]";
static const Double_t gCenterBin = 0.5;  // decide which value represents a bin

//////////////////////////////////////////////////////////////////
// 
// Constructor
// input is the same as for TH1F, but histogram will not be inmediately created
//
HdNdE::HdNdE(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup)
  : fHdNdESignal(NULL), fName(name), fTitle(title), fNFineBins(nbinsx), fFineLEMin(xlow), fFineLEMax(xup), fHdNdESignalChecked(kFALSE)
{
}


//////////////////////////////////////////////////////////////////
// 
// Destructor
//
HdNdE::~HdNdE()
{
  if(fHdNdESignal) delete fHdNdESignal;
}

//////////////////////////////////////////////////////////////////
// 
// To be called after the histogram is created
//
Int_t HdNdE::SetHistogramPars()
{
  if(!fHdNdESignal)
    {
      cout << "HdNdE::SetHistogramPars Warning: histogram does not exist" << endl;
      return 1;
    }
  
  fHdNdESignal->SetDirectory(0);
  if(strcmp(fName,"")==0) fHdNdESignal->SetName(gName);
  if(strcmp(fTitle,"")==0) fHdNdESignal->SetName(gTitle);
  fHdNdESignal->SetXTitle(gXTitle);
  fHdNdESignal->SetYTitle(gYTitle);
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Set directly the dN/dE histogram for signal (e.g. from Damasco)
// Replacement of existing histo is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t HdNdE::SetdNdESignal(TH1F* hdNdESignal)
{
  // pathologies
  if(!hdNdESignal)
    {
      cout << "HdNdE::SetdNdESignal Warning: input histo does not exist" << endl;
      return 1;
    }
  
  if(ResetdNdESignal())
    return 1;
      
  ReadAndInterpolate(hdNdESignal,fHdNdESignal);
  
  return 0;
}


//////////////////////////////////////////////////////////////////
// 
// Recreate a fresh new version of fHdNdESignal
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t HdNdE::ResetdNdESignal()
{
  // Delete existing fHdNdESignal and create empty one
  if(fHdNdESignal) 
    delete fHdNdESignal; 
  
  // Create histo
  fHdNdESignal = new TH1F(fName,fTitle,fNFineBins,fFineLEMin,fFineLEMax);

  // check it was successfully created
  if(!fHdNdESignal)
    {
      cout << "HdNdE::ResetdNdESignal Warning: problem creating fHdNdESignal histogram" << endl;
      return 1;
    }

  // set default parameters for this histogram
  if(SetHistogramPars())
    return 1;

  // Set unchecked status
  fHdNdESignalChecked=kFALSE;
  
  // exit
  return 0;
}


//////////////////////////////////////////////////////////////////
// 
// Add a new component to previously existing dN/dE histogram
// for signal using histogram from file <filename> and assuming
// a branching ratio <br>
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t HdNdE::AdddNdESignal(TString filename,Float_t br)
{
  // open file and look for histo
  TFile* dNdESignalFile  = new TFile(filename);
  TH1F*  hdNdESignal = (TH1F*) dNdESignalFile->Get("hdNdE");

  // pathologies
  if(!hdNdESignal)
    {
      cout << "HdNdE::AdddNdESignal Warning: input histo does not exist" << endl;
      return 1;
    }

  // transform the histogram to the Iact1dUnbinnedLkl internal format
  TH1F* fHdNdESignal_temp = new TH1F(*fHdNdESignal);
  if(ReadAndInterpolate(hdNdESignal,fHdNdESignal_temp))
    return 1;

  for(Int_t ibin=0;ibin<fHdNdESignal->GetNbinsX();ibin++)
    fHdNdESignal->SetBinContent(ibin+1,fHdNdESignal->GetBinContent(ibin+1)+br*(fHdNdESignal_temp->GetBinContent(ibin+1)));

  // clean and exit
  delete dNdESignalFile;
  fHdNdESignalChecked = kFALSE;
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Create and save a dN/dE histogram for signal according toa 
// given <function> of parameters <p0>, <p1>, <p2>... (add more in case of need...)
// and branching ration <br>
//
// (see HdNdE::AdddNdESignalFunction for details on available
// functions)
//
// If fdNdESignal exists it will be replaced
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t HdNdE::SetdNdESignalFunction(TString function,Float_t br,Float_t p0,Float_t p1,Float_t p2)
{  
  // Delete existing fHdNdESignal and create empty one
  ResetdNdESignal();

  // call to add the function
  if(AdddNdESignalFunction(function,br,p0,p1,p2))
    return 1;
  
  // exit
  return 0;
}
//////////////////////////////////////////////////////////////////
// 
// Add to a previously existing dN/dE histogram for signal according to 
// given <function> and parameters <p0>, <p1>, <p2> ... (add more in case of need...)
// and branching ration <br>
//
// If <function>=="line"
// <p0> = E0, the energy of the line
// <p1> = scale, i.e. the number of photons per annihilation (e.g 2 for XX->gamma gamma)
//
// If <function>=="box"
// <p0> = Emin, lower border of the box
// <p1> = Emax, upper border of the box
// <p2> = scale, i.e. the number of photons per annihilation (e.g 4 for XX->pi0 pi0)
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t HdNdE::AdddNdESignalFunction(TString function,Float_t br,Float_t p0,Float_t p1,Float_t p2)
{    
  // Check function
  enum function_t {gLine=0,gBox};
  Int_t functionType=-1;
  if(function.CompareTo("line",TString::kIgnoreCase)==0)
    functionType=gLine;
  else if(function.CompareTo("box",TString::kIgnoreCase)==0)
    functionType=gBox;
  else
    {
      cout << "HdNdE::AdddNdESignalFunction Error: no known function type (" << function << "), check documentation for allowed values" << endl;
      return 1;
    }

  // get  histogram binning parameters
  Float_t xMin  = fHdNdESignal->GetXaxis()->GetXmin();
  Float_t xMax  = fHdNdESignal->GetXaxis()->GetXmax();
  Int_t   nBins = fHdNdESignal->GetNbinsX();
  
  // Add to the histogram according to specified function and parameters
  if(functionType==gLine)
    {
      Float_t E0    = p0;
      Float_t scale = p1;
      
      Float_t log10E0   = TMath::Log10(E0);
 
      
      // check range
      if(xMin>log10E0 || xMax<log10E0)
	{
	  cout << "HdNdE::AdddNdESignalFunction Warning: trying to set a line at an energy E0="
	       << E0 << ", out of range (" << TMath::Power(10,xMin) << "," << TMath::Power(10,xMax) << ")" << endl;
	  return 1;
	}
         
      // fill histo
      Int_t ibin   = nBins*(log10E0-xMin)/(xMax-xMin);
      Float_t Emin = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibin+1));
      Float_t Emax = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibin+1)+fHdNdESignal->GetBinWidth(ibin+1));
      Float_t dE   = Emax-Emin;
      fHdNdESignal->SetBinContent(ibin+1,fHdNdESignal->GetBinContent(ibin+1)+br*(scale/dE));
    }
  else if(functionType==gBox)
    {      
      // check range
      Float_t Emin  = p0;
      Float_t Emax  = p1;
      Float_t scale = p2;

      if(Emin>Emax)
	{
	  cout << "HdNdE::AdddNdESignalFunction Warning: Emin (" << Emin << ") larger than Emax (" << Emax << ")" << endl;
	  return 1;
	}
      
      Float_t log10Emin = TMath::Log10(Emin);
      Float_t log10Emax = TMath::Log10(Emax);

         
      // fill histo
      Int_t ibinmin = nBins*(log10Emin-xMin)/(xMax-xMin);
      Int_t ibinmax = nBins*(log10Emax-xMin)/(xMax-xMin);
      
      Float_t realEmin  = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibinmin+1));
      Float_t realEmax  = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibinmax+1)+fHdNdESignal->GetBinWidth(ibinmax+1));
      Float_t dE        = realEmax-realEmin;
      for(Int_t ibin=ibinmin;ibin<=ibinmax;ibin++)
	fHdNdESignal->SetBinContent(ibin+1,fHdNdESignal->GetBinContent(ibin+1)+br*(scale/dE));
    }
  else
    {
      cout << "HdNdE::AdddNdESignalFunction Warning: this should never happen... " << endl;
      return 1;
    }
      
  // exit
  fHdNdESignalChecked = kFALSE;
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Create and save a dN/dE histogram for signal according to 
// given <function> (TF1 object, created and parameterized before
// passing it as an argument)
//
// Set only values for which energy is between <emin> and <emax>
//
// If a previous histogram exists it will be replaced
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t HdNdE::SetdNdESignalFunction(TF1* function,Float_t br,Float_t emin,Float_t emax)
{  
  // Delete existing fHdNdESignal and create empty one
  ResetdNdESignal();

  // call to add the function
  if(AdddNdESignalFunction(function,br,emin,emax))
    return 1;
  
  // exit
  return 0;
}
//////////////////////////////////////////////////////////////////
// 
// Add to a previously existing dN/dE histogram for signal according to 
// given <function> (TF1 defined and parameterized outside this method)
//
// Add only values for which energy is between <emin> and <emax>
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t HdNdE::AdddNdESignalFunction(TF1* function,Float_t br,Float_t emin,Float_t emax)
{
  // Check function
  if(!function)
    {
      cout << "HdNdE::AdddNdESignalFunction Error: function does not exist" << endl;
      return 1;
    }
  
  for(Int_t ibin=0;ibin<fHdNdESignal->GetNbinsX();ibin++)
    {
      Float_t etest = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibin+1)+gCenterBin*fHdNdESignal->GetBinWidth(ibin+1));
      if(etest>emin && etest<emax)
      fHdNdESignal->SetBinContent(ibin+1,fHdNdESignal->GetBinContent(ibin+1)+br*(function->Eval(etest)));
    }
  
  // exit
  fHdNdESignalChecked = kFALSE;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Read dN/dE for signal from file in the Segue Stereo input format
//
// Replacement of existing file is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t HdNdE::ReaddNdESignal(TString filename)
{
  // open file and look for histo
  TFile* dNdESignalFile  = new TFile(filename);
  TH1F*  hProvdNdESignal = (TH1F*) dNdESignalFile->Get("hdNdE");

  Int_t status = 0;
  if(!hProvdNdESignal)
    status = 1;
  else
    status = SetdNdESignal(hProvdNdESignal);

  // clean and exit
  dNdESignalFile->Close();
  delete dNdESignalFile;
  return status;
}

//////////////////////////////////////////////////////////////////
//
// Copy the contents of ih into oh by interpolating if necessary
// The x axis of ih and oh must be in log-scale and even binning
// The interpolation is done linearly in log(y)
// The units of the x-axis of ih are those of oh times 10^scale
// If isDiff=kTRUE the input histogram is differential (default)
// If isDiff=kFALSE the input histogram is bin-integrated
// The output histogram is ALWAYS differential
//
Int_t HdNdE::ReadAndInterpolate(TH1F* ih,TH1F* oh,Double_t scale,Bool_t isDiff)
{
  // check for minimum requirements
  if(!ih || !oh)
    {
      cout << " HdNdE::ReadAndInterpolate Error: either the input or output directories are NULL pointers" << endl;
      return 1;
    }
  
  // input histogram binning
  Double_t imine   = ih->GetXaxis()->GetXmin()+scale; // minimum log(E) in input histo
  Double_t imaxe   = ih->GetXaxis()->GetXmax()+scale; // maximum log(E) in input histo
  Int_t    inbinse = ih->GetNbinsX();          
  Double_t ide     = (imaxe-imine)/inbinse;

  // output histogram binning
  Double_t omine   = oh->GetXaxis()->GetXmin(); // minimum log(E) in output histo
  Double_t omaxe   = oh->GetXaxis()->GetXmax(); // maximum log(E) in output histo
  Int_t    onbinse = oh->GetNbinsX();          
  Double_t ode     = (omaxe-omine)/onbinse;
  
  // interpolate values
  Double_t x0,x1,dx0i,dx1i,y0,y1,etest,etestbin,y,rtest,ycontent0,ycontent1;
  for(Int_t ibin=0;ibin<onbinse;ibin++)
    {
      etest = omine+ode*(ibin+gCenterBin); 

      // log-log interpolation
      if(etest<imine+gCenterBin*ide) // extrapolation of values below the minimum input energy
	{
	  oh->SetBinContent(ibin+1,0);
	  continue;
	}
      else if(etest>imaxe-(1-gCenterBin)*ide) // extrapolation of values above the maximum input energy
	{
	  oh->SetBinContent(ibin+1,0);
	  continue;
	}
      else // interpolation of values in the range of provided energies
	{
	  etestbin  = Int_t((etest-imine-gCenterBin*ide)/ide); // corresponding bin in ih histo

	  ycontent0 = ih->GetBinContent(etestbin+1);
	  ycontent1 = ih->GetBinContent(etestbin+2);
	  
	  if(ycontent0<=0 || ycontent1<=0)
	    {
	      oh->SetBinContent(ibin+1,0);
	      continue;
	    }
	  
	  x0 = imine+(etestbin+gCenterBin)*ide;
	  x1 = imine+(etestbin+1+gCenterBin)*ide;
	  
	  if(!isDiff)
	    {
	      dx0i = TMath::Power(10,imine+(etestbin+1)*ide)-TMath::Power(10,imine+etestbin*ide);
	      dx1i = TMath::Power(10,imine+(etestbin+2)*ide)-TMath::Power(10,imine+(etestbin+1)*ide);
	      y0 = TMath::Log10(ycontent0/dx0i);
	      y1 = TMath::Log10(ycontent1/dx1i);	  
	    }
	  else
	    {
	      y0 = TMath::Log10(ycontent0);
	      y1 = TMath::Log10(ycontent1);	  
	    }
	}
	  
      y     = y0+(etest-x0)*(y1-y0)/(x1-x0);
      rtest = TMath::Power(10,y);
      oh->SetBinContent(ibin+1,rtest);
    }
  return 0;
}
