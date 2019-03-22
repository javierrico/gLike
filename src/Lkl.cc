/* ======================================================================== *\
!
!   Author: Jelena Aleksic      06/2012 <mailto:jelena@ifae.es>
!   Author: Javier Rico         01/2015 <mailto:jrico@ifae.es>
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
// Lkl
//
// This is an abstract class (i.e. must be inherited by others to be
// used) to perform a profile likelihood maximization (in reality a
// profile -2logL minimization; the expressions "likelihood", "profile
// likelihood" and "-2logL" are used as synonimous of "profile -2logL"
// in this code).
//
// The likelihood function (to be defined in the daughter class) has
// one free parameter only, called G or g (which you can think of as
// the number of signal events in a given data sample, but can be
// anything). It can also contain as many nuisance parameters as
// wanted, which will be profiled in the -2logL minimization process.
// Depending on the dauther class and the nuisance parameter, nuisance
// parameters can be optionally fixed during the fit (e.g. tau in
// Iact1dUnbinnedLkl).
//
// The data member fUnitsOfG (default=1, can be changed using
// SetUnitsOfG) defines the units in which g is expressed. That means
// that the results will be provided for vs g*fUnitsOfG, instead of g
// (internally the fits will be done with g as free parameter,
// though). 
//
// In addition, fUnitsOfG is used internally by JointLkl to relate
// the values of g for each of the added Lkl objects (see JointLkl
// documentation) and therefore MUST be defined for all those Lkl
// objects. In addition, g_i*fUnitsOfG_i MUST be equal for all
// added Lkl objects, i.e. proportional to a physical parameter
// common to all added samples.
//
// The uncertainty on fUnitsOfG (fDUnitsOfG) can be set using
// SetDUnitsOfG(arg1,arg2), for which several types of pdf are
// available, namely:
// * fUnitsOfG          distributes as a Gaussian with sigma = fDUnitsOfG
// * 1/fUnitsOfG        distributes as a Gaussian with sigma = fDUnitsOfG
// * log10(fUnitsOfG)   distributes as a Gaussian with sigma = fDUnitsOfG
// * 1/log10(fUnitsOfG) distributes as a Gaussian with sigma = fDUnitsOfG
// arg1 is the value of fDUnitsOfG and arg2 determines the type of pdf
// (see Lkl.h)
//
// If several Lkl objects share a common uncertainty in units of G,
// do NOT set fDUnitsOfG for each of them, instead, add them all to a
// JointLkl, for which you set the common value of fDUnitsOfG.
// 
// Lkl daughter classes may define methods to compute fUnitsOfG for
// different physical processes of interest, see e.g.
// Iact1dUnbinnedLkl::SetDMAnnihilationUnitsForG.
//
// Lkl class creates, configures and calls a TMinuit object,
// computing the profile -2logL vs g close to its minimum, and
// provides acces to the results. The profile over nuisance parameters
// is done using the Migrad minimization algorithm. For the
// computation of -2logL vs G use:
// - ComputeLklVsG: computes the curve -2logL vs g for an optimized
//                  range of g around the minimum of the curve. If
//                  fDUnitsOfG>0, it profiles the obtained curve for
//                  fUnitsOfG. The final curve can be obtained with
//                  the method GetLklVsG, and the value of g crossing
//                  -2logLmin+fErrorDef with GetGLklMinErr or
//                  GetGForLkl(fErrorDef).This is the default
//                  preferred method for getting the correct limits to
//                  g, in the proper units, and with all nuisance
//                  profiles properly profiled.
//
// Any daughter class inheriting from Lkl MUST construct its
// Lkl part with Lkl(npars), where npars is the total number of 
// free+nuisance parameters.
// Also, it MUST override the following methods:
// - MakeChecks: performes the checks needed before the minimization is
//               called. If successful, it should return 0, 1 otherwise.
//               Any method changing some of the material to be checked
//               should also call the method SetChecked(kFALSE).
// - SetFunctionAndPars: where fMinuit->SetFCN is called, thereby passing the 
//               function that computes
//               the likelihood (in the format defined by TMinuit). 
//               Also it has to call the method SetParameters to pass
//               the name, start value and precision of the
//               free+nuisance parameters.
// - SetMinuitLink: where simply a link is made between the fMinuit 
//                  object and a static TMinuit copy, due to some 
//                  subtelties of TMinuit.
// (see, e.g. Iact1dUnbinnedLkl for a practical implementation of these methods).
//
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"
#include "TH1.h"
#include "TPRegexp.h"

#include "Lkl.h"

ClassImp(Lkl);

using namespace std;

// static constants
static const Double_t gDUofGMargin  = 0.25;   // relative number of extra points taken at each side of the scanned G range
static const Double_t gNSigma       = 2;      // number of sigmas to be covered (at least) by the likelihood parabolas


////////////////////////////////////////////////////////////////
//
// Default constructor, all daughter classes MUST call it with 
// the relevant number of parameter (free+nuisance)
//
Lkl::Lkl(Int_t npars,TString inputString,TString name,TString title) :
  TNamed(name,title), fMinuit(NULL), fNPars(npars),
  fParName(NULL), fParStart(NULL), fParDelta(NULL),
  fParVal(NULL),  fParErr(NULL),   fIsParFixed(NULL),
  fErrorDef(1), fGLklVsG(NULL), fLklMin(0), 
  fUnitsOfG(1), fDUnitsOfG(0), fDUofGType(none), fDUofGEst(0),
  fGIsPositive(kFALSE), fGShift(0),
  fIsChecked(kFALSE), fStatus(-1)
{
  if(InterpretInputString(inputString))
    cout << "Lkl::Lkl Warning: there were problems interpreting the input string" << endl;

}

//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
// Basically check if the string has a specification for the
// error in the units of G and the type of distribution to be assumed
//
// DlinG=<val>:    takes fDUnitsOfG as <val>, interpreted as the relative error of fUnitsOfG
// DinvlinG=<val>: takes fDUnitsOfG as <val>, interpreted as the relative error of 1/fUnitsOfG
// DlogG=<val>:    takes fDUnitsOfG as <val>, interpreted as the error of log10(fUnitsOfG)
// DinvlogG=<val>: takes fDUnitsOfG as <val>, interpreted as the error of 1/log10(fUnitsOfG)
// DlogJ=<val>:    same as DinvlogG
//
Int_t Lkl::InterpretInputString(TString inputString)
{
  TPMERegexp re("\\s+");
  TPMERegexp fldre("=");
  TString inputfileName = "";
  
  // split the inputString into the different fields, search for the Lkl options and save values specified in each of them
  UInt_t nfields = re.Split(inputString);
  for(UInt_t ifield=0;ifield<nfields;ifield++)
    {
      UInt_t noptfields = fldre.Split(re[ifield]);
      if(noptfields!=2)
	{
	  cout << "Lkl::InterpretInputString Warning: ignoring option " << re[ifield] << " for bad format (should be \"<optname>=<optval>\")" << endl;
	  continue;
	}
      TString optname = fldre[0];
      if(optname.CompareTo("DlinG",TString::kIgnoreCase)==0)
	{
	  fDUnitsOfG=fldre[1].Atof();
	  fDUofGType=Lkl::lin;
	}
      else if(optname.CompareTo("DinvlinG",TString::kIgnoreCase)==0)
	{
	  fDUnitsOfG=fldre[1].Atof();
	  fDUofGType=Lkl::invlin;
	}
      else if(optname.CompareTo("DlogG",TString::kIgnoreCase)==0)
	{
	  fDUnitsOfG=fldre[1].Atof();
	  fDUofGType=Lkl::log;
	}
      else if(optname.CompareTo("DinvlogG",TString::kIgnoreCase)==0 || optname.CompareTo("DlogJ",TString::kIgnoreCase)==0)
	{
	  fDUnitsOfG=fldre[1].Atof();
	  fDUofGType=Lkl::invlog;
	}
    }
  
  return 0;
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
Lkl::~Lkl()
{
  if(fMinuit)     delete fMinuit;
  if(fParName)    delete [] fParName;
  if(fParStart)   delete [] fParStart;
  if(fParDelta)   delete [] fParDelta;
  if(fParVal)     delete [] fParVal;
  if(fParErr)     delete [] fParErr;
  if(fIsParFixed) delete [] fIsParFixed;
  if(fGLklVsG)    delete fGLklVsG;  
}


////////////////////////////////////////////////////////////////
// 
// Compute -2logL vs g around the minimum and store the results
// in fGLklVsG (TGraph*)
//
// The curve is computed in <npoints> points, in a range from
// gmin-gNSigma*Dg to gmin+gNSigma*Dg, where gmin is the value of g
// minimizing -2logL, and Dg is the error in g as computed by migrad,
// with the confidence interval as set by SetErroDef.
// These limits are expanded in case we do not reach the value
// -2logL_min+fErrorDef, if fDUnitsOfG=0 or
// -2logL_min+fErrorDef*(1+gExpandFactor*fDUnitsOfG), otherwise
//
// In case fDUnitsOfG>0, the method ApplyDUnitsOfG is called to profile
// the curve for the effect of the uncertainty in fUnitsOfG
//
// Before saving fGLklVsG, all values in the curve are subtracted the
// minimum in such a way that the value at gmin is 0
//
// if centerAtZero=kTRUE (default kFALSE) the curve minimum is set for g=0
//
// Return the lkl minimum 
//
Double_t Lkl::ComputeLklVsG(Bool_t centerAtZero,Int_t npoints,Double_t glow,Double_t gupp,Bool_t isVerbose)
{
  // delete previously existing curve
  ResetGLklVsG();
  
  Bool_t definedRange = ((glow!=0 || gupp!=0) && gupp>glow);
  Bool_t applyDUofG   = (fDUnitsOfG>0 && fDUofGType!=none);
  Double_t expCoeff   = 1;
  if(!definedRange)
    {
      // first find where the minimum is and estimate the error
      if(isVerbose)
	cout << "Lkl::ComputeLklVsG (" << GetName() << ") Message: Finding minimum of -2logL... " << endl;
      Lkl::MinimizeLkl();
      FindGLowAndGUpp(glow,gupp,centerAtZero);
    }

  npoints = 200;
  Double_t sigmaVlow = 1.e-28;
  Double_t sigmaVupp = 2.e-22;

  Double_t unitsOfG = GetUnitsOfG();
  Double_t glow_fixrange = sigmaVlow/unitsOfG;
  Double_t gupp_fixrange = sigmaVupp/unitsOfG;
  glow = glow_fixrange;
  gupp = gupp_fixrange;
 
  // expansion of glow and gupp from daughter likelihoods
  expCoeff = GetExpansionCoefficient();
  // expansion of glow and gupp from uncertainties in units of G
  ExpandGLowAndGUpp(glow,gupp,expCoeff);

  const Int_t nmaxcounts=10;
  Int_t ncounts=0;
  while(1) // break when Lkl(gupp)>Lkl_min+fErrorDef or ncounts>=nmaxcounts
    {
      ncounts++;
      
      // delete previously existing curve
      ResetGLklVsG();

      // if there are uncertainties in units of G, we need to compute beyond the limits, 
      // to compute the profile within the limits
      Int_t    realNPoints = npoints;
      Double_t newglow     = glow;
      Double_t newgupp     = gupp;

      if(applyDUofG)
	{
	  realNPoints = npoints+2*Int_t(npoints*gDUofGMargin);
	  newglow     = glow-(gupp-glow)*gDUofGMargin;
	  newgupp     = gupp+(gupp-glow)*gDUofGMargin;
	  if(fGIsPositive && newglow<0)
	    { 
	      Int_t removepoints = Int_t(-newglow*realNPoints/(newgupp-newglow)+.1);
	      realNPoints -= removepoints;
	      newglow      = 0;
	    }
	}
      
      // compute the value of -2logL for g between glow and gupp
      Double_t shiftg=0;
      
      // maybe there is something to be done first
      if(PrepareForLklScan(centerAtZero,npoints,newglow,newgupp,isVerbose))
	{
	  cout << "Lkl::ComputeLklVsG (" << GetName() << ") Warning: PrepareForLklScan failed" << endl;
	  return 0;
	}

      // write the range and number of points we will explore
      if(isVerbose)	
	cout << "Lkl::ComputeLklVsG (" << GetName() << ") Message: computing -2logL in " << realNPoints << " points between g=" << newglow << "(" << newglow*fUnitsOfG << "), and g=" << newgupp << "(" << newgupp*fUnitsOfG << "), this could take a while" << endl;
     
      // create the arrays to hold the parabola x and y values
      Double_t* xval = new Double_t[realNPoints];
      Double_t* yval = new Double_t[realNPoints];
  
      // and now compute the value of the full lkl for the g range
      Double_t gstep  = (newgupp-newglow)/(realNPoints-1.);
      Int_t    point0 = TMath::FloorNint(-newglow/gstep);

      for(Int_t ipoint=0;ipoint<realNPoints;ipoint++)
	{
	  // "progress bar"
	  if(isVerbose)
	    if((ipoint+1)%10==0)	
	      cout << "." << flush;

	  xval[ipoint] = (ipoint==point0? 0 :newglow+gstep*ipoint);      
	  yval[ipoint] = MinimizeLkl(xval[ipoint],kTRUE,kFALSE);

	}

      TString className = GetName();
      if(className.CompareTo("JointLkl_00")==0 && glow==glow_fixrange && gupp==gupp_fixrange)
        {
          ofstream myfile;
          myfile.open("./plots/Gloryduck_example_file.txt", std::ios_base::app);
          /*
          for(Int_t ipoint=0;ipoint<realNPoints;ipoint++)
            {
              myfile << xval[ipoint]*unitsOfG << " ";
            }
          myfile << "\n";
          */
          for(Int_t ipoint=0;ipoint<realNPoints;ipoint++)
            {
              myfile << yval[ipoint] << " ";
            }
          myfile << "\n";
          myfile.close();
        }
      if(isVerbose)
	cout << " Completed " << realNPoints << " points" << endl;
      
      // save it as TGraph
      fGLklVsG = new TGraph(realNPoints,xval,yval);

      // delete unneeded arrays
      delete [] xval;
      delete [] yval;  

      Double_t savegupp = gupp;
      Double_t saveglow  = glow;

      // move so that lklmin is at g=0 if requested
      if(centerAtZero && applyDUofG)
	//OJO
	{
	  Bool_t repeat = kFALSE;
	  
	  shiftg = CenterAtZero();
	  if(shiftg>newgupp-gupp)
	    {
	      repeat = kTRUE;
	      gupp+=shiftg;
	    }
	  if(shiftg<newglow-glow)
	    {
	      repeat = kTRUE;
	      glow+=shiftg;
	    }
	  if(repeat)
	    {
	      cout << "Lkl::ComputeLklVsG (" << GetName()
		   << ") Message: insufficient range for profiling, trying with glow,gupp = "
		   << glow << "," << gupp << " (before it was " << saveglow <<","<< savegupp << ")" <<  endl;
	      continue;
	    }
	}
      
      // profile for the error in units of G
      if(applyDUofG)
	ApplyDUnitsOfG(npoints,glow,gupp);

      // Put the minimum of LklVsG to 0, find point crossing fErrorDef...
      if(!definedRange)
	FixLklVsG(shiftg);

      // check if convergence condition is met (or maximum number of trials) and leave
      Double_t lklmax = fGLklVsG->GetY()[fGLklVsG->GetN()-1];
      if(definedRange || (lklmax>fErrorDef && lklmax<4*fErrorDef) || ncounts>nmaxcounts) break;

      // otherwise, modify the range of g
      if(lklmax<fErrorDef)
	gupp*=(gupp>0? 2 : -1);
      if(lklmax>4*fErrorDef)
	{
	  Double_t estval = GetGForLkl(fErrorDef,kFALSE);
	  if(estval>0)
	    gupp=estval*1.2;
	  else
	    gupp=estval*0.8;
	}
      
      Double_t estgmin = GetGForLkl(0,kFALSE);
            
      if(!fGIsPositive)
	{
	  if(glow==estgmin)
	    glow = -TMath::Abs(gupp-2*(gupp-estgmin));
	  else
	    glow = gupp-2*(gupp-estgmin);
	}
	  
      if(isVerbose)
	cout << "Lkl::ComputeLklVsG (" << GetName() << ") Message: -2logL(gupp) = " << lklmax 
	     << ", trying with glow,gupp = " << glow << "," << gupp << " (before it was " << saveglow
	     <<","<< savegupp << " and estimated g_min = " << estgmin << ")" << endl;	  
    }

  if(ncounts>nmaxcounts)
    cout << "Lkl::ComputeLklVsG (" << GetName() << ") Warning: maximum number of trials to find optimal glow-gupp range exeeded, results will be probably wrong/inacurate, check -2logL curves!" << endl;
  
  // Print results
  if(!definedRange && isVerbose)
    cout << "Lkl::ComputeLklVsG (" << GetName() << ") Message:"
	 << " g_min = " << fParVal[gGParIndex] << " +/- " << fParErr[gGParIndex]
	 << " (" << fParVal[gGParIndex]*fUnitsOfG << " +/- " << fParErr[gGParIndex]*fUnitsOfG
	 << "), -2logLmin = "<< fLklMin << endl;
  
  return fLklMin;
}

////////////////////////////////////////////////////////////////
// 
// Assuming -2logL has a parabolic shape, find the range of g
// values for which the lkl curve needs be computed so that the
// region close to the minimum of the -2logL curve is scanned,
// i.e. the curve is above fLklMin+fErrorDef for the minimum and
// maximum g values, and it contains the absolute minimum
//
void Lkl::FindGLowAndGUpp(Double_t& glow,Double_t& gupp,Bool_t centerAtZero)
{
  // guess the range of g to be explored assuming Migrad found a parabola right
  Double_t gmin    = fParVal[gGParIndex];
  Double_t gdlt    = gmin+fParErr[gGParIndex];
  Double_t fgmin   = fLklMin;
  Double_t fgdlt   = fLklMin+fErrorDef;
  Double_t targetL = fLklMin+gNSigma*fErrorDef;

  // Parabola passing by (gmin,fgmin) and (gdlt,fgdlt) with minimum at gmin
  Double_t a = (fgdlt-fgmin)/TMath::Power(gdlt-gmin,2);
  Double_t b = -2*a*gmin;
  Double_t c = fgmin+a*gmin*gmin;
  
  // solve a*g*g+b*g+c=targetL
  if(b*b-4*a*(c-targetL)>=0)
    {
      glow             = (-b-TMath::Sqrt(b*b-4*a*(c-targetL)))/(2.*a);
      gupp             = (-b+TMath::Sqrt(b*b-4*a*(c-targetL)))/(2.*a);
    }
  else
    {
      glow = gmin-1;
      gupp = gmin+1;
    }
  
  // if the range was too narrow, expand it
  Double_t lklval;

  // contract the lower end
  while((lklval=MinimizeLkl(glow,kTRUE,kFALSE,kTRUE))>fgmin+4*fErrorDef)
    {
      Double_t save  = glow;
      glow+=TMath::Abs(glow-gmin)/2.;
      cout << "Lkl::FindGLowAndGUpp (" << GetName() << ") Message: -2logL(glow) = " << lklval << ", glow rised from " << save << " to " << glow <<" (fLklMin = " << fgmin << ")"<< endl;
    }
  
  // expand the upper end
  while((lklval=MinimizeLkl(glow,kTRUE,kFALSE,kTRUE))<fgmin+fErrorDef || (centerAtZero && glow>0))
    {
      Double_t save  = glow;      
      glow-=(glow>0? 2*TMath::Abs(glow) : TMath::Abs(glow));
      cout << "Lkl::FindGLowAndGUpp (" << GetName() << ") Message: -2logL(gupp) = " << lklval << ", glow lowered from " << save << " to " << glow <<" (fLklMin = " << fgmin << ")"<< endl;
    }
  
  // if only positive values set glow to zero
  if(fGIsPositive && glow<0)
    {
      Double_t save  = glow;
      glow = 0;
      cout << "Lkl::FindGLowAndGUpp (" << GetName() << ") Message: glow value changed from " << save << " to " << glow << " (negative values excluded)" << endl;
      
      if(fParVal[gGParIndex]<0)
	{
	  fgmin = fLklMin = MinimizeLkl(0,kTRUE,kFALSE,kTRUE);
	  gmin    = 0;
	}
    }

  // contract the upper end
  while((lklval=MinimizeLkl(gupp,kTRUE,kFALSE,kTRUE))>fgmin+4*fErrorDef)
    {
      Double_t save  = gupp;      
      gupp-=TMath::Abs(gupp-gmin)/2.;
      cout << "Lkl::FindGLowAndGUpp (" << GetName() << ") Message: -2logL(gupp) = " << lklval << ", gupp lowered from " << save << " to " << gupp <<" (fLklMin = " << fgmin << ")"<< endl;
    }
  
  // expand the upper end
  while((lklval=MinimizeLkl(gupp,kTRUE,kFALSE,kTRUE))<fgmin+fErrorDef || (centerAtZero && gupp<0))
    {
      Double_t save  = gupp;      
      gupp+=(gupp<0? 2*TMath::Abs(gupp) : TMath::Abs(gupp));
      cout << "Lkl::FindGLowAndGUpp (" << GetName() << ") Message: -2logL(gupp) = " << lklval << ", gupp rised from " << save << " to " << gupp <<" (fLklMin = " << fgmin << ")"<< endl;
    }
  fLklMin=fgmin;
}

////////////////////////////////////////////////////////////////
// 
// Expand g range by factor <stretch> to cope with profile
// over UnitsOfG
//
void Lkl::ExpandGLowAndGUpp(Double_t& glow,Double_t& gupp,Double_t stretch)
{
  Double_t saveupp = gupp;
  Double_t savelow = glow;
  if(!fGIsPositive)
    glow-=TMath::Max(TMath::Abs(glow),TMath::Abs(gupp))*(stretch-1);
  gupp+=TMath::Max(TMath::Abs(glow),TMath::Abs(gupp))*(stretch-1);
  
  if(savelow!=glow)
    cout << "Lkl::ExpandGLowAndGUpp (" << GetName() << ") Message: glow value changed from " << savelow << " to " << glow << " to account for uncertainty in Units of G" << endl;
  
  if(saveupp!=gupp)
    cout << "Lkl::ExpandGLowAndGUpp (" << GetName() << ") Message: gupp value changed from " << saveupp << " to " << gupp << " to account for uncertainty in Units of G" << endl;
}
		
////////////////////////////////////////////////////////////////
//
// Call Migrad to find the minimum and width of -2logL, unless:
// If -2logL curve has already been computed, return its value at g
// (if gIsFixed) or at the minimum (if not fixed)
// If there are uncertainties in fUnitsOf G, and force=kFALSE, compute
// the profile -2logL curve. This is necessary because to profile the
// curve you need to previously know its value for all g.
// Other options:
//
// If gIsFixed  = kTRUE (default=kFALSE), g is fixed during the fit
//                (i.e. only the nuisance parameters are let free)
// If isVerbose = kTRUE (default=kTRUE), the values of the parameters
//                after the fit are reported
//
// return the value of -2logL (at g or minimum) or 0 if it fails.
//
Double_t Lkl::MinimizeLkl(Double_t g,Bool_t gIsFixed,Bool_t isVerbose,Bool_t force)
{
  // if parabola has already been computed, give back requested lkl value
  if(IsChecked() && fGLklVsG && !force) 
    {
      if(gIsFixed)
	{
	  fParVal[gGParIndex] = g;
	  FixPar(gGParIndex,gIsFixed);
	  // OJO the 1.01 factor is there to allow for linear interpolation close to the upper bound, this avoids some ununderstood weird shape of the parabole due probably to rounding problem	  
	  if(g>=fGLklVsG->GetX()[0] && g<=fGLklVsG->GetX()[fGLklVsG->GetN()-1]*1.01)
	    return fGLklVsG->Eval(g);
	}
      else return GetLklMin();
    }

  // configure input, TMinuit, input parameters etc, if needed
  if(MakeChecks())
    {
      cout << "Lkl::MinimizeLkl (" << GetName() << ") Warning: checks were not successfull!!!!" << endl;
      return 0;
    }
     
  // initialize minuit
  InitMinuit();
  SetMinuitLink();

  // assign value to g and fix it if requested
  fMinuit->DefineParameter(gGParIndex,fParName[gGParIndex],g+fGShift,fParDelta[gGParIndex],0,0);

  fMinuit->Release(gGParIndex);
  FixPar(gGParIndex,kFALSE);
  if(gIsFixed)
    {
      fMinuit->FixParameter(gGParIndex);
      FixPar(gGParIndex);
      if(GetNFreePars()==0)
	{
	  fParVal[gGParIndex] = g;
	  fStatus = 0;
	  return GetLklVal();
	}
    }

  // write results if requested
  if(isVerbose)
    cout << "Lkl::MinimizeLkl (" << GetName() << ") Message: minimizing -2logL"
	 << (gIsFixed? Form(" for fixed g=%.4e (%.4e)",g,g*fUnitsOfG) : "") 
	 << endl;

  // call minimization (with requested verbosity)
  Int_t strategy = (!gIsFixed)*2; // when it's fixed, we can relax the strategy since there is no interest in errors
  Double_t lkl = CallMinimization(g,isVerbose,strategy);

  // keep value of minimum lkl and associated g (with error)  
  for(Int_t ipar=0;ipar<fNPars;ipar++)
    {
      Double_t parval,parerr; 
      fMinuit->GetParameter(ipar,parval,parerr);
      fParVal[ipar] = parval;
      fParErr[ipar] = parerr;
    }
  fLklMin     = lkl;

  return lkl;
}

////////////////////////////////////////////////////////////////
//
// Initialize TMinuit object 
//
void Lkl::InitMinuit(Double_t ginit)
{
  // configure minuit
  if(fMinuit) delete fMinuit;

  fMinuit = new TMinuit(fNPars);
  fMinuit->SetPrintLevel(-1);
  fMinuit->SetErrorDef(fErrorDef);   // number of std deviations of error
  fMinuit->SetObjectFit(this);

  SetFunctionAndPars(ginit);
}

////////////////////////////////////////////////////////////////
// 
// Call minimization of fMinuit with iterations varying the step on g
// to increase the chances of convergence
//
// g is the initial or fixed value of g
// If isVerbose=kTRUE, dump the resulting values for free+nuisance parameters
//
// Returns the value of -2logL 
//
Double_t Lkl::CallMinimization(Double_t g,Bool_t isVerbose,Int_t strategy)
{
  Double_t *arglist = new Double_t [10];
  Int_t iflag=-1;

  // Compute limits to <sv> using migrad errors	(g is a free parameter)
  arglist[0] = strategy;
  fMinuit->mnexcm("SET STR", arglist, 1, iflag);
  arglist[0] = 10000;
  iflag = -1;
  Int_t counter = 0;
  const Int_t maxcounts = 10;

  // try until convergence is achieved
  while((iflag!=0 || TMath::IsNaN(GetParErr(0))) && counter<maxcounts) // try until the fit converges
    {
      fMinuit->mnexcm("MIGrad", arglist, 1, iflag);
      fStatus = iflag;
      
      // Under request, display values of parameters
      if(isVerbose)
	{
	  cout << "Lkl::CallMinimization (" << GetName() << ") Results: Trial #" << counter+1 <<", ";
	  Double_t parval,parerr; 		
	  for(Int_t ipar=0;ipar<fNPars;ipar++)
	    {
	      fMinuit->GetParameter(ipar,parval,parerr);
	      cout << fParName[ipar] << ": " << Form("%.2e +/- %.2e",parval,parerr)
		   << (ipar==gGParIndex? Form(" (%.2e +/- %.2e)",parval*fUnitsOfG,parerr*fUnitsOfG) :"") << "; ";
	    }
	  cout << "-2logL = " << GetLklVal() << "; iflag = " << iflag << (iflag==0? " (converged)" : " (not converged)") << endl;
	}

      // try changing precision
      if(iflag!=0 || TMath::IsNaN(GetParErr(gGParIndex)))
	fMinuit->DefineParameter(gGParIndex,fParName[gGParIndex],g,fParDelta[gGParIndex]/10.*TMath::Power(1.8,counter), 0, 0);
      
      counter++;
      if(counter==maxcounts)
	cout << "Lkl::CallMinimization (" << GetName() << ") Warning: No convergence reached for g = " << g << " after " << maxcounts << " trials: check your -2logL curves for features" << endl;
    }

  delete [] arglist;

  return (counter<maxcounts? GetLklVal() : 9e99);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Apply the error on the units of G (or 1/UofG or log10(UofG) or
// 1/log10(UofG), depending on the value of the variable
// fDUofGType:
//
// fDUofGType = none   -> fUnitsOfG has no error, similar to fDUnitsOfG=0 [default]
// fDUofGType = lin    -> fDUnitsOfG interpreted as the relative error of fUnitsOfG
// fDUofGType = invlin -> fDUnistOfG interpreted as the relative error of 1/fUnitsOfG
// fDUofGType = log    -> fDUnitsOfG interpreted as the error of log10(fUnitsOfG)
// fDUofGType = invlog -> fDUnitsOfG interpreted as the error of 1/log10(fUnitsOfG)
//
// In all the cases it is assumed that the uncertainty on UofG is
// uncorrelated with the uncertainty on G, and the pdf for fDUnitsOfG
// is a Gaussian on the parameter specified by fDUofGType
//
// OJO: only case invlog fully tested!
//
// Return 1 if problem is found, 0 otherwise
//
Int_t Lkl::ApplyDUnitsOfG(Int_t npoints,Double_t glow,Double_t gupp)
{
  if(fDUofGType==none) return 0;

  // get -2logL vs G
  TGraph* m2logLvsG = GetLklVsG(kFALSE);
  if(!m2logLvsG)
    {
      cout << "Lkl::ApplyDUnitsOfG (" << GetName() << ") Warning: curve -2logL vs g has not been computed" << endl;
      return 1;
    }
  cout << "Lkl::ApplyDUnitsOfG (" << GetName() << ") Message: profiling -2logL curve over uncertainty in units of G (" << fDUnitsOfG;
  if(fDUofGType==lin)
    cout << ", linear mode)";
  else if(fDUofGType==invlin)
    cout << ", inverse linear mode)";
  else if(fDUofGType==log)
    cout << ", log10 mode)";
  else if(fDUofGType==invlog)
    cout << ", inverse log10 mode)";
  cout << "... " << flush;
  
  // compute graph with Lkl vs G
  const Int_t nvalsI    = m2logLvsG->GetN();
  Double_t gval[npoints];
  Double_t gmin= m2logLvsG->GetX()[0];
  Double_t gmax= m2logLvsG->GetX()[nvalsI-1];
  for(Int_t ivalF=0;ivalF<npoints;ivalF++)
    gval[ivalF]   = glow+(gupp-glow)/(npoints-1)*ivalF;

  // compute the Lkl for g*fUnitsOfG
  Double_t nsigma = 6;
  Int_t    nlvals = 10000; // number of steps in the profiling 
  Double_t profileLklval[npoints];
  Double_t j,l;
  Int_t    nout   = 0; // number of times I need a value I cannot find
  for(Int_t ivalF=0;ivalF<npoints;ivalF++)
    {
      profileLklval[ivalF]=5e99;
      Double_t s =  gval[ivalF];
      Int_t lminval=-1;
      for(Int_t lval=0;lval<nlvals;lval++)
	{
	  Double_t g=0;
	  Double_t lLkl=0;	  
	  Double_t maxdev = TMath::Abs(nsigma*fDUnitsOfG);
	  if(fDUofGType==lin)
	    {
	      Double_t dl   = 2*TMath::Abs(maxdev)/(nlvals-1);
	      l             = 1-maxdev+lval*dl;

	      lLkl = -2*TMath::Log(TMath::Gaus(l,1,fDUnitsOfG,kTRUE));
	      g = s*l;
	    }
	  else if(fDUofGType==invlin)
	    {
	      Double_t dl   = 2*TMath::Abs(maxdev)/(nlvals-1);
	      l             = 1-maxdev+lval*dl;

	      lLkl = -2*TMath::Log(TMath::Gaus(l,1,fDUnitsOfG,kTRUE)*l*l);
	      g = s/l;
	    }
	  else if(fDUofGType==log)
	    {
	      Double_t dl   = 2*maxdev/(nlvals-1);
	      l             = -maxdev+lval*dl;
	      j             = TMath::Power(10,l);

	      lLkl = -2*TMath::Log(TMath::Gaus(l,0,fDUnitsOfG,kTRUE)/j/TMath::Log(10.));
	      g = s*j;
	    }
	  else if(fDUofGType==invlog)
	    {
	      Double_t dl   = 2*maxdev/(nlvals-1);
	      l             = -maxdev+lval*dl;
	      j             = TMath::Power(10,l);
	      lLkl = -2*TMath::Log(TMath::Gaus(l,0,fDUnitsOfG,kTRUE)/TMath::Log(10.));     // new Fermi prescription (correct)
	      g = s/j;
	    }
	  
	  if(g<gmin || g>gmax) continue;
	  Double_t gLkl   = m2logLvsG->Eval(g);
	  Double_t totLkl = lLkl+gLkl;
	  if(totLkl<profileLklval[ivalF])
	    {
	      lminval = lval;
	      profileLklval[ivalF] = totLkl;
	      switch(fDUofGType)
	      {
	      case lin:    fDUofGEst = l-1; break;
	      case invlin: fDUofGEst = l-1; break;
	      case log:    fDUofGEst = l; break;
	      case invlog: fDUofGEst = l; break;		
	      default:     fDUofGEst = 1;
	      }
	    }
	}
      if(lminval<=0 || lminval>=nlvals-1) nout++;
    }

  if(nout>0)
    cout << endl << "Lkl::ApplyDUnitsOfG (" << GetName() << ") Warning: " << nout << " out-or-range values were requested" << endl;
      
  
  // save it as TGraph
  if(fGLklVsG) delete fGLklVsG;  
  fGLklVsG = new TGraph(npoints,gval,profileLklval);
  
  delete m2logLvsG;
  cout << "Done" << endl;

  return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Set the names, initial values and precisions of the free+nuisance
// parameters. This function is thought to be called by
// SetFunctionAndPars, in the daughter class
//
void Lkl::SetParameters(const Char_t** parname, Double_t* parstart, Double_t* pardelta) 
{
  if(fParName)    delete [] fParName;
  if(fParStart)   delete [] fParStart;
  if(fParDelta)   delete [] fParDelta;
  if(fParVal)     delete [] fParVal;
  if(fParErr)     delete [] fParErr;
  if(fIsParFixed) delete [] fIsParFixed;

  fParName    = new TString[fNPars];
  fParStart   = new Double_t[fNPars];
  fParDelta   = new Double_t[fNPars];
  fParVal     = new Double_t[fNPars];
  fParErr     = new Double_t[fNPars];
  fIsParFixed = new Bool_t[fNPars];

  for(Int_t ipar = 0; ipar<fNPars; ipar++)
    {
      fParName[ipar]    = parname[ipar];
      fParStart[ipar]   = parstart[ipar];
      fParDelta[ipar]   = pardelta[ipar];
      fParVal[ipar]     = 0;
      fParErr[ipar]     = 0;
      fIsParFixed[ipar] = kFALSE;
      if(fMinuit)
	fMinuit->DefineParameter(ipar, fParName[ipar], fParStart[ipar], fParDelta[ipar], 0, 0);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Move the x-axis so that the minimum is at g=0 
// Usefull when computing sensitivity curves with units as nuisance
//
Double_t Lkl::CenterAtZero()
{
  // compute minimum
  Double_t minLogL = 9e99;
  Double_t gminLogL;
  for(Int_t ipoint=0;ipoint<fGLklVsG->GetN();ipoint++)
    if(fGLklVsG->GetY()[ipoint]<minLogL){minLogL = fGLklVsG->GetY()[ipoint];gminLogL=fGLklVsG->GetX()[ipoint];}

  // g of minimum Lkl
  for(Int_t ipoint=0;ipoint<fGLklVsG->GetN();ipoint++)
    fGLklVsG->GetX()[ipoint]-=gminLogL;

  cout << "Lkl::CenterAtZero (" << GetName() << ") Message: g values shifted by "
       << gminLogL << Form("(%.2e)",gminLogL*fUnitsOfG) << endl;

  fGShift = gminLogL;
  
  return gminLogL;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Find the minimum of GLklVsG, keep its value in fLklMin, the
// associated values of the free parameters, and subtract fLklMin to
// the whole curve. Find the point of curve crossing fErrorDef and
// keep it in as fParErr[gGParIndex].
//
void Lkl::FixLklVsG(Double_t shiftg)
{
  // compute minimum
  Double_t minLogL = 9e99;
  Double_t gminLogL;
  for(Int_t ipoint=0;ipoint<fGLklVsG->GetN();ipoint++)
    if(fGLklVsG->GetY()[ipoint]<minLogL){minLogL = fGLklVsG->GetY()[ipoint];gminLogL=fGLklVsG->GetX()[ipoint];}

  // subtract minimum Lkl
  for(Int_t ipoint=0;ipoint<fGLklVsG->GetN();ipoint++)
    fGLklVsG->GetY()[ipoint]-=minLogL;
  
  // compute and save the values of nuisance parameter for g at the minimum
  fLklMin = MinimizeLkl(gminLogL+shiftg,kTRUE,kFALSE,kTRUE);
  FixPar(gGParIndex,kFALSE);
  fParVal[0] = gminLogL;
  fParErr[0] = GetGForLkl(fErrorDef,kFALSE)-fParVal[0];
  if(fNPars>1)
    {
      Double_t parval,parerr; 
      for(Int_t ipar=1;ipar<fNPars;ipar++)
	{
	  fMinuit->GetParameter(ipar,parval,parerr);
	  fParVal[ipar] = parval;
	  fParErr[ipar] = parerr;
	}
    }
  SpreadFixLklVsG(gminLogL+shiftg);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Creates and saves fGLklVsG using input points
//
void Lkl::SetGLklVsG(Int_t npoints,Double_t* x,Double_t* y)
{
  if(fGLklVsG) delete fGLklVsG;
  fGLklVsG = new TGraph(npoints,x,y);
  SetChecked(kFALSE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Return the value of g for which the -2logL is the minimum (0) plus
// lkl Normally there are two (or more) solutions, take the larger one
// within the computed range.
// If units==kTRUE, return g*fUnitsOfG, g otherwise
// 
// Returns 0 if the -2logL has not been computed yet
//
Double_t Lkl::GetGForLkl(Double_t lkl,Bool_t units) const
{
  if(!fGLklVsG) return 0;

  const Int_t npts = fGLklVsG->GetN();
  Double_t xlow=0,xup=0,ylow=0,yup=0;
  Double_t factor = (units? fUnitsOfG : 1);

  // treat case when there is only one data point in the curve
  if(npts==0)
    {
      Double_t noval = fGLklVsG->GetX()[0];
      cout << "Lkl::GetGForLkl (" << GetName() << ") Warning: only one point in -2logL curve, returning it = " << noval*factor << " -2logL = " << fGLklVsG->GetY()[0] << endl;
      return noval*factor;
    }

  // if the right-most point is already below the searched value, return it
  if(fGLklVsG->GetY()[npts-1]<lkl)
    {
      Double_t noval = fGLklVsG->GetX()[npts-1];
      cout << "Lkl::GetGForLkl (" << GetName() << ") Warning: highest g below searched -2logL value, returning highest g = " << noval*factor << " -2logL = " << fGLklVsG->GetY()[npts-1] << endl;
      return noval*factor;
    } 

  // find where the searched for value is and interpolate from closest values
  for(Int_t ipt = 0;ipt<npts;ipt++)
    if(fGLklVsG->GetY()[npts-1-ipt]>=lkl && fGLklVsG->GetY()[npts-2-ipt]<=lkl)
      {
	xlow = fGLklVsG->GetX()[npts-2-ipt];
	xup  = fGLklVsG->GetX()[npts-1-ipt];
	ylow = fGLklVsG->GetY()[npts-2-ipt];
	yup  = fGLklVsG->GetY()[npts-1-ipt];
	break;
      }

  // compute interpolation of closest values and return it
  if((yup-ylow)!=0)
    return (xlow+(lkl-ylow)*(xup-xlow)/(yup-ylow))*factor;

  // something went wrong, at least warn about it
  Double_t noval = fGLklVsG->GetX()[npts-1];
  cout << "Lkl::GetGForLkl (" << GetName() << ") Warning: no good cutting value, returning highest g = " <<  noval*factor  << " -2logL = " << fGLklVsG->GetY()[npts-1] << endl;
  return noval*factor;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Compute the value of the -2logL function for the current values of parameters
// Returns 0 if no fit has been performed yet
//
Double_t Lkl::GetLklVal() const
{
  if(!fMinuit) return 0;

  Double_t parval[fNPars];   // To retrieve fitted values
  Double_t parerr[fNPars];   // To retrieve fit errors
  Double_t grad[fNPars];
  Double_t lklval;

  for(Int_t ipar=0;ipar<fNPars;ipar++)
    fMinuit->GetParameter(ipar, parval[ipar], parerr[ipar]);
  fMinuit->Eval(fNPars,grad,lklval,parval, 3);	

  return lklval;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Print information about the lkl object and the fit results
//
void Lkl::PrintOverview(Int_t level)
{
  // Lkl object name
  Margin(level); cout << "  Name              = " << GetName() << endl;

  // Convergence status of last Migrad call
  Margin(level); cout << "  Status            = " << fStatus << " (";
  if(fStatus==0)
    cout << "converged";
  else if(fStatus>0)
    cout << "not converged";
  else
    cout << "undefined";
  cout << ")" << endl;

  // Defined error level
  Margin(level); cout << "  Delta(2logL)      = " << fErrorDef << (fGIsPositive? " (G is positive only)" : "") << endl;

  // free and nuisance parameter values with errors and units for the case of G
  Margin(level); cout << "  # of parameters   = " << fNPars << " (" << GetNFreePars() << " free):" << endl;
  for(Int_t ipar=0;ipar<fNPars;ipar++)
    {
      Margin(level); cout << Form("      %13s = ",fParName[ipar].Data());
      if(fIsParFixed[ipar])
	{
	  cout << fParVal[ipar];
	  if(ipar==gGParIndex)
	    cout  << " (" << fParVal[ipar]*fUnitsOfG << ")";
	}
      else
	{
	  cout << fParVal[ipar] <<" +/- " << fParErr[ipar];
	  if(ipar==gGParIndex)
	    cout  << " (" << fParVal[ipar]*fUnitsOfG << " +/- " << fParErr[ipar]*fUnitsOfG << ")";
	}
      cout << (fIsParFixed[ipar]? " (fixed)" : "") 
	   << (fGIsPositive?      " (positive only)" : "") 
	   << endl;
    }

  // units of G and uncertainty
  Margin(level); cout << "  Units of G        = " << fUnitsOfG << endl;
  if(fDUnitsOfG>0)
    {  
      Margin(level); cout << "  Relative D(units) = " << fDUnitsOfG << " (";
      switch(fDUofGType)
	{
	case none   : cout << "None";                break;
	case lin    : cout << "Linear";              break;
	case invlin : cout << "Inverse-linear";      break;
	case log    : cout << "Logarithmic";         break;
	case invlog : cout << "Inverse-logarithmic"; break;
	default     : cout << "Undefined";           break;
	}
      cout << " PDF)" << endl;
      Margin(level); cout << "  Est      D(units) = " << fDUofGEst << endl;
    }

  // lkl value
  Margin(level); cout << "  -2logL_min        = " << fLklMin << endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Create and returns a TGraph with -2logL vs 
// G*fUnitsOfG if applyUnits=kTRUE (default), or
// G if applyUnits=kFALSE
//
// It is the responsability of the code calling this function to
// delete the created TGraph object.
//
TGraph* Lkl::GetLklVsG(Bool_t applyUnits) const
{
  // basic check
  if(!fGLklVsG) return NULL;

  // copy fGlklVsG with x-axis in units of fUnitsOfG
  Int_t ngvals    = fGLklVsG->GetN();
  Double_t* svval = new Double_t[ngvals];

  for(Int_t i=0;i<ngvals;i++)
    svval[i] = fGLklVsG->GetX()[i]*(applyUnits? fUnitsOfG :1);
  TGraph* thegraph = new TGraph(ngvals,svval,fGLklVsG->GetY());
  thegraph->SetName("fGLklVsG");

  delete [] svval;
  return thegraph;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Print gLike Banner
//
void Lkl::PrintGLikeBanner()
{
  std::cout << "***********************************************************************************" << std::endl;
  std::cout << "***********************************************************************************" << std::endl;
  std::cout << "***                                                                             ***" << std::endl;
  std::cout << "***                        88           88  88                                  ***" << std::endl;
  std::cout << "***                        88           \"\"  88                                  ***" << std::endl;
  std::cout << "***                        88               88                                  ***" << std::endl;
  std::cout << "***            ,adPPYb,d8  88           88  88   ,d8   ,adPPYba,                ***" << std::endl;
  std::cout << "***           a8\"    `Y88  88           88  88 ,a8\"   a8P_____88                ***" << std::endl;
  std::cout << "***           8b       88  88           88  8888[     8PP\"\"\"\"\"\"\"                ***" << std::endl;
  std::cout << "***           \"8a,   ,d88  88           88  88`\"Yba,  \"8b,   ,aa                ***" << std::endl;
  std::cout << "***            `\"YbbdP\"Y8  88888888888  88  88   `Y8a  `\"Ybbd8\"'                ***" << std::endl;
  std::cout << "***            aa,    ,88                                                       ***" << std::endl;
  std::cout << "***             \"Y8bbdP\"                                                        ***" << std::endl;
  std::cout << "***                                                                             ***" << std::endl;
  std::cout << "***                             V0.4 Oct 2017                                   ***" << std::endl;
  std::cout << "***                              J. Rico et al.                                 ***" << std::endl;
  std::cout << "***********************************************************************************" << std::endl;
  std::cout << "***********************************************************************************" << std::endl;
  std::cout << "***  IMPORTANT NOTE: THE USE OF THIS CODE TO PRODUCE PAPERS OF THE MAGIC        ***" << std::endl;
  std::cout << "***  AND/OR CTA COLLABORATIONS IS ALLOWED FOLLOWING THEIR RESPECTIVE            ***" << std::endl;
  std::cout << "***  PUBLICATION POLICIES FOR FULL-COLLABORATION PAPERS. FOR                    ***" << std::endl;
  std::cout << "***  PUBLICATIONS OUTSIDE THOSE FRAMEWORKS PLEASE CONTACT FIRST THE             ***" << std::endl;
  std::cout << "***  AUTHORS (Jelena Aleksic <jelena@ifae.es> AND Javier Rico                   ***" << std::endl;
  std::cout << "***  <mailto:jrico@ifae.es>), WHO COULD CLAIM AUTHORSHIP OVER THE               ***" << std::endl;
  std::cout << "***  RESULTING PAPER                                                            ***" << std::endl;
  std::cout << "***                                                                             ***" << std::endl;
  std::cout << "***  PLEASE CITE:                                                               ***" << std::endl;
  std::cout << "***  Aleksic, Rico & Martinez JCAP 10 (2012) 032                                ***" << std::endl;
  std::cout << "***********************************************************************************" << std::endl;
  std::cout << "***********************************************************************************" << std::endl;
}

////////////////////////////////////////////////////////////////
//
// integrate a histogram in logx between minimum and maximum values
// of log(x), but with respect to dx
//
Double_t Lkl::IntegrateLogE(const TH1F* h1,Double_t lmin,Double_t lmax) const
{
  // find the bins of integration limits
  const Int_t    nbins        = h1->GetNbinsX();
  const Double_t histolowedge = h1->GetBinLowEdge(1);
  const Double_t histouppedge = h1->GetBinLowEdge(nbins);
  if(lmin<histolowedge) lmin=histolowedge;
  if(lmax<histolowedge) lmax=histolowedge;
  if(lmin>histouppedge) lmin=histouppedge;
  if(lmax>histouppedge) lmax=histouppedge;

  Int_t lminbin = 1;
  Int_t lmaxbin = nbins;

  for(Int_t ibin = 0;ibin<lmaxbin;ibin++)
    {
      Double_t binlowedge = h1->GetBinLowEdge(ibin+1);
      Double_t binuppedge = h1->GetBinLowEdge(ibin+1)+h1->GetBinWidth(ibin+1);
      if(lmin>=binlowedge && lmin<binuppedge)
	lminbin = ibin+1;
      if(lmax>=binlowedge && lmax<binuppedge)
	lmaxbin = ibin+1;
    }
  if(lmax==h1->GetBinLowEdge(nbins)+h1->GetBinWidth(nbins))
     lmaxbin=nbins;
     
  // sanity checks
  if(lmin>lmax)
    {
      cout << "IntegrateLogE (" << GetName() << ") Error: bad integration limits (lmin = " << Form("%.6f",lmin) << ", lmax = " << Form("%.6f",lmax) << ", lminbin = " << lminbin << ", lmaxbin = " << lmaxbin << ")"   << endl;
      return -1;
    }
  

  // Compute integral in dE where x is logE
  Double_t integral = 0;
  Double_t ebinmin,ebinmax;
  Double_t eintmin = TMath::Power(10,lmin);
  Double_t eintmax = TMath::Power(10,lmax);

  if(lminbin == lmaxbin)
    integral+=h1->GetBinContent(lminbin)*(eintmax-eintmin);
  else
    {
      // first bin might not be completely integrated
      ebinmax = TMath::Power(10,h1->GetBinLowEdge(lminbin)+h1->GetBinWidth(lminbin)); 
  
      integral+=h1->GetBinContent(lminbin)*(ebinmax-eintmin);

      // loop over fully integrated bins
      for(Int_t ibin=lminbin+1;ibin<lmaxbin;ibin++)
	{
	  ebinmin = TMath::Power(10,h1->GetBinLowEdge(ibin));                       
	  ebinmax = TMath::Power(10,h1->GetBinLowEdge(ibin)+h1->GetBinWidth(ibin)); 

	  integral+=h1->GetBinContent(ibin)*(ebinmax-ebinmin);
	}
      
      // last bin might not be completely integrated
      ebinmin = TMath::Power(10,h1->GetBinLowEdge(lmaxbin));                       

      integral+=h1->GetBinContent(lmaxbin)*(eintmax-ebinmin);
    }
  
  return integral;    
}

