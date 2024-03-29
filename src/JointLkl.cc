/* ======================================================================== *\
!
!   Author(s): Jelena Aleksic      06/2012 <mailto:jelena@ifae.es>
!              Javier Rico         12/2014 <mailto:jrico@ifae.es>
!
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// JointLkl
//
// This class performs a joint likelihood maximization (minimization
// of -2logL) of a set of samples (each represented by an Lkl object)
// for which g_i share a common physical origin.
//
// The only free parameter is g_0, i.e. g for the *reference* sample
// (normally the first sample added to the JointLkl object, but not
// always, see note below).  The value of g in subsequent samples
// (g_1,..., g_n-1) depend on g_0. Each
// Lkl object can contain different nuisance parameters. Nuisance
// parameters from different Lkl instances will be considered
// independent from each other. If several Lkl objects share a common
// uncertainty in units of G, do NOT set fDUnitsOfG for each of them,
// since in such a case they will be considered as independent
// nuisance parameters. Instead, add them all to a JointLkl instance, for
// which you set the common value of fDUnitsOfG.
//
// VERY IMPORTANT: You need to define fUnitsOfG (with SetUnitsOfG) for 
// every Lkl object added to a given JointLKl object. This is so, because
// fUnitsOfG is used to compute the values of g_1,...,g_n given a free
// parameter g_0. fUnitsOfG must be such that for all added objects
// g_i*fUnitsOfG_i has the same value, i.e. a common physical quantity.
//
// IMPORTANT: the reference sample (that which g_0 refers to) may not
// necessarily be the one added first to the JointLkl object.
// JointLkl internally rearranges the Lkl objects in case of
// need (e.g. when the first introduced sample is not sensitive to a
// given signal model, see JointLkl::CheckSamples). Use the function 
// JointLkl::GetReferenceSample to access the reference sample after
// the likelihood minimization has finished.
//
// Usage example:
// --------------
// (for a fully working example see exec/jointLkl.cc)
//
// Iact1dUnbinnedLkl** fLkl = new Iact1dUnbinnedLkl*[nsample];
// JointLkl* jLkl = new JointLkl(jointInputString);
// 
// // Configure each Iact1dUnbinnedLkl sample
// for(Int_t isample=0;isample<nsample;isample++)
//   {
//      // Configure the samples
//      fLkl[isample] = new Iact1dUnbinnedLkl(inputString[isample]);
//      fLkl[isample]->AdddNdESignal(dNdESignalFileName);
//      fLkl[isample]->SetDMAnnihilationUnitsForG(mass);	
//
//      // add the sample to the JointLkl object
//      jLkl->AddSample(fLkl[isample]);
//   }
//
// // Configure the Joint Likelihood
// jLkl->SetErrorDef(deltaLogLkl);              // set the error correponding to the required CL
// jLkl->SetGIsPositive();                      // only G values are considered (Fermi prescription)
//
// // Compute and profile -2logL curve
// jLkl->ComputeLklVsG();                        // compute joint -2logL vs g around minimum
//
// // get the results
// Double_t gminval     = jLkl->GetGLklMin();    // value of g minimizing joint -2logL
// Double_t glimval     = jLkl->GetGForLkl(2.7); // value of g for which -2logL = minimum+2.7
// Double_t gminvalErr  = glimval-gminval;       // determination of 95% CL error 
//
// jLkl->GetLklVsG()->Draw("al");                // draw joint -2logL vs g (in specified units) 
//
//////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "TMath.h"

#include "JointLkl.h"

ClassImp(JointLkl);

using namespace std;

// static constants
static const TString  gName            = "JointLkl";
static const TString  gTitle           = "Joint Likelihood";
static const Int_t    gNPars           = 1;         // Number of free+nuisance parameters
static const TString  gInputString     = "DlinG=0";
  
// -2logL function for minuit
void jointLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;

////////////////////////////////////////////////////////////////
// 
// Default constructor
// if sampleArray is provided, it is used as initial list of
// samples, which can be later expanded by AddSample
//
JointLkl::JointLkl(TString inputString) : 
  Lkl(gNPars,(inputString==""?gInputString:inputString) ,gName,gTitle), fSampleArray(NULL)
{
  fSampleArray = new TObjArray();
}


////////////////////////////////////////////////////////////////
// 
// Destructor
//
JointLkl::~JointLkl()
{
  delete fSampleArray;
}
				   
////////////////////////////////////////////////////////////////
//
// Check that the samples are correct
//
// Return 0 in case of success, 1 otherwise
//
Int_t JointLkl::MakeChecks()
{
  if(ReorderSamples()) return 1;
  
  // Check all added objects
  TObjArrayIter* iter    = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    if(!sample->IsChecked())
	sample->MakeChecks();
  delete iter;
  
  // units of g are those of the reference sample
  SetUnitsOfG(GetReferenceSample()->GetUnitsOfG());

  SetChecked();

  return 0;
}

////////////////////////////////////////////////////////////////
//
// Check that the object has some daughter objects before
// calling the computation of Lkl vs G
//
// Check Lkl::ComputeLklVsG for explanation about the arguments
//
Double_t  JointLkl::ComputeLklVsG(Double_t gInitWithUnits,Bool_t centerAtZero,Int_t npoints,Double_t glow,Double_t gupp,Bool_t isVerbose)
{
  if(fSampleArray->GetEntries()<=0)
    {
      cout << "JointLkl::ComputeLklVsG Error: object " << GetName() << " is of type JointLkl, but no other object depend on it" << endl;
      return 1;
    }

  return Lkl::ComputeLklVsG(gInitWithUnits,centerAtZero,npoints,glow,gupp,isVerbose);
}
////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the free parameter
//
void  JointLkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(jointLkl);

  fMinuit->SetName(Form("%s_Minuit",GetName()));

  // Initialize minuit for added objects
  TObjArrayIter* iter    = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      if(GetGIsPositive())
	sample->SetGIsPositive();
      sample->InitMinuit();
      sample->SetMinuitLink();
    }
  delete iter;
  
  // take the parameter configuration from the reference sample
  const Char_t* parName[gNPars] = {GetReferenceSample()->GetParName()[0].Data()};    
  Double_t      pStart[gNPars]  = {ginit};   
  Double_t      pDelta[gNPars]  = {GetReferenceSample()->GetParDelta()[0]};   

  SetParameters(parName,pStart,pDelta);  
}  
	
////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of minuit subtleties)
//
void JointLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN fulllkl
}

////////////////////////////////////////////////////////////////
//
// If needed, reorder samples so that the first (or "reference")
// one has non-zero predicted signal events
//
// Return 0 even if no sample produce signal because this could be ok
// if this JointLkl is in turn part of another JointLkl list in which
// there are classes producing a signal
//
Int_t JointLkl::ReorderSamples()
{
  // check that integral for sample 0 is not 0 (used as denominator later)
  // if it is, move it to the end of the array and try with next
  TObjArrayIter* iter = (TObjArrayIter*)fSampleArray->MakeIterator();
  Lkl* sample0  = (Lkl*) iter->Next();
  Int_t nchecks=0;
  Int_t nmaxchecks = fSampleArray->GetEntries();
  while(sample0->GetUnitsOfG()==0)
    {
      TString refname = sample0->GetName();
      fSampleArray->Remove(sample0);
      fSampleArray->Add(sample0);
      sample0 = (Lkl*) iter->Next();
  
      cout << "JointLkl::ReorderSamples (" << GetName() << ") Message: sample " <<refname<<" has been sent to the end" << endl;
      if(++nchecks>=nmaxchecks)
	{
	  cout << "JointLkl::ReorderSamples (" << GetName() << ") Warning: none of the samples will produce any signal event" << endl;
	  break;
	}
    }
  fSampleArray->Compress();

  delete iter;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Reset the -2logL curve for all added samples
//
void JointLkl::ResetGLklVsG()
{
  Lkl::ResetGLklVsG();
  // reset GLklVsG curve in all added samples
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    sample->ResetGLklVsG();
  delete iter;

  SetChecked(kFALSE);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// compute the values of nuisance parameter of added samples
// This is an internal function called in Lkl::FixLklVsG
//
void JointLkl::SpreadFixLklVsG(Double_t g)
{
  // compute the values of nuisance parameter of added samples
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();

  // read weight for reference sample
  Lkl*          sample0 = (Lkl*) iter->Next();
  Double_t       w0      = sample0->GetUnitsOfG();
  iter->Reset();

  // loop over all samples and compute the -2logL curve in the corresponding range
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      Double_t   w = sample->GetUnitsOfG(); 
      if(w>0)
	sample->MinimizeLkl(g*w0/w,kTRUE,kFALSE,kTRUE);
    }
  delete iter;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Compute the -2logL curve for all added samples between glow and gupp
// This is an internal function called in Lkl::ComputeLklVsG
//
Int_t JointLkl::PrepareForLklScan(Bool_t centerAtZero,Int_t npoints,Double_t glow,Double_t gupp,Bool_t isVerbose)
{
  cout << "JointLkl::PrepareForLklScan (" << GetName() << ") Message: call ComputeLklVsG for " << GetSampleArray()->GetEntries()<< " samples:" << endl;
  // compute the -2logL curve in all added samples
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();

  // read weight for reference sample
  Lkl*          sample0 = (Lkl*) iter->Next();
  Double_t       w0      = sample0->GetUnitsOfG();
  iter->Reset();

  // loop over all samples and compute the -2logL curve in the corresponding range
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      Double_t   w = sample->GetUnitsOfG(); 
      if(w>0)
	sample->ComputeLklVsG(0.5*(glow+gupp)*w,centerAtZero,npoints,glow*w0/w,gupp*w0/w,isVerbose);
    }

  delete iter;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Return the expansion factor to take into account uncertainty in units of G 
// from added samples
//
Double_t JointLkl::GetExpansionCoefficient() const
{
  // compute the -2logL curve in all added samples
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();

  // loop over all samples and average fDUnitsOfG and take the highest value of fDUofGType
  Double_t avDUnitsOfG = 1;
  DUofGType_t type     = none;
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      Double_t    w     = sample->GetUnitsOfG(); 
      DUofGType_t stype = sample->GetDUofGType();
      if(w>0)
	{
	  avDUnitsOfG*=sample->GetDUnitsOfG();
	  if(stype>type) type=stype;
	}
    }
  delete iter;

  // compute average
  Int_t nsamples = GetSampleArray()->GetEntries();
  if(nsamples)
    avDUnitsOfG=TMath::Power(avDUnitsOfG,1./nsamples);

  // compute the expansion factor for the average error and highest type
  Double_t stretch = 1;
  if(avDUnitsOfG>0 && type!=none)
    {
      if(type==lin)
	stretch = 1+GetErrorDef()*avDUnitsOfG;
      else if(type==invlin)
	stretch = 1+1./(GetErrorDef()*avDUnitsOfG);
      else if(type==log)
	stretch = TMath::Power(10,GetErrorDef()*avDUnitsOfG);
      else if(type==invlog)		  
	stretch = 1./TMath::Power(10,-GetErrorDef()*avDUnitsOfG);
    }

  return stretch;  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Print information about all added samples and their fit results
//
void JointLkl::PrintOverview(Int_t level) 
{
  // plot info about the JointLkl
  Lkl::PrintOverview(level);
  Margin(level); cout << "  # of samples      = " << GetNSamples() << endl;

  // plot info about added samples
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  Margin(level); cout << "***************************************" << endl;
  while((sample=(Lkl*)iter->Next()))
    {
      sample->PrintOverview(level+1);
      Margin(level+1); cout << "***************************************" << endl;
    }
  delete iter;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Print information about all added samples configuration
//
void JointLkl::PrintData(Int_t level) 
{
  Lkl::PrintData(level);
  Margin(level); cout << "             # of samples = " << GetNSamples() << endl;

  // plot info about added samples
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  Margin(level); cout << "***************************************" << endl;
  while((sample=(Lkl*)iter->Next()))
    {
      sample->PrintData(level+1);
      Margin(level+1); cout << "***************************************" << endl;
    }
  delete iter;
}

////////////////////////////////////////////////////////////////////////
//
// Joint likelihood function (-2logL) for all samples
// To be minimized by TMinuit
// Free parameters:
// par[0] = g0 (total estimated number of signal events in On region of sample 0)
//
void jointLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;
  f = 0;

  minuit->GetObjectFit();
  
  // get internal object
  JointLkl*     mylkl   = dynamic_cast<JointLkl*>(minuit->GetObjectFit());

  // g for sample 0 (taken as reference) and corresponding weight
  TObjArrayIter* iter    = (TObjArrayIter*) mylkl->GetSampleArray()->MakeIterator();
  Lkl*          sample0 = (Lkl*) iter->Next();
  Double_t       w0      = sample0->GetUnitsOfG();
  Double_t       g0      = par[0];

  // loop over samples and compute the corresponding g, and compute the contribution to the likelihood from the sample
  iter->Reset();
  Lkl* sample;

  while((sample=(Lkl*)iter->Next()))
    {
      Double_t   w = sample->GetUnitsOfG();
      if(w>0)
	{
	  Double_t   g = g0*w0/w;      
	  f+=sample->MinimizeLkl(g,kTRUE,kFALSE);      
	}
    }   
  
  // recover the link to TMinuit that might be "stolen" if there are JointLkl samples
  mylkl->SetMinuitLink();

  delete iter;
} 
