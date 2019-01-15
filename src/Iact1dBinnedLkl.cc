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
// Iact1dBinnedLkl
//
// Class to perform BINNED likelihood maximization (minimization of -2logL)
// to estimate the presence of signal events following a certain spectral 
// shape (no physics origin is assumed and the class can be used as long as 
// one has an a priori knowledge of the signal spectral shape).
//
// The free parameter g is the total number of signal events in the On
// region. The Off/On normalization (tau, needs to be set when calling
// ReadDataSamplesSegueStereo or any other method to read the input
// data) is a nuisance parameter, and it's common to all bins. You can
// fix tau by making Dtau<=0. In each bin, the number of background
// events in the On region (b) is also a nuisance parameter internally
// by the class.
//
// Upon object creation, we specify the maximum and minimum E' of the
// events to be considered, and the number of bins that we want to use
// in the likelihood. If fMinBinContent>0 (default is 0), the On and Off
// histograms will be rebinned in such a
// way that every bin of both histograms have at least fMinBinContent
// events (can be set with SetMinBinContent). Bins for which the On or
// On histograms have less events are merged to subsequent bins until
// fMinBinContent are gathered. If after this process the number
// of events of the last bin of the On or Off distributions is still
// below the threshold, the last two bins are merged.
//
// Iact1dBinnedLkl inherits Iact1dUnbinnedLkl, so that all methods for the
// definition of the signal and background dN/dE' histograms, on and
// off samples, etc can be used. From the integral of dN/dE'_signal in
// the different bins we compute the corresponding g_i in each bin for
// a given value of g. This, together with Non and Noff values for the
// bins, and Tau, is used to defined an PoissonLkl object per
// bin. Those objects are associated to the Iact1dBinnedLkl through its
// JointLkl inheritance.
//
//
// Usage example:
// --------------
// (for a fully working example see macro jointLkl.C)
//
// Iact1dBinnedLkl* bLkl = new Iact1dBinnedLkl(Emin,Emax,nbins);
// bLkl->SetErrorDef(2.7);  // for 1-sided 95% CL
// bLkl->ReadAeffSegueStereo(aEffFileName);
// bLkl->ReadEResoAndBiasSegueStereo(energyRFileName);
// bLkl->ReaddNdEpBkgSegueStereo(bkgModFileName);
// bLkl->ReadDataSamplesSegueStereo(evtFileName,tau,dtau);
// bLkl->ReaddNdESignalSegueStereo(dNdESignalFileName);
//
// // Set units for DM interpretation of results
// bLkl->SetDMAnnihilationUnitsForG(Teff,mass,log10_J);
//
// // Compute the profile -2logL curve
// bLkl->ComputeLklVsG();
// 
// // Get the results
// Double_t gminval     = bLkl->GetGLklMin();    // value of g minimzing -2logL
// Double_t gminvalErr  = bLkl->GetGLklMinErr(); // 1-sided 95% CL error bar from Migrad
// Double_t glimval     = bLkl->GetGForLkl(2.7); // value of g for which -2logL = minimum+2.7
// Double_t gminvalErr2 = glimval-gminval;       // better determination of 95% CL error 
//
// // Plot the -2logL curve
// bLkl->GetLklVsG()->Draw("al");                // draw -2logL vs g 
//
//////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "TMath.h"
#include "TPRegexp.h"

#include "Iact1dBinnedLkl.h"
#include "PoissonLkl.h"

ClassImp(Iact1dBinnedLkl);

using namespace std;

// static constants
static const TString  gName            = "Iact1dBinnedLkl";
static const TString  gTitle           = "Binned Full Likelihood";
static const Int_t    gNPars            = 2;            // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars]  = {"g","tau"};  // Name of parameters

// -2logL function for minuit
void binnedIact1dUnbinnedLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;

void GetRebinning(TH1F*hOn,TH1F*hOff,UInt_t minNInBin,UInt_t& inewbin,Double_t* newbin);


////////////////////////////////////////////////////////////////
// 
// String constructor
//
//
Iact1dBinnedLkl::Iact1dBinnedLkl(TString inputString) :
  Lkl(gNPars,inputString,gName,gTitle), Iact1dUnbinnedLkl(inputString), JointLkl(inputString),
  fNBins(gDefNBins), fMinBinContent(gDefMinBinContent), fNRemovedBins(0), fKnownBackground(kFALSE), 
  fTauEDepFluct(kFALSE), fdNdEpBkgFromOff(kFALSE), fHNOn(NULL), fHNOff(NULL)
{
  if(InterpretInputString(inputString))
    cout << "Iact1dBinnedLkl::Iact1dBinnedLkl Warning: there were problems interpreting the input string" << endl;
}

//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
// The string has been already passed to Iact1dUnbinnedLkl::InterpretInputString in the constructor
// Here we search for the following options:
//
// <nbins>=val:            is the number of bins 
// <minbincontent>=val:    is the minimum number of events in each bin
// <tauEDepFluct>=val:     val can be TRUE or FALSE (default). If true, tau will be allowed to fluctuate independently from energy bin to energy bin
// <knownBackground>=val:  val can be TRUE or FALSE (default). If true, b_i (i=1,...,Nbins) are fixed parameters (as opposed to nuisance)
// <dndepbkgfromoff>=val:  val can be TRUE or FALSE (default). If true, the background distribution (fHdNdEpBkg, for simulations) will be constructed from Off data distribution
//
Int_t Iact1dBinnedLkl::InterpretInputString(TString inputString)
{
  // Prepeare to interpret inputString
  TPMERegexp re("\\s+");
  TPMERegexp fldre("=");
  
  // split the inputString into the different fields, and check the option and values specified in each of them
  UInt_t nfields = re.Split(inputString);
  for(UInt_t ifield=0;ifield<nfields;ifield++)
    {
      fldre.Split(re[ifield]);
      TString optname = fldre[0];
      if(optname.CompareTo("nbins",TString::kIgnoreCase)==0)
	fNBins=fldre[1].Atof();
      else if(optname.CompareTo("minbincontent",TString::kIgnoreCase)==0)
	fMinBinContent=fldre[1].Atoi();      
      else if(optname.CompareTo("knownBackground",TString::kIgnoreCase)==0)
	{
	  if (fldre[1].CompareTo("TRUE",TString::kIgnoreCase)==0)
	    fKnownBackground=kTRUE;
	}
      else if(optname.CompareTo("tauEDepFluct",TString::kIgnoreCase)==0)
	{
	  if (fldre[1].CompareTo("TRUE",TString::kIgnoreCase)==0)
	    fTauEDepFluct=kTRUE;
	}
      else if(optname.CompareTo("dNdEpBkgFromOff",TString::kIgnoreCase)==0)
	{
	  if (fldre[1].CompareTo("TRUE",TString::kIgnoreCase)==0)
	    {
	      fdNdEpBkgFromOff=kTRUE;
	      ResetHdNdEpBkg();
	    }
	}
    }

  return 0;
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
Iact1dBinnedLkl::~Iact1dBinnedLkl()
{
  if(fHNOn)  delete fHNOn;
  if(fHNOff) delete fHNOff;
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of TMinuit subtleties)
//
void Iact1dBinnedLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN fulllkl
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance),
//
void  Iact1dBinnedLkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(binnedIact1dUnbinnedLkl);
  fMinuit->SetName(Form("%s_Minuit",GetName()));

  // Initialize minuit for added objects
  TObjArrayIter* iter    = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      sample->InitMinuit();
      sample->SetMinuitLink();
    }
  delete iter;
  
  // setup parameters
  Double_t pStart[gNPars] = {ginit,GetTau()};
  Double_t gPDelta[gNPars]   = {TMath::Sqrt(GetNoff())/10.,GetDTau()};  // Precision of parameters during minimization

  SetParameters(gParName,pStart,gPDelta);

  // Fix tau if requested
  fMinuit->Release(Iact1dBinnedLkl::gTauIndex);
  FixPar(Iact1dBinnedLkl::gTauIndex,kFALSE);
  if(GetDTau()<=0 || fTauEDepFluct)
    { 
      fMinuit->FixParameter(Iact1dBinnedLkl::gTauIndex);
      FixPar(Iact1dBinnedLkl::gTauIndex);
    }
}

////////////////////////////////////////////////////////////////
//
// Check that all needed input is present and correct
//
// Return 0 in case of success, 1 otherwise
//
Int_t Iact1dBinnedLkl::MakeChecks()
{
  if(IsChecked()) return 0;

  // Check the Iact1dUnbinnedLkl part
  ///////////////////////////
  // check if all needed histos are there
  if(CheckHistograms(kFALSE))
    {
      cout << "Iact1dBinnedLkl::MakeChecks (" << GetName() << ") Warning: missing information, cannot perform fit, check your code!" << endl;
      return 1;
    }

  // normalize unnormalized histos
  if(GetdNdEpSignalIntegral()<=0)
    {
      cout << "Iact1dBinnedLkl::MakeChecks (" << GetName() << ") Warning: fHdNdEpSignal is not normalized!" << endl;
      return 1;
    }
  
  if(GetHdNdEpSignalOff())
    if(GetdNdEpSignalOffIntegral()<=0)
      {
	cout << "Iact1dBinnedLkl::MakeChecks (" << GetName() << ") Warning: fHdNdEpSignalOff is not normalized!" << endl;
	return 1;
      }
  
  if(GetdNdESignalIntegral()<=0) 
    {
      cout << "Iact1dBinnedLkl::MakeChecks (" << GetName() << ") Warning: fHdNdESignal is not normalized!" << endl;
      return 1;
    }

  if(GetHdNdEpFrg())
    if(GetdNdEpFrgIntegral()<=0) 
      {
	cout << "Iact1dBinnedLkl::MakeChecks (" << GetName() << ") Warning: fHdNdEpFrg is not normalized!" << endl;
	return 1;
      }
  
  // Check the Iact1dBinnedLkl specific part
  /////////////////////////////////////////
  // get the dN/dE' histograms for On and Off and check binning
  if(BuildAndBinOnOffHistos())
    {
      cout << "Iact1dBinnedLkl::MakeChecks (" << GetName() << ") Warning: problems building On and/or Off histos!" << endl;
      return 1;
    }

  // configure the JointLkl part
  if(ConfigureJointLkl())
    {
      cout << "Iact1dBinnedLkl::MakeChecks (" << GetName() << ") Warning: problems configuring JointLkl list!" << endl;
      return 1;
    }

  // Check the JointLkl part
  ////////////////////////////
  if(ReorderSamples())
    return 1;

  // Check all bins
  TObjArrayIter* iter    = (TObjArrayIter*) GetSampleArray()->MakeIterator();
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    if(!sample->IsChecked())
      sample->MakeChecks();
  delete iter;

  SetChecked();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// compute the values of nuisance parameter of added samples
// This is an internal function called in Lkl::FixLklVsG
//
void Iact1dBinnedLkl::SpreadFixLklVsG(Double_t g)
{
  // compute the values of nuisance parameter of added samples
  TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();

  // loop over all samples and compute the -2logL curve in the corresponding range
  Lkl* sample;
  while((sample=(Lkl*)iter->Next()))
    {
      Double_t   w = sample->GetUnitsOfG(); 
      if(w>0)
	sample->MinimizeLkl(g*w,kTRUE,kFALSE,kTRUE);
    }
  delete iter;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Compute the -2logL curve for all added samples between glow and gupp
// This is an internal function called in Lkl::ComputeLklVsG
//
Int_t Iact1dBinnedLkl::PrepareForLklScan(Bool_t centerAtZero,Int_t npoints,Double_t glow,Double_t gupp,Bool_t isVerbose)
{
  // cout << "JointLkl::PrepareForLklScan (" << GetName() << ") Message: call ComputeLklVsG for " << GetSampleArray()->GetEntries()<< " samples" << endl;
  // compute the -2logL curve in all added samples

  // if tau is not nuisance then we can precompute the bin(s) likelihoods
  if(IsParFixed(Iact1dBinnedLkl::gTauIndex))
    {
      TObjArrayIter* iter = (TObjArrayIter*) GetSampleArray()->MakeIterator();
      
      // loop over all samples and compute the -2logL curve in the corresponding range
      Lkl* sample;
      while((sample=(Lkl*)iter->Next()))
	{
	  Double_t   w = sample->GetUnitsOfG(); 
	  if(w>0)
	    sample->ComputeLklVsG(centerAtZero,npoints,glow*w,gupp*w,isVerbose*0);
	}
      delete iter;
    }

  return 0;
}
////////////////////////////////////////////////////////////////
//
// Take the list of On and Off events and make histograms out of 
// them with the specified binning.
// The minimum number of entries per bin in any of the histograms
// is fMinBinContent (default 10), which can be changed through
// SetMinBinContent.
// If fMinBinContent>0 (default 0),
// the binning is changed so that the condition on
// the minimum number of entries is fulfilled.
//
// Return 0 in case of success, 1 otherwise
//
Int_t Iact1dBinnedLkl::BuildAndBinOnOffHistos()
{
  if(!fHNOn || !fHNOff)
    {      
      // Get the E' distribution for On and Off events
      TH1F* provHNOn  = GetHdNdEpOn(kFALSE,fNBins);
      TH1F* provHNOff = GetHdNdEpOff(kFALSE,fNBins);
      
      if(!provHNOn || !provHNOff)
	{
	  cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: problems creating On and/or Off histograms, On = " << provHNOn << ", Off = " << provHNOff << "." << endl;
	  if(provHNOn)  delete provHNOn;
	  if(provHNOff) delete provHNOff;
	  return 1;
	}

      Bool_t done=kFALSE;
      fNRemovedBins=0;
      // Rebin if necessary and requested
      if(fMinBinContent>0)
	{
	  UInt_t nnewbins;
	  Double_t* newbin = new Double_t[fNBins+1];
	  GetRebinning(provHNOn,provHNOff,fMinBinContent,nnewbins,newbin);

	  if(nnewbins<fNBins)
	    {
	      fNRemovedBins = fNBins-nnewbins;
	      cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Message: fHNON/fHNOff original number of bins = " << fNBins << ", rebinned to " << nnewbins << " bins to keep minimum of " << fMinBinContent << " events per bin" << endl;

	      // Get the rebinned On/Off histograms
	      TH1F* hRebinOn  =  (TH1F*)  provHNOn->Rebin(nnewbins,"hRebinOn", newbin);
	      TH1F* hRebinOff =  (TH1F*) provHNOff->Rebin(nnewbins,"hRebinOff",newbin);
	      hRebinOn->SetDirectory(0);
	      hRebinOff->SetDirectory(0);
	      
	      // replace the fHNOn and fHNOff histograms by the rebinned ones
	      fHNOn  = new TH1I("fHNOn", "E' distribution of On events", nnewbins,newbin);
	      fHNOff = new TH1I("fHNOff","E' distribution of Off events",nnewbins,newbin);
	      fHNOn->SetDirectory(0);
	      fHNOff->SetDirectory(0);
	      for(Int_t ibin=0;ibin<nnewbins;ibin++)
		{
		  fHNOn->SetBinContent(ibin+1,Int_t(hRebinOn->GetBinContent(ibin+1)));
		  fHNOff->SetBinContent(ibin+1,Int_t(hRebinOff->GetBinContent(ibin+1)));
		}
	      
	      delete hRebinOn;
	      delete hRebinOff;
	      done = kTRUE;
	    }

	  delete [] newbin;
	}


      if(!done)
	{
	  fHNOn  = new TH1I("fHNOn", "E' distribution of On events", fNBins,provHNOn->GetXaxis()->GetXmin(),provHNOn->GetXaxis()->GetXmax());
	  fHNOff = new TH1I("fHNOff","E' distribution of Off events",fNBins,provHNOff->GetXaxis()->GetXmin(),provHNOff->GetXaxis()->GetXmax());
	  fHNOn->SetDirectory(0);
	  fHNOff->SetDirectory(0);

	  for(Int_t ibin=0;ibin<fNBins;ibin++)
	    {
	      fHNOn->SetBinContent(ibin+1,Int_t(provHNOn->GetBinContent(ibin+1)));
	      fHNOff->SetBinContent(ibin+1,Int_t(provHNOff->GetBinContent(ibin+1)));
	    }    
	}

      if(provHNOn)  delete provHNOn;
      if(provHNOff) delete provHNOff;
    }

  // check that everything went fine before leaving
  if(!fHNOn || !fHNOff)
    {
      cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: problems rebinning On and/or Off histograms" << endl;
      return 1;
    }

  // check that every bin has more than the minimum required number of events (fMinBinContent)
  Bool_t binsAreOk = kTRUE;
  for(Int_t ibin=0;ibin<fNBins-fNRemovedBins;ibin++)
    {
      Int_t onevts  = fHNOn->GetBinContent(ibin+1);
      Int_t offevts = fHNOff->GetBinContent(ibin+1);
      if(onevts<fMinBinContent)
	{
	  binsAreOk = kFALSE;
	  cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: number of On events in bin #" << ibin+1 
	       << " (" << onevts << ") is below the minimum requested (" << fMinBinContent << ")" << endl;
	}
      if(offevts<fMinBinContent)
	{
	  binsAreOk = kFALSE;
	  cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: number of Off events in bin #" << ibin+1 
	       << " (" << offevts << ") is below the minimum requested (" << fMinBinContent << ")" << endl;
	}
    }

  
  if(!binsAreOk) return 1;
  
  return 0;
}

////////////////////////////////////////////////////////////////
//
// Feed the JointLkl with one PoissonLkl per E' bin 
//
// Return 0 in case of success, 1 otherwise
//
Int_t Iact1dBinnedLkl::ConfigureJointLkl()
{
  if(!fHNOn || !fHNOff)
    {
      cout << "Iact1dBinnedLkl::ConfigureJointLkl (" << GetName() << ") Warning: On and/or Off histograms do no exist" << endl;
      return 1;
    }

  ClearSampleList(); // delete samples in the list

  // fill the Poisson Lkl classes for every bin in the list
  Double_t totalw = 0;
  for(Int_t ibin=0;ibin<fNBins-fNRemovedBins;ibin++)
    {
      PoissonLkl* pLkl = new PoissonLkl(fHNOn->GetBinContent(ibin+1),fHNOff->GetBinContent(ibin+1),GetTau(),(fTauEDepFluct? GetDTau() : 0),0,Form("%s_Bin_%d",GetName(),ibin));
      Double_t lemin  = fHNOn->GetBinLowEdge(ibin+1);
      Double_t lemax  = fHNOn->GetBinLowEdge(ibin+1)+fHNOn->GetBinWidth(ibin+1);
      
      Double_t weight = IntegrateLogE(GetHdNdEpSignal(),lemin,lemax);
      pLkl->SetUnitsOfG(weight>0? weight : 0);

      // is there signal leakeage in the Off region?
      if(GetHdNdEpSignalOff())
	{
	  Double_t gfractioninoff = IntegrateLogE(GetHdNdEpSignalOff(),lemin,lemax)*GetdNdEpSignalOffIntegral()/(weight*GetdNdEpSignalIntegral());
	  pLkl->SetGFractionInOff(gfractioninoff);
	}

      // are there gamma-ray foreground events in the signal region?
      if(GetHdNdEpFrg())
	{
	  Double_t nFrgEvts = IntegrateLogE(GetHdNdEpFrg(),lemin,lemax)*GetdNdEpFrgIntegral();
	  pLkl->SetFrgNEvents(nFrgEvts);
	}

      // is b a fixed (as opposed to nuisance) parameter? mostly for tests
      if(fKnownBackground) pLkl->SetKnownBackground();
	

      AddSample(pLkl);
      
      totalw+=weight;
    }
  SetOwner(); // so that all created pLkl are deleted with the JointLkl object

  // check that the sum of all weights is approximately equal to 1
  if(TMath::Abs(totalw-1)>0.01)
    {
      cout << "Iact1dBinnedLkl::ConfigureJointLkl (" << GetName() << ") Warning: total sum of weights should be 1! (It's " << totalw << "), too many bins maybe??" << endl;
      return 1;
    }

  return 0;
}

//////////////////////////////////////////////////////////////////
//
// Produce the E' distribution of On events and return the 
// corresponding histogram.
// If fHNOn does not exist, then fall back to base method
// defined in Iact1dUnbinnedLkl. Otherwise, return a copy of fHNOn
// If isDifferential=kTRUE (default is kTRUE), return the 
// dN/dE' distribution (i.e. number of entries in the bin divided
// by the size of the bin in E' unit). Otherwise the returned
// histogram is the number of events per bin of E'.
// The returned histogram will not be deleted by the destructor
// of Iact1dBinnedLkl.
//
TH1F* Iact1dBinnedLkl::GetHdNdEpOn(Bool_t isDifferential,Int_t nbins) const
{
  if(!fHNOn || (nbins>0 && nbins!=fNBins-fNRemovedBins))
    return Iact1dUnbinnedLkl::GetHdNdEpOn(isDifferential,nbins);

  // check if bin width is constant
  Bool_t binWidthIsConstant = kTRUE;
  for(Int_t ibin=0;ibin<fNBins-fNRemovedBins-1;ibin++)
    if(TMath::Abs(Float_t(fHNOn->GetBinWidth(ibin+1))-Float_t(fHNOn->GetBinWidth(ibin+2)))>1e-9)
      binWidthIsConstant = kFALSE;

  TH1F* h;
  if(binWidthIsConstant)
    h = new TH1F("dNdEpOn","dN/dE' for On events",fNBins-fNRemovedBins,fHNOn->GetXaxis()->GetXmin(),fHNOn->GetXaxis()->GetXmax());
  else
    h = new TH1F("dNdEpOn","dN/dE' for On events",fNBins-fNRemovedBins,fHNOn->GetXaxis()->GetXbins()->GetArray());

  h->SetDirectory(0);
  h->SetXTitle("log_{10}(E' [GeV])");
  h->SetYTitle("dN/dE' [GeV^{-1}]");
  for(Int_t ibin=0;ibin<fNBins-fNRemovedBins;ibin++)
    {
      h->SetBinContent(ibin+1,Float_t(fHNOn->GetBinContent(ibin+1)));
      h->SetBinError(ibin+1,TMath::Sqrt(fHNOn->GetBinContent(ibin+1)));
    }

   // divide by bin width
  if(isDifferential)
    for(Int_t ibin=0;ibin<fNBins-fNRemovedBins;ibin++)
      {
	Double_t leminbin = h->GetBinLowEdge(ibin+1);
	Double_t lemaxbin = leminbin+h->GetBinWidth(ibin+1);
	Double_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);
	h->SetBinError(ibin+1,h->GetBinError(ibin+1)/deltaE);
	h->SetBinContent(ibin+1,h->GetBinContent(ibin+1)/deltaE);
      }

  return h;
}


//////////////////////////////////////////////////////////////////
//
// Produce the E' distribution of Off events and return the 
// corresponding histogram.
// If fHNOff does not exist, then fall back to base method
// defined in Iact1dUnbinnedLkl. Otherwise, return a copy of fHNOff
// If isDifferential=kTRUE (default is kTRUE), return the 
// dN/dE' distribution (i.e. number of entries in the bin divided
// by the size of the bin in E' unit). Otherwise the returned
// histogram is the number of events per bin of E'.
// The returned histogram will not be deleted by the destructor
// of Iact1dBinnedLkl.
//
TH1F* Iact1dBinnedLkl::GetHdNdEpOff(Bool_t isDifferential,Int_t nbins) const
{
  if(!fHNOff || (nbins>0 && nbins!=fNBins-fNRemovedBins))
    return Iact1dUnbinnedLkl::GetHdNdEpOff(isDifferential,nbins);

  // check if bin width is constant
  Bool_t binWidthIsConstant = kTRUE;
  for(Int_t ibin=0;ibin<fNBins-fNRemovedBins-1;ibin++)
    if(TMath::Abs(Float_t(fHNOff->GetBinWidth(ibin+1))-Float_t(fHNOff->GetBinWidth(ibin+2)))>1e-9)
      binWidthIsConstant = kFALSE;

  TH1F* h;
  if(binWidthIsConstant)
    h = new TH1F("dNdEpOff","dN/dE' for Off events",fNBins-fNRemovedBins,fHNOff->GetXaxis()->GetXmin(),fHNOff->GetXaxis()->GetXmax());
  else
    h = new TH1F("dNdEpOff","dN/dE' for Off events",fNBins-fNRemovedBins,fHNOff->GetXaxis()->GetXbins()->GetArray());
  h->SetDirectory(0);
  h->SetXTitle("log_{10}(E' [GeV])");
  h->SetYTitle("dN/dE' [GeV^{-1}]");
  for(Int_t ibin=0;ibin<fNBins-fNRemovedBins;ibin++)
    {
      h->SetBinContent(ibin+1,Float_t(fHNOff->GetBinContent(ibin+1)));
      h->SetBinError(ibin+1,TMath::Sqrt(fHNOff->GetBinContent(ibin+1)));
    }

   // divide by bin width
  if(isDifferential)
    for(Int_t ibin=0;ibin<fNBins-fNRemovedBins;ibin++)
      {
	Double_t leminbin = h->GetBinLowEdge(ibin+1);
	Double_t lemaxbin = leminbin+h->GetBinWidth(ibin+1);
	Double_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);
	h->SetBinError(ibin+1,h->GetBinError(ibin+1)/deltaE);
	h->SetBinContent(ibin+1,h->GetBinContent(ibin+1)/deltaE);
      }

  return h;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Simulate list of On and Off events
// If no background model has been provided, the use the distribution
// of Off events as model
// see also Iact1dUnbinnedLkl::SimulateDataSamples for more details
//
// Returns: 0 in case of success
//          1 otherwise
//
Int_t Iact1dBinnedLkl::SimulateDataSamples(UInt_t seed,Float_t meanG)
{
  // basic check
  if(fTauEDepFluct && !fdNdEpBkgFromOff)
    {
      cout << "Iact1dBinnedLkl::SimulateDataSamples Error: requested energy dependent uncertainty in tau but continuous background model has been provided" << endl;
      return 1;
    }

  //
  TH1F* provHNOff = NULL;
  if(!GetHdNdEpBkg())
    {
      BuildAndBinOnOffHistos();
      provHNOff = GetHdNdEpOff(kTRUE);
      if(!provHNOff->GetEntries())
	{
	  cout << "Iact1dBinnedLkl::SimulateDataSamples Error: missing Off dataset, cannot simulate anything" << endl;
	  return 1;
	}
      provHNOff->Scale(1./GetDataObsTime()/GetDataTau());

      TransformAndSavedNdEpBkg(provHNOff,kFALSE);
    }

  if(fHNOn)  delete fHNOn;
  if(fHNOff) delete fHNOff;
  fHNOn = fHNOff = NULL;
  
  return Iact1dUnbinnedLkl::SimulateDataSamples(seed,meanG);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Histograms hdNdEpBkg and hdNdEpSignalOff will contain entries obtained from those from fHdNdEpBkg and fHdNdEpSignalOff,
// respectively, each fluctuating independently according to the uncertainty in tau
// Note that fHdNdEpBkg is the expected distribution of background events in the On region
// and fHdNdEpSignalOff the expected distribution of signal events in the total Off region
// (total meaning that if tau=3 the effective area to consder is that of the three subregions)
//
Int_t Iact1dBinnedLkl::GetRealBkgAndGoffHistos(TRandom3* rdm,TH1F*& hdNdEpBkg,TH1F*& hdNdEpSignalOff) 
{
  // if fHdNdBkg was provided as input, use global normalization
  if(!fdNdEpBkgFromOff) return Iact1dUnbinnedLkl::GetRealBkgAndGoffHistos(rdm,hdNdEpBkg,hdNdEpSignalOff);

  // create the new histos with contents of the existing ones
  if(GetHdNdEpBkg()) hdNdEpBkg = new TH1F(*GetHdNdEpBkg());
  else
    {
      cout << "Iact1dBinnedLkl::GetRealBkgAndGoffHistos Warning: fHdNdEpBkg histo does not exist" << endl;
      return 1;
    }
  if(GetHdNdEpSignalOff()) hdNdEpSignalOff = new TH1F(*GetHdNdEpSignalOff());
  else hdNdEpSignalOff = NULL;

  // if no uncertainty in tau, that's all we need to do
  if(GetDTau()<=0) return 0;
  
  //  check the histograms have the same number of events
  if(hdNdEpSignalOff)
    if(hdNdEpBkg->GetNbinsX()!=hdNdEpSignalOff->GetNbinsX())
      {
	cout << "Iact1dBinnedLkl::GetRealBkgAndGoffHistos Error: realHdNdEpBkg and realHdNdEpSignalOff should have equal number of bins" << endl;
	return 1;
      }

  // unnormalize histos
  hdNdEpBkg->Scale(hdNdEpBkg->GetBinContent(0));
  hdNdEpBkg->SetBinContent(0,0);
  if(hdNdEpSignalOff)
    {
      hdNdEpSignalOff->Scale(hdNdEpSignalOff->GetBinContent(0));
      hdNdEpSignalOff->SetBinContent(0,0);
    }

  // make each bin fluctuate according to the value of fDTau
  SetTrueTau(rdm->Gaus(GetTau(),GetDTau()));
  Float_t binval = hdNdEpBkg->GetBinContent(1);
  for(Int_t ibin=0;ibin<hdNdEpBkg->GetNbinsX();ibin++)    
    {
      if(fTauEDepFluct)
	if(binval!=hdNdEpBkg->GetBinContent(ibin+1))
	  {
	    SetTrueTau(rdm->Gaus(GetTau(),GetDTau()));
	    binval = hdNdEpBkg->GetBinContent(ibin+1);
	  }

      if(GetTrueTau()>0)
	{
	  hdNdEpBkg->SetBinContent(ibin+1,hdNdEpBkg->GetBinContent(ibin+1)*GetTrueTau()/GetTau());
	  if(hdNdEpSignalOff)
	    hdNdEpSignalOff->SetBinContent(ibin+1,hdNdEpSignalOff->GetBinContent(ibin+1)*GetTrueTau()/GetTau());
	}
      else
	{
	  cout << "Iact1dBinnedLkl::GetRealBkgAndGoffHistos Error: negative or null tau value, we cannot work like that!" << endl;
	  return 1;
	}
    }

    
  NormalizedNdEHisto(hdNdEpBkg);
  if(hdNdEpSignalOff)
    NormalizedNdEHisto(hdNdEpSignalOff);
  
  return 0;
}

//////////////////////////////////////////////////////////////////
//
// Function called by Lkl::PrintOverview
// Print info about the experimental data
//
void Iact1dBinnedLkl::PrintData(Int_t level) 
{
  Iact1dUnbinnedLkl::PrintData(level);
  Margin(level); cout << "                # of bins = " << fNBins-fNRemovedBins << (fNRemovedBins? Form("(%d-%d)",fNBins,fNRemovedBins) : "") <<endl;
  Margin(level); cout << "                Rebinned? = " << (fMinBinContent>0? "YES" : "NO") 
		      << (fMinBinContent>0? Form(" (to %d min entries)",fMinBinContent) : "") <<endl;
}

//////////////////////////////////////////////////////////////////
// Fill newbin with optimal binning which contains at least minnevts
// in all bins of hOn and hOff, and fill inewbin with the number of
// bins
void GetRebinning(TH1F* hOn,TH1F* hOff,UInt_t minnevts,UInt_t& inewbin,Double_t* newbin)
{
  Int_t nibins=hOn->GetNbinsX();
  inewbin = 0;

  newbin[0] = hOn->GetBinLowEdge(1);

  for(Int_t ibin=0;ibin<nibins;ibin++,inewbin++)
    {
      Int_t non  = hOn->GetBinContent(ibin+1);
      Int_t noff = hOff->GetBinContent(ibin+1);
      while((non<minnevts || noff<minnevts) && ibin<nibins-1)
	{
	  ibin++;
	  non  += hOn->GetBinContent(ibin+1);
	  noff += hOff->GetBinContent(ibin+1);
	}
      
      newbin[inewbin+1] = hOn->GetBinLowEdge(ibin+1)+hOn->GetBinWidth(ibin+1);
    
      if((non<minnevts || noff<minnevts) && inewbin>1) // last bin does not comply with minimal statistics condition
	{
	  newbin[inewbin] = newbin[inewbin+1];
	  inewbin--;
	}      
    }
}


		
////////////////////////////////////////////////////////////////////////
// joint likelihood function (-2logL) for all bins
// To be minimized by TMinuit
// Free parameters:
// par[0] = g (total estimated number of signal events in all bins)
// par[1] = global tau for all bins (Off/On normalization)
//
void binnedIact1dUnbinnedLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;
  f = 0;

  // get internal object
  Iact1dBinnedLkl* mylkl          = dynamic_cast<Iact1dBinnedLkl*>(minuit->GetObjectFit());
  Bool_t          tauIsFixed    = mylkl->IsParFixed(Iact1dBinnedLkl::gTauIndex);
  Float_t         tau           = mylkl->GetTau();
  Float_t         dTau          = mylkl->GetDTau();
  
  // total g is the sum of g for all bins
  Double_t       g       = par[0];
  Double_t       trytau  = par[1];
  TObjArrayIter* iter    = (TObjArrayIter*) mylkl->GetSampleArray()->MakeIterator();
  PoissonLkl* sample;
  
  while((sample=dynamic_cast<PoissonLkl*>(iter->Next())))
    {
      sample->SetTau(trytau);

      Double_t   w_i = sample->GetUnitsOfG();
      Double_t   g_i = g*w_i;
     
      if(w_i>0)
	f+=sample->MinimizeLkl(g_i,kTRUE,kFALSE);
    }
  if(!tauIsFixed)
    f+= -2*TMath::Log(TMath::Gaus(trytau, tau, dTau, kTRUE));
  
  delete iter;
}
