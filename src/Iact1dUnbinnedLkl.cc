/* ======================================================================== *\
!
!
!   Author(s): Jelena Aleksic      06/2012 <mailto:jelena@ifae.es>
!              Javier Rico         12/2014 <mailto:jrico@ifae.es>
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
// Iact1dUnbinnedLkl
//
// Class to perform full likelihood maximization (minimization of -2logL)
// to estimate the presence of signal events following a certain spectral 
// shape. No physics origin is assumed and the class can be used as long as 
// one has an a priori knowledge of the signal spectral shape.
// 
// The class is based in the method described and characterize in
// Aleksic, Rico & Martinez JCAP 10 (2012) 032 [arXiv:1209.5589]
//
// Nomenclature, units and conventions:
// ------------------------------------
// E          = true energy
// E' (or Ep) = estimated energy
// Aeff       = effective area (vs E) for signal events in the ON region
// AeffOff    = effective area (vs E) for signal events in the OFF region (by default is zero but can be finite for extended sources)
//              IMPORTANT NOTE: is Aeff in the whole Off region
//              (e.g. if there are three off regions is the Aeff of the three considered together)
// g          = number of signal ray events in On region (free parameter during fit)
// b          = number of background events in On region (nuisance parameter)
// tau        = the ratio between Off and On exposures (e.g. number of Off regions)
//              (nuisance parameter, but can be fixed if provided error is 0) 
// Non        = number of selected On events after ALL cuts
// Noff       = number of selected Off events after ALL cuts
//
// E, E' are always in GeV
// Spatial dimentions always in cm
// Time always in s
// Therefore, e.g., Aeff is in cm^2 vs GeV, differential flux in GeV-1 cm-2 s-1, etc
//
// Basic Input:
// ------------
// fEpmin/fEpmax = lower/upper E' cuts, must be provided to the constructor
// fOnSample/fOffSample = list of On/Off event E' values, provided through 
//                        ReadDataSamplesSegueStereo or similar functions 
//                        (to be implemented for different input formats)
// fTau/fDTau = Off/On normalization and error, provided through 
//              ReadDataSamplesSegueStereo or similar functions 
//              (to be implemented for different input formats)
// fHAeff = histogram with Aeff vs logE for signal events (g), provided through
//          ReadAeffSegueStereo or similar functions (to be implemented for 
//          different input formats, e.g. ReadCTAIRF).
// fHAeffOff = histogram with AeffOff vs logE for signal events in the Off region, provided through
//             ReadAeffOffSegueStereo or similar functions (to be implemented for 
//             different input formats, e.g. ReadCTAIRF).
// fGEreso/fGEbias = graph with relative energy resolution/bias vs logE,  
//                   provided by ReadEResoAndBiasSegueStereo or similar functions 
//                   (to be implemented for different input formats). They are 
//                   overriden by fMigMatrix if the latter exists.
// fMigMatrix = migration matrix, normalized for fixed E. log10(E') is in X-axis,
//              log10(E) in Y-axis. It overrides fGEreso and fGEbias. 
//              provided by ReadCTAIRF or similar functions (to be implemented 
//              for different input formats). 
// fHdNdEpBkg = histogram with dN/dE'dt vs logE' for background events, 
//              provided through ReaddNdEpBkgSegueStereo or similar functions 
//              (to be implemented for different input formats)
//              The histogram is normalized and the integral stored in bin #0, which 
//              can be accessed through function GetdNdEpBkgIntegral. 
// fHdNdEpFrg = histogram with dN/dE'dt vs logE' for foreground events (i.e. those
//              produced by a nuisance gamma-ray source in the On region),
//              provided through ReaddNdEpFrgSegueStereo or similar functions 
//              (to be implemented for different input formats)
//              The histogram is normalized and the integral stored in bin #0, which 
//              can be accessed through function GetdNdEpFrgIntegral. 
// fHdNdEpSignal = histogram with convolution of dN/dE_signal*Aeff with energy dispersion 
//                 function provided through ReaddNdEpSignalSegueStereo or similar functions 
//                 (to be implemented for different input formats). Alternatively
//                 if fHdNdESignal is provided, it will be computed using fHAeff, plus
//                 fGEreso and fGEbias or fMigMatrix.
//                 The histogram is normalized and the integral stored in bin #0, which can 
//                 be accessed through function GetdNdEpSignalIntegral.
// fHdNdEpSignalOff = histogram with convolution of dN/dE_signal*AeffOff with energy dispersion 
//                    function provided through ReaddNdEpSignalOffSegueStereo or similar functions 
//                    (to be implemented for different input formats). Alternatively
//                    if fHdNdESignal is provided, it will be computed using fHAeffOff, plus
//                    fGEreso and fGEbias or fMigMatrix.
//                    The histogram is normalized and the integral stored in bin #0, which can 
//                    be accessed through function GetdNdEpSignalOffIntegral.
// fHdNdESignal = histogram with dN/dE vs logE for signal events
//                provided through ReaddNdESignalSegueStereo or similar functions 
//                (to be implemented for different input formats). 
//                The histogram is normalized and the integral stored in bin #0, which can 
//                be accessed through function GetdNdESignalIntegral. 
// 
// The previous histograms are used as functions (likelihood is
// unbinned). Therefore a large number of bins is needed. The number
// of bins, lower and upper logE/logE' bounds have default values:
// gNFineBins, gFineLEMin and gFineLEMax. Those can be configured by
// hand (before reading the input data) with SetNFineBins,
// SetFineLEMin and SetFineLEMax, but should cover at least the range
// of interest in E' (those of the selected events) and more (because
// the energy resolution and bias are finite). Do not change them
// unless you know what you're doing.
//
// The input (background rate, effectiva area, etc) can take different
// formats.  For each of those, the corresponding methods to transform
// them into the internal format must be provided. If you can't find
// support for your input format, implement it!
// 
// If we combine a set of Iact1dUnbinnedLkl objects representing observations
// of a given steady source into a JointLkl object (to minimize the 
// joint likelihood defined by the product of all added Iact1dUnbinnedLkl)
// make sure to call SetUnitsOfG for each sample, in such a way that
// g_i*fUnitsOfG_i is constant for all Iact1dUnbinnedLkl objects (i.e. fUnitsOfG
// is the factor to transform g into some common physical quantity).
//
// There is a built-in method (SetDMAnnihilationUnitsForG) to set
// fUnitsOfG in such a way that g*fUnitsOfG=<sv> (i.e. the thermal
// cross section for DM annihilation). For other physical constants of
// your interest, it is recommeded that you add a new method, so that
// the user cannot make mistakes in setting fUnitsOfG.
//
//
// Usage example:
// --------------
// (for a fully working example see macro jointLkl.C)
//
// Iact1dUnbinnedLkl* fLkl = new Iact1dUnbinnedLkl(Emin,Emax);
// fLkl->SetErrorDef(2.7);  // for 1-sided 95% CL
// fLkl->ReadAeffSegueStereo(aEffFileName);
// fLkl->ReadEResoAndBiasSegueStereo(energyRFileName);
// fLkl->ReaddNdEpBkgSegueStereo(bkgModFileName);
// fLkl->ReadDataSamplesSegueStereo(evtFileName,tau,dtau);
// fLkl->ReaddNdESignalSegueStereo(dNdESignalFileName);
//
// // Set units for DM interpretation of results
// fLkl->SetDMAnnihilationUnitsForG(Teff,mass,log10_J);
//
// // Compute the profile -2logL curve
// fLkl->ComputeLklVsG();
// 
// // Get the results
// Double_t gminval     = fLkl->GetGLklMin();    // value of g minimzing -2logL
// Double_t gminvalErr  = fLkl->GetGLklMinErr(); // 1-sided 95% CL error bar from Migrad
// Double_t glimval     = fLkl->GetGForLkl(2.7); // value of g for which -2logL = minimum+2.7
// Double_t gminvalErr2 = glimval-gminval;       // better determination of 95% CL error 
//
// // Plot the -2logL curve
// fLkl->GetLklVsG()->Draw("al");                // draw -2logL vs g 
//
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "TMath.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPRegexp.h"

#include "Iact1dUnbinnedLkl.h"
#include "IactEventListIrf.h"

ClassImp(Iact1dUnbinnedLkl);

using namespace std;

// static constants
static const TString  gName            = "Iact1dUnbinnedLkl";
static const TString  gTitle           = "Unbinned Full Likelihood";
static const Int_t    gNPars           = 3;                      // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"g","b","tau"};        // Name of parameters
static const Int_t    gNBins           = 100;                    // default number of histograms for dN/dE plots
static const Int_t    gNFineBins       = 5000;                   // default number of fine bins for internal histos
static const Double_t gFineLEMin       = TMath::Log10(10);       // default minimum log(energy[GeV]) for internal histos
static const Double_t gFineLEMax       = TMath::Log10(1000000);  // default maximum log(energy[GeV]) for internal histos
static const Double_t gCenterBin       = 0.5;                    // decide which value represents bin in histogram (= 0 for lower bin edge, 0.5 for the middle, 1 for the right edge)

// static functions (for internal processing of input data)
static Int_t  SmearHistogram(TH1F* sp,TH1F* smsp,TGraph* grreso,TGraph* grbias);
static Int_t  SmearHistogram(TH1F* sp,TH1F* smsp,TH2F* mm);
static void  readAndInterpolate(TH1F* ih,TH1F* oh,Double_t scale=0,Bool_t isDiff=kTRUE);
static Int_t copyBinByBin(TH1F* ih,TH1F* oh,Double_t scale=0,Bool_t isDiff=kTRUE);
static TH1F* GetResidualsHisto(TH1F* hModel,TH1F* hData);
static void  NormalizeMatrix(TH2F* im,TH2F* om,Double_t scale=0);
//static void  makeLogAndInterpolate(TH1D* ih,TH1F* oh);

// -2logL function for minuit
void fullLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;


//////////////////////////////////////////////////////////////////////////////
//
// String constructor
// The string contains the elements for the constructor 
//
Iact1dUnbinnedLkl::Iact1dUnbinnedLkl(TString inputString) : 
  Lkl(gNPars,inputString,gName,gTitle), fEpmin(gEpmin), fEpmax(gEpmax), 
  fNon(0), fNoff(0), fOnSample(NULL), fOffSample(NULL),
  fDataTau(1), fDataDTau(0), fDataObsTime(0), fTau(1), fDTau(0), fTrueTau(1), fObsTime(0),
  fLogJ(0), fIsOffAsOn(kFALSE),
  fNFineBins(gNFineBins), fFineLEMin(gFineLEMin), fFineLEMax(gFineLEMax),
  fHdNdESignal(NULL), fHAeff(NULL), fHAeffOff(NULL), fGEreso(NULL), fGEbias(NULL), fMigMatrix(NULL),
  fHdNdEpBkg(NULL), fHdNdEpFrg(NULL), fHdNdEpSignal(NULL), fHdNdEpSignalOff(NULL)
{
  if(InterpretInputString(inputString))
    cout << "Iact1dUnbinnedLkl::Iact1dUnbinnedLkl Warning: there were problems interpreting the input string" << endl;      
}

//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
// The string has been passed to Lkl::InterpretInputString in constructor
// here it is searched for the following options:
//
// logJ=<val>:          mean value of the log10 of the J factor, in units of GeV^2 cm^-5 (annihilation) or GeV cm^-3 (decay)
// path=<val>:          path of the input file (will be appended to inputFileName)
// inputfile=<val>:     name of the input file, which contain the IRFs and/or data, which is subsequently interpreted
// obsTime=<val>:       assumed observation time, in units of hour. Be aware observation time is stored in the data input file,
//                      but it can be useful to override it for tests and simulations
// tau=<val>            assumed tau value.  Be aware tau value is stored in the data input file,
//                      but it can be useful to override it for tests and simulations
// DtauStat=<val>       statistical uncertainty in tau (the Off/On normalization factor).  Be aware DtauStat value is stored 
//                      in the data input file, but it can be useful to override it for tests and simulations
// DtauSyst=<val>:      relative systematic uncertainty in tau (the Off/On normalization factor),
//                      e.g. if DtauSyst=0.015 an uncertaintyof 0.015*tau will be added in quadrature to the statistical
//                      uncertainty (normally stored in data file)
// fineEmin:            if provieded, overrides the value of fFineLEMin with TMath::Log10(fineEmin)
// fineEmax:            if provieded, overrides the value of fFineLEMax with TMath::Log10(fineEmax)
//
Int_t Iact1dUnbinnedLkl::InterpretInputString(TString inputString)
{
  // Prepeare to interpret inputString
  TPMERegexp re("\\s+");
  TPMERegexp fldre("=");
  TString path          = "";
  TString inputfileName = " (No file has been specified) ";
  Float_t relDTauSyst = 0;
  Float_t obstime = -1;
  Float_t tau=-1;
  Float_t dTauStat=-1;

  // split the inputString into the different fields, and check the option and values specified in each of them
  UInt_t nfields = re.Split(inputString);
  for(UInt_t ifield=0;ifield<nfields;ifield++)
    {
      fldre.Split(re[ifield]);
      TString optname = fldre[0];
      if(optname.CompareTo("logJ",TString::kIgnoreCase)==0)
	fLogJ=fldre[1].Atof();
      else if(optname.CompareTo("path",TString::kIgnoreCase)==0)
	path=fldre[1];      
      else if(optname.CompareTo("inputfile",TString::kIgnoreCase)==0)
	inputfileName=fldre[1];      
      else if(optname.CompareTo("obsTime",TString::kIgnoreCase)==0)
	obstime=fldre[1].Atof()*60*60;
      else if(optname.CompareTo("tau",TString::kIgnoreCase)==0)
	tau=fldre[1].Atof();
      else if(optname.CompareTo("DtauStat",TString::kIgnoreCase)==0)
	dTauStat=fldre[1].Atof();
      else if(optname.CompareTo("DtauSyst",TString::kIgnoreCase)==0)
	relDTauSyst=fldre[1].Atof();
      else if(optname.CompareTo("isoffason",TString::kIgnoreCase)==0)
	{
	  if (fldre[1].CompareTo("TRUE",TString::kIgnoreCase)==0)
	    fIsOffAsOn=kTRUE;
	}
      else if(optname.CompareTo("fineEmin",TString::kIgnoreCase)==0)
	fFineLEMin=TMath::Log10(fldre[1].Atof());
      else if(optname.CompareTo("fineEmax",TString::kIgnoreCase)==0)
	fFineLEMax=TMath::Log10(fldre[1].Atof());
    }
  
  // open and read input files with data and IRFs
  IactEventListIrf* dataSet = new IactEventListIrf("dataSet", "", path+(path==""?"":"/")+inputfileName);

  // extract info from file 
  fEpmin       = dataSet->GetEpmin();
  fEpmax       = dataSet->GetEpmax();
  fTrueTau     = fDataTau     = fTau      = dataSet->GetTau();
  fDataDTau    = fDTau     = TMath::Sqrt(TMath::Power(dataSet->GetDTau(),2)+TMath::Power(relDTauSyst*fTau,2));
  fDataObsTime = fObsTime  = dataSet->GetObsTime();
  fNon         = dataSet->GetOnSample()->GetEntries();
  fNoff        = dataSet->GetOffSample()->GetEntries();
  
  // extract data
  Double_t eventOnE,eventOffE;
  dataSet->SetOnBranchAddress("E",&eventOnE);
  dataSet->SetOffBranchAddress("E",&eventOffE);
 
  fOnSample  = new Float_t[fNon];
  fOffSample = new Float_t[fNoff];

  for(Int_t i=0;i<fNon;i++)
    {
      dataSet->GetOnEntry(i);
      fOnSample[i] = TMath::Log10(eventOnE);
    }
  for(Int_t i=0;i<fNoff;i++)
    {
      dataSet->GetOffEntry(i);
      fOffSample[i] = TMath::Log10(eventOffE);
    }
  
  // create IRFs object only if their correspondent in IactEventListIrf is non empty
  if(dataSet->GetHAeff()->GetEntries())
    if(SetAeff(dataSet->GetHAeff()))
      cout << "Iact1dUnbinnedLkl::Iact1dUnbinnedLkl Warning: problems setting AEff histogram" << endl;
  if(dataSet->GetHAeffOff()->GetEntries())
    if(SetAeffOff(dataSet->GetHAeffOff()))
      cout << "Iact1dUnbinnedLkl::Iact1dUnbinnedLkl Warning: problems setting AEffOff histogram" << endl;      
  if(dataSet->GetGEreso()->GetN() && dataSet->GetGEbias()->GetN())
    if(SetEResoAndBias(dataSet->GetGEreso(),dataSet->GetGEbias()))
      cout << "Iact1dUnbinnedLkl::Iact1dUnbinnedLkl Warning: problems setting Ereso and Ebias graphs" << endl;
  if(dataSet->GetMigMatrix()->GetEntries())
    if(SetMigMatrix(dataSet->GetMigMatrix()))
      cout << "Iact1dUnbinnedLkl::Iact1dUnbinnedLkl Warning: problems setting MigMatrix histogram" << endl;   
  if(dataSet->GetHdNdEpBkg()->GetEntries())
    if(SetdNdEpBkg(dataSet->GetHdNdEpBkg()))
      cout << "Iact1dUnbinnedLkl::Iact1dUnbinnedLkl Warning: problems setting dNdEpBkg histogram" << endl;   
  if(dataSet->GetHdNdEpFrg()->GetEntries())
    if(SetdNdEpFrg(dataSet->GetHdNdEpFrg()))
      cout << "Iact1dUnbinnedLkl::Iact1dUnbinnedLkl Warning: problems setting dNdEpFrg histogram" << endl;   

  // if observation time was specified, override the value in the input data file (for simulation of different obs time)
  if(obstime>0)
    fObsTime = obstime;
  
  // if tau was specified, override the value in the input data file (for simulation of different tau)
  if(tau>0)
    fTrueTau = fTau = tau;
  
  // if dTauStat was specified, override the value in the input data file (for simulation of different dTauStat)
  if(dTauStat>=0)
    fDTau  = TMath::Sqrt(TMath::Power(dTauStat,2)+TMath::Power(relDTauSyst*fTau,2));
  
  // sanity checks
  if(CheckEnergyLimits())
    cout << "Iact1dUnbinnedLkl::Iact1dUnbinnedLkl Warning: energy limits out of allowed bounds!" << endl;
  
  delete dataSet;  
  return 0;
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
Iact1dUnbinnedLkl::~Iact1dUnbinnedLkl()
{
  if(fOnSample)        delete [] fOnSample;
  if(fOffSample)       delete [] fOffSample;
  if(fHdNdESignal)     delete fHdNdESignal;
  if(fHAeff)           delete fHAeff;
  if(fGEreso)          delete fGEreso;
  if(fGEbias)          delete fGEbias;
  if(fMigMatrix)       delete fMigMatrix;
  if(fHdNdEpBkg)       delete fHdNdEpBkg;
  if(fHdNdEpFrg)       delete fHdNdEpFrg;
  if(fHdNdEpSignal)    delete fHdNdEpSignal;
  if(fHdNdEpSignalOff) delete fHdNdEpSignalOff;
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of TMinuit subtleties)
//
void Iact1dUnbinnedLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN fulllkl
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance),
// fix tau if requested
//
void Iact1dUnbinnedLkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(fullLkl);
  fMinuit->SetName(Form("%s_Minuit",GetName()));
  Double_t pStart[gNPars] = {ginit, fNoff/fTau, fTau};
  Double_t pDelta[gNPars] = {TMath::Sqrt(fNoff)/10.,TMath::Sqrt(fNoff)/10.,fDTau/10.};    // Precision of parameters during minimization

  SetParameters(gParName,pStart,pDelta);

  // Fix tau if requested (both in minuit and in lkl)
  fMinuit->Release(Iact1dUnbinnedLkl::gTauIndex);
  FixPar(Iact1dUnbinnedLkl::gTauIndex,kFALSE);
  if(GetDTau()<=0)
    { 
      fMinuit->FixParameter(Iact1dUnbinnedLkl::gTauIndex);
      FixPar(Iact1dUnbinnedLkl::gTauIndex);
    }
}		
		

////////////////////////////////////////////////////////////////
//
// Check that all elements needed for the fit are present, otherwise
// try to create/compute them
//
// Return 0 in case checks are Ok, 1 otherwise
//
Int_t Iact1dUnbinnedLkl::MakeChecks()
{
  if(IsChecked()) return 0;

  // check if all needed histos are there
  if(CheckHistograms())
    {
      cout << "Iact1dUnbinnedLkl::MakeChecks Warning: missing information, cannot perform fit, check your code!" << endl;
      return 1;
    }

  SetChecked();
  return 0;
}	      
			   
////////////////////////////////////////////////////////////////
//
// Check that the needed histograms are present and that we are ready 
// for calling minimization 
//
// If needed, convolute dNdESignal*Aeff with Eres and Ebias
// to get dNdEpSignal
//
// If needed, convolute dNdESignal*AeffOff with Eres and Ebias
// to get dNdEpSignalOff
//
// Return 0 in case of success, 1 otherwise
//
Int_t Iact1dUnbinnedLkl::CheckHistograms(Bool_t checkdNdEpBkg)
{
  // if fHdNdEpSignal is missing, try to construct it from fHdNdESignal, fHAeff fGEreso and fGEbias
  if(!fHdNdEpSignal && (fHdNdESignal && fHAeff && ((fGEreso && fGEbias) || fMigMatrix)))
    {
      if(fMigMatrix)
	cout << "Iact1dUnbinnedLkl::CheckHistograms Message: will create fHdNdEpSignal from fHdNdESignal, fHAeff & fMigMatrix... " << flush;
      else
	cout << "Iact1dUnbinnedLkl::CheckHistograms Message: will create fHdNdEpSignal from fHdNdESignal, fHAeff, fGEreso & fGEbias... " << flush;
      
      // multiply dNdESignal times Aeff
      TH1F* hdNdESignalAeff = new TH1F("hdNdESignalAeff","Product of dN/dE for Signal and Aeff",fNFineBins,fFineLEMin,fFineLEMax);
      hdNdESignalAeff->SetDirectory(0);      
      hdNdESignalAeff->Multiply(fHAeff,fHdNdESignal);

      // create fHdNdEpSignal   
      fHdNdEpSignal         = new TH1F("fHdNdEpSignal","dN/dE' for Signal",fNFineBins,fFineLEMin,fFineLEMax);
      fHdNdEpSignal->SetDirectory(0);

      // smear hdNdESignalAeff
      if(fMigMatrix)
	{
	  if(SmearHistogram(hdNdESignalAeff,fHdNdEpSignal,fMigMatrix))
	    return 1;
	}
      else
	{
	  if(SmearHistogram(hdNdESignalAeff,fHdNdEpSignal,fGEreso,fGEbias))
	    return 1;
	}	    

      cout << "Done! " << endl;
      // clean
      delete hdNdESignalAeff;
    }

  // if fHdNdEpSignalOff is missing but fHAeffOff exists, try to construct it from fHdNdESignal, fHAeffOff, fGEreso and fGEbias
  if(!fHdNdEpSignalOff && (fHdNdESignal && fHAeffOff && ((fGEreso && fGEbias) || fMigMatrix)))
    {
      if(fMigMatrix)
	cout << "Iact1dUnbinnedLkl::CheckHistograms Message: will create fHdNdEpSignalOff from fHdNdESignal, fHAeffOff & fMigMatrix" << endl;
      else
	cout << "Iact1dUnbinnedLkl::CheckHistograms Message: will create fHdNdEpSignalOff from fHdNdESignal, fHAeffOff, fGEreso & fGEbias" << endl;
      
      // multiply dNdESignal times AeffOff
      TH1F* hdNdESignalAeffOff = new TH1F("hdNdESignalAeffOff","Product of dN/dE for Signal and AeffOff",fNFineBins,fFineLEMin,fFineLEMax);
      hdNdESignalAeffOff->SetDirectory(0);
      hdNdESignalAeffOff->Multiply(fHAeffOff,fHdNdESignal);

      // create fHdNdEpSignalOff   
      fHdNdEpSignalOff      = new TH1F("fHdNdEpSignalOff","dN/dE' for Signal in Off region",fNFineBins,fFineLEMin,fFineLEMax);
      fHdNdEpSignalOff->SetDirectory(0);
      fHdNdEpSignalOff->SetLineColor(2);

      // smear hdNdESignalAeffOff
      if(fMigMatrix)
	{
	  if(SmearHistogram(hdNdESignalAeffOff,fHdNdEpSignalOff,fMigMatrix))
	    return 1;
	}
      else
	{
	  if(SmearHistogram(hdNdESignalAeffOff,fHdNdEpSignalOff,fGEreso,fGEbias))
	    return 1;
	}

      // clean
      delete hdNdESignalAeffOff;
    }
  
  // normalize unnormalized histos
  NormalizedNdEHisto(fHdNdESignal);
  NormalizedNdEHisto(fHdNdEpSignal);

  NormalizedNdEHisto(fHdNdEpSignalOff);

    
  if(checkdNdEpBkg)
    NormalizedNdEHisto(fHdNdEpBkg);
    
  if(fHdNdEpFrg)
    NormalizedNdEHisto(fHdNdEpFrg);
  
  
  // if there are the dNdE' histograms for signal and background + data we're ready to go
  if(checkdNdEpBkg)
    {
      if(fHdNdEpBkg && fHdNdEpSignal && fOnSample && fOffSample)
	return 0;
    }
  else
    {
      if(fHdNdEpSignal && fOnSample && fOffSample)
	return 0;
    }

  if(checkdNdEpBkg && !fHdNdEpBkg)
    cout << "Iact1dUnbinnedLkl::CheckHistograms Warning: fHdNdEpBkg histogram missing!!" << endl;
  if(!fHdNdEpSignal)
    cout << "Iact1dUnbinnedLkl::CheckHistograms Warning: fHdNdEpSignal histogram missing!!" << endl;
  if(!fOnSample)
    cout << "Iact1dUnbinnedLkl::CheckHistograms Warning: fOnSample histogram missing!!" << endl;
  if(!fOffSample)
    cout << "Iact1dUnbinnedLkl::CheckHistograms Warning: fOffSample histogram missing!!" << endl;

  return 1;
}

////////////////////////////////////////////////////////////////
//
// Normalize dN/dE' histos for signal and background and dN/dE for
// signal. 
// Save integral in bin 0.
// If bin 0 contains a non-zero value, do not normalize
//
// Return 0 in case of success
//        1 otherwise
//
Int_t Iact1dUnbinnedLkl::NormalizedNdEHisto(TH1F* histo)
{
// basic check
  if(!histo) return 1;
  
  // do not normalize if it's already normalized
  if(histo->GetBinContent(0)>0) return 0;

  // normalize and keep normalization in bin 0
  Double_t intSignal = IntegrateLogE(histo,TMath::Log10(fEpmin),TMath::Log10(fEpmax));
  histo->Scale(1./intSignal);
  histo->SetBinContent(0,intSignal);

  return 0;
}

////////////////////////////////////////////////////////////////
//
// Check energy limits are safe returning 0 (1) if they are (not)
//
Int_t Iact1dUnbinnedLkl::CheckEnergyLimits() const
{
  if(fEpmin<TMath::Power(10,fFineLEMin)) return 1;
  if(fEpmax<TMath::Power(10,fFineLEMin)) return 1;
  if(fEpmin>TMath::Power(10,fFineLEMax)) return 1;
  if(fEpmax>TMath::Power(10,fFineLEMax)) return 1;
  if(fEpmax<fEpmin) return 1;
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Set directly the dN/dE histogram for signal (e.g. from Damasco)
// Replacement of existing histo is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::SetdNdESignal(TH1F* hdNdESignal)
{
  // pathologies
  if(!hdNdESignal)
    {
      cout << "Iact1dUnbinnedLkl::SetdNdESignal Warning: input histo does not exist" << endl;
      return 1;
    }

  // replace existing fHdNdESignal and delete fHdNdEpSignal and fHdNdEpSignalOff
  if(fHdNdESignal) 
    delete fHdNdESignal; 
    
  if(fHdNdEpSignal)
    {
      delete fHdNdEpSignal;
      fHdNdEpSignal=NULL;
    }

  if(fHdNdEpSignalOff)
    {
      delete fHdNdEpSignalOff;
      fHdNdEpSignalOff=NULL;
    }
   
  // transform it to the Iact1dUnbinnedLkl format
  fHdNdESignal = new TH1F("fHdNdESignal","dN/dE for signal events",fNFineBins,fFineLEMin,fFineLEMax);
  fHdNdESignal->SetDirectory(0);
  fHdNdESignal->SetXTitle("log_{10}(E [GeV])");
  fHdNdESignal->SetYTitle("dN/dE [GeV^{-1}]");
  readAndInterpolate(hdNdESignal,fHdNdESignal);
  
  // clean and exit
  SetChecked(kFALSE);
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Add to a previously existing dN/dE histogram for signal from file
// in the Segue Stereo input format produced by Jelena
// No replacement of existing histo is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::AdddNdESignal(TString filename,Float_t br)
{
  // open file and look for histo
  TFile* dNdESignalFile  = new TFile(filename);
  TH1F*  hdNdESignal = (TH1F*) dNdESignalFile->Get("hdNdE");

  // pathologies
  if(!hdNdESignal)
    {
      cout << "Iact1dUnbinnedLkl::AdddNdESignal Warning: input histo does not exist" << endl;
      return 1;
    }

  // transform it to the Iact1dUnbinnedLkl format
  TH1F* fHdNdESignal_temp = new TH1F("fHdNdESignal_temp","dN/dE for signal events",fNFineBins,fFineLEMin,fFineLEMax);
  readAndInterpolate(hdNdESignal,fHdNdESignal_temp);

  for(Int_t ibin=0;ibin<fHdNdESignal->GetNbinsX();ibin++)
      fHdNdESignal->SetBinContent(ibin+1,fHdNdESignal->GetBinContent(ibin+1)+br*(fHdNdESignal_temp->GetBinContent(ibin+1)));

  // clean and exit
  delete dNdESignalFile;
  SetChecked(kFALSE);
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Recreate a fresh new version of fHdNdESignal
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::ResetdNdESignal()
{  
  // Delete existing fHdNdESignal and create empty one
  if(fHdNdESignal) 
    delete fHdNdESignal; 
  
  // Create histo
  fHdNdESignal = new TH1F("fHdNdESignal","dN/dE for signal events",fNFineBins,fFineLEMin,fFineLEMax);
  fHdNdESignal->SetDirectory(0);
  fHdNdESignal->SetXTitle("log_{10}(E [GeV])");
  fHdNdESignal->SetYTitle("dN/dE [GeV^{-1}]");
  
 // Delete existing fHdNdEpSignal and fHdNdEpSignalOff 
  if(fHdNdEpSignal)
    {
      delete fHdNdEpSignal;
      fHdNdEpSignal=NULL;
    }
  
  if(fHdNdEpSignalOff)
    {
      delete fHdNdEpSignalOff;
      fHdNdEpSignalOff=NULL;
    }
  SetChecked(kFALSE);
  
  // exit
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Create and save a dN/dE histogram for signal according to 
// given <function> and parameters <p0>, <p1>, ...
// (see Iact1dUnbinnedLkl::AdddNdESignalFunction for details on available
// functions)
//
// If fdNdESignal exists it will be replaced
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::SetdNdESignalFunction(TString function,Float_t p0,Float_t p1,Float_t p2,Float_t p3,Float_t p4,Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{  
  // Delete existing fHdNdESignal and create empty one
  ResetdNdESignal();

  // call to add the function
  AdddNdESignalFunction(function,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);
  
  // exit
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Add to a previously existing dN/dE histogram for signal according to 
// given <function> and parameters <p0>, <p1>, ...
//
// If <function>=="line"
// <p0> = E0, the energy of the line
// <p1> = scale, i.e. the number of photons per annihilation (e.g 2 for XX->gamma gamma)
// <p2> = br, the branching ratio
//
// If <function>=="box"
// <p0> = Emin, lower border of the box
// <p1> = Emax, upper border of the box
// <p2> = scale, i.e. the number of photons per annihilation (e.g 4 for XX->pi0 pi0)
// <p3> = br, the branching ratio
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::AdddNdESignalFunction(TString function,Float_t p0,Float_t p1,Float_t p2,Float_t p3,Float_t p4,Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{    
  // Check that fHdNdESignal exists
  if(!fHdNdESignal)
    {
      cout << "Iact1dUnbinnedLkl::AdddNdESignalFunction Error: fHdNdESignal does not exist" << endl;
      return 1;
    }

  // Check function
  enum function_t {gLine=0,gBox};
  Int_t functionType=-1;
  if(function.CompareTo("line",TString::kIgnoreCase)==0)
    functionType=gLine;
  else if(function.CompareTo("box",TString::kIgnoreCase)==0)
    functionType=gBox;
  else
    {
      cout << "Iact1dUnbinnedLkl::AdddNdESignalFunction Error: no known function type (" << function << "), check documentation for allowed values" << endl;
      return 1;
    }
    
  // Add to the histogram according to specified function and parameters
  if(functionType==gLine)
    {
      Float_t E0    = p0;
      Float_t scale = p1;
      Float_t br    = p2;
      
      Float_t log10E0 = TMath::Log10(E0);
      
      // check range
      if(fFineLEMin>log10E0 || fFineLEMax<log10E0)
	{
	  cout << "Iact1dUnbinnedLkl::AdddNdESignal Warning: trying to set a line at an energy E0="
	       << E0 << ", out of range (" << TMath::Power(10,fFineLEMin) << "," << TMath::Power(10,fFineLEMax) << ")" << endl;
	  return 1;
	}
         
      // fill histo
      Int_t ibin   = fNFineBins*(log10E0-fFineLEMin)/(fFineLEMax-fFineLEMin);
      Float_t Emin = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibin+1));
      Float_t Emax = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibin+1)+fHdNdESignal->GetBinWidth(ibin+1));
      Float_t dE   = Emax-Emin;
      fHdNdESignal->SetBinContent(ibin+1,fHdNdESignal->GetBinContent(ibin+1)+br*(scale/dE));
    }
  else if(functionType==gBox)
    {
      Float_t Emin  = p0;
      Float_t Emax  = p1;
      Float_t scale = p2;
      Float_t br    = p3;
      
      Float_t log10Emin = TMath::Log10(Emin);
      Float_t log10Emax = TMath::Log10(Emax);
      
      // check range
      if(Emin>Emax)
	{
	  cout << "Iact1dUnbinnedLkl::AdddNdESignal Warning: Emin (" << Emin << ") larger than Emax (" << Emax << ")" << endl;
	  return 1;
	}
         
      // fill histo
      Int_t ibinmin = fNFineBins*(log10Emin-fFineLEMin)/(fFineLEMax-fFineLEMin);
      Int_t ibinmax = fNFineBins*(log10Emax-fFineLEMin)/(fFineLEMax-fFineLEMin);
      
      Float_t realEmin  = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibinmin+1));
      Float_t realEmax  = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibinmax+1)+fHdNdESignal->GetBinWidth(ibinmax+1));
      Float_t dE        = realEmax-realEmin;
      for(Int_t ibin=ibinmin;ibin<=ibinmax;ibin++)
	fHdNdESignal->SetBinContent(ibin+1,fHdNdESignal->GetBinContent(ibin+1)+br*(scale/dE));
    }
  else
    {
      cout << "Iact1dUnbinnedLkl::AdddNdESignalFunction Warning: this should never happen... " << endl;
      return 1;
    }
      
  // exit
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
// If fdNdESignal exists it will be replaced
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::SetdNdESignalFunction(TF1* function,Float_t emin,Float_t emax)
{  
  // Delete existing fHdNdESignal and create empty one
  ResetdNdESignal();

  // call to add the function
  AdddNdESignalFunction(function,emin,emax,1.0);
  
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
Int_t Iact1dUnbinnedLkl::AdddNdESignalFunction(TF1* function,Float_t emin,Float_t emax,Float_t br)
{
  // Check that fHdNdESignal exists
  if(!fHdNdESignal)
    {
      cout << "Iact1dUnbinnedLkl::AdddNdESignalFunction Error: fHdNdESignal does not exist" << endl;
      return 1;
    }

  // Check function
  if(!function)
    {
      cout << "Iact1dUnbinnedLkl::AdddNdESignalFunction Error: function does not exist" << endl;
      return 1;
    }
  
  for(Int_t ibin=0;ibin<fHdNdESignal->GetNbinsX();ibin++)
    {
      Float_t etest = TMath::Power(10,fHdNdESignal->GetBinLowEdge(ibin+1)+gCenterBin*fHdNdESignal->GetBinWidth(ibin+1));
      if(etest>emin && etest<emax)
      fHdNdESignal->SetBinContent(ibin+1,fHdNdESignal->GetBinContent(ibin+1)+br*(function->Eval(etest)));
    }
  
  // exit
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Read dN/dE for signal from file in the Segue Stereo input format
// produced by Jelena
// Replacement of existing file is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::ReaddNdESignal(TString filename)
{
  // open file and look for histo
  TFile* dNdESignalFile  = new TFile(filename);
  TH1F*  hProvdNdESignal = (TH1F*) dNdESignalFile->Get("hdNdE");

  Int_t status = 0;
  if(!hProvdNdESignal)
    status = 1;
  else
    status = SetdNdESignal(hProvdNdESignal);
  delete dNdESignalFile;
  return status;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Read dN/dE' for signal from file in the Segue Stereo input format
// produced by Jelena
// Replacement of existing file is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::ReaddNdEpSignal(TString filename)
{
  // open file and look for histo
  TFile* dNdEpSignalFile = new TFile(filename);
  TH1F*  hNdEpSignalProv =    (TH1F*) dNdEpSignalFile->Get("fHdNdEpSignal");
  if(!hNdEpSignalProv) {delete dNdEpSignalFile; return 1;}

  // copy histo 
  if(fHdNdEpSignal) delete fHdNdEpSignal; 
  fHdNdEpSignal = new TH1F(*hNdEpSignalProv);
  fHdNdEpSignal->SetXTitle("log_{10}(E' [GeV])");
  fHdNdEpSignal->SetYTitle("cm^{2} dN/dE' [GeV^{-1}]");
  fHdNdEpSignal->SetDirectory(0);

  // check what we've read is good
  Int_t status = 0;
  Int_t nbins = fHdNdEpSignal->GetNbinsX();
  Float_t xmin = fHdNdEpSignal->GetXaxis()->GetXmin();
  Float_t xmax = fHdNdEpSignal->GetXaxis()->GetXmax();
  if(nbins!=fNFineBins || xmin!=fFineLEMin || xmax!=fFineLEMax)
    {
      cout << "Iact1dUnbinnedLkl::ReaddNdEpSignal, histogram read from file " << filename << " has " << nbins << " (should be " << fNFineBins << "), xmin = " << xmin << " (should be " << fFineLEMin << ") and xmax = " << xmax << " (should be " << fFineLEMax << ")" << endl;
      status = 1;
    }
     
  // clean and exit
  SetChecked(kFALSE);
  delete dNdEpSignalFile;
  return status;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Read dN/dE' for signal in the off region from file in the Segue Stereo input format
// produced by Jelena
// Replacement of existing file is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::ReaddNdEpSignalOff(TString filename)
{
  // open file and look for histo
  TFile* dNdEpSignalOffFile = new TFile(filename);
  TH1F*  hNdEpSignalOffProv =    (TH1F*) dNdEpSignalOffFile->Get("fHdNdEpSignalOff");
  if(!hNdEpSignalOffProv) {delete dNdEpSignalOffFile; return 1;}

  // copy histo 
  if(fHdNdEpSignalOff) delete fHdNdEpSignalOff; 
  fHdNdEpSignalOff = new TH1F(*hNdEpSignalOffProv);
  fHdNdEpSignalOff->SetXTitle("log_{10}(E' [GeV])");
  fHdNdEpSignalOff->SetYTitle("dN/dE' [cm^{2} GeV^{-1}]");
  fHdNdEpSignalOff->SetDirectory(0);

    // check what we've read is good
  Int_t status = 0;
  Int_t nbins = fHdNdEpSignalOff->GetNbinsX();
  Float_t xmin = fHdNdEpSignalOff->GetXaxis()->GetXmin();
  Float_t xmax = fHdNdEpSignalOff->GetXaxis()->GetXmax();
  if(nbins!=fNFineBins || xmin!=fFineLEMin || xmax!=fFineLEMax)
    {
      cout << "Iact1dUnbinnedLkl::ReaddNdEpSignalOff, histogram read from file " << filename << " has " << nbins << " (should be " << fNFineBins << "), xmin = " << xmin << " (should be " << fFineLEMin << ") and xmax = " << xmax << " (should be " << fFineLEMax << ")" << endl;
      status=1;
    }

  // clean and exit
  SetChecked(kFALSE);
  delete dNdEpSignalOffFile;
  return status;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Save dNdHpBkg histogram in gLike format
// Return 0 in case of success
//        1 otherwise
//
Int_t Iact1dUnbinnedLkl::TransformAndSavedNdEpBkg(TH1F* hProvdNdEBkg,Bool_t interpolate,Double_t scale,Bool_t isDiff)
{
  if(!hProvdNdEBkg) return 1;
  
  // transform it to the Iact1dUnbinnedLkl format
  if(fHdNdEpBkg) delete fHdNdEpBkg;
  fHdNdEpBkg = new TH1F("fHdNdEpBkg","dN/dE'dt for background events",fNFineBins,fFineLEMin,fFineLEMax);
  fHdNdEpBkg->SetDirectory(0);

  if(interpolate)
    readAndInterpolate(hProvdNdEBkg,fHdNdEpBkg,scale,isDiff);
  else
    {
      if(copyBinByBin(hProvdNdEBkg,fHdNdEpBkg,scale,isDiff))
	{
	  cout << "Iact1dUnbinnedLkl::TransformAndSavedNdEpBkg error copying the histogram bin by bin" << endl;
	  return 1;
	}
    }

  // labels
  fHdNdEpBkg->SetXTitle("log_{10}(E' [GeV])");
  fHdNdEpBkg->SetYTitle("dN/dE'dt [GeV^{-1} s^{-1}]");
  fHdNdEpBkg->SetStats(0);

  NormalizedNdEHisto(fHdNdEpBkg);
  SetChecked(kFALSE);

  return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Save dNdHpFrg histogram in gLike format
// Return 0 in case of success
//        1 otherwise
//
Int_t Iact1dUnbinnedLkl::TransformAndSavedNdEpFrg(TH1F* hProvdNdEFrg,Bool_t interpolate,Double_t scale,Bool_t isDiff)
{
  if(!hProvdNdEFrg) return 1;
  
  // transform it to the Iact1dUnbinnedLkl format
  if(fHdNdEpFrg) delete fHdNdEpFrg;
  fHdNdEpFrg = new TH1F("fHdNdEpFrg","dN/dE'dt for foreground events",fNFineBins,fFineLEMin,fFineLEMax);
  fHdNdEpFrg->SetDirectory(0);
  if(interpolate)
    readAndInterpolate(hProvdNdEFrg,fHdNdEpFrg,scale,isDiff);
  else
    if(copyBinByBin(hProvdNdEFrg,fHdNdEpFrg,scale,isDiff))
      {
	cout << "Iact1dUnbinnedLkl::TransformAndSavedNdEpFrg error copying the histogram bin by bin" << endl;
	return 1;
      }
		 
  // labels
  fHdNdEpFrg->SetXTitle("log_{10}(E' [GeV])");
  fHdNdEpFrg->SetYTitle("dN/dE'dt [GeV^{-1} s^{-1}]");
  fHdNdEpFrg->SetStats(0);
    
  NormalizedNdEHisto(fHdNdEpFrg);
  SetChecked(kFALSE);

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Transform and Set hProvAeff as the Aeff histogram 
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::SetAeff(TH1F* hProvAeff)
{
  if(fHAeff) delete fHAeff;
  fHAeff = new TH1F("fHAeff","Effective area",fNFineBins,fFineLEMin,fFineLEMax);
  fHAeff->SetDirectory(0);
  readAndInterpolate(hProvAeff,fHAeff);

  // configure
  fHAeff->SetMinimum(1e4);
  fHAeff->SetXTitle("log_{10}(E [GeV])");
  fHAeff->SetYTitle("Aeff [cm^{2}]");
  fHAeff->SetStats(0);

  // clean and exit
  SetChecked(kFALSE);
  return 0;
}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Transform and Set hProvAeff as the AeffOff histogram 
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::SetAeffOff(TH1F* hProvAeff)
{  
  if(fHAeffOff) delete fHAeffOff;
  fHAeffOff = new TH1F("fHAeffOff","Effective area for signal events in the Off region",fNFineBins,fFineLEMin,fFineLEMax);
  fHAeffOff->SetDirectory(0);
  readAndInterpolate(hProvAeff,fHAeffOff);
  
  // configure
  fHAeffOff->SetMinimum(1e4);
  fHAeffOff->SetXTitle("log_{10}(E [GeV])");
  fHAeffOff->SetYTitle("cm^{2}");
  fHAeffOff->SetLineColor(2);
  fHAeffOff->SetStats(0);

  // clean and exit
  SetChecked(kFALSE);
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Set energy resolution and bias graphs
// Return 0 in case of success
//        1 otherwise
//
Int_t Iact1dUnbinnedLkl::SetEResoAndBias(TGraph* ereso,TGraph* ebias)
{
  if(!ereso || !ebias)
    {
      cout << "Iact1dUnbinnedLkl::SetEResoAndBias error: trying to set a null graph" << endl;
      return 1;
    }

  if(fGEreso) delete fGEreso;
  if(fGEbias) delete fGEbias;

  fGEreso = new TGraph(ereso->GetN(),ereso->GetX(),ereso->GetY());
  fGEbias = new TGraph(ebias->GetN(),ebias->GetX(),ebias->GetY());

  fGEreso->SetName(ereso->GetName());
  fGEbias->SetName(ebias->GetName());
  
  fGEbias->SetLineStyle(2);

  SetChecked(kFALSE);
  return 0;
}
  
  
///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Read Aeff, migration matrix and dN/dE'_Bkg from IFAE CTA prod 2 file
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t Iact1dUnbinnedLkl::ReadCTAIRF(TString filename)
{
  const Float_t scale    = 3;
  fEpmin = 10;
  fEpmax = 2e5;
  // open file and look for histos
  TFile* ctaIRFFile    = new TFile(filename);
  TH2F*  provMM        = (TH2F*) ctaIRFFile->Get("MigMatrix");
  TH1F*  hProvdNdEpBkg = (TH1F*) ctaIRFFile->Get("BGRate");
  TH1F*  hProvAeff     = (TH1F*) ctaIRFFile->Get("EffectiveAreaEtrue");

  // check if all histos are found
  Bool_t isOk = kTRUE;
  if(!provMM) 
    {
      cout << "Iact1dUnbinnedLkl::ReadCTAIRF Warning: no Migration Matrix found" << endl;
      isOk = kFALSE;
    }
  else
    {
      // Save migration matrix in Iact1dUnbinnedLkl format
      ///////////////////////////////////////////
      
      // get binning and transform it to log(E [GeV])
      const Int_t   nbinslep = provMM->GetNbinsX(); // number of bins in logE'
      const Float_t lepmin   = provMM->GetXaxis()->GetXmin()+scale;
      const Float_t lepmax   = provMM->GetXaxis()->GetXmax()+scale;
      
      const Int_t   nbinsle  = provMM->GetNbinsY(); // number of bins in logE
      const Float_t lemin    = provMM->GetYaxis()->GetXmin()+scale;
      const Float_t lemax    = provMM->GetYaxis()->GetXmax()+scale;
      
      // normalize and save
      TH2F* provMigMatrix = new TH2F("provMigMatrix","Migration Matrix",nbinslep,lepmin,lepmax,nbinsle,lemin,lemax);
      provMigMatrix->SetDirectory(0);
      NormalizeMatrix(provMM,provMigMatrix,scale);
      if(SetMigMatrix(provMigMatrix))
	isOk = kFALSE;
    }
  if(TransformAndSavedNdEpBkg(hProvdNdEpBkg,kFALSE,scale,kFALSE)) 
    {
      cout << "Iact1dUnbinnedLkl::ReadCTAIRF Warning: no Background Rate histo found" << endl;
      isOk = kFALSE;
    }
  
  if(!hProvAeff) 
    {
      cout << "Iact1dUnbinnedLkl::ReadCTAIRF Warning: no Effectiva Area histo found" << endl;
      isOk = kFALSE;
    }
  else
    {
      // Save effective area in Iact1dUnbinnedLkl format
      ///////////////////////////////////////////
      // transform it to the Iact1dUnbinnedLkl format
      hProvAeff->Scale(1e4); // convert it to cm^2
      if(fHAeff) delete fHAeff;
      fHAeff = new TH1F("fHAeff","Effective area",fNFineBins,fFineLEMin,fFineLEMax);
      fHAeff->SetDirectory(0);
      copyBinByBin(hProvAeff,fHAeff,3);
      
      // configure
      fHAeff->SetMinimum(1e4);
      fHAeff->SetXTitle("log_{10}(E [GeV])");
      fHAeff->SetYTitle("cm^{2}");
      fHAeff->SetStats(0);
    }

  SetChecked(kFALSE);
  
  // Clean and exit
  //////////////////
  delete ctaIRFFile;
  return !isOk;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Set migration matrix
//
// Return 0 in case of success
//        1 otherwise
//
Int_t Iact1dUnbinnedLkl::SetMigMatrix(TH2F* provMM)
{
  // sanity check
  if(!provMM) return 1;
  
  // create matrix
  if(fMigMatrix) delete fMigMatrix;
  fMigMatrix = new TH2F(*provMM);
  fMigMatrix->SetDirectory(0);
					       		
  // configure histo
  fMigMatrix->SetName("fMigMatrix");
  fMigMatrix->SetXTitle("log_{10}(E' [GeV])");
  fMigMatrix->SetYTitle("log_{10}(E [GeV])");
  fMigMatrix->SetStats(0);
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Simulate list of On and Off events
// according to the total pdf described by fHdNdEpBkg and fHdNdEpSignal
// and the observation time in fObsTime and tau taking randomly from a
// gaussian of mean fTau and mean fDTau
//
// meanGwithUnits (default 0)    = assumed value of g (in units of fUnitsOfG)
// rdm            (default NULL) = random generator (if NULL take gRandom)
//
// IF meanGwithUnits<0, do not simulate independent ON events, use OFF sample
// also as ON
//
// Return 0 in case of success
//        1 otherwise
//
Int_t Iact1dUnbinnedLkl::SimulateDataSamples(Float_t meanGwithUnits,TRandom* rdm)
{
  if(meanGwithUnits<0) meanGwithUnits=0;
  
  // Sanity checks
  if(!fHdNdEpBkg)
    {
      cout << endl << "Iact1dUnbinnedLkl::SimulateDataSamples Error: fHdNdEpBkg missing, cannot simulate events" << endl;
      return 1;
    }
  if(meanGwithUnits>0 && !fHdNdEpSignal)
    {
      cout << endl << "Iact1dUnbinnedLkl::SimulateDataSamples Error: you want to simulate signal events but fHdNdEpSignal is missing, cannot simulate events" << endl;
      return 1;
    }
  if(fIsOffAsOn && (fTau!=1 || fDTau!=0))
    {
      cout << endl << "Iact1dUnbinnedLkl::SimulateDataSamples Error: the OnAsOff test only makes sense for fTau=1 && fDTau=0, but values are fTau="<<fTau<<", fDTau="<<fDTau<<", no events generated" << endl;
      return 1;
    }

  // compute weights for different pdf components
  if(!rdm) rdm = gRandom;
  TRandom*  saverdm = gRandom;
  gRandom = rdm;

  TH1F* realHdNdEpBkg       = NULL;
  TH1F* realHdNdEpSignalOff = NULL;

  if(GetRealBkgAndGoffHistos(rdm,realHdNdEpBkg,realHdNdEpSignalOff)) return 1;


  Float_t meanG    = meanGwithUnits/GetUnitsOfG();   // remove units to compute expected number of signal events
  Float_t meanB    = GetdNdEpBkgIntegral()*fObsTime;  
  Float_t meanBoff = realHdNdEpBkg->GetBinContent(0)*fObsTime*fTau;
  Float_t meanF    = GetdNdEpFrgIntegral()*fObsTime;
  Float_t meanGoff = ((fHdNdEpSignal && realHdNdEpSignalOff)? meanG*realHdNdEpSignalOff->GetBinContent(0)/GetdNdEpSignalIntegral() : 0); 
  Float_t meanNon  = meanB+meanF+meanG;
  Float_t meanNoff = meanBoff+meanGoff;

  // setup histogram to build pdfs
  TH1F* hBkgBinIntegrated     = new TH1F("hBkgBinIntegrated",     "Histogram for Off background event generation",          fNFineBins,fFineLEMin,fFineLEMax);
  TH1F* hRealBkgBinIntegrated = new TH1F("hRealBkgBinIntegrated", "Histogram for  On background event generation",          fNFineBins,fFineLEMin,fFineLEMax);
  TH1F* hFrgBinIntegrated     = new TH1F("hFrgBinIntegrated",     "Histogram for foreground event generation",              fNFineBins,fFineLEMin,fFineLEMax);
  TH1F* hSigBinIntegrated     = new TH1F("hSigBinIntegrated",     "Histogram for signal event generation",                  fNFineBins,fFineLEMin,fFineLEMax);
  TH1F* hSigOffBinIntegrated  = new TH1F("hSigOffBinIntegrated",  "Histogram for signal event generation in the Off region",fNFineBins,fFineLEMin,fFineLEMax);
  TH1F* hOnBinIntegrated      = new TH1F("hOnBinIntegrated",      "Histogram for  On event generation",                     fNFineBins,fFineLEMin,fFineLEMax);
  TH1F* hOffBinIntegrated     = new TH1F("hOffBinIntegrated",     "Histogram for Off event generation",                     fNFineBins,fFineLEMin,fFineLEMax);
  
  hBkgBinIntegrated->SetDirectory(0);
  hRealBkgBinIntegrated->SetDirectory(0);
  hFrgBinIntegrated->SetDirectory(0);
  hSigBinIntegrated->SetDirectory(0);
  hSigOffBinIntegrated->SetDirectory(0);
  hOnBinIntegrated->SetDirectory(0);
  hOffBinIntegrated->SetDirectory(0);

  // build pdf
  for(Int_t ibin = 0;ibin<fNFineBins;ibin++)
    {
      Double_t leminbin = hBkgBinIntegrated->GetBinLowEdge(ibin+1);
      Double_t lemaxbin = leminbin+hBkgBinIntegrated->GetBinWidth(ibin+1);
      if(TMath::Power(10,leminbin) < fEpmax && TMath::Power(10,lemaxbin) > fEpmin)
	{
	  Float_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin); 
	  hBkgBinIntegrated->SetBinContent(ibin+1,fHdNdEpBkg->GetBinContent(ibin+1)*deltaE);
	  hRealBkgBinIntegrated->SetBinContent(ibin+1,realHdNdEpBkg->GetBinContent(ibin+1)*deltaE);
	  if(fHdNdEpSignal)
	    hSigBinIntegrated->SetBinContent(ibin+1,fHdNdEpSignal->GetBinContent(ibin+1)*deltaE);
	  else
	    hSigBinIntegrated->SetBinContent(ibin+1,0);
	  if(fHdNdEpSignalOff)
	    hSigOffBinIntegrated->SetBinContent(ibin+1,realHdNdEpSignalOff->GetBinContent(ibin+1)*deltaE);
	  else
	    hSigOffBinIntegrated->SetBinContent(ibin+1,0);
	  if(fHdNdEpFrg)
	    hFrgBinIntegrated->SetBinContent(ibin+1,fHdNdEpFrg->GetBinContent(ibin+1)*deltaE);
	  else
	    hFrgBinIntegrated->SetBinContent(ibin+1,0);
	}
      else
	{
	  hBkgBinIntegrated->SetBinContent(ibin+1,0);
	  hRealBkgBinIntegrated->SetBinContent(ibin+1,0);
	  hSigBinIntegrated->SetBinContent(ibin+1,0);
	  hSigOffBinIntegrated->SetBinContent(ibin+1,0);
	  hFrgBinIntegrated->SetBinContent(ibin+1,0);
	}
    }

  // sum histos with their proper weight
  // on histogram
  hOnBinIntegrated->Add(hBkgBinIntegrated,meanB/hBkgBinIntegrated->Integral());
  if(fHdNdEpSignal)
    hOnBinIntegrated->Add(hOnBinIntegrated,hSigBinIntegrated,1,meanG/hSigBinIntegrated->Integral());
  if(fHdNdEpFrg)
    hOnBinIntegrated->Add(hOnBinIntegrated,hFrgBinIntegrated,1,meanF/hFrgBinIntegrated->Integral());

  // off histogram
  hOffBinIntegrated->Add(hRealBkgBinIntegrated,meanBoff/hRealBkgBinIntegrated->Integral());
  if(fHdNdEpSignalOff)
    hOffBinIntegrated->Add(hOffBinIntegrated,hSigOffBinIntegrated,1,meanGoff/hSigOffBinIntegrated->Integral());

  // simulate events    
  UInt_t noff = rdm->Poisson(meanNoff);
  UInt_t non  = (fIsOffAsOn? noff :rdm->Poisson(meanNon));
  vector<Float_t> offEnergy;
  vector<Float_t> onEnergy;
  offEnergy.reserve(noff);
  onEnergy.reserve(non);
  
  for(UInt_t ioff=0;ioff<noff;ioff++)
    offEnergy.push_back(hOffBinIntegrated->GetRandom());
  for(UInt_t ion=0;ion<non;ion++)
    if(fIsOffAsOn)
      onEnergy.push_back(offEnergy[ion]);
    else
      onEnergy.push_back(hOnBinIntegrated->GetRandom());
      
  cout << "Iact1dUnbinnedLkl::SimulateDataSamples Message: Generated " << non+noff << " events: [" << non << " On and " << noff 
       << " Off, TrueTau = " << fTrueTau << " (" << fTau << "+/-" << fDTau  << ")]" << endl;

  // save data as data members
  fNon  = onEnergy.size();
  fNoff = offEnergy.size();
  if(fOnSample)  delete [] fOnSample;
  if(fOffSample) delete [] fOffSample;
  fOnSample = new Float_t[fNon];
  fOffSample = new Float_t[fNoff];
  memcpy(fOnSample,onEnergy.data(),sizeof(Float_t)*fNon);
  memcpy(fOffSample,offEnergy.data(),sizeof(Float_t)*fNoff);
  
  SetChecked(kFALSE);
  
  // clean and exit
  gRandom = saverdm;

  delete hBkgBinIntegrated;
  delete hRealBkgBinIntegrated;
  delete hFrgBinIntegrated;
  delete hSigBinIntegrated;
  delete hSigOffBinIntegrated;
  delete hOnBinIntegrated;
  delete hOffBinIntegrated;
  
  if(realHdNdEpBkg)
    delete realHdNdEpBkg;
  if(realHdNdEpSignalOff)
    delete realHdNdEpSignalOff;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Histograms hdNdEpBkg and hdNdEpSignalOff will contain same entries as fHdNdEpBkg and fHdNdEpSignalOff,
// respectively, but with normalization fluctuating according to the uncertainty in tau
// Note that fHdNdEpBkg is the expected distribution of background events in the On region
// and fHdNdEpSignalOff the expected distribution of signal events in the total Off region
// (total meaning that if tau=3 the effective area to consder is that of the three subregions)
//
Int_t Iact1dUnbinnedLkl::GetRealBkgAndGoffHistos(TRandom* rdm,TH1F*& hdNdEpBkg,TH1F*& hdNdEpSignalOff) 
{
  // create new histos with contents of the existing ones
  if(fHdNdEpBkg) hdNdEpBkg = new TH1F(*fHdNdEpBkg);
  else
    {
      cout << "Iact1dUnbinnedLkl::GetRealBkgAndGoffHistos Warning: fHdNdEpBkg histo does not exist" << endl;
      return 1;
    }
  if(fHdNdEpSignalOff) hdNdEpSignalOff = new TH1F(*fHdNdEpSignalOff);

  // if no uncertainty in tau, that's all we need to do
  if(fDTau<=0) return 0;

  // chose the true tau for this simulated sample according to the pdf
  fTrueTau = rdm->Gaus(fTau,fDTau);
  
  if(fTrueTau>0)
    {
      hdNdEpBkg->SetBinContent(0,hdNdEpBkg->GetBinContent(0)*fTrueTau/fTau);
      if(hdNdEpSignalOff)
	hdNdEpSignalOff->SetBinContent(0,hdNdEpSignalOff->GetBinContent(0)*fTrueTau/fTau);
    }
  else
    {
      cout << "Iact1dUnbinnedLkl::GetRealBkgAndGoffHistos Error: negative or null tau value, we cannot work like that!" << endl;
      return 1;
    }
  
  return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Plot at the input canvas all material used to minimize the -2logL:
// - Aeff
// - Ereso and Ebias
// - dN/dE  for signal
// - dN/dE' for background compared to the On and Off distributions with residuals
// - dN/dE' for signal
// 
//
void Iact1dUnbinnedLkl::PlotHistosAndData(TCanvas* canvas) 
{
  MakeChecks();
  SetChecked(kFALSE);
  
  // create and divide canvas
  if(!canvas)
    canvas = new TCanvas("histosAndDataCanvas","Iact1dUnbinnedLkl histos and data used to minimize -2logL", 1000, 1500);
  canvas->Divide(2,3);

  // draw plots
  
  ///////////////////////////
  // CD(1) EFFECTIVE AREA
  ///////////////////////////
  canvas->cd(1);
  if(fHAeff)     fHAeff->DrawCopy();
  if(fHAeffOff)  fHAeffOff->DrawCopy("same");
  gPad->SetLogy();
  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();    

  ///////////////////////////
  // CD(2)  dN/dE' background vs data
  ///////////////////////////
  canvas->cd(2);

  // dN/dE' for bkg compared to Non and Noff distributions
  TH1F* hdNdEpBkg = NULL;
  if(fHdNdEpBkg)
    {
      hdNdEpBkg = new TH1F(*fHdNdEpBkg);
      hdNdEpBkg->SetDirectory(0);
    }
  TH1F* hOn  = GetHdNdEpOn();
  TH1F* hOff = GetHdNdEpOff();
  hOn->SetDirectory(0);  
  hOff->SetDirectory(0);
  hOff->Scale(1./fTau);

  if(hdNdEpBkg && fNoff>1)
    hdNdEpBkg->Scale(fNoff/fTau);

  // Foreground contribution (if any)
  TH1F* hdNdEpFrg = NULL;
  Float_t dNdEpFrgNorm = 1;
  if(fHdNdEpFrg)
    {
      hdNdEpFrg = new TH1F(*fHdNdEpFrg);
      hdNdEpFrg->SetDirectory(0);
      dNdEpFrgNorm = GetdNdEpFrgIntegral()*fObsTime;
      hdNdEpFrg->Scale(dNdEpFrgNorm);
    }

  // set the framework plot
  TH1I *dummya = new TH1I("dummya", "dN/dE' bkg model vs On and Off distributions",1,TMath::Log10(fEpmin),TMath::Log10(fEpmax));
  dummya->SetStats(0);
  if(fNon>1)
    {
      dummya->SetMinimum(hOn->GetMinimum(0)/2.);
      dummya->SetMaximum(hOn->GetMaximum()*2);
    }
  else if(fNoff>1)
    {
      dummya->SetMinimum(hOff->GetMinimum(0)/2.);
      dummya->SetMaximum(hOff->GetMaximum()*2);
    }
  else if(hdNdEpBkg)
    {
      dummya->SetMinimum(hdNdEpBkg->GetMinimum());
      dummya->SetMaximum(hdNdEpBkg->GetMaximum());
    }
    
  if(!hdNdEpBkg) dummya->SetTitle("dN/dE' distributions for On and Off event samples");
  dummya->SetXTitle("log_{10}(E' [GeV])");
  dummya->SetYTitle("dN/dE' [GeV^{-1}]");
  dummya->DrawCopy();

  // configure and plot the different histograms
  if(hdNdEpBkg)
    {
      hdNdEpBkg->SetMarkerStyle(1);
      hdNdEpBkg->SetLineColor(2);
      hdNdEpBkg->DrawCopy("hist same");
    }
  if(hdNdEpFrg)
    {
      hdNdEpBkg->SetMarkerStyle(2);
      hdNdEpBkg->SetLineColor(2);
      hdNdEpFrg->DrawCopy("hist same");
    }
  hOn->SetLineColor(4);
  hOn->SetMarkerColor(4);
  hOn->SetMarkerStyle(8);
  hOn->SetMarkerSize(0.5);
  hOn->DrawCopy("esame");
  hOff->SetLineColor(2);
  hOff->SetMarkerColor(2);
  hOff->SetMarkerStyle(8);
  hOff->SetMarkerSize(0.5);
  hOff->DrawCopy("esame");
 
  gPad->SetLogy();
  gPad->SetGrid();

  // legend
  TLegend* hleg = new TLegend(0.6, 0.65, 0.92, 0.92);
  hleg->SetFillColor(0);
  hleg->SetMargin(0.40);
  hleg->SetBorderSize(0);
  hleg->AddEntry(hOn,"On events","P");
  hleg->AddEntry(hOff,"Off events","P");
  if(hdNdEpBkg)
    hleg->AddEntry(hdNdEpBkg,"Background model","L");
  if(hdNdEpFrg)
    hleg->AddEntry(hdNdEpFrg,"Foreground model","L");    
  hleg->Draw();

  gPad->Modified();
  gPad->Update();    

  ///////////////////////////
  // CD(3)  energy dispersion
  ///////////////////////////
  canvas->cd(3);
  
  if(fMigMatrix)
    fMigMatrix->DrawCopy("colz");
  else
    {
      TH1I *dummye;
      if(fGEreso)
	dummye = new TH1I("dummye", "Energy resolution and bias",1,fGEreso->GetX()[0],fGEreso->GetX()[fGEreso->GetN()-1]);
      else
	dummye = new TH1I("dummye", "Energy resolution and bias",1,2,4);
      dummye->SetStats(0);
      dummye->SetMinimum(-0.1);
      dummye->SetMaximum(0.4);
      dummye->SetXTitle("log_{10}(E [GeV])");
      dummye->SetYTitle("Energy resolution and bias");
      dummye->DrawCopy();
      if(fGEreso) fGEreso->Draw("l");
      if(fGEbias) fGEbias->Draw("l");
      
      TLegend* hleg2 = new TLegend(0.30, 0.69, 0.55, 0.89);
      hleg2->SetFillColor(0);
      hleg2->SetMargin(0.40);
      hleg2->SetBorderSize(0);
      if(fGEreso) hleg2->AddEntry(fGEreso,"Resolution","L");
      if(fGEbias) hleg2->AddEntry(fGEbias,"Bias","LP");
      hleg2->Draw();
      delete dummye;
    }

  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();

  /////////////////////////////////////////////////
  // CD(4) On and Off residuals wrt background model
  /////////////////////////////////////////////////
  canvas->cd(4);
  TH1F* hResidualsOn  = NULL;
  TH1F* hResidualsOff = NULL;
  if(hdNdEpBkg && hOn && hOff)
    {
      hResidualsOn  =  GetResidualsHisto(hdNdEpBkg,hOn);
      hResidualsOff =  GetResidualsHisto(hdNdEpBkg,hOff);
    }
  
  dummya->SetMinimum(TMath::Min(hResidualsOn->GetMinimum(),-3.));
  dummya->SetMaximum(TMath::Max(hResidualsOn->GetMaximum(),3.));
  dummya->SetTitle("Residuals");
  dummya->SetXTitle("log_{10}(E' [GeV])");
  if(hdNdEpBkg)
    dummya->SetYTitle("(Data-dN/dE')/ #Delta Data");
  else
    dummya->SetTitle("(On-Off)/ #Delta On");
  dummya->DrawCopy();
  if(hdNdEpBkg && hOn && hOff)
    {
      hResidualsOn->SetLineColor(4);
      hResidualsOn->SetMarkerColor(4);
      hResidualsOn->SetMarkerStyle(8);
      hResidualsOn->SetMarkerSize(0.5);  
      hResidualsOn->DrawCopy("esame");
      hResidualsOff->SetLineColor(2);
      hResidualsOff->SetMarkerColor(2);
      hResidualsOff->SetMarkerStyle(8);
      hResidualsOff->SetMarkerSize(0.5);  
      hResidualsOff->DrawCopy("esame");
    }
  else if(hOff && hOn)
    {
      hResidualsOn  =  GetResidualsHisto(hOff,hOn);
      hResidualsOn->DrawCopy("esame");
    }
  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();

 
  /////////////////////////
  // CD(5) dN/dE for signal
  /////////////////////////
  canvas->cd(5);
  dummya->SetMinimum(1e-7);
  dummya->SetMaximum(1e0);
  dummya->SetTitle("dN/dE for signal events");
  dummya->SetXTitle("log_{10}(E [GeV])");
  dummya->SetYTitle("dN/dE [GeV^{-1}]");
  dummya->DrawCopy();
  TH1F* hdNdESignal = NULL;

  if(fHdNdESignal)
    {
      hdNdESignal = new TH1F(*fHdNdESignal);
      hdNdESignal->SetDirectory(0);
      
      Double_t scale = GetdNdESignalIntegral();     
      hdNdESignal->Scale(scale);
      hdNdESignal->SetLineWidth(1);
      hdNdESignal->SetLineStyle(1);
      hdNdESignal->SetLineColor(1);
      hdNdESignal->DrawCopy("hist same");
      // old way
      //Double_t scale = GetdNdESignalIntegral();
      //cout<<"scale normalization = "<<scale<<endl;
      //fHdNdESignal->Scale(scale);
      //fHdNdESignal->DrawCopy("same");
      //fHdNdESignal->Scale(1./scale);      
    }
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->Modified();
  gPad->Update();

  //////////////////////////
  // CD(6) dN/dE' for signal
  //////////////////////////
  canvas->cd(6);
  dummya->SetMinimum(1e2);
  dummya->SetMaximum(1e8);
  dummya->SetTitle("dN/dE' (#times Aeff) for signal events");
  dummya->SetXTitle("log_{10}(E' [GeV])");
  dummya->SetYTitle("dN/dE'(#times A_{eff}) [GeV^{-1}cm^{2}]");
  dummya->DrawCopy();

  TH1F* hdNdEpSignal = NULL;
  if(fHdNdEpSignal)
    {
      hdNdEpSignal = new TH1F(*fHdNdEpSignal);
      hdNdEpSignal->SetDirectory(0);

      Double_t scale = GetdNdEpSignalIntegral();
      hdNdEpSignal->Scale(scale);      
      hdNdEpSignal->DrawCopy("hist same");
    }

  TH1F* hdNdEpSignalOff = NULL;
  if(fHdNdEpSignalOff)
    {
      hdNdEpSignalOff = new TH1F(*fHdNdEpSignalOff);
      hdNdEpSignalOff->SetDirectory(0);
      
      Double_t scale = GetdNdEpSignalOffIntegral();
      hdNdEpSignalOff->Scale(scale);
      hdNdEpSignalOff->SetLineStyle(2);
      hdNdEpSignalOff->Draw("hist same");
    }
  
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->Modified();
  gPad->Update();

  // clean and exit
  if(hdNdEpBkg) delete hdNdEpBkg;
  if(hdNdEpFrg) delete hdNdEpFrg;
  if(hdNdESignal) delete hdNdESignal;
  if(hdNdEpSignal) delete hdNdEpSignal;
  if(hdNdEpSignalOff) delete hdNdEpSignalOff;
  delete hOn;
  delete hOff;
  if(hResidualsOn)  delete hResidualsOn;
  if(hResidualsOff) delete hResidualsOff;
  delete dummya;
}

//////////////////////////////////////////////////////////////////
//
// Print the values of the parameters of the Lkl-based object
//
void Iact1dUnbinnedLkl::PrintData(Int_t level) 
{
  Lkl::PrintData(level);
  
  Margin(level); cout << "                E' range  = [" << fEpmin << "," << fEpmax << "] GeV" <<endl;
  Margin(level); cout << "                 Non/Noff = " << fNon << "/" << fNoff << endl;  
  Margin(level); cout << "             Measured tau = " << fTau << " +/- " << fDTau << endl;
  Margin(level); cout << "         Observation time = " << fObsTime/3600. << " h" << endl;
  Margin(level); cout << "                   log(J) = ("<< fLogJ << " +/- "<< (GetDUofGType()==Lkl::invlog? GetDUnitsOfG() : -9999) << ") GeV^2 cm^-5  or  GeV cm^-2" << endl;

  Margin(level); cout << "               Aeff histo : " << (fHAeff?    "YES" : "NO") << endl;
  Margin(level); cout << "           Off Aeff histo : " << (fHAeffOff? "YES" : "NO") << endl;
  Margin(level); cout << "       E Reso/Bias graphs : " << ((fGEreso && fGEbias)? "YES" : "NO") << endl;
  Margin(level); cout << "         Migration matrix : " << (fMigMatrix? "YES" : "NO") << endl;
  Margin(level); cout << "    Background rate histo : " << (fHdNdEpBkg? "YES" : "NO") << endl;
  Margin(level); cout << "    Foreground rate histo : " << (fHdNdEpFrg? "YES" : "NO") << endl;

  Margin(level); cout << "       Signal dN/dE histo : " << (fHdNdESignal? "YES" : "NO") << endl;
  Margin(level); cout << "      Signal dN/dE' histo : " << (fHdNdEpSignal? "YES" : "NO") << endl;
  Margin(level); cout << "  Signal dN/dE' Off histo : " << (fHdNdEpSignalOff? "YES" : "NO") << endl;
  if(fHdNdEpSignalOff && fHdNdEpSignal)
    {Margin(level); cout << "            Signal in Off = " << GetdNdEpSignalOffIntegral()/GetdNdEpSignalIntegral()*100  << "% of that in On" << endl;}
}

//////////////////////////////////////////////////////////////////
//
// Produce the E' distribution of On events and return the 
// corresponding histogram.
// if isDifferential=kTRUE (default is kTRUE), return the 
// dN/dE' distribution (i.e. number of entries in the bin divided
// by the size of the bin in E' unit). Otherwise the returned
// histogram is the number of events per bin of E'.
// The returned histogram will not be deleted by the destructor
// of Iact1dUnbinnedLkl.
//
TH1F* Iact1dUnbinnedLkl::GetHdNdEpOn(Bool_t isDifferential,Int_t nbins) const
{
  // we need a positive number of bins
  if(nbins<=0) nbins = gNBins;
  
  // create histo
  TH1F* h = new TH1F("dNdEpOn","dN/dE' for On events",nbins,TMath::Log10(fEpmin),TMath::Log10(fEpmax));
  h->SetDirectory(0);
  h->SetXTitle("log_{10}(E' [GeV])");
  h->SetYTitle("dN/dE' [GeV^{-1}]");

  // fill histo
  for(Int_t i=0;i<fNon;i++)
    h->Fill(fOnSample[i]);

  // divide by bin width
  if(h->GetEntries()>0)
    if(isDifferential)
      for(Int_t ibin=0;ibin<nbins;ibin++)
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
// if isDifferential=kTRUE (default is kTRUE), return the 
// dN/dE' distribution (i.e. number of entries in the bin divided
// by the size of the bin in E' unit). Otherwise the returned
// histogram is the number of events per bin of E'.
// The returned histogram will not be deleted by the destructor
// of Iact1dUnbinnedLkl.
//
TH1F* Iact1dUnbinnedLkl::GetHdNdEpOff(Bool_t isDifferential,Int_t nbins) const
{
  // we need a positive number of bins
  if(nbins<=0) nbins = gNBins;

  // create histo
  TH1F* h = new TH1F("dNdEpOff","dN/dE' for Off events",nbins,TMath::Log10(fEpmin),TMath::Log10(fEpmax));
  h->SetDirectory(0);
  h->SetXTitle("log_{10}(E' [GeV])");
  h->SetYTitle("dN/dE' [GeV^{-1}]");

  // fill histo
  for(Int_t i=0;i<fNoff;i++)
    h->Fill(fOffSample[i]);
  
  // divide by bin width
  if(h->GetEntries()>0)
    if(isDifferential)
      
      for(Int_t ibin=0;ibin<nbins;ibin++)
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
// Copy the contents of ih into oh by interpolating if necessary
// the x axis of ih and oh must be in log-scale and even binning
// the interpolation is done linearly in log(y)
// the units of the x-axis of ih are those of oh times 10^scale
// if isDiff=kTRUE the input histogram is differential (default)
// if isDiff=kFALSE the input histogram is bin-integrated
// the output histogram is ALWAYS differential
//
void readAndInterpolate(TH1F* ih,TH1F* oh,Double_t scale,Bool_t isDiff)
{  
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
	  // ycontent0 = ih->GetBinContent(1);
	  // ycontent1 = ih->GetBinContent(2);
	  
	  // if(ycontent0<=0 || ycontent1<=0)
	  //   {
	  //     oh->SetBinContent(ibin+1,(ycontent0+ycontent1)/2.);
	  //     continue;
	  //   }

	  // x0 = imine+gCenterBin*ide;
	  // x1 = imine+(1.+gCenterBin)*ide;

	  // if(!isDiff)
	  //   {
	  //     dx0i = TMath::Power(10,imine+ide)-TMath::Power(10,imine);
	  //     dx1i = TMath::Power(10,imine+2*ide)-TMath::Power(10,imine+ide);

	  //     y0  = TMath::Log10(ycontent0/dx0i);
	  //     y1  = TMath::Log10(ycontent1/dx1i);	  
	  //   }
	  // else
	  //   {
	  //     y0 = TMath::Log10(ycontent0);
	  //     y1 = TMath::Log10(ycontent1);	  
	  //   }
	  oh->SetBinContent(ibin+1,0);
	  continue;
	}
      else if(etest>imaxe-(1-gCenterBin)*ide) // extrapolation of values above the maximum input energy
	{
	  // ycontent0 = ih->GetBinContent(onbinse-1);
	  // ycontent1 = ih->GetBinContent(onbinse);
	  
	  // if(ycontent0<=0 || ycontent1<=0)
	  //   {
	  //     oh->SetBinContent(ibin+1,(ycontent0+ycontent1)/2.);
	  //     continue;
	  //   }

	  // x0 = imaxe-(2.-gCenterBin)*ide;
	  // x1 = imaxe-(1.-gCenterBin)*ide;
	  // if(!isDiff)
	  //   {
	  //     dx0i = TMath::Power(10,imine+(inbinse-1)*ide)-TMath::Power(10,imine+(inbinse-2)*ide);
	  //     dx1i = TMath::Power(10,imine+inbinse*ide)-TMath::Power(10,imine+(inbinse-1)*ide);
	  //     y0  = TMath::Log10(ycontent0/dx0i);
	  //     y1  = TMath::Log10(ycontent1/dx1i);	    
	  //   }
	  // else
	  //   {
	  //     y0 = TMath::Log10(ycontent0);
	  //     y1 = TMath::Log10(ycontent1);	    
	  //   }
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
      // cout << ibin+1 << ": ycontent0 = " << ycontent0 << ", ycontent1 = " << ycontent1 << "y0 = " << y0 << ", y1 = " << y1 << ", x0 = " << x0 << ", x1 = " << x1 << ", rtest = " << rtest << endl;
    }
}
//////////////////////////////////////////////////////////////////
//
// Copy the contents of ih into oh bin by bin
// the x axis of ih and oh must be in log-scale and oh in even binning (not necessary for ih)
// Coarse bins in ih are copied into finner bins in oh, but
// the shape (i.e. the discontinuities between bin limits) is preserved 
// the units of the x-axis of ih are those of oh times 10^scale
// if isDiff=kTRUE the input histogram is differential (default)
// if isDiff=kFALSE the input histogram is bin-integrated
// the output histogram is ALWAYS differential
Int_t copyBinByBin(TH1F* ih,TH1F* oh,Double_t scale,Bool_t isDiff)
{
  // input histogram binning
  Double_t imine   = ih->GetXaxis()->GetXmin()+scale; // minimum log(E) in input histo
  Double_t imaxe   = ih->GetXaxis()->GetXmax()+scale; // maximum log(E) in input histo
  Int_t    inbinse = ih->GetNbinsX();          
  

  // output histogram binning
  Double_t omine   = oh->GetXaxis()->GetXmin(); // minimum log(E) in output histo
  Double_t omaxe   = oh->GetXaxis()->GetXmax(); // maximum log(E) in output histo
  Int_t    onbinse = oh->GetNbinsX();          
  Double_t ode     = (omaxe-omine)/onbinse;

  // copy values
  for(Int_t ibin=0;ibin<onbinse;ibin++)
    {
      Double_t etest = omine+ode*(ibin+0.5); 
      
      // copy bin by bin, zero outside limits
      if(etest<imine || etest>imaxe)	
	oh->SetBinContent(ibin+1,0);
      else
	{
	  // corresponding bin in ih histo
	  Int_t etestbin;
	  for(etestbin=0;etestbin<inbinse;etestbin++)
	    if(ih->GetBinLowEdge(etestbin+1)+scale<etest && ih->GetBinLowEdge(etestbin+1)+ih->GetBinWidth(etestbin+1)+scale>etest)
	      break;
	  if(etestbin>=inbinse)
	    {
	      cout << "copyBinByBin Warning: the two histograms must have very different energy ranges" << endl;
	      return 1;
	    }
	  Float_t dE = 1;
	  if(!isDiff)
	    {
	      Double_t leminbin = ih->GetBinLowEdge(etestbin+1)+scale;
	      Double_t lemaxbin = ih->GetBinLowEdge(etestbin+1)+ih->GetBinWidth(etestbin+1)+scale;
	      dE                = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);

	    }
	  oh->SetBinContent(ibin+1,ih->GetBinContent(etestbin+1)/dE);    
	}      
    }
  return 0;
}


////////////////////////////////////////////////////////////////////////
// 
// Given a histogram hModel representing a model
// and another one hData representing data (with properly computed errors)
// return the histogram of residuals (hData-hModel)
//
TH1F* GetResidualsHisto(TH1F* hModel,TH1F* hData)
{
  // basic check
  if(!hModel || !hData) return NULL;

  // get number of data bins
  UInt_t nbins = hData->GetNbinsX();

  // check if bin width is constant
  Bool_t binWidthIsConstant = kTRUE;
  for(Int_t ibin=0;ibin<nbins-1;ibin++)
    if(TMath::Abs(Float_t(hData->GetBinWidth(ibin+1))-Float_t(hData->GetBinWidth(ibin+2)))>1e-9)
      binWidthIsConstant = kFALSE;

  // create reiduals histo
  TH1F* hResiduals;
  if(binWidthIsConstant)
    hResiduals = new TH1F("hResiduals","Residuals data-Model",nbins,hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  else
    hResiduals = new TH1F("hResiduals","Residuals data-Model",nbins,hData->GetXaxis()->GetXbins()->GetArray());
  hResiduals->SetDirectory(0);
  
  // fill histo
  for(Int_t ibin=0;ibin<nbins;ibin++)
    {
      Double_t letest = hResiduals->GetBinLowEdge(ibin+1)+hResiduals->GetBinWidth(ibin+1)/2.;
      Int_t modbin = hModel->FindBin(letest);
      Double_t errbar = hData->GetBinError(ibin+1);
      if(errbar>0)
	{
	  hResiduals->SetBinContent(ibin+1,(hData->GetBinContent(ibin+1)-hModel->GetBinContent(modbin))/errbar);
	  hResiduals->SetBinError(ibin+1,1);
	}
      else
	{
	  hResiduals->SetBinContent(ibin+1,0);
	  hResiduals->SetBinError(ibin+1,0);
	}      
    }
  return hResiduals;
}

////////////////////////////////////////////////////////////////
// Normalize the rows of im (ie. for a fixed Etrue value)
// in such a way that the integral vs dE' for a given line (fixed E)
// is 1
// the units of the y-axis of im are those of om times 10^scale
// and store values in om
// im and om must have same binning in both axes
void NormalizeMatrix(TH2F* im,TH2F* om,Double_t scale)
{
  Int_t nbinste   = im->GetNbinsY();  // number of bins in Etrue
  Int_t nbinsme   = im->GetNbinsX();  // number of bins in Emeas
  TAxis* emaxis   = im->GetXaxis();
	    
  // normalize for each Etrue value
  for(Int_t ibinte=0;ibinte<nbinste;ibinte++)
    {
      // compute total number of events with same Etrue
      Double_t tot=0;
      for(Int_t ibinme=0;ibinme<nbinsme;ibinme++)
	tot+=im->GetBinContent(ibinme+1,ibinte+1);

      // normalize matrix lines (fixed Etrue)
      if(tot>0)
       	for(Int_t ibinme=0;ibinme<nbinsme;ibinme++)
	  {
	    Double_t minbinme   = TMath::Power(10,emaxis->GetBinLowEdge(ibinme+1));
	    Double_t maxbinme   = TMath::Power(10,emaxis->GetBinUpEdge(ibinme+1));
	    Double_t dem        = (maxbinme-minbinme)*TMath::Power(10,scale);  // true energy bin size [GeV]
      
	    om->SetBinContent(ibinme+1,ibinte+1,im->GetBinContent(ibinme+1,ibinte+1)/tot/dem);
	  }
    }
}

////////////////////////////////////////////////////////////////
// Smear a histogram <sp> (in log of true energy) using 
// the energy dispersion function described by grreso and grbias
// which are energy resolution and bias vs log(E[GeV])
// and put the result in histogram <smsp> (in measured energy)
//
// Return 0 in case of success
//        1 otherwise
//
Int_t SmearHistogram(TH1F* sp,TH1F* smsp,TGraph* grreso,TGraph* grbias)
{  
  // checks
  if(!sp || !smsp || !grreso || !grbias)
    {
      cout << "SmearHistogram Warning: missing histos" << endl;
      return 1;
    }

  Int_t nbinste = sp->GetNbinsX();              // number of bins in input histo
  Int_t nbinsme = smsp->GetNbinsX();            // number of bins in output histo
  
  // do the convolution of sp with energy resolution and store result in smsp
  for(Int_t ibinte=0;ibinte<nbinste;ibinte++)
    {
      Double_t let  = sp->GetBinCenter(ibinte+1); // log true energy
      Double_t et   = TMath::Power(10,let);       // true energy
      Double_t reso = grreso->Eval(let);          // energy resolution
      Double_t bias = grbias->Eval(let);          // energy bias
      
      // bin size in true enery
      Double_t minbinte = TMath::Power(10,sp->GetBinLowEdge(ibinte+1));
      Double_t maxbinte = TMath::Power(10,sp->GetBinLowEdge(ibinte+1)+sp->GetBinWidth(ibinte+1));
      Double_t det      = maxbinte-minbinte;  // bin size [GeV]
      
      for(Int_t ibinme=0;ibinme<nbinsme;ibinme++)
	{	  	  
	  // compute smearing factor
	  Double_t minbinme = TMath::Power(10,smsp->GetBinLowEdge(ibinme+1));
	  Double_t maxbinme = TMath::Power(10,smsp->GetBinLowEdge(ibinme+1)+smsp->GetBinWidth(ibinme+1));
	  Double_t dem      = maxbinme-minbinme;  // bin size [GeV]

	  Double_t mingausme    = (minbinme-et*(1+bias))/(et*reso);
	  Double_t maxgausme    = (maxbinme-et*(1+bias))/(et*reso);
	  Double_t gausintegral = (TMath::Erf(maxgausme/TMath::Sqrt(2))-TMath::Erf(mingausme/TMath::Sqrt(2)))/2.;	  
	  Double_t smfactor     = gausintegral/dem;
	  smsp->SetBinContent(ibinme+1,smsp->GetBinContent(ibinme+1)+sp->GetBinContent(ibinte+1)*det*smfactor);
	}
    }
  return 0;
}

////////////////////////////////////////////////////////////////
// smear a spectral shape <sp> (in true energy) using 
// the energy dispersion function from migration matrix <mm>
// and put the result in histogram <smsp> (in measured energy)
// <mm> is computed as N_ij/(N_j*DeltaE_j)
// with N_ij number of events passing all analysis cuts, and
// with true energy in DeltaE_j and recontructed energy in DeltaE_i;
// N_j total number of events passing all analysis cuts with true energy in DeltaE_j;
// and DeltaE_j the size [GeV] of the DeltaE_j energy bin
// 
Int_t SmearHistogram(TH1F* sp,TH1F* smsp,TH2F* mm)
{
  // checks
  if(!sp || !smsp || !mm)
    {
      cout << "SmearHistogram Warning: missing histos" << endl;
      return 1;
    }

  Int_t nbinste  = sp->GetNbinsX();     // number of bins in input histo
  TH1F* provsmsp = new TH1F("provsmsp","Provisonal smeared histo",mm->GetXaxis()->GetNbins(),mm->GetXaxis()->GetXmin(),mm->GetXaxis()->GetXmax());
  provsmsp->SetDirectory(0);
  
  // do the convolution of sp with mm and store result in smsp
  for(Int_t ibinte=0;ibinte<nbinste;ibinte++)
    {
      Double_t et = sp->GetBinCenter(ibinte+1); // log of true energy

      // bin size in true energy
      Double_t minbinte = TMath::Power(10,sp->GetBinLowEdge(ibinte+1));
      Double_t maxbinte = TMath::Power(10,sp->GetBinLowEdge(ibinte+1)+sp->GetBinWidth(ibinte+1));
      Double_t det      = maxbinte-minbinte;  // bin size

      // do the convolution of sp with the matrix row j=ibinte and fill
      // a histogram with same binning as matrix column i=ibinme
      for(Int_t iprovbinme=0;iprovbinme<mm->GetXaxis()->GetNbins();iprovbinme++)
	{
	  Double_t em       = provsmsp->GetBinCenter(iprovbinme+1);
	  Int_t    gidbin   = mm->FindBin(em,et);
	  Double_t smfactor = mm->GetBinContent(gidbin);
	  provsmsp->SetBinContent(iprovbinme+1,provsmsp->GetBinContent(iprovbinme+1)+sp->GetBinContent(ibinte+1)*det*smfactor);	  
	}
    }

  // transfer provsmsp to smsp (with different binning, in general...)
  if(copyBinByBin(provsmsp,smsp))
    return 1;
  
  delete provsmsp;
  return 0;
}


////////////////////////////////////////////////////////////////////////
// full likelihood function (-2logL) 
// To be minimized by TMinuit
// Free parameters:
// par[0] = g (total estimated number of signal events in On region)
// par[1] = b (total estimated number of background events in On region)
// par[2] = estimated value of tau
//
void fullLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;

  // get internal object, histos, values, etc
  Iact1dUnbinnedLkl*        mylkl           = dynamic_cast<Iact1dUnbinnedLkl*>(minuit->GetObjectFit());
  const TH1F*      hdNdEpSignal    = mylkl->GetHdNdEpSignal();
  const TH1F*      hdNdEpSignalOff = mylkl->GetHdNdEpSignalOff();
  const TH1F*      hdNdEpBkg       = mylkl->GetHdNdEpBkg();
  const TH1F*      hdNdEpFrg       = mylkl->GetHdNdEpFrg();
  const Float_t*   onSample        = mylkl->GetOnSample();
  const Float_t*   offSample       = mylkl->GetOffSample();
  UInt_t           Non             = mylkl->GetNon();
  UInt_t           Noff            = mylkl->GetNoff();
  Float_t          tau             = mylkl->GetTau();
  Float_t          dTau            = mylkl->GetDTau();
  Float_t          Tobs            = mylkl->GetObsTime();
 
  const Int_t      nbins           = hdNdEpBkg->GetNbinsX();
  const Double_t   xmin            = hdNdEpBkg->GetXaxis()->GetXmin();
  const Double_t   xmax            = hdNdEpBkg->GetXaxis()->GetXmax();
  
  // Estimated number of background events in signal and background regions
  Double_t g       = par[0];
  Double_t b       = par[1];
  Double_t tauest  = par[2];
  Double_t frg     = (hdNdEpFrg? mylkl->GetdNdEpFrgIntegral()*Tobs : 0);
  Double_t boff    = b*tauest;
  Double_t goff    = (hdNdEpSignalOff? tauest*g*mylkl->GetdNdEpSignalOffIntegral()/mylkl->GetdNdEpSignalIntegral() : 0);
  Double_t fnorm   = g+b+frg+boff+goff;
  
  
  // sum signal and background (and maybe foreground) contributions and normalize resulting pdf (On + Off)
  TH1F* hdNdEpOn  = new TH1F("hdNdEpOn", "On  event rate vs E'",nbins,xmin,xmax);
  hdNdEpOn->Reset();

  hdNdEpOn->Add(hdNdEpSignal,hdNdEpBkg,g,b);
  if(hdNdEpFrg)
    hdNdEpOn->Add(hdNdEpOn,hdNdEpFrg,1,frg);

  // normalize
  if(fnorm>0)
    hdNdEpOn->Scale(1./fnorm);
  else
    mylkl->NormalizedNdEHisto(hdNdEpOn);

  TH1F* hdNdEpOff = new TH1F("hdNdEpOff","Off event rate vs E'", nbins,xmin,xmax);
  hdNdEpOff->Reset();
  if(hdNdEpSignalOff)
    hdNdEpOff->Add(hdNdEpSignalOff,hdNdEpBkg,goff,boff); 
  else
    hdNdEpOff->Add(hdNdEpBkg,boff); 

  // normalize
  if(fnorm>0)  
    hdNdEpOff->Scale(1./fnorm);
  else
    mylkl->NormalizedNdEHisto(hdNdEpOff);
 
  // -2 log-likelihood
  f = 0;

  // On events
  for(ULong_t ievent=0; ievent<Non; ievent++)
    {
      Float_t val = hdNdEpOn->GetBinContent(hdNdEpOn->FindBin(onSample[ievent]));
      if(val>0)
	f += -2*TMath::Log(val);
      else
	f += gLklValVeryHigh;
    }
  
  // Off events
  for(ULong_t ievent=0; ievent<Noff; ievent++)
    {
      Float_t val = hdNdEpOff->GetBinContent(hdNdEpOff->FindBin(offSample[ievent]));
      if(val>0)
      f += -2*TMath::Log(val);
        
      else
	f += gLklValVeryHigh;
    }
    
  // nuisance tau
  if(dTau>0)
    f+=-2*TMath::Log(TMath::Gaus(tauest, tau, dTau, kTRUE));

  // tot Nevts and nuisance Noff
  if(g+b+frg>0)
    f += -2*TMath::Log(TMath::Poisson(Non,g+b+frg));
  else
    f += gLklValVeryHigh;

  if(goff+boff>0)
    f += -2*TMath::Log(TMath::Poisson(Noff,goff+boff));
  else
    f += gLklValVeryHigh;  

  delete hdNdEpOn;
  delete hdNdEpOff;
}		
