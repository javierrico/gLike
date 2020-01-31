/* ======================================================================== *\
!
!   Author: Javier Rico         03/2017 <mailto:jrico@ifae.es>
!   Author: Joaquim Palacio     03/2017 <mailto:jpalacio@ifae.es>
!
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//						
// IactEventListIrf
//
// Class to hold the data and IRFs from IACTs to be used as input by
// Iact1dUnbinnedLkl class 
//
//////////////////////////////////////////////////////////////////////////////

#include "IactEventListIrf.h"

ClassImp(IactEventListIrf);

using namespace std;

static const Double_t gEpmin = 1e00; // [GeV] default value of minimum E_est
static const Double_t gEpmax = 1e06; // [GeV] default value of maximum E_est
static const Int_t   gBuffSize = 512000; // buffer size for Ntuples

// static functions, helpers to load data from the FITS Header Data Units (HDUs)
static void FillEventListFromHDU(TFITSHDU* hdu, TNtupleD* onSample);
static TH1F* GetHAeffFromHDU(TFITSHDU* hdu);
static TGraph* GetGEbiasFromHDU(TFITSHDU* hdu);
static TGraph* GetGEresoFromHDU(TFITSHDU* hdu);

const Double_t IactEventListIrf::gDefEVal      = 0.;    // default value when energy is not provided
const Double_t IactEventListIrf::gDefRADECVal  = 9999.; // default value when dRA and dDEC are not provided
const Double_t IactEventListIrf::gDefTVal      = -1.;    // default value when time is not provided
const Double_t IactEventListIrf::gDefHadVal    = -1.;    // default value when hadronness is not provided

////////////////////////////////////////////////////////////////
//
// Default constructor, just create empty ntuples for data samples
//
IactEventListIrf::IactEventListIrf(TString name,TString title) :
  TNamed(name,title), 
  fOnSample(NULL), fOffSample(NULL),
  fEpmin(gEpmin), fEpmax(gEpmax), fTau(1), fDTau(0),
  fObsTime(0), fHAeff(NULL), fHAeffOff(NULL), fGEreso(NULL),
  fGEbias(NULL), fMigMatrix(NULL), fHdNdEpBkg(NULL), fHdNdEpFrg(NULL)
{
  // create the event lists
  fOnSample  = new TNtupleD("fOnSample", "On data set", "E:pointRA:pointDEC:dRA:dDEC:t:had",gBuffSize);
  fOffSample = new TNtupleD("fOffSample","Off data set","E:pointRA:pointDEC:dRA:dDEC:t:had",gBuffSize);
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
IactEventListIrf::~IactEventListIrf()
{
  if(fOnSample)     delete fOnSample;
  if(fOffSample)    delete fOffSample;
  
  if(fHAeff)        delete fHAeff;    
  if(fHAeffOff)     delete fHAeffOff;
  if(fGEreso)       delete fGEreso;
  if(fGEbias)       delete fGEbias;
  if(fMigMatrix)    delete fMigMatrix;

  if(fHdNdEpBkg)    delete fHdNdEpBkg;
  if(fHdNdEpFrg)    delete fHdNdEpFrg;
}

////////////////////////////////////////////////////////////////
// 
// Load the IactEventListIrf from a FITS input file 
// TODO : generalize this function to a Load() member able to 
// select between ROOT and FITS inputs
// 
void IactEventListIrf::LoadFITS(TString inputFileName){
  cout << "IactEventListIrf::LoadFITS Info: loading the dataset in: " << inputFileName.Data() << endl;
  
  // open all the Header Data Units (HDUs)
  TFITSHDU* hduHeader = new TFITSHDU(inputFileName, 0);
  TFITSHDU* hduOn = new TFITSHDU(inputFileName, 1);
  TFITSHDU* hduOff = new TFITSHDU(inputFileName, 2);
  TFITSHDU* hduAeff = new TFITSHDU(inputFileName, 3);
  TFITSHDU* hduEdisp = new TFITSHDU(inputFileName, 4);
  
  // fetch the effective time and tau from the HEADER 
  Double_t obsTime = hduHeader->GetKeywordValue("TEFF").Atof();
  Double_t acceptanceOn = hduHeader->GetKeywordValue("ACC").Atof();
  Double_t acceptanceOff = hduHeader->GetKeywordValue("ACC_OFF").Atof();
  Double_t tau = acceptanceOff / acceptanceOn;
  
  // set all the attributes
  fObsTime = obsTime;
  fTau = tau;
  FillEventListFromHDU(hduOn, fOnSample);
  FillEventListFromHDU(hduOff, fOffSample);
  fHAeff = GetHAeffFromHDU(hduAeff);
  fGEbias = GetGEbiasFromHDU(hduEdisp);
  fGEreso = GetGEresoFromHDU(hduEdisp);

  delete hduHeader;
  delete hduOn; 
  delete hduOff; 
  delete hduAeff; 
  delete hduEdisp;
}

void IactEventListIrf::Print(Option_t* o) const
{
  cout << "Energy Range: " <<  fEpmin << " - " << fEpmax << " GeV " << endl;
  cout << "ON/OFF Norm (tau): " <<  fTau << " +- " << fDTau << " (Prob.: " << fTauPValue << ")" << endl;
  cout << "Observation Time: " <<  fObsTime  << " s "  << endl;
}

////////////////////////////////////////////////////////////////
// 
// fill the ON / OFF event list from the Header Data Unit of the FITS
// NOTE : energies in the FITS files are given in TeV
//
void FillEventListFromHDU(TFITSHDU* hdu, TNtupleD* eventList){
  TVectorD* energies = hdu->GetTabRealVectorColumn("ENERGY");
  TVectorD* times = hdu->GetTabRealVectorColumn("TIME");
  for(Int_t i=0; i<energies->GetNrows(); ++i){
    eventList->Fill(energies[0][i] * 1e3, 0., 0., 0., 0., times[0][i], 0.);
  }
}

////////////////////////////////////////////////////////////////
// 
// get the effective area from the Header Data Unit of the FITS
// NOTE : energies in the FITS files are given in TeV
// 
TH1F* GetHAeffFromHDU(TFITSHDU* hdu){
  TVectorD* lowEtrueEdges = hdu->GetTabRealVectorColumn("ENERG_LO");
  TVectorD* highEtrueEdges = hdu->GetTabRealVectorColumn("ENERG_HI");
  TVectorD* aeffValues = hdu->GetTabRealVectorColumn("AEFF");
  Int_t nBinsAeff = aeffValues->GetNrows();
  Double_t log10EtrueEdges[nBinsAeff+1];
  for (Int_t i=0; i<nBinsAeff; ++i){
    log10EtrueEdges[i] = TMath::Log10(1e3 * lowEtrueEdges[0][i]);
  } 
  // the last bin edge is the last element of the highEtrueEdges
  log10EtrueEdges[nBinsAeff] = TMath::Log10(1e3 * highEtrueEdges[0][nBinsAeff-1]);
  // define the TH1F to be returned, log10(E / GeV) on the X axis
  TH1F* hAeff = new TH1F();
  hAeff->SetBins(nBinsAeff, log10EtrueEdges);
  for (Int_t i=0; i<nBinsAeff; ++i){
    hAeff->SetBinContent(i, aeffValues[0][i]);
  }
  return hAeff;
}

////////////////////////////////////////////////////////////////
// 
// get the bias from the Header Data Unit of the FITS
// NOTE : energies in the FITS files are given in TeV
// 
TGraph* GetGEbiasFromHDU(TFITSHDU* hdu){
  TVectorD* eTrue = hdu->GetTabRealVectorColumn("E_TRUE");
  TVectorD* bias = hdu->GetTabRealVectorColumn("BIAS");
  Int_t nBinsEdisp = bias->GetNrows();
  // TGraph to be returned, log10(E / GeV) on the X axis
  TGraph* biasTGraph = new TGraph();
  for (Int_t i=0; i<nBinsEdisp; ++i){
    biasTGraph->SetPoint(i, TMath::Log10(1e3 * eTrue[0][i]), bias[0][i]);
  }
  return biasTGraph;
}

////////////////////////////////////////////////////////////////
// 
// get the resolution from the Header Data Unit of the FITS
// NOTE : energies in the FITS files are given in TeV
// 
TGraph* GetGEresoFromHDU(TFITSHDU* hdu){
  TVectorD* eTrue = hdu->GetTabRealVectorColumn("E_TRUE");
  TVectorD* resolution = hdu->GetTabRealVectorColumn("RES");
  Int_t nBinsEdisp = resolution->GetNrows();
  // TGraph to be returned, log10(E / GeV) on the X axis
  TGraph* resolutionTGraph = new TGraph();
  for (Int_t i=0; i<nBinsEdisp; ++i){
    resolutionTGraph->SetPoint(i, TMath::Log10(1e3 * eTrue[0][i]), resolution[0][i]);
  }
  return resolutionTGraph;
}