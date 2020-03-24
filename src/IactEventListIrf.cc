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

#include <stdlib.h>
#include "TObject.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "IactEventListIrf.h"

ClassImp(IactEventListIrf);

using namespace std;

static const Double_t gEpmin = 1e00; // [GeV] default value of minimum E_est
static const Double_t gEpmax = 1e06; // [GeV] default value of maximum E_est
static const Int_t   gBuffSize = 512000; // buffer size for Ntuples

// static functions, helpers to load data from the FITS Header Data Units (HDUs)
static void FillEventListFromHDU(TFITSHDU* hdu, TNtupleD* eventList);
static TH1F* GetHAeffFromHDU(TFITSHDU* hdu);
static TH2F* GetMigMatrixFromHDU(TFITSHDU* hduMatrix, TFITSHDU* hduEbounds);

const Double_t IactEventListIrf::gDefEVal      = 0.;    // default value when energy is not provided
const Double_t IactEventListIrf::gDefRADECVal  = 9999.; // default value when dRA and dDEC are not provided
const Double_t IactEventListIrf::gDefTVal      = -1.;    // default value when time is not provided
const Double_t IactEventListIrf::gDefHadVal    = -1.;    // default value when hadronness is not provided

////////////////////////////////////////////////////////////////
//
// Default constructor, just create empty ntuples for data samples
//
IactEventListIrf::IactEventListIrf(TString name, TString title) : TNamed(name, title)
{
  _initialize_me();
  // create the event lists
  fOnSample  = new TNtupleD("fOnSample", "On data set", "E:pointRA:pointDEC:dRA:dDEC:t:had", gBuffSize);
  fOffSample = new TNtupleD("fOffSample","Off data set","E:pointRA:pointDEC:dRA:dDEC:t:had", gBuffSize);
}

////////////////////////////////////////////////////////////////
//
// constructor with input file
//
IactEventListIrf::IactEventListIrf(TString name, TString title, TString fileName) : TNamed(name, title)
{
  _initialize_me();
  // create the event lists
  fOnSample  = new TNtupleD("fOnSample", "On data set", "E:pointRA:pointDEC:dRA:dDEC:t:had", gBuffSize);
  fOffSample = new TNtupleD("fOffSample","Off data set","E:pointRA:pointDEC:dRA:dDEC:t:had", gBuffSize);
  // fill with the FITS file
  if (fileName.EndsWith(".fits")) LoadFITSFile(fileName);
}

////////////////////////////////////////////////////////////////
//
// initialise empty elements, function to be called by both constructors
//  
void IactEventListIrf::_initialize_me()
{
  fOnSample = NULL;
  fOffSample = NULL;
  fEpmin = gEpmin; 
  fEpmax = gEpmax; 
  fTau = 1; 
  fDTau = 0; 
  fObsTime = 0; 
  fHAeff = NULL; 
  fHAeffOff = NULL; 
  fGEreso = NULL;
  fGEbias = NULL; 
  fMigMatrix = NULL; 
  fHdNdEpBkg = NULL; 
  fHdNdEpFrg = NULL;
}

////////////////////////////////////////////////////////////////
// 
// Load the IactEventListIrf from a FITS input file 
// TODO : generalize this function to a Load() member able to 
// select between ROOT and FITS inputs
// 
void IactEventListIrf::LoadFITSFile(TString inputFileName)
{
  Info("LoadFITS", "loading the dataset in: %s", inputFileName.Data());
  
  // open all the Header Data Units (HDUs)
  TFITSHDU* hduOn = new TFITSHDU(inputFileName, 1);
  TFITSHDU* hduOff = new TFITSHDU(inputFileName, 2);
  TFITSHDU* hduAeff = new TFITSHDU(inputFileName, 3);
  TFITSHDU* hduMatrix = new TFITSHDU(inputFileName, 4);
  TFITSHDU* hduEbounds = new TFITSHDU(inputFileName, 5);
  
  // fetch the effective time and tau from the HEADER 
  Double_t obsTime = hduOn->GetKeywordValue("LIVETIME").Atof();
  Double_t acceptanceOn = hduOn->GetKeywordValue("ACC").Atof();
  Double_t acceptanceOff = hduOff->GetKeywordValue("ACC").Atof();
  Double_t tau = acceptanceOff / acceptanceOn;
  
  // set attributes
  SetObsTime(obsTime);
  SetTau(tau);
  // fill ON event list
  TVectorD* onEnergies = hduOn->GetTabRealVectorColumn(3);
  TVectorD* onTimes = hduOn->GetTabRealVectorColumn(0);
  for(Int_t i = 0; i < onEnergies->GetNrows(); ++i)
    FillOnEvent(onEnergies[0][i] * 1e3, 0., 0., 0., 0., onTimes[0][i], 0.);
  // fill OFF event list
  TVectorD* offEnergies = hduOff->GetTabRealVectorColumn(3);
  TVectorD* offTimes = hduOff->GetTabRealVectorColumn(0);
  for(Int_t i = 0; i < offEnergies->GetNrows(); ++i)
    FillOffEvent(offEnergies[0][i] * 1e3, 0., 0., 0., 0., offTimes[0][i], 0.);
  // set IRFs
  SetHAeff(GetHAeffFromHDU(hduAeff));
  SetMigMatrix(GetMigMatrixFromHDU(hduMatrix, hduEbounds));

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

void IactEventListIrf::Print(Option_t* o) const
{
  cout << "Energy Range: " <<  fEpmin << " - " << fEpmax << " GeV " << endl;
  cout << "ON/OFF Norm (tau): " <<  fTau << " +- " << fDTau << " (Prob.: " << fTauPValue << ")" << endl;
  cout << "Observation Time: " <<  fObsTime  << " s "  << endl;
}

////////////////////////////////////////////////////////////////
// 
// Plot an overview of the histograms of the data
//
void IactEventListIrf::PlotOverview(Bool_t logY)
{
  // histograms of the ON and OFF events
  Double_t onEnergy;
  Double_t offEnergy;
  SetOnBranchAddress("E", &onEnergy);
  SetOffBranchAddress("E", &offEnergy);
  Int_t nOnEntries = GetOnSample()->GetEntries();
  Int_t nOffEntries = GetOffSample()->GetEntries();
  
  // fill ON histogram 
  // take the estimated energy binning from the migration matrix
  Int_t nBinsEest = GetMigMatrix()->GetYaxis()->GetNbins();
  Double_t minLog10Eest = GetMigMatrix()->GetYaxis()->GetXmin();
  Double_t maxLog10Eest = GetMigMatrix()->GetYaxis()->GetXmax();
  
  // fill estimated energy historgams with ON and OFF counts
  TH1D *histoEestOn = new TH1D("hEestOn", "", nBinsEest, minLog10Eest, maxLog10Eest);
  TH1D *histoEestOff = new TH1D("hEestOff", "", nBinsEest, minLog10Eest, maxLog10Eest);
  for (Int_t i = 0; i < nOnEntries; ++i) {
    GetOnEntry(i);
    histoEestOn->Fill(TMath::Log10(onEnergy));
  }
  for (Int_t i = 0; i < nOffEntries; ++i) {
    GetOffEntry(i);
    histoEestOff->Fill(TMath::Log10(offEnergy));
  }
  histoEestOff->Scale(1 / GetTau());

  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 500); 
  c1->Divide(3, 1);
  TPad *p1 = (TPad *) c1->cd(1);
  if (logY) p1->SetLogy();
  histoEestOn->SetStats(0);
  histoEestOn->SetLineColor(46);
  histoEestOn->SetLineWidth(2);
  histoEestOff->SetStats(0);
  histoEestOff->SetLineColor(4);
  histoEestOff->SetLineWidth(2);
  histoEestOn->GetXaxis()->SetTitle("log_{10}(E' / GeV)");
  TLegend *legend = new TLegend(0.65, 0.65, 0.9, 0.9);
  legend->AddEntry(histoEestOn, "ON events", "le");
  legend->AddEntry(histoEestOff, "#tau x OFF events", "le");
  legend->SetTextSize(0.03);
  histoEestOn->Draw("E");
  histoEestOff->Draw("ESAME");
  legend->Draw("SAME");
  
  TPad *p2 = (TPad *) c1->cd(2);
  if (logY) p2->SetLogy();
  TH1F *hAeff = GetHAeff();
  hAeff->SetStats(0);
  hAeff->SetLineWidth(2);
  hAeff->SetLineColor(4);
  hAeff->GetXaxis()->SetTitle("log_{10}(E / GeV)");
  hAeff->GetYaxis()->SetTitle("A_{eff} / m^{2}");
  hAeff->Draw("E");
    
  TPad *p3 = (TPad *) c1->cd(3);
  TH2F *migMatrix = GetMigMatrix();
  migMatrix->SetStats(0);
  migMatrix->GetXaxis()->SetTitle("log_{10}(E' / GeV)");
  migMatrix->GetYaxis()->SetTitle("log_{10}(E / GeV)");
  migMatrix->Draw("COLZ");
}

////////////////////////////////////////////////////////////////
// 
// get the effective area from the Header Data Unit of the FITS
// NOTE : energies in the FITS files are given in TeV
// 
TH1F* GetHAeffFromHDU(TFITSHDU* hdu)
{
  TVectorD* lowEtrueEdges = hdu->GetTabRealVectorColumn("ENERG_LO");
  TVectorD* highEtrueEdges = hdu->GetTabRealVectorColumn("ENERG_HI");
  TVectorD* aeffValues = hdu->GetTabRealVectorColumn("SPECRESP");
  const Int_t nBinsAeff = aeffValues->GetNrows();
  Double_t log10EtrueEdges[nBinsAeff+1];
  for (Int_t i=0; i < nBinsAeff; ++i)
    log10EtrueEdges[i] = TMath::Log10(1e3 * lowEtrueEdges[0][i]);
  // the last bin edge is the last element of the highEtrueEdges
  log10EtrueEdges[nBinsAeff] = TMath::Log10(1e3 * highEtrueEdges[0][nBinsAeff-1]);
  // define the TH1F to be returned, log10(E / GeV) on the X axis
  TH1F* hAeff = new TH1F();
  hAeff->SetBins(nBinsAeff, log10EtrueEdges);
  for (Int_t i = 0; i < nBinsAeff; ++i)
    hAeff->SetBinContent(i, aeffValues[0][i]);
  return hAeff;
}

////////////////////////////////////////////////////////////////
// 
// get the migration matrix from the Header Data Unit of the FITS
// NOTE : energies in the FITS files are given in TeV
// 
TH2F* GetMigMatrixFromHDU(TFITSHDU* hduMatrix, TFITSHDU* hduEbounds)
{
  // note energies in FITS files are in TeV
  TVectorD* lowEestEdges = hduEbounds->GetTabRealVectorColumn("E_MIN");
  TVectorD* highEestEdges = hduEbounds->GetTabRealVectorColumn("E_MAX");
  TVectorD* lowEtrueEdges = hduMatrix->GetTabRealVectorColumn("ENERG_LO"); 
  TVectorD* highEtrueEdges = hduMatrix->GetTabRealVectorColumn("ENERG_HI");

  //-> set migration matrix bins and edges
  
  // number of bins
  const Int_t nBinsEtrue = lowEtrueEdges->GetNrows();
  const Int_t nBinsEest = lowEestEdges->GetNrows();
  
  // true energy edges
  Double_t log10EtrueEdges[nBinsEtrue+1];
  for (Int_t i = 0; i < nBinsEtrue; ++i)
    log10EtrueEdges[i] = TMath::Log10(1e3 * lowEtrueEdges[0][i]);
  log10EtrueEdges[nBinsEtrue] = TMath::Log10(1e3 * highEtrueEdges[0][nBinsEtrue-1]);
  
  // estimated energy edges
  Double_t log10EestEdges[nBinsEest+1];
  for (Int_t i = 0;i < nBinsEest; ++i)
    log10EestEdges[i] = TMath::Log10(1e3 * lowEestEdges[0][i]);
  log10EestEdges[nBinsEest] = TMath::Log10(1e3 * highEestEdges[0][nBinsEest-1]);
  
  // set the matrix content
  TH2F *migMatrix = new TH2F();
  migMatrix->SetBins(nBinsEtrue, log10EestEdges, nBinsEtrue, log10EtrueEdges);
  
  // fill the matrix looping through the rows
  for (Int_t rownum = 0; rownum < nBinsEtrue; rownum++) {
    TArrayD *fChanArray = hduMatrix->GetTabVarLengthVectorCell(rownum, "F_CHAN");
    TArrayD *nChanArray = hduMatrix->GetTabVarLengthVectorCell(rownum, "N_CHAN"); 
    TArrayD *matrix = hduMatrix->GetTabVarLengthVectorCell(rownum, "MATRIX"); 
    // loop in the variable length arrays
    Int_t lastIndex = 0;
    for (Int_t i = 0; i < fChanArray->GetSize(); ++i) {
      // starting column index in which the migration is not null
      Int_t colnum = (Int_t) fChanArray->At(i); 
      // column span, after the starting column index, of non-null values
      Int_t arrayIndex = (Int_t) nChanArray->At(i);
      for (Int_t j = lastIndex; j < arrayIndex + lastIndex; ++j) {
        migMatrix->SetBinContent(rownum, colnum + j, matrix->At(j));
      }
      lastIndex += arrayIndex;
    }
  }
  return migMatrix;
}
