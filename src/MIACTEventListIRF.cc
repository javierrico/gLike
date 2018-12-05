/* ======================================================================== *\
!
!   Author: Javier Rico         03/2017 <mailto:jrico@ifae.es>
!   Author: Joaquim Palacio     03/2017 <mailto:jpalacio@ifae.es>
!
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//						
// MIACTEventListIRF
//
// Class to hold the data and IRFs from IACTs to be used as input by
// MFullLkl class 
//
//////////////////////////////////////////////////////////////////////////////

#include "MIACTEventListIRF.h"

ClassImp(MIACTEventListIRF);


using namespace std;

static const Float_t gEpmin = 1e00; // [GeV] default value of minimum E_est
static const Float_t gEpmax = 1e06; // [GeV] default value of maximum E_est

const Float_t MIACTEventListIRF::gDefEVal      = 0.;    // default value when energy is not provided
const Float_t MIACTEventListIRF::gDefRADECVal  = 9999.; // default value when dRA and dDEC are not provided
const Float_t MIACTEventListIRF::gDefTVal      = -1.;    // default value when time is not provided
const Float_t MIACTEventListIRF::gDefHadVal    = -1.;    // default value when hadronness is not provided

////////////////////////////////////////////////////////////////
//
// Default constructor, just create empty ntuples for data samples
//
MIACTEventListIRF::MIACTEventListIRF(TString name,TString title) :
  TNamed(name,title), 
  fOnSample(NULL), fOffSample(NULL),
  fEpmin(gEpmin), fEpmax(gEpmax), fTau(1), fDTau(0),
  fObsTime(0), fHAeff(NULL), fHAeffOff(NULL), fGEreso(NULL),
  fGEbias(NULL), fMigMatrix(NULL), fHdNdEpBkg(NULL), fHdNdEpFrg(NULL)
{
  // create the event lists
  fOnSample  = new TNtuple("fOnSample", "On data set", "E:pointRA:pointDEC:dRA:dDEC:t:had");
  fOffSample = new TNtuple("fOffSample","Off data set","E:pointRA:pointDEC:dRA:dDEC:t:had");
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
MIACTEventListIRF::~MIACTEventListIRF()
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

void MIACTEventListIRF::Print(Option_t* o) const
{
  cout << "Energy Range: " <<  fEpmin << " - " << fEpmax << " GeV " << endl;
  cout << "ON/OFF Norm (tau): " <<  fTau << " +- " << fDTau << " (Prob.: " << fTauPValue << ")" << endl;
  cout << "Observation Time: " <<  fObsTime  << " s "  << endl;
}
