/////////////////////////////////////////////////////////////////////////////
// Class to hold the data and IRFs from IACTs to be used as input by Iact1dUnbinnedLkl class 
/////////////////////////////////////////////////////////////////////////////

#ifndef IACTEVENTLISTIRF
#define IACTEVENTLISTIRF

#include <iostream>

#include "TNamed.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TNtupleD.h"
#include "TVectorT.h"


// structure to read the events from the NTuple
typedef struct {
  Double_t  E;        // [GeV] measured energy of event
  Double_t  pointRA;  // [deg] RA of telescope pointing direction
  Double_t  pointDEC; // [deg] DEC of telescope pointing direction
  Double_t  dRA;      // [deg] distance to pointing direction in the RA axis
  Double_t  dDEC;     // [deg] distance to pointing direction in the DEC axis
  Double_t  t;        // [MJD] arrival time
  Double_t  had;      // [1]   hadronness
} IactEvent_t;

// class to hold the list of events and their corresponding IRFs
class IactEventListIrf : public TNamed
{
 public:  

  // static constants
  static const Double_t gDefEVal      ;    //! default value when energy is not provided
  static const Double_t gDefRADECVal  ;    //! default value when dRA and dDEC are not provided
  static const Double_t gDefTVal      ;    //! default value when time is not provided
  static const Double_t gDefHadVal    ;    //! default value when hadronness is not provided

  // constructor
  void Initialize(); // initialize empty resources, used by both constructors
  IactEventListIrf(TString name="IactEventListIrf", TString title="");
  IactEventListIrf(TString fileName, TString name, TString title);

  // destructor
  virtual ~IactEventListIrf();

  // load ON, OFF event lists and IRFs from a ROOT or FITS input file
  void LoadFITSFile(TString inputFileName);
  
  // fill event with input values
  void FillOnEvent(Double_t E=gDefEVal, Double_t pointRA=gDefRADECVal, Double_t pointDEC=gDefRADECVal, Double_t dRA=gDefRADECVal, Double_t dDEC=gDefRADECVal, Double_t t=gDefTVal, Double_t had=gDefHadVal)
  {fOnSample->Fill(E,pointRA,pointDEC,dRA,dDEC,t,had);}
  
  void FillOffEvent(Double_t E=gDefEVal, Double_t pointRA=gDefRADECVal, Double_t pointDEC=gDefRADECVal, Double_t dRA=gDefRADECVal, Double_t dDEC=gDefRADECVal, Double_t t=gDefTVal, Double_t had=gDefHadVal)
  {fOffSample->Fill(E,pointRA,pointDEC,dRA,dDEC,t,had);}

  // plot a resume of the input
  void PlotOverview(Bool_t logY = kTRUE);
    
  // setters
  void SetTau(Double_t tau, Double_t dTau=0,Double_t dPValue=-1.)  {fTau=tau; fDTau=dTau; fTauPValue=dPValue;}
  void SetEnergyRange(Double_t Epmin, Double_t Epmax)              {fEpmin=Epmin; fEpmax=Epmax;}
  void SetEpmin(Double_t Epmin)                                    {fEpmin=Epmin;}
  void SetEpmax(Double_t Epmax)                                    {fEpmax=Epmax;}
  void SetObsTime(Double_t obsTime)                                {fObsTime  = obsTime;}
  void SetHAeff(const TH1F* hAeff)                                 {if(fHAeff)     delete fHAeff;     fHAeff     = new TH1F(*hAeff);     fHAeff->SetName("hAeff");}
  void SetHAeffOff(const TH1F* hAeffOff)                           {if(fHAeffOff)  delete fHAeffOff;  fHAeffOff  = new TH1F(*hAeffOff);  fHAeffOff->SetName("hAeffOff");}
  void SetGEResoAndBias(const TGraph* gEreso,const TGraph* gEbias) {if(fGEreso)    delete fGEreso;    fGEreso    = new TGraph(*gEreso);  fGEreso->SetName("gEreso");
                                                                    if(fGEbias)    delete fGEbias;    fGEbias    = new TGraph(*gEbias);  fGEbias->SetName("gEbias");}     
  void SetMigMatrix(const TH2F* migMatrix)                         {if(fMigMatrix) delete fMigMatrix; fMigMatrix = new TH2F(*migMatrix); fMigMatrix->SetName("migMatrix");}
  void SetHdNdEpBkg(const TH1F* hdNdEpBkg)                         {if(fHdNdEpBkg) delete fHdNdEpBkg; fHdNdEpBkg = new TH1F(*hdNdEpBkg); fHdNdEpBkg->SetName("hdNdEpBkg");}
  void SetHdNdEpFrg(const TH1F* hdNdEpFrg)                         {if(fHdNdEpFrg) delete fHdNdEpFrg; fHdNdEpFrg = new TH1F(*hdNdEpFrg); fHdNdEpFrg->SetName("hdNdEpFrg");}

  void SetOnBranchAddress(const char* bname, Double_t* add) {fOnSample->SetBranchAddress(bname,add);}
  void SetOffBranchAddress(const char* bname, Double_t* add) {fOffSample->SetBranchAddress(bname,add);}

  void Print(Option_t* o="") const;
  
  // getters
  TNtupleD* GetOnSample()          {return fOnSample;}
  TNtupleD* GetOffSample()         {return fOffSample;}
  Int_t     GetOnEntry(Int_t iev)  {return fOnSample->GetEntry(iev);}
  Int_t     GetOffEntry(Int_t iev) {return fOffSample->GetEntry(iev);}

  Double_t  GetEpmin()     {return fEpmin;}
  Double_t  GetEpmax()     {return fEpmax;}
  Double_t  GetTau()       {return fTau;}
  Double_t  GetDTau()      {return fDTau;}
  Double_t  GetTauPValue() {return fTauPValue;}
  Double_t  GetObsTime()   {return fObsTime;}

  TH1F*     GetHAeff()     {return fHAeff;}
  TH1F*     GetHAeffOff()  {return fHAeffOff;}
  TGraph*   GetGEreso()    {return fGEreso;}
  TGraph*   GetGEbias()    {return fGEbias;}
  TH2F*     GetMigMatrix() {return fMigMatrix;}

  TH1F*     GetHdNdEpBkg() {return fHdNdEpBkg;}
  TH1F*     GetHdNdEpFrg() {return fHdNdEpFrg;}
  
 private:
  
  TNtupleD* fOnSample;        //-> ON  event list ("Energy:Point_RA:Point_DEC:Delta_RA:Delta_DEC:Arr_time:Hadronness")
  TNtupleD* fOffSample;       //-> OFF event list ("Energy:Point_RA:Point_DEC:Delta_RA:Delta_DEC:Arr_time:Hadronness")

  Double_t  fEpmin;           // [GeV] minimum estimated energy
  Double_t  fEpmax;           // [GeV] maximum estimated energy
  Double_t  fTau;             // normalization Noff/Non (e.g # of off regions)
  Double_t  fDTau;            // statistical error in fTau
  Double_t  fTauPValue;       // Probability value of agreement between On/Off
  Double_t  fObsTime;         // [s] observation time

  TH1F*     fHAeff       ;    //-> Effective area vs E histogram for signal events
  TH1F*     fHAeffOff    ;    //-> Effective area vs E histogram for signal events IN THE OFF REGION
  TGraph*   fGEreso      ;    //-> Graph with energy resolution vs energy values
  TGraph*   fGEbias      ;    //-> Graph with energy bias vs energy values
  TH2F*     fMigMatrix   ;    //-> Migration matrix (overrides fGEreso and fGEbias)

  TH1F*     fHdNdEpBkg   ;    //-> dN/dE'dt vs E' for background events (normalized)
  TH1F*     fHdNdEpFrg   ;    //-> dN/dE'dt vs E' for foreground (gamma-events from a different source) events (normalized)


  ClassDef(IactEventListIrf,2) // Class to hold the data and IRFs from IACTs to be used as input by Iact1dUnbinnedLkl class
};

#endif
