/////////////////////////////////////////////////////////////////////////////
// Class to hold the data and IRFs from IACTs to be used as input by MFullLkl class 
/////////////////////////////////////////////////////////////////////////////

#ifndef M_IACTEVENTLISTIRF
#define M_IACTEVENTLISTIRF

#include <iostream>

#include "TNamed.h"

#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TNtuple.h"

// structure to read the events from the NTuple
typedef struct {
  Float_t  E;        // [GeV] measured energy of event
  //  Float_t  pointRA;  // [deg] RA of telescope pointing direction
  //  Float_t  pointDEC; // [deg] DEC of telescope pointing direction
  Float_t  dRA;      // [deg] distance to pointing direction in the RA axis 
  Float_t  dDEC;     // [deg] distance to pointing direction in the DEC axis
  Float_t  t;        // [MJD] arrival time
  Float_t  had;      // [1]   hadronness
} MEvent_t;

// class to hold the list of events and their corresponding IRFs
class MIACTEventListIRF : public TNamed
{
 public:  

  // static constants
  static const Float_t gDefEVal      ;    //! default value when energy is not provided
  static const Float_t gDefRADECVal  ;    //! default value when dRA and dDEC are not provided
  static const Float_t gDefTVal      ;    //! default value when time is not provided
  static const Float_t gDefHadVal    ;    //! default value when hadronness is not provided

  // constructor
  MIACTEventListIRF(TString name="MIACTEventListIRF",TString title="");
  
  // destructor
  virtual ~MIACTEventListIRF();

  // fill event with input values
  void FillOnEvent( Float_t E=gDefEVal, Float_t pointRA=gDefRADECVal, Float_t pointDEC=gDefRADECVal, Float_t dRA=gDefRADECVal, Float_t dDEC=gDefRADECVal, Float_t t=gDefTVal, Float_t had=gDefHadVal)
  {fOnSample->Fill(E,pointRA,pointDEC,dRA,dDEC,t,had);}
  
  void FillOffEvent(Float_t E=gDefEVal, Float_t pointRA=gDefRADECVal, Float_t pointDEC=gDefRADECVal, Float_t dRA=gDefRADECVal, Float_t dDEC=gDefRADECVal, Float_t t=gDefTVal, Float_t had=gDefHadVal)
  {fOffSample->Fill(E,pointRA,pointDEC,dRA,dDEC,t,had);}
    
  // setters
  void SetTau(Float_t tau, Float_t dTau=0,Float_t dPValue=-1.)     {fTau=tau; fDTau=dTau; fTauPValue=dPValue;}
  void SetEnergyRange(Float_t Epmin, Float_t Epmax)                {fEpmin=Epmin; fEpmax=Epmax;}
  void SetEpmin(Float_t Epmin)                                     {fEpmin=Epmin;}
  void SetEpmax(Float_t Epmax)                                     {fEpmax=Epmax;}
  void SetObsTime(Float_t obsTime)                                 {fObsTime  = obsTime;}
  void SetHAeff(const TH1F* hAeff)                                 {if(fHAeff)     delete fHAeff;     fHAeff     = new TH1F(*hAeff);     fHAeff->SetName("hAeff");}
  void SetHAeffOff(const TH1F* hAeffOff)                           {if(fHAeffOff)  delete fHAeffOff;  fHAeffOff  = new TH1F(*hAeffOff);  fHAeffOff->SetName("hAeffOff");}
  void SetGEResoAndBias(const TGraph* gEreso,const TGraph* gEbias) {if(fGEreso)    delete fGEreso;    fGEreso    = new TGraph(*gEreso);  fGEreso->SetName("gEreso");
                                                                    if(fGEbias)    delete fGEbias;    fGEbias    = new TGraph(*gEbias);  fGEbias->SetName("gEbias");}     
  void SetMigMatrix(const TH2F* migMatrix)                         {if(fMigMatrix) delete fMigMatrix; fMigMatrix = new TH2F(*migMatrix); fMigMatrix->SetName("migMatrix");}
  void SetHdNdEpBkg(const TH1F* hdNdEpBkg)                         {if(fHdNdEpBkg) delete fHdNdEpBkg; fHdNdEpBkg = new TH1F(*hdNdEpBkg); fHdNdEpBkg->SetName("hdNdEpBkg");}
  void SetHdNdEpFrg(const TH1F* hdNdEpFrg)                         {if(fHdNdEpFrg) delete fHdNdEpFrg; fHdNdEpFrg = new TH1F(*hdNdEpFrg); fHdNdEpFrg->SetName("hdNdEpFrg");}

  void SetOnBranchAddress(const char* bname, Float_t* add) {fOnSample->SetBranchAddress(bname,add);}
  void SetOffBranchAddress(const char* bname, Float_t* add) {fOffSample->SetBranchAddress(bname,add);}

  void Print(Option_t* o="") const;
  
  // getters
  TNtuple* GetOnSample()          {return fOnSample;}
  TNtuple* GetOffSample()         {return fOffSample;}
  Int_t    GetOnEntry(Int_t iev)  {return fOnSample->GetEntry(iev);}
  Int_t    GetOffEntry(Int_t iev) {return fOffSample->GetEntry(iev);}    

  Float_t  GetEpmin()     {return fEpmin;}
  Float_t  GetEpmax()     {return fEpmax;}
  Float_t  GetTau()       {return fTau;}        
  Float_t  GetDTau()      {return fDTau;}
  Float_t  GetTauPValue() {return fTauPValue;}
  Float_t  GetObsTime()   {return fObsTime;}

  TH1F*    GetHAeff()     {return fHAeff;}  
  TH1F*    GetHAeffOff()  {return fHAeffOff;}
  TGraph*  GetGEreso()    {return fGEreso;}  
  TGraph*  GetGEbias()    {return fGEbias;}   
  TH2F*    GetMigMatrix() {return fMigMatrix;}

  TH1F*    GetHdNdEpBkg() {return fHdNdEpBkg;}
  TH1F*    GetHdNdEpFrg() {return fHdNdEpFrg;}

  
 private:
  
  TNtuple* fOnSample;        //-> ON  event list ("Energy:Point_RA:Point_DEC:Delta_RA:Delta_DEC:Arr_time:Hadronness")
  TNtuple* fOffSample;       //-> OFF event list ("Energy:Point_RA:Point_DEC:Delta_RA:Delta_DEC:Arr_time:Hadronness")

  Float_t  fEpmin;           // [GeV] minimum estimated energy
  Float_t  fEpmax;           // [GeV] maximum estimated energy
  Float_t  fTau;             // normalization Noff/Non (e.g # of off regions)
  Float_t  fDTau;            // statistical error in fTau
  Float_t  fTauPValue;       // Probability value of agreement between On/Off
  Float_t  fObsTime;         // [s] observation time

  TH1F*    fHAeff       ;    //-> Effective area vs E histogram for signal events
  TH1F*    fHAeffOff    ;    //-> Effective area vs E histogram for signal events IN THE OFF REGION
  TGraph*  fGEreso      ;    //-> Graph with energy resolution vs energy values
  TGraph*  fGEbias      ;    //-> Graph with energy bias vs energy values
  TH2F*    fMigMatrix   ;    //-> Migration matrix (overrides fGEreso and fGEbias)

  TH1F*    fHdNdEpBkg   ;    //-> dN/dE'dt vs E' for background events (normalized)
  TH1F*    fHdNdEpFrg   ;    //-> dN/dE'dt vs E' for foreground (gamma-events from a different source) events (normalized)

  ClassDef(MIACTEventListIRF,4) // Class to hold the data and IRFs from IACTs to be used as input by MFullLkl class
};

#endif
