//////////////////////////////////////////////////////////////////////
// Unbinned full likelihood
//////////////////////////////////////////////////////////////////////

#ifndef FULLLKL
#define FULLLKL

#include "TH2.h"
#include "TString.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"

#include "Lkl.h"

static const Float_t  gEpmin           = 1e01;                   // [GeV] default value of minimum E_est
static const Float_t  gEpmax           = 1e06;                   // [GeV] default value of maximum E_est

class Iact1dUnbinnedLkl : public virtual Lkl
{
 public:
  // enumerations
  enum parIndex_t {gGIndex,gBIndex,gTauIndex};                    // Indeces of parameters
  
  // constructors
  Iact1dUnbinnedLkl(TString inputString="");
  
  // destructor
  virtual ~Iact1dUnbinnedLkl();
    
  // specify obs time, mass and log10 of Jfactor and its 1-sigma error (nuisance parameter)
  virtual void SetDMAnnihilationUnitsForG(Double_t mass)
  {
    Double_t J     = TMath::Power(10,fLogJ);
    Double_t den   = J*fObsTime*GetdNdEpSignalIntegral()/(8*TMath::Pi()*mass*mass);
    Double_t units = (den>0? 1./den: 0);
    SetUnitsOfG(units);
  }
  
  virtual void SetDMDecayUnitsForG(Double_t mass)
  {
    Double_t J     = TMath::Power(10,fLogJ);
    Double_t den   = J*fObsTime*GetdNdEpSignalIntegral()/(4*TMath::Pi()*mass);
    Double_t units = (den>0? 1./den: 0);
    SetUnitsOfG(units);
  }

  // configure internal histos binning
  inline void SetNFineBins(Int_t nfinebins)    {fNFineBins = nfinebins;}
  inline void SetFineLEMin(Double_t finelemin) {fFineLEMin = finelemin;}
  inline void SetFineLEMax(Double_t finelemax) {fFineLEMax = finelemax;}

  // getters
  inline       Float_t  GetEmin()             const {return fEpmin;} 
  inline       Float_t  GetEmax()             const {return fEpmax;}
  inline       Float_t  GetNFineBins()        const {return fNFineBins;}
  inline       Float_t  GetFineLEMin()        const {return fFineLEMin;}
  inline       Float_t  GetFineLEMax()        const {return fFineLEMax;}
  inline       UInt_t   GetNon()              const {return fNon;}
  inline       UInt_t   GetNoff()             const {return fNoff;}
  inline       Float_t* GetOnSample()         const {return fOnSample;}
  inline       Float_t* GetOffSample()        const {return fOffSample;}
  inline       Float_t  GetTau()              const {return fTau;}
  inline       Float_t  GetDTau()             const {return fDTau;}
  inline       Float_t  GetTrueTau()          const {return fTrueTau;}
  inline       Float_t  GetObsTime()          const {return fObsTime;}
  inline       Float_t  GetDataTau()          const {return fDataTau;}
  inline       Float_t  GetDataDTau()         const {return fDataDTau;}
  inline       Float_t  GetDataObsTime()      const {return fDataObsTime;}
  inline       TH1F*    GetHAeff()            const {return fHAeff;}
  inline       TH1F*    GetHAeffOff()         const {return fHAeffOff;}
  virtual inline TH1F*  GetHdNdEpBkg()        const {return fHdNdEpBkg;}
  inline       TH1F*    GetHdNdEpFrg()        const {return fHdNdEpFrg;}
  virtual inline TH1F*  GetHdNdEpSignal()     const {return fHdNdEpSignal;}
  inline       TH1F*    GetHdNdEpSignalOff()  const {return fHdNdEpSignalOff;}
  inline       TH1F*    GetHdNdESignal()      const {return fHdNdESignal;}
  inline       TGraph*  GetGEreso()           const {return fGEreso;}
  inline       TGraph*  GetGEbias()           const {return fGEbias;}
  inline       TH2F*    GetMigMatrix()        const {return fMigMatrix;}
  inline       Float_t  GetLogJ()             const {return fLogJ;}

  inline Double_t GetdNdEpBkgIntegral()       {if(!fHdNdEpBkg) return 0; NormalizedNdEHisto(fHdNdEpBkg); return fHdNdEpBkg->GetBinContent(0);}
  inline Double_t GetdNdEpFrgIntegral()       {if(!fHdNdEpFrg) return 0; NormalizedNdEHisto(fHdNdEpFrg); return fHdNdEpFrg->GetBinContent(0);}
  inline Double_t GetdNdEpSignalIntegral()    {CheckHistograms(kFALSE); if(!fHdNdEpSignal) return 0; return fHdNdEpSignal->GetBinContent(0);}
  inline Double_t GetdNdEpSignalOffIntegral() {if(!fHdNdEpSignalOff) return 0; return fHdNdEpSignalOff->GetBinContent(0);}
  inline Double_t GetdNdESignalIntegral()     {if(!fHdNdESignal) return 0; NormalizedNdEHisto(fHdNdESignal); return fHdNdESignal->GetBinContent(0);}
  virtual TH1F*   GetHdNdEpOn(Bool_t isDifferential=kTRUE,Int_t nbins=0)  const;
  virtual TH1F*   GetHdNdEpOff(Bool_t isDifferential=kTRUE,Int_t nbins=0) const;
  
  // Read input dN/dE files and related functions
  virtual Int_t ResetdNdESignal();
  Int_t SetdNdESignal(TH1F* hdNdESignal);
  Int_t AdddNdESignal(TString filename,Float_t br=1.0);
  virtual Int_t SetdNdESignalFunction(TString function,Float_t p0=0,Float_t p1=0,Float_t p2=0,Float_t p3=0,Float_t p4=0,Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  Int_t AdddNdESignalFunction(TString function,Float_t p0=0,Float_t p1=0,Float_t p2=0,Float_t p3=0,Float_t p4=0,Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  virtual Int_t SetdNdESignalFunction(TF1* function,Float_t emin=0,Float_t emax=1e9);
  Int_t AdddNdESignalFunction(TF1* function,Float_t emin=0,Float_t emax=1e9,Float_t br=1.0);
  Int_t SetTrueTau(Float_t truetau) {fTrueTau=truetau; return 0;}
  Int_t ReaddNdESignal(TString filename);
  virtual Int_t ReaddNdEpSignal(TString filename);
  Int_t ReaddNdEpSignalOff(TString filename);
  Int_t ReadCTAIRF(TString filename);
  Int_t TransformAndSavedNdEpBkg(TH1F* provHNOff,Bool_t interpolate=kTRUE,Double_t scale=0,Bool_t isDiff=kTRUE);
  Int_t TransformAndSavedNdEpFrg(TH1F* provHNOff,Bool_t interpolate=kTRUE,Double_t scale=0,Bool_t isDiff=kTRUE);
  Int_t ResetHdNdEpBkg() {if(fHdNdEpBkg) delete fHdNdEpBkg; fHdNdEpBkg=NULL; return 0;}

  virtual Int_t SimulateDataSamples(UInt_t seed=0,Float_t meanG=0);

  // print data in the overview
  virtual void PrintData(Int_t level=0);
  
  // Plots
  virtual TCanvas* PlotHistosAndData();

  Int_t    NormalizedNdEHisto(TH1F* histo);

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual Int_t    GetRealBkgAndGoffHistos(TRandom3* rdm,TH1F*& hdNdEpBkg,TH1F*& hdNdEpSignalOff);
  
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();

  inline Int_t    SetEpmin(Float_t emin) {fEpmin   = emin;  return 0;}
  inline Int_t    SetEpmax(Float_t emax) {fEpmax   = emax;  return 0;}
  inline Int_t    SetObsTime(Float_t t)  {fObsTime = t;     return 0;}  
  inline Int_t    SetLogJ(Float_t logJ)  {fLogJ    = logJ;  return 0;}
  
  Int_t    SetAeff(TH1F* hProvAeff);
  Int_t    SetAeffOff(TH1F* hProvAeff);
  Int_t    SetEResoAndBias(TGraph* ereso,TGraph* ebias);
  Int_t    SetMigMatrix(TH2F* provMigMatrix);
  Int_t    SetdNdEpBkg(TH1F* hdNdEpBkg) {return TransformAndSavedNdEpBkg(hdNdEpBkg,kFALSE);}
  Int_t    SetdNdEpFrg(TH1F* hdNdEpFrg) {return TransformAndSavedNdEpFrg(hdNdEpFrg,kFALSE);}
  Int_t    CheckHistograms(Bool_t checkdNdEpBkg=kTRUE);
  Int_t    CheckEnergyLimits() const;
  
 private:  

  // Elements to be supplied by user
  Float_t  fEpmin;           // [GeV] Minimum measured energy of considered events
  Float_t  fEpmax;           // [GeV] Maximum measured energy of considered events

  UInt_t   fNon;             // Number of measured events in On region
  UInt_t   fNoff;            // Number of measured events in Off region
  Float_t* fOnSample;        //-> array of measured energies for On events 
  Float_t* fOffSample;       //-> array of measured energies for Off events 

  Float_t  fDataTau;         // Data normalization Noff/Non (e.g # of off regions). Equal to fTau except for simulations and tests
  Float_t  fDataDTau;        // Data statistical error in fTau (used for simualtions and in the likelihood). Equal to fDTau except for simulations and tests
  Float_t  fDataObsTime;     // Data [s] observation time. Equal to fObsTime except for simulations and tests

  Float_t  fTau;             // Assumed normalization Noff/Non (e.g # of off regions). Equal to fDataTau except for simulations and tests
  Float_t  fDTau;            // Assumed statistical error in fTau (used for simualtions and in the likelihood). Equal to fDataDTau except for simulations and tests
  Float_t  fTrueTau;         // True value of Tau used in simulations (equal to fTau if fDTau=0, if no simulations are performed its value is 0)
  Float_t  fObsTime;         // Assumed [s] observation time. Equal to fDataObsTime except for simulations and tests
  Float_t  fLogJ;            // log[GeV^2 cm^-5] or log[GeV cm^2] (ann/dec respectively), log of J factor
  Bool_t   fIsOffAsOn;       // if kTRUE, use the Off event sample as Off and as On
  
  Int_t    fNFineBins;       // number of fine bins for internal histos
  Double_t fFineLEMin;       // minimum log(energy[GeV]) for internal histos
  Double_t fFineLEMax;       // maximum log(energy[GeV]) for internal histos

  TH1F*    fHdNdESignal ;    //-> dN/dE vs E histogram for signal events
  TH1F*    fHAeff       ;    //-> Effective area vs E histogram for signal events
  TH1F*    fHAeffOff    ;    //-> Effective area vs E histogram for signal events IN THE OFF REGION
  TGraph*  fGEreso      ;    //-> Graph with energy resolution vs energy values
  TGraph*  fGEbias      ;    //-> Graph with energy bias vs energy values
  TH2F*    fMigMatrix   ;    //-> Migration matrix (overrides fGEreso and fGEbias)

  TH1F*    fHdNdEpBkg   ;    //-> dN/dE'dt vs E' for background events (normalized)
  TH1F*    fHdNdEpFrg   ;    //-> dN/dE'   vs E' for foreground (gamma-events from a different source) events (normalized)
  TH1F*    fHdNdEpSignal;    //-> dN/dE' vs E' histogram for signal events (normalized)
  TH1F*    fHdNdEpSignalOff; //-> dN/dE' vs E' histogram for signal events in the off region (normalized)

  // Flags internally used by the class

  ClassDef(Iact1dUnbinnedLkl,1) // Full Likelihood (unbinned)
};

#endif
