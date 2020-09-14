//////////////////////////////////////////////////////////////////////
// Likelihood for line search
//////////////////////////////////////////////////////////////////////

#ifndef LINESEARCHLKL
#define LINESEARCHLKL

#include "TCanvas.h"

#include "Iact1dUnbinnedLkl.h"

class LineSearchLkl : public Iact1dUnbinnedLkl
{
 public:
  
  // constructors
  LineSearchLkl(TString inputString="");
  
  // destructor
  virtual ~LineSearchLkl();

  inline  TH1F*    GetHdNdEpBkg()             const {return fHdNdEpBkg;}
  inline  TH1F*    GetHdNdEpSignal()          const {return fHdNdEpSignal;}

  inline  Double_t GetdNdEpSignalIntegral()   {CheckHistograms(kFALSE); if(!fHdNdEpSignal) return 0; return fHdNdEpSignal->GetBinContent(0);}
  Int_t ReaddNdEpSignal(TString filename);
  virtual TH1F* GetHdNdEpModelBkg(Bool_t isDifferential=kTRUE,Int_t nbins=0) const;
  virtual TH1F* GetHdNdEpOn(Bool_t isDifferential=kTRUE,Int_t nbins=0)  const; 

  // Plots
  TCanvas* PlotHistosAndData();

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();

  Int_t ComputeBkgModelFromOnHisto();

  Int_t CheckHistograms(Bool_t checkdNdEpBkg=kTRUE);
  Int_t ResetdNdESignal();
  Int_t SetdNdESignalFunction(TString function,Float_t p0=0,Float_t p1=0,Float_t p2=0,Float_t p3=0,Float_t p4=0,Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  Int_t SetdNdESignalFunction(TF1* function,Float_t emin=0,Float_t emax=1e9);

 private:
  TH1F*    fHdNdEpBkg;              //-> dN/dE'dt vs E' for background events (normalized)
  TH1F*    fHdNdEpSignal;           //-> dN/dE' vs E' histogram for signal events (normalized)

  ClassDef(LineSearchLkl,1) // Line Search Likelihood
};

#endif
