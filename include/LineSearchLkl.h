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

  inline  Float_t  GetRelativePeakIntensity() const {return fRelativePeakIntensity;};
  inline  Float_t  GetBkgRegionWigth()        const {return fBkgRegionWidth;};
  inline  TH1F*    GetHdNdEpBkg()             const {return fHdNdEpBkg;}
  inline  TH1F*    GetHdNdEpSignal()          const {return fHdNdEpSignal;}
  
  //new implemented 04/09/2020
  //inline  TF1*    GetfFEventbkg()          const {return fFEventbkg;}

  inline  Double_t GetdNdEpSignalIntegral()   {CheckHistograms(kFALSE); if(!fHdNdEpSignal) return 0; return fHdNdEpSignal->GetBinContent(0);}
  Int_t ReaddNdEpSignal(TString filename);
  virtual TH1F* GetHdNdEpModelBkg(Bool_t isDifferential=kTRUE,Int_t nbins=0) const;
  virtual TH1F*   GetHdNdEpOn(Bool_t isDifferential=kTRUE,Int_t nbins=0)  const; // Tomo
  inline  void  SetRelativePeakIntensity(Float_t relativeIntensity) {fRelativePeakIntensity = relativeIntensity;};
  inline  void  SetBkgRegionWidth(Float_t bkgWidth) {fBkgRegionWidth = bkgWidth;};

  //newly implemeted on 27th July
  Int_t SetEnergyWindow(Int_t DMmass, Float_t windowsize);
  // Plots
  TCanvas* PlotHistosAndData();

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0.);
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
  //TF1*     fFEventbkg;              //-> fitting function to event distribution (part of fOnsample) inside slinding window  
  Float_t  fRelativePeakIntensity;  // Relative peak intensity (compare to signal maximum) to define regions where the backrground is estimated. By default it is 10%.
  Float_t  fBkgRegionWidth;         // Relative peak intensity (compare to signal maximum) to define regions where the backrground is estimated. By default it is 10%.

  ClassDef(LineSearchLkl,1) // Line Search Likelihood
};

#endif
