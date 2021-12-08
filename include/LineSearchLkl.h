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

  // Plots
  TCanvas* PlotHistosAndData();

  // getters
  inline Int_t GetEventsInEnergyWindow()             const {return fEventsInEnergyWindow;}
  inline Float_t  GetEminWindow()                    const {return fEpmin_window;}
  inline Float_t  GetEmaxWindow()                    const {return fEpmax_window;}
  inline Float_t  GetEwindowWidth()                  const {return fEwindowWidth;}
  inline Float_t  GetEwindowThreshold()              const {return fEwindowThreshold;}

  virtual Int_t SimulateDataSamples(Float_t meanGwithUnits=0,TRandom* rdm=NULL);

  inline Int_t     SetEpminWindow(Float_t emin) {fEpmin_window   = emin; SetChecked(kFALSE); return 0;}
  inline Int_t     SetEpmaxWindow(Float_t emax) {fEpmax_window   = emax; SetChecked(kFALSE); return 0;}

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();

  Int_t ComputeBkgModelFromOnHisto();

 private:

  Int_t fEventsInEnergyWindow;

  Float_t  fEpmin_window;           // [GeV] Minimum measured energy of considered events
  Float_t  fEpmax_window;           // [GeV] Maximum measured energy of considered events
  Float_t  fEwindowWidth;           // [sigma of the energy resolution] width of the window, ie the signal region, around the mass for which a line is searched
  Float_t  fEwindowThreshold;       // [GeV] minimum energy of the window (used to avoid thresold effect at low energy)

  ClassDef(LineSearchLkl,1) // Line Search Likelihood
};

#endif
