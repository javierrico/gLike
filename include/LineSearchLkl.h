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

  virtual Int_t SimulateDataSamples(Float_t meanGwithUnits=0,TRandom* rdm=NULL);

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();

  Int_t ComputeBkgModelFromOnHisto();

 private:

  Int_t fEventsInEnergyWindow;

  ClassDef(LineSearchLkl,1) // Line Search Likelihood
};

#endif
