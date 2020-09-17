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

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();

  Int_t ComputeBkgModelFromOnHisto();

 private:

  ClassDef(LineSearchLkl,1) // Line Search Likelihood
};

#endif
