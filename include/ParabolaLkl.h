//////////////////////////////////////////////////////////////////////
// -2logL parabola
//////////////////////////////////////////////////////////////////////

#ifndef PARABOLALKL
#define PARABOLALKL

#include "Lkl.h"

class ParabolaLkl : public Lkl
{
 public:

  // constructor
  ParabolaLkl(Int_t n,Double_t* x,Double_t* y,TString name="",TString title="");
  ParabolaLkl(TString fileName,TString name="",TString title="");
  
  // destructor
  virtual ~ParabolaLkl();
  
  // getters
  Int_t   GetN()                           const {return GetLklVsG()->GetN();}
  TGraph* GetParabola(Bool_t units= kTRUE) const {return GetLklVsG(units);}
  

 protected:
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();
  
 private:  

  ClassDef(ParabolaLkl,1) // -2logL vs G
};

#endif
