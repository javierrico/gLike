#include "TMath.h"
#include "PoissonLkl.h"


Bool_t assertPoissonLklMin(Int_t nOn, Int_t nOff, Double_t tau, Double_t errorDef, Double_t unitsOfG, Double_t gMinExpected){
  PoissonLkl* poissLkl = new PoissonLkl(nOn, nOff, tau);  
  poissLkl->SetErrorDef(errorDef);
  poissLkl->SetUnitsOfG(unitsOfG);
  poissLkl->ComputeLklVsG();
  Double_t gMin = poissLkl->GetGLklMin();
  Bool_t assertClose = TMath::AreEqualRel(gMin, gMinExpected, 1e-7);
  delete poissLkl;
  return assertClose;
}


Bool_t assertPoissonLklErr(Int_t nOn, Int_t nOff, Double_t tau, Double_t errorDef, Double_t unitsOfG, Double_t gErrExpected){
  PoissonLkl* poissLkl = new PoissonLkl(nOn, nOff, tau);  
  poissLkl->SetErrorDef(errorDef);
  poissLkl->SetUnitsOfG(unitsOfG);
  poissLkl->ComputeLklVsG();
  Double_t gErr = poissLkl->GetGLklMinErr();
  Bool_t assertClose = TMath::AreEqualRel(gErr, gErrExpected, 1e-7);
  delete poissLkl;
  return assertClose;
}



