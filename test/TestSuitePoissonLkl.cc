// collection of functions (suite) to test the properties of the the
// PoissonLkl.h class
#include "PoissonLkl.h"
#include "TMath.h"
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
//
// check if the minimum and the error of g returned by the PoissonLkl are
// consistent with the expected values gMinExpected and gErrExpected
//
void checkPoissonLklMinAndErr(Int_t nOn, Int_t nOff, Double_t tau,
                              Double_t errorDef, Double_t unitsOfG,
                              Double_t gMinExpected, Double_t gErrExpected,
                              Bool_t &assertCloseMin, Bool_t &assertCloseErr) {
  PoissonLkl *poissLkl = new PoissonLkl(nOn, nOff, tau);
  poissLkl->SetErrorDef(errorDef);
  poissLkl->SetUnitsOfG(unitsOfG);
  poissLkl->ComputeLklVsG();
  Double_t gMin = poissLkl->GetGLklMin();
  Double_t gErr = poissLkl->GetGLklMinErr();
  assertCloseMin = TMath::AreEqualRel(gMin, gMinExpected, 1e-7);
  assertCloseErr = TMath::AreEqualRel(gErr, gErrExpected, 1e-7);
  delete poissLkl;
}

////////////////////////////////////////////////////////////////////////////////
//
// assert if the minimum and the error on g of the PoissonLkl are the same as
// the value returned by gLike v00.08.00, obtained as:
//
// const Int_t nOn = 130;
// const Int_t nOff = 90;
// const Double_t tau = 1.;
// const Double_t errorDef = 4;
// const Double_t unitsOfG = 2;
// PoissonLkl* poissLkl = new PoissonLkl(nOn, nOff, tau);
// poissLkl->SetErrorDef(errorDef);
// poissLkl->SetUnitsOfG(unitsOfG);
// poissLkl->ComputeLklVsG();
// poissLkl->GetGLklMin();
// poissLkl->GetGLklMinErr();
//
void testPoissonLkl() {
  TString testName = "TEST::testPoissonLkl ";
  cout << testName + "Info: running test on PoissonLkl" << endl;
  // input values
  Int_t nOn = 130;
  Int_t nOff = 90;
  Double_t tau = 1.;
  Double_t errorDef = 4;
  Double_t unitsOfG = 2;
  // expected output values
  Double_t gMinExpected = 79.583416;
  Double_t gErrExpected = 60.483922;
  // check consistency
  Bool_t assertCloseMin = kFALSE;
  Bool_t assertCloseErr = kFALSE;
  checkPoissonLklMinAndErr(nOn, nOff, tau, errorDef, unitsOfG, gMinExpected,
                           gErrExpected, assertCloseMin, assertCloseErr);

  cout << testName + "Info: asserting the minimum of g" << endl;
  assert(assertCloseMin);
  cout << testName + "Info: correct assertion of the minimum of g!" << endl;

  cout << testName + "Info: asserting the error of g" << endl;
  assert(assertCloseErr);
  cout << testName + "Info: correct assertion of the error of g!" << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// TODO: add more tests on PoissonLkl.h here
//