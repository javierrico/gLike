// collection of functions (suite) to test the properties of the the JointLkl.h
// class
#include "JointLkl.h"
#include "PoissonLkl.h"
#include "TMath.h"
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
//
// check if the minimum and the error of g returned by two joint PoissonLkl are
// consistent with the expected values gMinExpected and gErrExpected
//
void checkJointPoissonLklMinAndErr(Int_t nSamples, Int_t *nOn, Int_t *nOff,
                                   Double_t *tau, Double_t errorDef,
                                   Double_t unitsOfG, Double_t gMinExpected,
                                   Double_t gErrExpected,
                                   Bool_t &assertCloseMin,
                                   Bool_t &assertCloseErr) {
  const Int_t nTerms = nSamples;
  PoissonLkl *poissLkl[nTerms];
  JointLkl *jointLkl = new JointLkl();
  jointLkl->SetErrorDef(errorDef);
  jointLkl->SetUnitsOfG(unitsOfG);
  for (Int_t i = 0; i < nTerms; ++i) {
    poissLkl[i] = new PoissonLkl(nOn[i], nOff[i], tau[i]);
    jointLkl->AddSample(poissLkl[i]);
  }
  jointLkl->ComputeLklVsG();
  Double_t gMin = jointLkl->GetGLklMin();
  Double_t gErr = jointLkl->GetGLklMinErr();
  assertCloseMin = TMath::AreEqualRel(gMin, gMinExpected, 1e-7);
  assertCloseErr = TMath::AreEqualRel(gErr, gErrExpected, 1e-7);
  delete jointLkl;
}

////////////////////////////////////////////////////////////////////////////////
//
// assert if the minimum and the error on g of two joint PoissonLkl are the same
// as the value returned by gLike v00.08.00, obtained as:
//
// const Int_t nSamples = 2;
// Int_t nOn[nSamples] {130, 90};
// Int_t nOff[nSamples] = {105, 110};
// Double_t tau[nSamples] = {1., 1.};
// Double_t errorDef = 4;
// Double_t unitsOfG = 2;
// PoissonLkl *poissLkl[nSamples];
// JointLkl *jointLkl = new JointLkl();
// jointLkl->SetErrorDef(errorDef);
// jointLkl->SetUnitsOfG(unitsOfG);
// for (Int_t i = 0; i < nSamples; ++i) {
//   poissLkl[i] = new PoissonLkl(nOn[i], nOff[i], tau[i]);
//   jointLkl->AddSample(p[i]);
// }
// jointLkl->ComputeLklVsG();
// jointLkl->GetGLklMin();
// jointLkl->GetGLklMinErr();
//
void testJointPoissonLkl() {
  TString testName = "TEST::testJointPoissonLkl ";
  cout << testName + "Info: running test on JointLkl using two PoissonLkl"
       << endl;
  // input values
  const Int_t nSamples = 2;
  Int_t nOn[nSamples] = {130, 90};
  Int_t nOff[nSamples] = {105, 110};
  Double_t tau[nSamples] = {1., 1.};
  Double_t errorDef = 4;
  Double_t unitsOfG = 2;
  // expected output values
  Double_t gMinExpected = 0.84559242;
  Double_t gErrExpected = 20.825172;
  // check consistency
  Bool_t assertCloseMin = kFALSE;
  Bool_t assertCloseErr = kFALSE;
  checkJointPoissonLklMinAndErr(nSamples, nOn, nOff, tau, errorDef, unitsOfG,
                                gMinExpected, gErrExpected, assertCloseMin,
                                assertCloseErr);

  cout << testName + "Info: asserting the minimum of g" << endl;
  assert(assertCloseMin);
  cout << testName + "Info: correct assertion of the minimum of g!" << endl;

  cout << testName + "Info: asserting the error of g" << endl;
  assert(assertCloseErr);
  cout << testName + "Info: correct assertion of the error of g!" << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// TODO: add more tests on JointLkl.h here
//