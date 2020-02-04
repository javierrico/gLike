#include "testSuitePoissonLkl.cpp"
#include <assert.h>

void testPoissonLkl(){
  cout << "Running test on PoissonLkl" << endl;
  const Int_t nOn = 130;
  const Int_t nOff = 90;
  const Double_t tau = 1.;
  const Double_t errorDef = 4;
  const Double_t unitsOfG = 2;
  const Double_t gMinExpected = 79.583416;
  const Double_t gErrExpected = 60.483922;
  cout << "asserting the minimum of the likelihood:" << endl;
  assert(assertPoissonLklMin(nOn, nOff, tau, errorDef, unitsOfG, gMinExpected));
  cout << "asserting the error on the likelihood:" << endl;
  assert(assertPoissonLklErr(nOn, nOff, tau, errorDef, unitsOfG, gErrExpected));
}

int main(){
  // run all the tests
  testPoissonLkl(); 

}
