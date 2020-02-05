// collection of functions (suite) to test the properties of the the
// IactEventListIrf.h class
#include "IactEventListIrf.h"
#include "TFile.h"
#include "TMath.h"

////////////////////////////////////////////////////////////////////////////////
//
// check if the expected number of entries in one of the IactEventList members
// is consistent with what expected.
// The function is overloaded for different cases e.g. TNtuples return Long_t
// when GetEntries() is used, while TH1* and TH2* return Dpuble_t
//
Bool_t checkMemberEntries(Int_t nEntries, Int_t nEntriesExpected) {
  return (nEntries == nEntriesExpected);
}
Bool_t checkMemberEntries(Long_t nEntries, Long_t nEntriesExpected) {
  return (nEntries == nEntriesExpected);
}
Bool_t checkMemberEntries(Double_t nEntries, Double_t nEntriesExpected) {
  // when used on TH1 and TH2 GetEntries() returns a double
  return TMath::AreEqualRel(nEntries, nEntriesExpected, 1e-7);
}

////////////////////////////////////////////////////////////////////////////////
//
// assert if the IactEventListIrf members have been properly intialised
// i.e. none of them corresponds to a null pointer
//
void assertMembersExistence(IactEventListIrf *iactDataSet) {
  TString testName = "TEST::assertMembersExistence ";
  cout << testName + "Info: asserting the individual members existence" << endl;

  cout << testName + "Info: asserting the On sample existence" << endl;
  assert(iactDataSet->GetOnSample() != NULL);
  cout << testName + "Info: correct assertion of the On sample existence!"
       << endl;

  cout << testName + "Info: asserting the Off sample existence" << endl;
  assert(iactDataSet->GetOffSample() != NULL);
  cout << testName + "Info: correct assertion of the Off sample existence!"
       << endl;

  cout << testName + "Info: asserting the effective area existence" << endl;
  assert(iactDataSet->GetHAeff() != NULL);
  cout << testName + "Info: correct assertion of the effective area existence!"
       << endl;

  cout << testName + "Info: asserting the OFF effective area existence" << endl;
  assert(iactDataSet->GetHAeffOff() != NULL);
  cout << testName + "Info: correct assertion of the effective area existence!"
       << endl;

  cout << testName + "Info: asserting the resolution existence" << endl;
  assert(iactDataSet->GetGEreso() != NULL);
  cout << testName + "Info: correct assertion of the resolution existence!"
       << endl;

  cout << testName + "Info: asserting the bias existence" << endl;
  assert(iactDataSet->GetGEbias() != NULL);
  cout << testName + "Info: correct assertion of the bias existence!" << endl;

  cout << testName + "Info: asserting the migration matrix existence" << endl;
  assert(iactDataSet->GetMigMatrix() != NULL);
  cout << testName +
              "Info: correct assertion of the migration matrix existence!"
       << endl;

  cout << testName + "Info: asserting the background dNdE existence" << endl;
  assert(iactDataSet->GetHdNdEpBkg() != NULL);
  cout << testName + "Info: correct assertion of the background dNdE existence!"
       << endl;

  cout << testName + "Info: asserting the foreground dNdE existence" << endl;
  assert(iactDataSet->GetHdNdEpFrg() != NULL);
  cout << testName + "Info: correct assertion of the foreground dNdE existence!"
       << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// assert if the IactEventListIrf members contain the same number of entries as
// we expect, only for the GEreso and GEbias we specify a number of points,
// as they are TGraphs and do not have entries
//
void assertMembersEntries(IactEventListIrf *iactDataSet,
                          Long_t onEntriesExpected, Long_t offEntriesExpected,
                          Double_t hAeffEntriesExpected,
                          Double_t hAeffOffEntriesExpected,
                          Int_t nGEresoExpected, Int_t nGEbiasExpected,
                          Double_t migMatrixEntriesExpected,
                          Double_t hdNdEpBkgEntriesExpected,
                          Double_t hdNdEpFrgEntriesExpected) {
  TString testName = "TEST::assertMembersEntries ";
  cout << testName + "Info: asserting the individual members entries";
  // entries read from the file
  Long_t onEntries = iactDataSet->GetOnSample()->GetEntries();
  Long_t offEntries = iactDataSet->GetOffSample()->GetEntries();
  Double_t hAeffEntries = iactDataSet->GetHAeff()->GetEntries();
  Double_t hAeffOffEntries = iactDataSet->GetHAeffOff()->GetEntries();
  Int_t nGEreso = iactDataSet->GetGEreso()->GetN();
  Int_t nGEbias = iactDataSet->GetGEbias()->GetN();
  Double_t migMatrixEntries = iactDataSet->GetMigMatrix()->GetEntries();
  Double_t hdNdEpBkgEntries = iactDataSet->GetHdNdEpBkg()->GetEntries();
  Double_t hdNdEpFrgEntries = iactDataSet->GetHdNdEpFrg()->GetEntries();
  // run all the assertions
  cout << testName + "Info: asserting the On Sample entries" << endl;
  assert(checkMemberEntries(onEntries, onEntriesExpected));
  cout << testName + "Info: correct assertion of the On Sample entries!"
       << endl;

  cout << testName + "Info: asserting the Off Sample entries" << endl;
  assert(checkMemberEntries(offEntries, offEntriesExpected));
  cout << testName + "Info: correct assertion of the Off Sample entries!"
       << endl;

  cout << testName + "Info: asserting the effective area entries" << endl;
  assert(checkMemberEntries(hAeffEntries, hAeffEntriesExpected));
  cout << testName + "Info: correct assertion of the effective area entries!"
       << endl;

  cout << testName + "Info: asserting the OFF effective area entries" << endl;
  assert(checkMemberEntries(hAeffOffEntries, hAeffOffEntriesExpected));
  cout << testName +
              "Info: correct assertion of the OFF effective area entries!"
       << endl;

  cout << testName + "Info: asserting the N points of the resolution TGraph"
       << endl;
  assert(checkMemberEntries(nGEreso, nGEresoExpected));
  cout
      << testName +
             "Info: correct assertion of the N points of the resolution TGraph!"
      << endl;

  cout << testName + "Info: asserting the N points of the bias TGraph" << endl;
  assert(checkMemberEntries(nGEbias, nGEbiasExpected));
  cout << testName +
              "Info: correct assertion of the N points of the bias TGraph!"
       << endl;

  cout << testName + "Info: asserting the entries of the migration matrix"
       << endl;
  assert(checkMemberEntries(migMatrixEntries, migMatrixEntriesExpected));
  cout << testName +
              "Info: correct assertion of the entries of the migration matrix!"
       << endl;

  cout << testName + "Info: asserting the entries of the background dNdE"
       << endl;
  assert(checkMemberEntries(hdNdEpBkgEntries, hdNdEpBkgEntriesExpected));
  cout << testName +
              "Info: correct assertion of the entries of the background dNdE!"
       << endl;

  cout << testName + "Info: asserting the entries of the foreground dNdE"
       << endl;
  assert(checkMemberEntries(hdNdEpFrgEntries, hdNdEpFrgEntriesExpected));
  cout << testName +
              "Info: correct assertion of the entries of the foreground dNdE!"
       << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// assert if the members of IactEventListIrf loaded from
// "data/genericIact_dataIRF_01.root" exists and if they contain the same
// number of entries as expected (the latters are obtained opening the file via
// ROOT)
//
void testRootFile() {
  TString testName = "TEST::testRootFile ";
  TString inputFileName = "data/genericIact_dataIRF_01.root";
  cout << testName + "Info: running the test on " + inputFileName << endl;
  TFile *rootFile = new TFile(inputFileName, "READ");
  IactEventListIrf *iactDataSet = new IactEventListIrf();
  iactDataSet = (IactEventListIrf *)rootFile->Get("IactEventListIrf");
  // assert the existence of all the members
  assertMembersExistence(iactDataSet);
  // asssert the expected entries for all the members
  Long_t onEntriesExpected = 11572;
  Long_t offEntriesExpected = 11337;
  Double_t hAeffEntriesExpected = 10000.000;
  Double_t hAeffOffEntriesExpected = 0.0000000;
  Int_t nGEresoExpected = 27;
  Int_t nGEbiasExpected = 27;
  Double_t migMatrixEntriesExpected = 0.0000000;
  Double_t hdNdEpBkgEntriesExpected = 10002.000;
  Double_t hdNdEpFrgEntriesExpected = 0.0000000;
  assertMembersEntries(iactDataSet, onEntriesExpected, offEntriesExpected,
                       hAeffEntriesExpected, hAeffOffEntriesExpected,
                       nGEresoExpected, nGEbiasExpected,
                       migMatrixEntriesExpected, hdNdEpBkgEntriesExpected,
                       hdNdEpFrgEntriesExpected);

  delete rootFile;
  delete iactDataSet;
}
