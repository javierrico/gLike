// to make the assert() function and crash the execution you have to run
// .x test/main.cpp++g
#include "TestSuiteIactEventListIrf.cc"
#include "TestSuiteJointLkl.cc"
#include "TestSuitePoissonLkl.cc"

int main() {
  testPoissonLkl();
  testJointPoissonLkl();
  testRootFile();
}
