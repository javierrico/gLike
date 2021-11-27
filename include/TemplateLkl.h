//////////////////////////////////////////////////////////////////////
// Template for new Lkl class
//////////////////////////////////////////////////////////////////////

#ifndef TEMPLATELKL
#define TEMPLATELKL

#include "Lkl.h"

class TemplateLkl : public Lkl
{
 public:
  
  // constructors
  TemplateLkl(TString inputString="");
  
  // destructor
  virtual ~TemplateLkl();

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();

  
 private:  

  ClassDef(TemplateLkl,1) // Template for new Lkl class
};

#endif
