//////////////////////////////////////////////////////////////////////
// On vs Off Poisson likelihood
//////////////////////////////////////////////////////////////////////

#ifndef POISSONLKL
#define POISSONLKL

#include "Lkl.h"

class PoissonLkl : public virtual Lkl
{
 public:
  // enumerations
  enum         parIndex_t {gGIndex,gBIndex,gTauIndex,gEffIndex};        // Indeces of parametersstatic
  
  // constructors
  PoissonLkl(UInt_t non,UInt_t noff,Float_t tau=1,Float_t dTau=0,TString name="",TString title="");
  PoissonLkl(TString fileName,TString name="",TString title="");
  
  // destructor
  virtual ~PoissonLkl();

  // configure minuit
  inline  void     SetGFractionInOff(Float_t f)       {fGFractionInOff  = f;}
  inline  void     SetFrgNEvents(Float_t frgn)        {fFrgNEvents      = frgn;}
  inline  void     SetTau(Float_t tau)                {fTau             = tau;}
  inline  void     SetDEff(Float_t deff)              {fDEff            = deff;}
  inline  void     SetKnownBackground(Bool_t k=kTRUE) {fKnownBackground = k;}
  
  // getters
  inline  UInt_t   GetNon()            const    {return fNon;}
  inline  UInt_t   GetNoff()           const    {return fNoff;}
  inline  Float_t  GetTau()            const    {return fTau;}
  inline  Float_t  GetDTau()           const    {return fDTau;}
  inline  Float_t  GetDEff()           const    {return fDEff;}
  inline  Float_t  GetGFractionInOff() const    {return fGFractionInOff;}
  inline  Float_t  GetFrgNEvents()     const    {return fFrgNEvents;}
  
 protected:
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();
  
 private:  
  UInt_t   fNon;          // Number of measured events in On region
  UInt_t   fNoff;         // Number of measured events in Off region
  Float_t  fTau;          // normalization Noff/Non (e.g # of off regions)
  Float_t  fDTau;         // statistical error in fTau
  Float_t  fDEff;         // uncertainty (Gaussian) in gamma-ray efficiency (following Rolke et al. 2005)

  // Flags internally used by the class
  Float_t  fGFractionInOff;  // what fraction of G is expected in the Off region (0 by default but can be finite for extended sources)
  Float_t  fFrgNEvents;      // number of foreground events in the On region (0 by default)
  Bool_t   fKnownBackground; // if kTRUE, the parameter b is considered fixed (for tests mostly)
  
  ClassDef(PoissonLkl,1)   // On vs Off Poisson likelihood
};

#endif
