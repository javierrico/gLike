//////////////////////////////////////////////////////////////////////
// multi-sample joint lilelihood
//////////////////////////////////////////////////////////////////////

#ifndef JOINTLKL
#define JOINTLKL

#include "TObjArray.h"

#include "Lkl.h"

class JointLkl : public virtual Lkl
{
 public:
  // constructor
  JointLkl(TString inputString="");
  
  // destructor
  virtual ~JointLkl();

  // make the TObjArray own its contents
  void SetOwner(Bool_t set=kTRUE)            {fSampleArray->SetOwner(set);}
  
  // add samples
  inline void AddSample(Lkl* sample)        {fSampleArray->Add(sample); SetChecked(kFALSE);}

  // getters
  inline TObjArray* GetSampleArray()     const {return fSampleArray;}
  inline Int_t      GetNSamples()        const {return fSampleArray->GetEntries();}
  inline Lkl*       GetReferenceSample() const {return (Lkl*)fSampleArray->At(0);}
  inline Lkl*       GetSample(Int_t i)   const {return (Lkl*)fSampleArray->At(i);}

  // print results
  virtual void      PrintOverview(Int_t level=0);
  virtual void      PrintData(Int_t level=0);

  virtual void      ResetGLklVsG();
  virtual void      SetMinuitLink();
  
  // virtual Double_t  ComputeLklVsG(Bool_t centerAtZero=kFALSE,Int_t npoints=200,Double_t glow=0,Double_t gupp=0,Bool_t isVerbose=kTRUE);
  //reimplemented on 27th July byt Tomo
  virtual      Double_t   ComputeLklVsG(Bool_t centerAtZero=kFALSE,Int_t npoints=2000,Double_t glow=0,Double_t gupp=0,Bool_t isVerbose=kTRUE);
 
 protected:
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  inline  void     ClearSampleList()  {fSampleArray->Clear();}
  virtual Int_t    ReorderSamples();
  virtual Int_t    PrepareForLklScan(Bool_t centerAtZero,Int_t npoints,Double_t glow,Double_t gupp,Bool_t isVerbose);
  virtual Double_t GetExpansionCoefficient() const;
  virtual void     SpreadFixLklVsG(Double_t g);

 private:  
  TObjArray*  fSampleArray;  //-> Array of Lkl objects (one for each sample)

  ClassDef(JointLkl,1) // Joint likelihood for multiple samples
};

#endif
