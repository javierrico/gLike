//////////////////////////////////////////////////////////////////////
// Binned flux likelihood
//////////////////////////////////////////////////////////////////////

#ifndef GLORYDUCKLKL
#define GLORYDUCKLKL

#include "TObjArray.h"
#include "TCanvas.h"

#include "Lkl.h"

class GloryDuckTables2019Lkl : public Lkl
{
 public:
  // constructors
  GloryDuckTables2019Lkl(TString inputString="");
  
  // destructor
  virtual ~GloryDuckTables2019Lkl();
 
  // make the TObjArray own its contents
  void SetOwner(Bool_t set=kTRUE)            {fSampleArray->SetOwner(set);}

  // add samples
  inline void AddSample(Lkl* sample)        {fSampleArray->Add(sample); SetChecked(kFALSE);}

  // getters
  inline const TObjArray* GetSampleArray()     const {return fSampleArray;}
  inline       Int_t      GetNSamples()        const {return fSampleArray->GetEntries();}
  inline       Lkl*       GetReferenceSample() const {return (Lkl*)fSampleArray->At(0);}
  inline       Lkl*       GetSample(Int_t i)   const {return (Lkl*)fSampleArray->At(i);}
  inline       Double_t   GetLogJ()            const {return fLogJ;}
  inline       UInt_t     GetNMasses()         const {return fNMasses;}
  inline       UInt_t     GetActiveMass()      const {return fActiveMass;}

  // unable SetUnitsOfG, dPhi/dE_signal must have proper units
  virtual void SetUnitsOfG(Double_t unit);

  // set mass index to the current active one
  void SetActiveMass(Double_t mass);

  // print data in the overview
  virtual void PrintData(Int_t level=0);

  // plots
  TCanvas* PlotInputData();

  // data input
  Int_t ReadGloryDuckInputData(TString inputfilename);

 protected:
          Int_t InterpretInputString(TString inputString);
  virtual void  SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t MakeChecks();
  virtual void  SetMinuitLink();
  inline  void  ClearSampleList()  {fSampleArray->Clear();}

  // data input
  Int_t CreateAndAddNewParabola(Double_t mass,Double_t svmin,Double_t svmax,Int_t npoints,Double_t* sv,Double_t* logL);

 private:
  UInt_t      fNMasses;       //   Number of masses
  UInt_t      fActiveMass;    //   Position in the array fMass of the current active mass
  TObjArray*  fSampleArray;   //-> Array of Lkl objects (one for each mass) 
  Double_t*   fsvMin;         //-> [fNMasses] array with <sv> lower values
  Double_t*   fsvMax;         //-> [fNMasses] array with <sv> upper values
  Double_t*   fMass;          //-> [GeV] array with dark matter particle mass
  Double_t    fLogJ;          //   log10 [GeV^2/cm^5] of estimated J-factor

  ClassDef(GloryDuckTables2019Lkl,1) // Likelihood vs <sv>
};

#endif
