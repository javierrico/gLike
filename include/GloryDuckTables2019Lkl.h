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
 
  // add samples
  inline void AddSample(Lkl* sample)        {fSampleArray->Add(sample); SetChecked(kFALSE);}

  // getters
  inline const TObjArray* GetSampleArray()     const {return fSampleArray;}
  inline       Int_t      GetNSamples()        const {return fSampleArray->GetEntries();}
  inline       Lkl*       GetSample(Int_t i)   const {return (Lkl*)fSampleArray->At(i);}
  inline       UInt_t     GetNMasses()         const {return fNMasses;}
  inline       UInt_t     GetActiveMass()      const {return fActiveMass;}

  // unable SetUnitsOfG, dPhi/dE_signal must have proper units
  virtual void SetUnitsOfG(Double_t unit);

  // set mass index to the current active one
  Int_t SetActiveMass(Double_t mass); // using the mass
  Int_t SetActiveMass(Int_t index);   // using the index

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
  Int_t CreateAndAddNewParabola(Double_t mass,Int_t npoints,Double_t* sv,Double_t* logL);

 private:
  UInt_t      fNMasses;       //   Number of masses
  UInt_t      fActiveMass;    //   Position in the array fMass of the current active mass
  Double_t*   fMass;          //-> [GeV] array with dark matter particle mass
  TObjArray*  fSampleArray;   //-> Array of Lkl objects (one for each mass) 
  UInt_t      fNsvVals;       //   Number of <sv> values
  Double_t*   fsvVals;        //-> Array with <sv> values

  ClassDef(GloryDuckTables2019Lkl,1) // Likelihood vs <sv>
};

#endif
