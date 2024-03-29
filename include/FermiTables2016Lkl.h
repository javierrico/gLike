//////////////////////////////////////////////////////////////////////
// Likelihood vs <sv> in E-bins from Fermi from their publication
// Ackermann, M. et al. Phys. Rev. Lett., 115 (2015) 231301
// and available at http://www-glast.stanford.edu/pub_data/1048/
//////////////////////////////////////////////////////////////////////

#ifndef BINNEDFLUXLKL
#define BINNEDFLUXLKL

#include "JointLkl.h"
#include "TCanvas.h"
#include "TMath.h"
#include "HdNdE.h"

class FermiTables2016Lkl :  public JointLkl, public HdNdE
{
 public:
  // constructors
  FermiTables2016Lkl(TString inputString="");
  
  // destructor
  virtual ~FermiTables2016Lkl();
  
  // getters
  inline UInt_t    GetNBins()         const {return fNBins;}
  inline Double_t* GetEFluxInt()      const {return fEFluxInt;}
  inline Double_t  GetLogJ()          const {return fLogJ;}

  // unable SetUnitsOfG, dPhi/dE_signal must have proper units
  virtual void SetUnitsOfG(Double_t unit);

  // set histogram for signal differential flux
  virtual Int_t ResetdNdESignal();
  Int_t ReaddNdESignalFromFermi(TString filename);

  // specify the mass and log10 of Jfactor and its 1-sigma error (nuisance parameter)
  virtual void SetDMMass(Double_t mass)
  {
    fMass  = mass; 
    SetChecked(kFALSE);
  }

  // print data in the overview
  virtual void PrintOverview(Int_t level=0) {Lkl::PrintOverview(level);}
  virtual void PrintData(Int_t level=0);
  virtual void ResetGLklVsG() {Lkl::ResetGLklVsG();}

  // Plots
  TCanvas* PlotInputData();
  
 protected:
  // include check of the dNdE histogram
  virtual Bool_t   IsChecked() const                   {return (Lkl::IsChecked() & HdNdE::IsHdNdESignalChecked());}
  virtual void     SetChecked(Bool_t status=kTRUE)     {Lkl::SetChecked(status); HdNdE::SetHdNdESignalChecked(status);}

  
          Int_t InterpretInputString(TString inputString);
  virtual void  SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t MakeChecks();
  virtual void  SetMinuitLink();

  // data input
  Int_t ReadFermiInputData(TString inputfilename);
  Int_t CreateAndAddNewParabola(Double_t emin,Double_t emax,Int_t npoints,Double_t* flux,Double_t* logL);
  
  // make data ready for minimization
  Int_t ComputeEFluxIntegrals();
  virtual Int_t PrepareForLklScan(Bool_t centerAtZero,Int_t npoints,Double_t glow,Double_t gupp,Bool_t isVerbose) {centerAtZero*=1;npoints*=1;glow*=1; gupp*=1; isVerbose*=1;return 0;}
  virtual void SpreadFixLklVsG(Double_t g);

 private:  
  UInt_t    fNBins;         //   Number of bins
  Double_t* fBinEMin;       //-> [fNBins] array with E bin lower bound
  Double_t* fBinEMax;       //-> [fNBins] array with E bin upper bound
  Double_t* fEFluxInt;      //-> [fNBins] array with E flux integrals in bins

  Double_t  fMass;          //   [GeV] dark matter particle mass
  Double_t  fLogJ;          //   log10 [GeV^2/cm^5] of estimated J-factor

  ClassDef(FermiTables2016Lkl,1) // Likelihood vs <sv> in E-bins from Fermi 2015 paper
};

#endif
