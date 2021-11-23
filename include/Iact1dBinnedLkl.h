//////////////////////////////////////////////////////////////////////
// Binned full likelihood
//////////////////////////////////////////////////////////////////////

#ifndef BINNEDFULLLKL
#define BINNEDFULLLKL

#include "Iact1dUnbinnedLkl.h"
#include "JointLkl.h"

class Iact1dBinnedLkl : public Iact1dUnbinnedLkl, public JointLkl
{
 public:
  // enumerations
  enum parIndex_t {gGIndex,gTauIndex};                    // Indeces of parameters

  // constants
  static const UInt_t  gDefNBins = 10;          // default number of bins
  static const UInt_t  gDefMinBinContent = 1;   // default minimum number of entries per bin (if != 0 allow rebinning)
  
  // constructors
  Iact1dBinnedLkl(TString inputFileName);
  
  // destructor
  virtual ~Iact1dBinnedLkl();

  // configure
  inline void SetMinBinContent(UInt_t minbincontent) {fMinBinContent=minbincontent; SetChecked(kFALSE);}

  inline void SetKnownBackground(Bool_t fix=kTRUE)   {fKnownBackground = fix;}
  
  inline UInt_t GetNBins()         const {return fNBins;}
  inline UInt_t GetMinBinContent() const {return fMinBinContent;}
  
  inline TH1I*  GetHNOn()          const {return fHNOn;}
  inline TH1I*  GetHNOff()         const {return fHNOff;}
  inline Bool_t GetTauEDepFluct()  const {return fTauEDepFluct;}

  virtual TH1F*   GetHdNdEpOn(Bool_t isDifferential=kTRUE,Int_t nbins=0)  const;
  virtual TH1F*   GetHdNdEpOff(Bool_t isDifferential=kTRUE,Int_t nbins=0) const;

  // print data in the overview
  virtual void PrintOverview(Int_t level=0)  {Lkl::PrintOverview(level);}
  virtual void PrintData(Int_t level=0);
  virtual Int_t  SimulateDataSamples(Float_t meanGwithUnits=0,TRandom* rdm=NULL);

 protected:
          Int_t  InterpretInputString(TString inputString);
  virtual Int_t  GetRealBkgAndGoffHistos(TRandom* rdm,TH1F*& hdNdEpBkg,TH1F*& hdNdEpSignalOff);
  
  virtual void  SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t MakeChecks();
  virtual void  SetMinuitLink();

  // prepare material
  Int_t          BuildAndBinOnOffHistos();
  Int_t          ConfigureJointLkl();
  virtual Int_t  PrepareForLklScan(Bool_t centerAtZero,Int_t npoints,Double_t glow,Double_t gupp,Bool_t isVerbose);
  virtual void   SpreadFixLklVsG(Double_t g);

 private:  
  UInt_t fNBins;            // Number of bins
  UInt_t fMinBinContent;    // Minimum allowed number of bins in On and Off samples (if !=0 allow rebinning)
  UInt_t fNRemovedBins;     // Number of removed bins by rebinning
  Bool_t fKnownBackground;  // if kTRUE, the parameter b of the double poisson is considered fixed (for tests mostly)
  Bool_t fTauEDepFluct;     // if kTRUE (default kFALSE) tau is allowed to fluctuate (according to fTau and fDTau independently in each bin)
  Bool_t fdNdEpBkgFromOff;  // if kTRUE (default kFALSE), the background distribution (fHdNdEpBkg, for simulations) will be constructed from Off data distribution (even if fHdNdepBkg provided in input file)
    
  
  TH1I* fHNOn;           //-> histogram for On  events vs E'
  TH1I* fHNOff;          //-> histogram for Off events vs E'

  ClassDef(Iact1dBinnedLkl,1) // Binned full likelihood
};

#endif
