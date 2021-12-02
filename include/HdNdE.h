//////////////////////////////////////////////////////////////////////
// Abstract class for Lkl classes using a dNdE Histogram
//////////////////////////////////////////////////////////////////////

#ifndef HDNDE
#define HDNDE

class TH1F;

class HdNdE
{
 public:
  
  // constructor
  HdNdE(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup);
  
  // destructor
  virtual ~HdNdE();

  // public functions
  inline  TH1F*    GetHdNdESignal()  const {return fHdNdESignal;}

  virtual Int_t    ResetdNdESignal();
  virtual Int_t    SetdNdESignal(TH1F* hdNdESignal);
  virtual Int_t    AdddNdESignal(TString filename,Float_t br);
  virtual Int_t    SetdNdESignalFunction(TString function,Float_t br=1.,Float_t p0=0,Float_t p1=0,Float_t p2=0);
  virtual Int_t    AdddNdESignalFunction(TString function,Float_t br=1.,Float_t p0=0,Float_t p1=0,Float_t p2=0);
  virtual Int_t    SetdNdESignalFunction(TF1* function,Float_t br=1.,Float_t emin=0,Float_t emax=1e9);
  virtual Int_t    AdddNdESignalFunction(TF1* function,Float_t br=1.,Float_t emin=0,Float_t emax=1e9);
  virtual Int_t    ReaddNdESignal(TString filename);
    
 protected:
  virtual Bool_t   IsHdNdESignalChecked() const {return fHdNdESignalChecked;}
  virtual void     SetHdNdESignalChecked(Bool_t status=kTRUE) {fHdNdESignalChecked=status;}
  static Int_t     ReadAndInterpolate(TH1F* ih,TH1F* oh,Double_t scale=0,Bool_t isDiff=kTRUE);
  Int_t            SetHistogramPars();
    
  // leave data members as protected for easy access by parent classes
  TH1F*    fHdNdESignal;        //-> dN/dE vs E histogram for signal events
  TString  fName;               //-> name of the histogram
  TString  fTitle;              //-> name of the histogram
  Int_t    fNFineBins;          // number of fine bins for the histogram
  Double_t fFineLEMin;          // minimum log(energy[GeV]) for the histogram
  Double_t fFineLEMax;          // maximum log(energy[GeV]) for the histogram
  Bool_t   fHdNdESignalChecked; // flag it up when applying any changes to the histogram
  
  ClassDef(HdNdE,1) // Abstract class for Lkl classes using a dNdE Histogram
};

#endif
