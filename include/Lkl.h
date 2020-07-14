/////////////////////////////////////////////////////////////////////////////
// abstract class for (profile) likelihood maximization (-2logL minimization)
/////////////////////////////////////////////////////////////////////////////

#ifndef LKL
#define LKL

#include "TNamed.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TH1.h"
#include "TString.h"
#include <string>
#include <iostream>

class Lkl : public TNamed
{
 public:
  // determines how to interpret fDUnitsOfG
  enum DUofGType_t {none=0,  // fUnitsOfG has no error, similar to fDUnitsOfG=0 [default]
		    lin,     // fDUnitsOfG interpreted as the relative error of fUnitsOfG
		    invlin,  // fDUnistOfG interpreted as the relative error of 1/fUnitsOfG
		    log,     // fDUnitsOfG interpreted as the error of log10(fUnitsOfG)
		    invlog,  // fDUnitsOfG interpreted as the error of 1/log10(fUnitsOfG)
		    tottypes};

  static const Int_t gGParIndex   = 0;      // index of g parameter in fMinuit
  static const Int_t gLevelMargin = 39;     // number of blanck spaces for PrintOverview

  // constructor
  Lkl(Int_t npars,TString inputString="",TString name="",TString title="");
  
  // destructor
  virtual ~Lkl();

  // call minimization of the  Likelihood function
  virtual Double_t ComputeLklVsG(Bool_t centerAtZero=kFALSE,Int_t npoints=200,Double_t glow=0,Double_t gupp=0,Bool_t isVerbose=kTRUE);

  static void PrintGLikeBanner();
  
  // Retrieve results from fit
  virtual      void      PrintOverview(Int_t level=0);
  virtual      void      PrintData(Int_t level=0)
  {
    Margin(level); std::cout << "              Object Name : " << GetName() << std::endl;
  }
  inline  TMinuit*  GetMinuit()                       const {return fMinuit;}
  TGraph*           GetLklVsG(Bool_t units=kTRUE)     const; 
  inline  Double_t  GetLklMin()                       const {return fLklMin;}
  inline  Double_t  GetGLklMin(Bool_t units=kTRUE)    const 
  {
    return (units? fParVal[gGParIndex]*fUnitsOfG: fParVal[gGParIndex]);
  }
  inline       Double_t  GetGLklMinErr(Bool_t units=kTRUE) const 
  {
    return (units? fParErr[gGParIndex]*fUnitsOfG:fParErr[gGParIndex]);
  }
  inline Double_t GetParVal(Int_t ipar)  const {return (ipar<fNPars?  fParVal[ipar] : 0);}
  inline Double_t GetParErr(Int_t ipar)  const {return (ipar<fNPars?  fParErr[ipar] : 0);}
  inline Bool_t   IsParFixed(Int_t ipar) const {return (ipar<fNPars?  fIsParFixed[ipar] : 0);}

  Double_t            GetGForLkl(Double_t lkl,Bool_t units=kTRUE) const;
  inline  Double_t    GetErrorDef()     const {return fErrorDef;}
  inline  Double_t*   GetParDelta()     const {return fParDelta;}
  inline  TString*    GetParName()      const {return fParName;}
  inline  Double_t*   GetParStart()     const {return fParStart;}
  virtual Double_t    GetUnitsOfG()     const {return fUnitsOfG;}
  virtual Double_t    GetDUnitsOfG()    const {return fDUnitsOfG;}
  virtual DUofGType_t GetDUofGType()    const {return fDUofGType;}
  inline  Bool_t      GetGIsPositive()  const {return fGIsPositive;}
  inline  Int_t       GetNFreePars()    const {Int_t n=0; for(Int_t i=0;i<fNPars;i++) n+=!IsParFixed(i); return n;}

  // set/get a global normalization factor (units of g if you will), compulsory if 
  // object is going to be added to an JointLkl instance
  inline void  SetErrorDef(Double_t errorDef)                 {fErrorDef = errorDef;}
  virtual void SetUnitsOfG(Double_t unit)                     {fUnitsOfG=unit;   SetChecked(kFALSE);}
  virtual void ResetUnitsOfG()                                {fUnitsOfG=1;      SetChecked(kFALSE);}
  virtual void SetDUnitsOfG(Double_t dUofG,DUofGType_t dtype) {fDUnitsOfG=dUofG; fDUofGType=dtype;}
  inline  void SetGIsPositive(Bool_t p=kTRUE)                 {fGIsPositive=p;}
  virtual void ResetGLklVsG()                                 {if(fGLklVsG){delete fGLklVsG; fGLklVsG=NULL;fGShift=0;}}

  void             InitMinuit(Double_t ginit=0);
  virtual Double_t MinimizeLkl(Double_t g,Bool_t isFixed=kFALSE,Bool_t isVerbose=kTRUE,Bool_t force=kFALSE);
  virtual Double_t MinimizeLkl() 
  { 
    if(!fParStart) return MinimizeLkl(0,kFALSE,kTRUE,kTRUE); 
    return MinimizeLkl(fParStart[gGParIndex],kFALSE,kTRUE,kTRUE);
  }
  virtual Int_t    ApplyDUnitsOfG(Int_t npoints,Double_t glow,Double_t gupp);
  inline Bool_t    IsChecked() const            {return fIsChecked;}
  virtual Int_t    MakeChecks()         = 0;  // NEEDS TO BE OVERRIDEN BY DAUGHTERS
  virtual void     SetMinuitLink()      = 0;  // NEEDS TO BE OVERRIDEN BY DAUGHTERS

 protected:
  Int_t            InterpretInputString(TString inputString);
  
  virtual void     SetFunctionAndPars(Double_t ginit=0) = 0;  // NEEDS TO BE OVERRIDEN BY DAUGHTERS

  void     SetParameters(const Char_t** parname, Double_t* pstart, Double_t* pdelta);
  Double_t CallMinimization(Double_t g=0,Bool_t isVerbose=kFALSE,Int_t strategy=2);
  
  inline void   FixPar(Int_t ipar,Bool_t fix=kTRUE) {if(ipar<fNPars) fIsParFixed[ipar] = fix;}
  inline void   SetChecked(Bool_t status=kTRUE)     {fIsChecked=status;}
  inline void   SetGLklVsG(TGraph* graph)           {if(fGLklVsG) delete fGLklVsG; fGLklVsG=graph;}

  Double_t      CenterAtZero();
  void          FixLklVsG(Double_t shiftg=0);
  virtual void  SpreadFixLklVsG(Double_t g) {g*=1;}

  Double_t         GetLklVal() const;
  inline  void     Margin(Int_t level) const {std::cout << " " ;for(Int_t i=0;i<level*gLevelMargin;i++) std::cout << " "; std::cout << "*";}
  void             FindGLowAndGUpp(Double_t& glow,Double_t& gupp,Bool_t center=kFALSE);
  void             ExpandGLowAndGUpp(Double_t& glow,Double_t& gupp,Double_t stretch=1);

  virtual Int_t    PrepareForLklScan(Bool_t centerAtZero,Int_t npoints,Double_t glow,Double_t gupp,Bool_t isVerbose) {Double_t dummy = centerAtZero*npoints*glow*gupp*isVerbose;return dummy*0;}
  
  virtual Double_t GetExpansionCoefficient() const {return 1;}

  // tools
  virtual Double_t IntegrateLogE(const TH1F* h1,Double_t lmin,Double_t lmax) const;
  virtual void     SetGLklVsG(Int_t npoints,Double_t* x,Double_t* y);

  // data members
  TMinuit*    fMinuit;       //-> TMinuit for minization
 private:  
  Int_t       fNPars;        //   Number of free+nuisance parameters 
  TString*    fParName;      //-> Names of free+nuisance parameters
  Double_t*   fParStart;     //-> Initial values of free+nuisance parameters
  Double_t*   fParDelta;     //-> Initial precision of free+nuisance parameters
  Double_t*   fParVal;       //-> Array of par values after fit
  Double_t*   fParErr;       //-> Array of par errors after fit
  Bool_t*     fIsParFixed;   //-> Array of fixed par flags
  Double_t    fErrorDef;     //   Delta(Lkl) defining the error returned by TMinuit (default = 1; 1-sided 95%CL = 2.70; 1-sided 5sigma CL = 4.495)
  TGraph*     fGLklVsG;      //-> Profile likelihood parabola (-2logL vs g)
  Double_t    fLklMin;       //   Minimum value of likelihood profile
  Double_t    fUnitsOfG;     //   the value you need to multiply g by to get it in your favorite units
  Double_t    fDUnitsOfG;    //   relative uncertainty in Units of G, 1/G, uncertainy on log10(G) or 1/log10(G) (see fDUofGType)
  DUofGType_t fDUofGType;    //   meaning of fDUnitsOfG
  Double_t    fDUofGEst;     //   estimated real value of real fDUnitsOfG (from profile)
  Bool_t      fGIsPositive;  //   restrict g to positive values only
  Double_t    fGShift;       //   applied shift (when centerning the -2logL curve so that min is at G=0)

  Bool_t      fIsChecked;    //   flag to be set once preliminary checks have been performed and we're ready to call minuit
  Int_t       fStatus;       //   status of fMinuit after last call (0=convergence)

  ClassDef(Lkl,1) // Abstract class for (profile) likelihood maximization (-2logL minimization)
};

#endif
