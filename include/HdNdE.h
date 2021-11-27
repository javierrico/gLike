//////////////////////////////////////////////////////////////////////
// Histogram containing dNdE information
//////////////////////////////////////////////////////////////////////

#ifndef HDNDE
#define HDNDE

class HdNdE : public TH1F
{
 public:
  
  // constructors
  HdNdE(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup);
  HdNdE(const char *name, const char *title, Int_t nbinsx, const Double_t *xbins);
  HdNdE(const char *name, const char *title, Int_t nbinsx, const Float_t *xbins);
  HdNdE(const TH1F &h1f);
  
  // destructor
  virtual ~HdNdE();

 protected:
  
 private:

  ClassDef(HdNdE,1) // Full Likelihood (unbinned)
};

#endif
