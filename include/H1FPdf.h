//////////////////////////////////////////////////////////////////////
// Histogram containing dNdE information
//////////////////////////////////////////////////////////////////////

#ifndef HDNDE
#define HDNDE

class H1FPdf : public TH1F
{
 public:
  
  // constructors
  H1FPdf(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup);
  H1FPdf(const char *name, const char *title, Int_t nbinsx, const Double_t *xbins);
  H1FPdf(const char *name, const char *title, Int_t nbinsx, const Float_t *xbins);
  H1FPdf(const TH1F &h1f);
  
  // destructor
  virtual ~H1FPdf();

 protected:
  
 private:

  ClassDef(H1FPdf,1) // Full Likelihood (unbinned)
};

#endif
