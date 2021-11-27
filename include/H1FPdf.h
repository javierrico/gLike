//////////////////////////////////////////////////////////////////////
// Histogram containing dNdE information
//////////////////////////////////////////////////////////////////////

#ifndef HDNDE
#define HDNDE

class H1FPdf : public TH1F
{
 public:
  
  // constructors
  H1FPdf(TString inputString="");
  
  // destructor
  virtual ~H1FPdf();

 protected:
  
 private:

  ClassDef(H1FPdf,1) // Full Likelihood (unbinned)
};

#endif
