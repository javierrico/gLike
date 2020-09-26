/* ======================================================================== *\
!
!   Author(s): Daniel Kerszberg  08/2019 <mailto:dkerszberg@ifae.es>
!              Tomohiro Inada    08/2019 <mailto:tomohiro@icrr.u-tokyo.ac.jp>
! 
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// LineSearchLkl
//
// Add documentation!
//
//////////////////////////////////////////////////////////////////////////////

// include c++ needed classes

// include Root needed classes
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"

// include gLike needed classes
#include "LineSearchLkl.h"

ClassImp(LineSearchLkl);

using namespace std;

// class name and title
static const TString  gName            = "LineSearchLkl";
static const TString  gTitle           = "Line Search Likelihood";
static const Double_t gCenterBin       = 0.5;                    // decide which value represents bin in histogram (= 0 for lower bin edge, 0.5 for the middle, 1 for the right edge
static const Int_t    gNBins           = 100;                    // default number of histograms for dN/dE plots

// List of free parameters.
static const Int_t    gNPars           = 3;                      // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"g","b","tau"};        // Name of parameters

// static functions (for internal processing of input data)
static TH1F* GetResidualsHisto(TH1F* hModel,TH1F* hData);
void differentiate(TH1F* h,Int_t error_flag);

// -2logL function for minuit
void lineSearchLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;


//////////////////////////////////////////////////////////////////////////////
//
// String constructor
// The string contains the elements for the constructor
// which you can read from e.g. an rc input file
//
LineSearchLkl::LineSearchLkl(TString inputString) :
  Lkl(gNPars,inputString,gName,gTitle), Iact1dUnbinnedLkl(inputString),
 fEventsInEnergyWindow(0)
{
  if(InterpretInputString(inputString))
    cout << "LineSearchLkl::LineSearchLkl Warning: there were problems interpreting the input string" << endl;      
}

//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
// The string has been passed to Lkl::InterpretInputString in constructor
// normally read from an input rc file
//
// Return 0 in case checks are Ok, 1 otherwise
//
Int_t LineSearchLkl::InterpretInputString(TString inputString)
{
  // decode your inputString and configure the class accordingly

  return 0;
}

////////////////////////////////////////////////////////////////
// 
// Destructor
// make sure you do not leave any memory leaks
//
LineSearchLkl::~LineSearchLkl()
{
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of TMinuit subtleties)
// Leave this unchanged
//
void LineSearchLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance)
// in case of need, some of the nuisance parameters can be fixed here
//
void LineSearchLkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(lineSearchLkl);  
  fMinuit->SetName(Form("%s_Minuit",GetName()));

  // set and pars initial and step values
  Double_t pStart[gNPars] = {ginit, fEventsInEnergyWindow,GetTau()};
  Double_t pDelta[gNPars] = {TMath::Sqrt(fEventsInEnergyWindow)/10.,TMath::Sqrt(fEventsInEnergyWindow)/10.,GetDTau()/10.};    // Precision of parameters during minimization

  // initialize the free (and nuisance) parameters
  SetParameters(gParName,pStart,pDelta);

  // Fix tau if requested (both in minuit and in lkl)
  fMinuit->Release(gTauIndex);
  FixPar(gTauIndex,kFALSE);
  if(GetDTau()<=0)
    {
      fMinuit->FixParameter(gTauIndex);
      FixPar(gTauIndex);
    }
}

////////////////////////////////////////////////////////////////
//
// Check that all elements needed for the fit are present, otherwise
// try to create/compute them
//
// Return 0 in case checks are Ok, 1 otherwise
//
Int_t LineSearchLkl::MakeChecks()
{

  if(IsChecked()) return 0;

  // compute background model
  ComputeBkgModelFromOnHisto();

  // check if all needed histos are there
  if(CheckHistograms())
    {
      cout << "LineSearchLkl::MakeChecks Warning: missing information, cannot perform fit, check your code!" << endl;
      return 1;
    }

  SetChecked();
  return 0;
}	      

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute background model from On data
//
Int_t LineSearchLkl::ComputeBkgModelFromOnHisto()
{
  const Float_t* onSample = GetOnSample();
  UInt_t              Non = GetNon();

  TH1F* hdNdEpOn = new TH1F("hdNdEpOn","dN/dE' for On data",GetNFineBins(),GetFineLEMin(),GetFineLEMax());

  // filling and counting number of events inside energy wondow
  int count=0;
  for(ULong_t ievent=0; ievent<Non; ievent++)
    {
      if(onSample[ievent] > TMath::Log10(GetEmin()) && onSample[ievent] < TMath::Log10(GetEmax())){
        hdNdEpOn->Fill(onSample[ievent]);
        count++;
      }
    }
  fEventsInEnergyWindow=count;

  differentiate(hdNdEpOn,1);

  TH1F* hdNdEpBkg = new TH1F("HdNdEpBkg","dN/dE' for background model",GetNFineBins(),GetFineLEMin(),GetFineLEMax());

  TF1* f1 = new TF1("f1","expo",GetFineLEMin(),GetFineLEMax());

  hdNdEpOn->Fit("f1","0","",TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));

  for(Int_t ibin=1; ibin < GetNFineBins()-1; ibin++)
    if (hdNdEpOn->GetBinContent(ibin)>0)
      hdNdEpBkg->SetBinContent(ibin,f1->Eval(hdNdEpBkg->GetBinCenter(ibin)));

  // replace fHdNdEpBkg with newly computed background histogram
  SetdNdEpBkg(hdNdEpBkg);

  delete hdNdEpOn;
  delete f1;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Plot a canvas with all material to minimize the -2logL:
// - Aeff
// - Ereso and Ebias
// - dN/dE  for signal
// - dN/dE' for background compared to the On and Off distributions with residuals
// - dN/dE' for signal
// 
// return a pointer to the canvas, in case you want to modify it
//
TCanvas* LineSearchLkl::PlotHistosAndData()
{
  MakeChecks();
  SetChecked(kFALSE);

  // create and divide canvas
  TCanvas* canvas = new TCanvas("histosAndDataCanvas","LineSearchLkl histos and data used to minimize -2logL", 1000, 1500);
  canvas->Divide(2,3);

  // draw plots

  // effective area
  canvas->cd(1);
  if(GetHAeff()) GetHAeff()->DrawCopy();
  gPad->SetLogy();
  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();

  // dN/dE' background vs data
  canvas->cd(2);
  TH1F* hOn  = GetHdNdEpOn();
  hOn->SetDirectory(0);
  TH1F* hdNdEpBkg = GetHdNdEpBkg();

  TH1I *dummya = new TH1I("dummya", "dN/dE' bkg model vs On distribution",1,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
  dummya->SetStats(0);
  if(GetNon()>1)
    {
      dummya->SetMinimum(hOn->GetMinimum(0)/2.);
      dummya->SetMaximum(hOn->GetMaximum()*2);
    }
  else if(hdNdEpBkg)
    {
      dummya->SetMinimum(hdNdEpBkg->GetMinimum());
      dummya->SetMaximum(hdNdEpBkg->GetMaximum());
    }

  if(!hdNdEpBkg) dummya->SetTitle("dN/dE' distributions for On event samples");
  dummya->SetXTitle("log_{10}(E' [GeV])");
  dummya->SetYTitle("dN/dE' [GeV^{-1}]");
  dummya->DrawCopy();

  hdNdEpBkg->SetDirectory(0);
  if(hdNdEpBkg)
    {
      hdNdEpBkg->SetLineColor(2);
      hdNdEpBkg->SetMarkerColor(2);
      hdNdEpBkg->SetMarkerStyle(2);
      hdNdEpBkg->DrawCopy("esame");
    }
  hOn->DrawCopy("esame");

  gPad->SetLogy();
  gPad->SetGrid();

  TLegend* hleg = new TLegend(0.65, 0.65, 0.85, 0.85);
  hleg->SetFillColor(0);
  hleg->SetMargin(0.40);
  hleg->SetBorderSize(0);
  hleg->AddEntry(hOn,"On events","LP");
  hleg->AddEntry(hdNdEpBkg,"Bkg model","LP");
  hleg->Draw();

  gPad->Modified();
  gPad->Update();

  // energy resolution and bias
  canvas->cd(3);
    if(GetMigMatrix()){
      GetMigMatrix()->DrawCopy("colz");
    }
    else{
  TH1I *dummye;
  if(GetGEreso())
    dummye = new TH1I("dummye", "Energy resolution and bias",1,GetGEreso()->GetX()[0],GetGEreso()->GetX()[GetGEreso()->GetN()-1]);
  else
    dummye = new TH1I("dummye", "Energy resolution and bias",1,2,4);
  dummye->SetStats(0);
  dummye->SetMinimum(-0.1);
  dummye->SetMaximum(0.4);
  dummye->SetXTitle("log_{10}(E [GeV])");
  dummye->SetYTitle("Energy resolution and bias [%]");
  dummye->DrawCopy();
  if(GetGEreso()) GetGEreso()->Draw("l");
  if(GetGEbias()) GetGEbias()->Draw("l");

  TLegend* hleg2 = new TLegend(0.65, 0.69, 0.90, 0.89);
  hleg2->SetFillColor(0);
  hleg2->SetMargin(0.40);
  hleg2->SetBorderSize(0);
  if(GetGEreso()) hleg2->AddEntry(GetGEreso(),"Resolution","L");
  if(GetGEbias()) hleg2->AddEntry(GetGEbias(),"Bias","LP");
  hleg2->Draw();
  delete dummye;
    }
  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();

  // residuals
  canvas->cd(4);
  dummya->SetMinimum(-10);
  dummya->SetMaximum(10);
  dummya->SetTitle("Residuals");
  dummya->SetXTitle("log_{10}(E' [GeV])");
  if(hdNdEpBkg)
    dummya->SetYTitle("(Data-dN/dE')/ #Delta Data");
  else
    dummya->SetTitle("(On-Off)/ #Delta On");
  dummya->DrawCopy();
  TH1F* hResidualsOn  = NULL;

  if(hdNdEpBkg && hOn)
    {
      hResidualsOn = GetResidualsHisto(hdNdEpBkg,hOn);
      hResidualsOn->DrawCopy("esame");
    }

  TLegend* hleg3 = new TLegend(0.65, 0.65, 0.90, 0.90);
  hleg3->SetFillColor(0);
  hleg3->SetMargin(0.40);
  hleg3->SetBorderSize(0);
  hleg3->AddEntry(hResidualsOn,"On events","LP");
  hleg3->Draw();

  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();

  TH1I *dummyabis = new TH1I("dummyabis", "dN/dE' bkg model vs On distribution",1,TMath::Log10(10.),TMath::Log10(100000.));
  dummyabis->SetStats(0);
  // dN/dE for signal
  canvas->cd(5);
  dummyabis->SetMinimum(1e-7);
  dummyabis->SetMaximum(1e1);
  dummyabis->SetTitle("dN/dE for signal events");
  dummyabis->SetXTitle("log_{10}(E [GeV])");
  dummyabis->SetYTitle("dN/dE [GeV^{-1}]");
  dummyabis->DrawCopy();
  if(GetHdNdESignal())
    {
      Double_t scale = GetdNdESignalIntegral();
      GetHdNdESignal()->Scale(scale);
      GetHdNdESignal()->DrawCopy("same");
      GetHdNdESignal()->Scale(1./scale);
    }
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->Modified();
  gPad->Update();

  canvas->cd(6);
  dummyabis->SetMinimum(1e2);
  dummyabis->SetMaximum(1e8);
  dummyabis->SetTitle("dN/dE' (#times Aeff) for signal events");
  dummyabis->SetXTitle("log_{10}(E' [GeV])");
  dummyabis->SetYTitle("dN/dE'(#times A_{eff}) [GeV^{-1}cm^{2}]");
  dummyabis->DrawCopy();
  TH1F* hdNdEpSignal = GetHdNdEpSignal();
  if(hdNdEpSignal)
    {
      Double_t scale = GetdNdEpSignalIntegral();
      hdNdEpSignal->Scale(scale);
      hdNdEpSignal->DrawCopy("same");
      hdNdEpSignal->Scale(1./scale);
    }
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->Modified();
  gPad->Update();

  delete hOn;
  if(hResidualsOn)  delete hResidualsOn;
  delete dummya;
  delete dummyabis;
  return canvas;
}

/////////////////////////////////////////////////////////////////
//
// Conversion histogram of event count into dNdE histogram
void differentiate(TH1F* h, Int_t error_flag=0)
{
    Int_t nbins = h->GetXaxis()->GetNbins();
    // divide by bin width
    if(h->GetEntries()>0)
        for(Int_t ibin=0;ibin<nbins;ibin++)
      {
        Double_t leminbin = h->GetBinLowEdge(ibin+1);
        Double_t lemaxbin = leminbin+h->GetBinWidth(ibin+1);
        Double_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);
          if(error_flag==0){
              h->SetBinError(ibin+1,h->GetBinError(ibin+1)/deltaE);
          }
          h->SetBinContent(ibin+1,h->GetBinContent(ibin+1)/deltaE);
      }
}

////////////////////////////////////////////////////////////////////////
// 
// Given a histogram hModel representing a model
// and another one hData representing data (with properly computed errors)
// return the histogram of residuals (hData-hModel)
//
TH1F* GetResidualsHisto(TH1F* hModel,TH1F* hData)
{
  // basic check
  if(!hModel || !hData) return NULL;

  // get number of data bins
  UInt_t nbins = hData->GetNbinsX();
 
  // check if bin width is constant
  Bool_t binWidthIsConstant = kTRUE;
  for(Int_t ibin=0;ibin<nbins-1;ibin++)
    if(TMath::Abs(Float_t(hData->GetBinWidth(ibin+1))-Float_t(hData->GetBinWidth(ibin+2)))>1e-9)
      binWidthIsConstant = kFALSE;

  // create reiduals histo
  TH1F* hResiduals;
  if(binWidthIsConstant)
    hResiduals = new TH1F("hResiduals","Residuals data-Model",nbins,hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  else
    hResiduals = new TH1F("hResiduals","Residuals data-Model",nbins,hData->GetXaxis()->GetXbins()->GetArray());
  hResiduals->SetDirectory(0);

  // fill histo
  for(Int_t ibin=0;ibin<nbins;ibin++)
    {
      Double_t letest = hResiduals->GetBinLowEdge(ibin+1)+hResiduals->GetBinWidth(ibin+1)/2.;
      Int_t modbin = hModel->FindBin(letest);
      Double_t errbar = hData->GetBinError(ibin+1);
      if(errbar>0)
        {
          hResiduals->SetBinContent(ibin+1,(hData->GetBinContent(ibin+1)-hModel->GetBinContent(modbin))/errbar);
          hResiduals->SetBinError(ibin+1,1);
        }
      else
        {
          hResiduals->SetBinContent(ibin+1,0);
          hResiduals->SetBinError(ibin+1,0);
        }
    }
  return hResiduals;
}

////////////////////////////////////////////////////////////////////////
// Line search likelihood function (-2logL) 
// To be minimized by TMinuit
// Free parameters:
// par[0] = g (total estimated number of signal events in On region)
// par[1] = b (total estimated number of background events in On region)
// par[2] = estimated value of tau
//
void lineSearchLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;

  // get internal object, histos, values, etc
  LineSearchLkl*   mylkl           = dynamic_cast<LineSearchLkl*>(minuit->GetObjectFit());
  const TH1F*      hdNdEpSignal    = mylkl->GetHdNdEpSignal();
  const TH1F*      hdNdEpBkg       = mylkl->GetHdNdEpBkg();
  const Float_t*   onSample        = mylkl->GetOnSample();
  UInt_t           Non             = mylkl->GetNon();
  Float_t          Tobs            = mylkl->GetObsTime();

  const Int_t      nbins           = hdNdEpBkg->GetNbinsX();
  const Double_t   xmin            = hdNdEpBkg->GetXaxis()->GetXmin();
  const Double_t   xmax            = hdNdEpBkg->GetXaxis()->GetXmax();
  
  //get energy window size
  Float_t lowE = mylkl->GetEmin();
  Float_t highE = mylkl->GetEmax();
 
  Float_t tau = mylkl->GetTau();
  Float_t dTau = mylkl->GetDTau();

  // Estimated number of signal and background events in signal region
  Double_t g = par[0];
  Double_t b = par[1];

  //new parameter for background normalization factor
  Double_t tauest = par[2];
    
  //Double_t fnorm   = g+b; // original for line search
  Double_t boff = b*tauest;
  Double_t fnorm   = g+boff;

  // sum signal and background (and maybe foreground) contributions and normalize resulting pdf (only On)
  TH1F* hdNdEpOn  = new TH1F("hdNdEpOn", "On  event rate vs E'",nbins,xmin,xmax);
  hdNdEpOn->Reset();
  hdNdEpOn->Add(hdNdEpSignal,hdNdEpBkg,g,boff);

  // normalize
  if(fnorm>0)
    hdNdEpOn->Scale(1./fnorm);
  else
    mylkl->NormalizedNdEHisto(hdNdEpOn);

  // -2 log-likelihood
  f = 0;

  // On events
  for(ULong_t ievent=0; ievent<Non; ievent++)
    {
      Float_t val = hdNdEpOn->GetBinContent(hdNdEpOn->FindBin(onSample[ievent]));
      if(onSample[ievent] > TMath::Log10(lowE) && onSample[ievent] < TMath::Log10(highE))
        {
          if(val>0)
            f += -2*TMath::Log(val);
          else
            f += 1e99;
        }
      else
        f += 0;
    }

  // nuisance tau
    if(dTau>0)
      f+=-2*TMath::Log(TMath::Gaus(tauest, tau, dTau, kTRUE));
    
  // tot Nevts
  if(g+b>0)
    f += -2*TMath::Log(TMath::Poisson(mylkl->GetEventsInEnergyWindow(),g+boff));
  else
    f += 1e99;

  delete hdNdEpOn;
}
