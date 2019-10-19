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
// Here we consider the minimal case with just the "g" free parameter
// As many nuisance parameters as wanted can be added to the list
// (check also LineSearchLkl::SetFunctionAndPars)

static const Int_t    gNPars           = 2;                      // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"g","b"};              // Name of parameters

// static functions (for internal processing of input data)
static Int_t SmearHistogram(TH1F* sp,TH1F* smsp,TGraph* grreso,TGraph* grbias);
static TH1F* GetResidualsHisto(TH1F* hModel,TH1F* hData);

// -2logL function for minuit
// you can change its name but the format is required by ROOT, you should not change it
// if you change the name, do it consistently in all this file
void lineSearchLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
// also, do not change this
static TMinuit* minuit = NULL;


//////////////////////////////////////////////////////////////////////////////
//
// String constructor
// The string contains the elements for the constructor
// which you can read from e.g. an rc input file
//
LineSearchLkl::LineSearchLkl(TString inputString) :
  Lkl(gNPars,inputString,gName,gTitle), Iact1dUnbinnedLkl(inputString),
  fHdNdEpSignal(NULL), fHdNdEpBkg(NULL), fRelativePeakIntensity(0.2), fBkgRegionWidth(0.1)
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
  if(fHdNdEpBkg)       delete fHdNdEpBkg;
  if(fHdNdEpSignal)    delete fHdNdEpSignal;
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
  Double_t pStart[gNPars] = {ginit, (Double_t)Iact1dUnbinnedLkl::GetNon()};
  Double_t pDelta[gNPars] = {TMath::Sqrt(Iact1dUnbinnedLkl::GetNon())/10.,TMath::Sqrt(Iact1dUnbinnedLkl::GetNon())/10.}; // Precision of parameters during minimization

  // initialize the free (and nuisance) parameters
  SetParameters(gParName,pStart,pDelta);

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
  const Float_t* onSample = Iact1dUnbinnedLkl::GetOnSample();
  UInt_t           Non             = Iact1dUnbinnedLkl::GetNon();

  TH1F* hdNdEpOn = new TH1F("hdNdEpOn","dN/dE' for On data",Iact1dUnbinnedLkl::GetNFineBins(),Iact1dUnbinnedLkl::GetFineLEMin(),Iact1dUnbinnedLkl::GetFineLEMax());

  for(ULong_t ievent=0; ievent<Non; ievent++)
    {
      hdNdEpOn->Fill(onSample[ievent]);
    }

  Float_t PeakMaximum = GetHdNdEpSignal()->GetMaximum();
  Int_t binLowBkgLow = 0;
  Bool_t binLowBkgLowBool = kFALSE;
  Int_t binLowBkgHigh = 0;
  Bool_t binLowBkgHighBool = kFALSE;
  Int_t binHighBkgLow = 0;
  Bool_t binHighBkgLowBool = kFALSE;
  Int_t binHighBkgHigh = 0;
  Bool_t binHighBkgHighBool = kFALSE;

  for(Int_t i=1; i< GetHdNdEpSignal()->GetNbinsX()-1; i++)
    {
      if (GetHdNdEpSignal()->GetBinContent(i) > PeakMaximum*(fRelativePeakIntensity-fBkgRegionWidth) && !binLowBkgLowBool)
        {
          binLowBkgLow=i;
          binLowBkgLowBool = kTRUE;
        }

      if (GetHdNdEpSignal()->GetBinContent(i) < PeakMaximum*(fRelativePeakIntensity+fBkgRegionWidth) && !binLowBkgHighBool)
        {
          binLowBkgHigh=i;
        }
      else binLowBkgHighBool = kTRUE;
    }

  for(Int_t i=binLowBkgHigh+1; i< GetHdNdEpSignal()->GetNbinsX()-1; i++)
    {
      if (GetHdNdEpSignal()->GetBinContent(i) > PeakMaximum*(fRelativePeakIntensity+fBkgRegionWidth) && !binHighBkgLowBool)
        {
          binHighBkgLow=i;
        }
      else binHighBkgLowBool = kTRUE;

      if (GetHdNdEpSignal()->GetBinContent(i) < PeakMaximum*(fRelativePeakIntensity-fBkgRegionWidth) && !binHighBkgHighBool)
        {
          binHighBkgHigh=i;
          binHighBkgHighBool = kTRUE;
        }
    }

  if(binLowBkgHigh < 100) binLowBkgLow = 1;
  else binLowBkgLow = binLowBkgHigh -100;
  if(binHighBkgLow > 4900) binHighBkgHigh = 4999;
  else binHighBkgHigh = binHighBkgLow + 100;

  fHdNdEpBkg = new TH1F("HdNdEpBkg","dN/dE' for background model",Iact1dUnbinnedLkl::GetNFineBins(),Iact1dUnbinnedLkl::GetFineLEMin(),Iact1dUnbinnedLkl::GetFineLEMax());

  TH1D* h1 = new TH1D("h1","h1",Iact1dUnbinnedLkl::GetNFineBins(),Iact1dUnbinnedLkl::GetFineLEMin(),Iact1dUnbinnedLkl::GetFineLEMax());

  for(Int_t ibin=1; ibin < Iact1dUnbinnedLkl::GetNFineBins()-1; ibin++)
    {
      if(ibin >= binLowBkgLow && ibin <= binLowBkgHigh) {h1->SetBinContent(ibin,hdNdEpOn->GetBinContent(ibin));}
      if(ibin >= binHighBkgLow && ibin <= binHighBkgHigh) {h1->SetBinContent(ibin,hdNdEpOn->GetBinContent(ibin));}
    }

  TF1* f1 = new TF1("f1","pol1",Iact1dUnbinnedLkl::GetFineLEMin(),Iact1dUnbinnedLkl::GetFineLEMax());
  h1->Fit("f1","0","",Iact1dUnbinnedLkl::GetFineLEMin(),Iact1dUnbinnedLkl::GetFineLEMax());

  for(Int_t ibin=1; ibin < Iact1dUnbinnedLkl::GetNFineBins()-1; ibin++)
    {
      if(binLowBkgHigh < ibin && ibin < binHighBkgLow)
        fHdNdEpBkg->SetBinContent(ibin,f1->Eval(fHdNdEpBkg->GetBinCenter(ibin)));
      else
        fHdNdEpBkg->SetBinContent(ibin,hdNdEpOn->GetBinContent(ibin));
    }

  Iact1dUnbinnedLkl::NormalizedNdEHisto(fHdNdEpBkg);

  delete hdNdEpOn;
  delete h1;
  delete f1;

  return 0;
}

////////////////////////////////////////////////////////////////
//
// Check that the needed histograms are present and that we are ready 
// for calling minimization 
//
// If needed, convolute dNdESignal*Aeff with Eres and Ebias
// to get dNdEpSignal
//
// Return 0 in case of success, 1 otherwise
//
Int_t LineSearchLkl::CheckHistograms(Bool_t checkdNdEpBkg)
{

  // if fHdNdEpSignal is missing, try to construct it from fHdNdESignal, fHAeff fGEreso and fGEbias
  if(!fHdNdEpSignal && (Iact1dUnbinnedLkl::GetHdNdESignal() && Iact1dUnbinnedLkl::GetHAeff() && ((Iact1dUnbinnedLkl::GetGEreso() && Iact1dUnbinnedLkl::GetGEbias()) || Iact1dUnbinnedLkl::GetMigMatrix())))
    {
      cout << "LineSearchLkl::CheckHistograms Message: will create fHdNdEpSignal from fHdNdESignal, fHAeff, fGEreso & fGEbias... " << flush;

      // multiply dNdESignal times Aeff
      TH1F* hdNdESignalAeff = new TH1F("hdNdESignalAeff","Product of dN/dE for Signal and Aeff",Iact1dUnbinnedLkl::GetNFineBins(),Iact1dUnbinnedLkl::GetFineLEMin(),Iact1dUnbinnedLkl::GetFineLEMax());
      hdNdESignalAeff->SetDirectory(0);
      hdNdESignalAeff->Multiply(Iact1dUnbinnedLkl::GetHAeff(), Iact1dUnbinnedLkl::GetHdNdESignal());

      // create fHdNdEpSignal   
      fHdNdEpSignal = new TH1F("fHdNdEpSignal","dN/dE' for Signal",Iact1dUnbinnedLkl::GetNFineBins(),Iact1dUnbinnedLkl::GetFineLEMin(),Iact1dUnbinnedLkl::GetFineLEMax());
      fHdNdEpSignal->SetDirectory(0);

      // smear hdNdESignalAeff
      if(SmearHistogram(hdNdESignalAeff,fHdNdEpSignal,Iact1dUnbinnedLkl::GetGEreso(),Iact1dUnbinnedLkl::GetGEbias()))
        return 1;

      cout << "Done! " << endl;
      // clean
      delete hdNdESignalAeff;
    }

  // normalize unnormalized histos
  Iact1dUnbinnedLkl::NormalizedNdEHisto(fHdNdEpSignal);

  ComputeBkgModelFromOnHisto();

  if(checkdNdEpBkg)
    Iact1dUnbinnedLkl::NormalizedNdEHisto(fHdNdEpBkg);

  // if there are the dNdE' histograms for signal and background + data we're ready to go
  if(checkdNdEpBkg)
    {
      if(fHdNdEpBkg && fHdNdEpSignal && Iact1dUnbinnedLkl::GetOnSample())
        return 0;
    }
  else
    {
      if(fHdNdEpSignal && Iact1dUnbinnedLkl::GetOnSample())
        return 0;
    }

  if(checkdNdEpBkg && !fHdNdEpBkg)
    cout << "LineSearchLkl::CheckHistograms Warning: fHdNdEpBkg histogram missing!!" << endl;
  if(!fHdNdEpSignal)
    cout << "LineSearchLkl::CheckHistograms Warning: fHdNdEpSignal histogram missing!!" << endl;
  if(!Iact1dUnbinnedLkl::GetOnSample())
    cout << "LineSearchLkl::CheckHistograms Warning: fOnSample histogram missing!!" << endl;

    return 1;
}

////////////////////////////////////////////////////////////////
// Smear a histogram <sp> (in log of true energy) using 
// the energy dispersion function described by grreso and grbias
// which are energy resolution and bias vs log(E[GeV])
// and put the result in histogram <smsp> (in measured energy)
//
// Return 0 in case of success
//        1 otherwise
//
Int_t SmearHistogram(TH1F* sp,TH1F* smsp,TGraph* grreso,TGraph* grbias)
{
  // checks
  if(!sp || !smsp || !grreso || !grbias)
    {
      cout << "SmearHistogram Warning: missing histos" << endl;
      return 1;
    }

  Int_t nbinste = sp->GetNbinsX();              // number of bins in input histo
  Int_t nbinsme = smsp->GetNbinsX();            // number of bins in output histo

  // do the convolution of sp with energy resolution and store result in smsp
  for(Int_t ibinte=0;ibinte<nbinste;ibinte++)
    {
      Double_t let  = sp->GetBinCenter(ibinte+1); // log true energy
      Double_t et   = TMath::Power(10,let);       // true energy
      Double_t reso = grreso->Eval(let);          // energy resolution
      Double_t bias = grbias->Eval(let);          // energy bias

      // bin size in true enery
      Double_t minbinte = TMath::Power(10,sp->GetBinLowEdge(ibinte+1));
      Double_t maxbinte = TMath::Power(10,sp->GetBinLowEdge(ibinte+1)+sp->GetBinWidth(ibinte+1));
      Double_t det      = maxbinte-minbinte;  // bin size [GeV]

      for(Int_t ibinme=0;ibinme<nbinsme;ibinme++)
        {
          // compute smearing factor
          Double_t minbinme = TMath::Power(10,smsp->GetBinLowEdge(ibinme+1));
          Double_t maxbinme = TMath::Power(10,smsp->GetBinLowEdge(ibinme+1)+smsp->GetBinWidth(ibinme+1));
          Double_t dem      = maxbinme-minbinme;  // bin size [GeV]

          Double_t mingausme    = (minbinme-et*(1+bias))/(et*reso);
          Double_t maxgausme    = (maxbinme-et*(1+bias))/(et*reso);
          Double_t gausintegral = (TMath::Erf(maxgausme/TMath::Sqrt(2))-TMath::Erf(mingausme/TMath::Sqrt(2)))/2.;
          Double_t smfactor     = gausintegral/dem;
          smsp->SetBinContent(ibinme+1,smsp->GetBinContent(ibinme+1)+sp->GetBinContent(ibinte+1)*det*smfactor);
        }
    }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Read dN/dE' for signal from file in the Segue Stereo input format
// produced by Jelena
// Replacement of existing file is allowed
// Return 0 in case of success
//        1 if file is not found
//
Int_t LineSearchLkl::ReaddNdEpSignal(TString filename)
{
  // open file and look for histo
  TFile* dNdEpSignalFile = new TFile(filename);
  TH1F*  hNdEpSignalProv =    (TH1F*) dNdEpSignalFile->Get("fHdNdEpSignal");
  if(!hNdEpSignalProv) {delete dNdEpSignalFile; return 1;}

  // copy histo 
  if(fHdNdEpSignal) delete fHdNdEpSignal;
  fHdNdEpSignal = new TH1F(*hNdEpSignalProv);
  fHdNdEpSignal->SetXTitle("log_{10}(E' [GeV])");
  fHdNdEpSignal->SetYTitle("cm^{2} dN/dE' [GeV^{-1}]");
  fHdNdEpSignal->SetDirectory(0);

  // check what we've read is good
  Int_t status = 0;
  Int_t nbins = fHdNdEpSignal->GetNbinsX();
  Float_t xmin = fHdNdEpSignal->GetXaxis()->GetXmin();
  Float_t xmax = fHdNdEpSignal->GetXaxis()->GetXmax();
  if(nbins!=GetNFineBins() || xmin!=GetFineLEMin() || xmax!=GetFineLEMax())
    {
      cout << "LineSearchLkl::ReaddNdEpSignal, histogram read from file " << filename << " has " << nbins << " (should be " << GetNFineBins() << "), xmin = " << xmin << " (should be " << GetFineLEMin() << ") and xmax = " << xmax << " (should be " << GetFineLEMax() << ")" << endl;
      status = 1;
    }

  // clean and exit
  SetChecked(kFALSE);
  delete dNdEpSignalFile;
  return status;
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
  if(Iact1dUnbinnedLkl::GetHAeff()) Iact1dUnbinnedLkl::GetHAeff()->DrawCopy();
  gPad->SetLogy();
  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();

  // dN/dE' background vs data
  canvas->cd(2);
  TH1F* hOn  = Iact1dUnbinnedLkl::GetHdNdEpOn();
  hOn->SetDirectory(0);

  TH1I *dummya = new TH1I("dummya", "dN/dE' bkg model vs On distribution",1,TMath::Log10(Iact1dUnbinnedLkl::GetEmin()),TMath::Log10(Iact1dUnbinnedLkl::GetEmax()));
  dummya->SetStats(0);
  if(Iact1dUnbinnedLkl::GetNon()>1)
    {
      dummya->SetMinimum(hOn->GetMinimum(0)/2.);
      dummya->SetMaximum(hOn->GetMaximum()*2);
    }
  else if(fHdNdEpBkg)
    {
      dummya->SetMinimum(fHdNdEpBkg->GetMinimum());
      dummya->SetMaximum(fHdNdEpBkg->GetMaximum());
    }

  if(!fHdNdEpBkg) dummya->SetTitle("dN/dE' distributions for On event samples");
  dummya->SetXTitle("log_{10}(E' [GeV])");
  dummya->SetYTitle("dN/dE' [GeV^{-1}]");
  dummya->DrawCopy();

  TH1F* hBkg  = GetHdNdEpModelBkg();
  hBkg->SetDirectory(0);
  if(fHdNdEpBkg)
    {
      Double_t scale = hBkg->GetBinContent(0);
      Double_t scale2 = hOn->GetBinContent(0);
      hBkg->SetLineColor(2);
      hBkg->SetMarkerColor(2);
      hBkg->SetMarkerStyle(2);
      hBkg->DrawCopy("esame");
    }
  hOn->DrawCopy("esame");

  gPad->SetLogy();
  gPad->SetGrid();

  TLegend* hleg = new TLegend(0.65, 0.65, 0.85, 0.85);
  hleg->SetFillColor(0);
  hleg->SetMargin(0.40);
  hleg->SetBorderSize(0);
  hleg->AddEntry(hOn,"On events","LP");
  hleg->AddEntry(hBkg,"Bkg model","LP");
  hleg->Draw();

  gPad->Modified();
  gPad->Update();

  // energy resolution and bias
  canvas->cd(3);

  TH1I *dummye;
  if(Iact1dUnbinnedLkl::GetGEreso())
    dummye = new TH1I("dummye", "Energy resolution and bias",1,Iact1dUnbinnedLkl::GetGEreso()->GetX()[0],Iact1dUnbinnedLkl::GetGEreso()->GetX()[Iact1dUnbinnedLkl::GetGEreso()->GetN()-1]);
  else
    dummye = new TH1I("dummye", "Energy resolution and bias",1,2,4);
  dummye->SetStats(0);
  dummye->SetMinimum(-0.1);
  dummye->SetMaximum(0.4);
  dummye->SetXTitle("log_{10}(E [GeV])");
  dummye->SetYTitle("Energy resolution and bias [%]");
  dummye->DrawCopy();
  if(Iact1dUnbinnedLkl::GetGEreso()) Iact1dUnbinnedLkl::GetGEreso()->Draw("l");
  if(Iact1dUnbinnedLkl::GetGEbias()) Iact1dUnbinnedLkl::GetGEbias()->Draw("l");

  TLegend* hleg2 = new TLegend(0.65, 0.69, 0.90, 0.89);
  hleg2->SetFillColor(0);
  hleg2->SetMargin(0.40);
  hleg2->SetBorderSize(0);
  if(Iact1dUnbinnedLkl::GetGEreso()) hleg2->AddEntry(Iact1dUnbinnedLkl::GetGEreso(),"Resolution","L");
  if(Iact1dUnbinnedLkl::GetGEbias()) hleg2->AddEntry(Iact1dUnbinnedLkl::GetGEbias(),"Bias","LP");
  hleg2->Draw();
  delete dummye;

  gPad->SetGrid();
  gPad->Modified();
  gPad->Update();

  // residuals
  canvas->cd(4);
  dummya->SetMinimum(-5);
  dummya->SetMaximum(5);
  dummya->SetTitle("Residuals");
  dummya->SetXTitle("log_{10}(E' [GeV])");
  if(fHdNdEpBkg)
    dummya->SetYTitle("(Data-dN/dE')/ #Delta Data");
  else
    dummya->SetTitle("(On-Off)/ #Delta On");
  dummya->DrawCopy();
  TH1F* hResidualsOn  = NULL;
  if(fHdNdEpBkg && hOn)
    {
      hResidualsOn  =  GetResidualsHisto(hBkg,hOn);
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

  // dN/dE for signal
  canvas->cd(5);
  dummya->SetMinimum(1e-7);
  dummya->SetMaximum(1e1);
  dummya->SetTitle("dN/dE for signal events");
  dummya->SetXTitle("log_{10}(E [GeV])");
  dummya->SetYTitle("dN/dE [GeV^{-1}]");
  dummya->DrawCopy();
  if(Iact1dUnbinnedLkl::GetHdNdESignal())
    {
      Double_t scale = Iact1dUnbinnedLkl::GetdNdESignalIntegral();
      Iact1dUnbinnedLkl::GetHdNdESignal()->Scale(scale);
      Iact1dUnbinnedLkl::GetHdNdESignal()->DrawCopy("same");
      Iact1dUnbinnedLkl::GetHdNdESignal()->Scale(1./scale);
    }
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->Modified();
  gPad->Update();

  canvas->cd(6);
  dummya->SetMinimum(1e2);
  dummya->SetMaximum(1e8);
  dummya->SetTitle("dN/dE' (#times Aeff) for signal events");
  dummya->SetXTitle("log_{10}(E' [GeV])");
  dummya->SetYTitle("dN/dE'(#times A_{eff}) [GeV^{-1}cm^{2}]");
  dummya->DrawCopy();
  if(fHdNdEpSignal)
    {
      Double_t scale = GetdNdEpSignalIntegral();
      fHdNdEpSignal->Scale(scale);
      fHdNdEpSignal->DrawCopy("same");
      fHdNdEpSignal->Scale(1./scale);
    }
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->Modified();
  gPad->Update();

  delete hOn;
  if(hResidualsOn)  delete hResidualsOn;
  delete dummya;
  return canvas;
}

//////////////////////////////////////////////////////////////////
//
// Produce the E' distribution of On events and return the 
// corresponding histogram.
// if isDifferential=kTRUE (default is kTRUE), return the 
// dN/dE' distribution (i.e. number of entries in the bin divided
// by the size of the bin in E' unit). Otherwise the returned
// histogram is the number of events per bin of E'.
// The returned histogram will not be deleted by the destructor
// of Iact1dUnbinnedLkl.
//
TH1F* LineSearchLkl::GetHdNdEpModelBkg(Bool_t isDifferential,Int_t nbins) const
{
  // we need a positive number of bins
  if(nbins<=0) nbins = gNBins;

  fHdNdEpBkg->Scale(fHdNdEpBkg->GetBinContent(0));

  // create histo
  TH1F* h = new TH1F("dNdEpBkg","dN/dE' for Bkg events",nbins,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
  h->SetDirectory(0);
  h->SetXTitle("log_{10}(E' [GeV])");
  h->SetYTitle("dN/dE' [GeV^{-1}]");

  // fill histo
  for(Int_t i=1;i<fHdNdEpBkg->GetNbinsX()-1;i++)
    {
      for(Int_t j=1;j<fHdNdEpBkg->GetBinContent(i)+0.5;j++)
        {
          h->Fill(fHdNdEpBkg->GetBinCenter(i));
        }
    }

  // divide by bin width
  if(h->GetEntries()>0)
    if(isDifferential)
      for(Int_t ibin=0;ibin<nbins;ibin++)
        {
          Double_t leminbin = h->GetBinLowEdge(ibin+1);
          Double_t lemaxbin = leminbin+h->GetBinWidth(ibin+1);
          Double_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);
          h->SetBinError(ibin+1,h->GetBinError(ibin+1)/deltaE);
          h->SetBinContent(ibin+1,h->GetBinContent(ibin+1)/deltaE);
        }

  fHdNdEpBkg->Scale(1./fHdNdEpBkg->GetBinContent(0));

  return h;
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

//////////////////////////////////////////////////////////////////
// 
// Recreate a fresh new version of fHdNdESignal
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t LineSearchLkl::ResetdNdESignal()
{
  // Delete existing fHdNdEpSignal 
  if(fHdNdEpSignal)
    {
      delete fHdNdEpSignal;
      fHdNdEpSignal=NULL;
    }

  SetChecked(kFALSE);

  // exit
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Create and save a dN/dE histogram for signal according to 
// given <function> and parameters <p0>, <p1>, ...
// (see Iact1dUnbinnedLkl::AdddNdESignalFunction for details on available
// functions)
//
// If fdNdESignal exists it will be replaced
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t LineSearchLkl::SetdNdESignalFunction(TString function,Float_t p0,Float_t p1,Float_t p2,Float_t p3,Float_t p4,Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{

  // Delete existing fHdNdESignal and create empty one
  ResetdNdESignal();

  // call to add the function
  Iact1dUnbinnedLkl::ResetdNdESignal();
  Iact1dUnbinnedLkl::AdddNdESignalFunction(function,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);

  // exit
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Create and save a dN/dE histogram for signal according to 
// given <function> (TF1 object, created and parameterized before
// passing it as an argument)
//
// Set only values for which energy is between <emin> and <emax>
//
// If fdNdESignal exists it will be replaced
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t LineSearchLkl::SetdNdESignalFunction(TF1* function,Float_t emin,Float_t emax)
{
  // Delete existing fHdNdESignal and create empty one
  ResetdNdESignal();

  // call to add the function
  Iact1dUnbinnedLkl::ResetdNdESignal();
  Iact1dUnbinnedLkl::AdddNdESignalFunction(function,emin,emax);

  // exit
  return 0;
}

////////////////////////////////////////////////////////////////////////
// Line search likelihood function (-2logL) 
// To be minimized by TMinuit
// Free parameters:
// par[0] = g (total estimated number of signal events in On region)
// par[1] = b (total estimated number of background events in On region)
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

  // Estimated number of signal and background events in signal region
  Double_t g       = par[0];
  Double_t b       = par[1];
  Double_t fnorm   = g+b;

  // sum signal and background (and maybe foreground) contributions and normalize resulting pdf (only On)
  TH1F* hdNdEpOn  = new TH1F("hdNdEpOn", "On  event rate vs E'",nbins,xmin,xmax);
  hdNdEpOn->Reset();
  
  hdNdEpOn->Add(hdNdEpSignal,hdNdEpBkg,g,b);

  // normalize
  if(fnorm>0)
    hdNdEpOn->Scale(1./fnorm);
  else
    mylkl->Iact1dUnbinnedLkl::NormalizedNdEHisto(hdNdEpOn);

  // -2 log-likelihood
  f = 0;

  // On events
  for(ULong_t ievent=0; ievent<Non; ievent++)
    {
      Float_t val = hdNdEpOn->GetBinContent(hdNdEpOn->FindBin(onSample[ievent]));
      if(val>0)
        f += -2*TMath::Log(val);
      else
        f += 1e99;
    }

  // tot Nevts
  if(g+b>0)
    f += -2*TMath::Log(TMath::Poisson(Non,g+b));
  else
    f += 1e99;

  delete hdNdEpOn;
}		
