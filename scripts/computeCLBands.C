//######################################################################
//##
//## computeCLBands.C
//##
//## Compute confidence level (CL) bands to Dark Matter annihilation
//## cross-section (<sv>) or decay lifetime (tau) using files obtained
//## with the jointLklDM.C script.
//##
//## You can specify two arguments:
//## 1. Name of the rc file (a TString): it should be the same as the 
//##    one used for jointLklDM.C
//## 2. The total number of simulation files (a Int_t) you computed 
//##    using jointLklDM.C. If you ran X simulations, it will assume
//##    that the files have seeds in the range [0,X[
//##
//## The script will automaticaaly stop if more than 5% of the simualtion
//## files are not found.
//##
//## The script produces (and saves) two files:
//## - One .root file containing the median curve, the right anf left 
//##   limits of the 68% CL intervall and the right and left limits of 
//##   the 95% CL intervall.
//## - One canavs limcanvas with the limits on <sv> or tau vs mass for
//##   the considered DM channel, range of masses and data. Overlaid on 
//##   this limit will be the 68% anf 95% CL bands as well as the median
//##   of the simulations.
//##
//## You need to run the compiled version of the script
//## by typing ".x computeCLBands.C+" in the ROOT command interpreter.
//##
//######################################################################

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "TPRegexp.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"

#include "Lkl.h"

using namespace std;

void setDefaultStyle();
Int_t GetNSkippedMasses(Int_t nm,const Double_t* vm,Double_t minm);

// number of bins in histograms
const Int_t nbins = 1000;

// show plots?
const Bool_t  makePlots   = kTRUE;

void computeCLBands(TString configFileName="$GLIKESYS/rcfiles/jointLklDM.rc",Int_t nSimuFiles=300)
{
  setDefaultStyle();

  // Look for configuration file
  gSystem->ExpandPathName(configFileName);
  TPMERegexp re("\\s+");
  
  Lkl::PrintGLikeBanner();
  
  // Print-out configuration info (part 1)
  cout << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***                          RUNNING computeCLBands.C                          ***" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***" << endl;
  cout << "*** CONFIGURATION FILE  : " << configFileName << endl;
  if (gSystem->AccessPathName(configFileName, kFileExists))
    {
      cout << endl << "    Oops! problems reading config file file " << configFileName << " <---------------- FATAL ERROR!!!"<< endl;
      return;
    }

  // Read configuration file
  TEnv*  env = new TEnv(configFileName);
  TString  label             = env->GetValue("jointLklDM.Label","");
  TString  channel           = env->GetValue("jointLklDM.Channel","bb");
  TString  process           = env->GetValue("jointLklDM.Process","ann");
  TString  fInputDataPath    = env->GetValue("jointLklDM.path","");
  TString  fPlotsDir         = fInputDataPath+"/"+env->GetValue("jointLklDM.plotsDir","")+"/";
  TString  massList          = env->GetValue("jointLklDM.MassList","");
 
  // fill the list of masses to be studied
  UInt_t  nmass0  = re.Split(massList);
  Double_t* massval0 = new Double_t[nmass0];
  for(UInt_t imass=0;imass<nmass0;imass++)
    massval0[imass] = re[imass].Atof();

  // Set some flags
  Bool_t  isDecay            = !process.CompareTo("dec",TString::kIgnoreCase);
  
  // define some labels according to input data
  TString strprocess         = (isDecay?       "Decay" : "Annihilation");
  
  // annihilation/decay channel string
  TString  strchannel;
  Double_t minmass = 0;
  if     (!channel.CompareTo("bb",TString::kIgnoreCase))         {strchannel = "b #bar{b}";         minmass=5;}
  else if(!channel.CompareTo("tautau",TString::kIgnoreCase))     {strchannel = "#tau^{+} #tau^{-}"; minmass=1.8;}
  else if(!channel.CompareTo("mumu",TString::kIgnoreCase))       {strchannel = "#mu^{+} #mu^{-}";   minmass=0.106;}
  else if(!channel.CompareTo("WW",TString::kIgnoreCase))         {strchannel = "W^{+} W^{-}";       minmass=80.3;}
  else if(!channel.CompareTo("gammagamma",TString::kIgnoreCase)) {strchannel = "#gamma#gamma";      minmass=0;}
  else if(!channel.CompareTo("pi0pi0",TString::kIgnoreCase))     {strchannel = "#pi^{0}#pi^{0}";    minmass=0.135;}
  else if(!channel.CompareTo("gammapi0",TString::kIgnoreCase))   {strchannel = "#pi^{0}#gamma";     minmass=0.135/2.;}
  else if(!channel.CompareTo("pi0gamma",TString::kIgnoreCase))   {strchannel = "#pi^{0}#gamma";     minmass=0.135/2.;}
  else strchannel = "";
  if(isDecay) minmass*=2;
  
  // Remove mass values below kinematical threshold
  Int_t           nskippedmass = GetNSkippedMasses(nmass0,massval0,minmass);
  Int_t           nmass        = nmass0-nskippedmass; 
  const Double_t* massval      = massval0+nskippedmass;

  // Print-out configuration info (part 2)
  cout << "*** I/O PATH            : " << fInputDataPath << endl;
  cout << "***" << endl;
  cout << "*** LABEL               : " << label <<  endl;
  cout << "*** PROCESS             : " << strprocess << endl;
  cout << "*** CHANNEL             : " << channel << endl;
      
  // how many and which mass values are we considering?
  cout << "*** NUMBER OF MASSES    : " << nmass << endl;
  cout << " ** Mass values         : "; 
  for(Int_t imass=0;imass<nmass;imass++)
    cout << massval[imass] << ((imass<nmass-1)? ", " : "");
  cout << " GeV" << endl;
  
  cout << "***" << endl;
  cout << "***********************************************************************************" << endl;
  cout << "***********************************************************************************" << endl;
  cout << endl;
  // End of print-out configuration info

  // where the root files for the simulations are
  const TString simusPath     = fPlotsDir + "MC/root/";

  // where the root files for the data are
  const TString dataPath      = fPlotsDir + "Data/root/";

  // where the output files will be
  const TString resultsPath   = fInputDataPath + "/results/";
  gSystem->Exec(Form("mkdir -p %s",resultsPath.Data()));
  gSystem->Exec(Form("mkdir -p %s/root",resultsPath.Data()));
  gSystem->Exec(Form("mkdir -p %s/pdf",resultsPath.Data()));

  // Arrays for results
  Double_t sv1sigmaL[nmass];
  Double_t sv1sigmaR[nmass];
  Double_t sv0sigma[nmass];
  Double_t sv2sigmaL[nmass];
  Double_t sv2sigmaR[nmass];

  // vectors to read input data
  vector<Double_t>* svLimDist = new vector<Double_t>[nmass];
  for(Int_t imass=0;imass<nmass;imass++)
    svLimDist[imass].reserve(nSimuFiles);

  // reading input files
  cout << "Reading " << nSimuFiles << " input files" << endl;
  Int_t nFilesRead=0;
  for(Int_t ifile=0;ifile<nSimuFiles;ifile++)
    {
      // file name
      TString seedtag  = Form("_%05d",ifile);
      TString filename = simusPath+label+"_MC_limits"+seedtag+".root";
      cout << "Opening file " << filename << endl;
      Int_t exist = gSystem->Exec("ls " + filename + " > /dev/null");
      if ( exist != 0 ){
        cout << "cannot find " << filename << endl << " ... skipping" << endl;
        continue;
      }
      else
        nFilesRead++;

      // open file and get the limits
      TFile*   infile  = TFile::Open(filename,"READ");
      TCanvas* canvas  = (TCanvas*) infile->Get("limcanvas");
      TGraph*  graph   = (TGraph*)  canvas->FindObject(isDecay? "grtaulim": "grsvlim");

      Char_t line[256];
      for(Int_t ival=0;ival<graph->GetN();ival++)
        {
          Double_t mass = graph->GetX()[ival];
          for(Int_t imass=0;imass<nmass;imass++)
            if(mass==massval[imass])
              {
                Double_t lim  = graph->GetY()[ival];
                svLimDist[imass].push_back(lim);
                break;
              }
        }
      delete infile;
    }
  cout << "Could read " << nFilesRead << "/" << nSimuFiles << " of the input files." << endl << endl;

  // if we miss 5% or more of the simulation files
  if((Double_t)nFilesRead/(Double_t)nSimuFiles <0.95)
    {
      cout << "Too few MC files could be read compare to the expected number of files. Did you really run " << nSimuFiles << " simulations?" << endl;
      return;
    }

  cout << "Computing the intervals for all masses..." << endl;
  // sort elements in the limit vector and put values in histos
  Int_t shift=0;
  TH1D* hsvLimDist[nmass];
  for(Int_t imass=0;imass<nmass;imass++)
    {
      hsvLimDist[imass]=NULL;
      // skip mass values with no associated data
      if(!svLimDist[imass].size())
        {
          shift = imass+1;
          continue;
        }

      // sort elements and fill them into histogram
      sort(svLimDist[imass].begin(),svLimDist[imass].end());
      hsvLimDist[imass] = new TH1D(Form("hsvLimDist%02d",imass),Form("hsvLimDist%02d",imass),nbins,svLimDist[imass][0]*0.99,svLimDist[imass][svLimDist[imass].size()-1]*1.01);
      hsvLimDist[imass]->SetStats(0);

      for(UInt_t ival=0;ival<svLimDist[imass].size();ival++)
        hsvLimDist[imass]->Fill(svLimDist[imass][ival]);

      // find the percentiles
      Double_t sv2sigmaLDist = 9e99;
      Double_t sv1sigmaLDist = 9e99;
      Double_t sv0sigmaDist  = 9e99;
      Double_t sv1sigmaRDist = 9e99;
      Double_t sv2sigmaRDist = 9e99;
      Double_t norm = hsvLimDist[imass]->Integral("width");
      for(Int_t ibin=0;ibin<nbins;ibin++)
        {
          Double_t perc = hsvLimDist[imass]->Integral(1,ibin+1,"width")/norm;
          if(TMath::Abs(perc-0.02275)<sv2sigmaLDist) {sv2sigmaLDist=TMath::Abs(perc-0.02275);sv2sigmaL[imass]=hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
          if(TMath::Abs(perc-0.1587)< sv1sigmaLDist) {sv1sigmaLDist=TMath::Abs(perc-0.1587); sv1sigmaL[imass]=hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
          if(TMath::Abs(perc-0.5)<    sv0sigmaDist)  {sv0sigmaDist= TMath::Abs(perc-0.5);    sv0sigma[imass]= hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
          if(TMath::Abs(perc-0.8414)< sv1sigmaRDist) {sv1sigmaRDist=TMath::Abs(perc-0.8414); sv1sigmaR[imass]=hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
          if(TMath::Abs(perc-0.97725)<sv2sigmaRDist) {sv2sigmaRDist=TMath::Abs(perc-0.97725);sv2sigmaR[imass]=hsvLimDist[imass]->GetBinLowEdge(ibin+1)+hsvLimDist[imass]->GetBinWidth(ibin+1);}
        }
      cout << "Mass = " << massval[imass] << ", read " << svLimDist[imass].size() << " values, limits: " <<sv2sigmaL[imass] << " " << sv1sigmaL[imass]<< " "<< sv0sigma[imass]<< " " <<  sv1sigmaR[imass] << " " << sv2sigmaR[imass] << endl;
    }

  // write results in a better format for using them later
  cout << endl;
  cout << "READY-TO-COPY results: " << endl;
  cout << endl;

  cout << "Double_t sv0sigma_"+channel+"[nmass]  = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv0sigma[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  cout << "Double_t sv1sigmaL_"+channel+"[nmass] = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv1sigmaL[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  cout << "Double_t sv2sigmaL_"+channel+"[nmass] = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv2sigmaL[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  cout << "Double_t sv1sigmaR_"+channel+"[nmass] = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv1sigmaR[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  cout << "Double_t sv2sigmaR_"+channel+"[nmass] = {";
  for(Int_t imass=shift;imass<nmass;imass++)
    cout << sv2sigmaR[imass] << (imass<nmass-1? "," : "");
  cout << "};" << endl;

  // Save results as graphs in a root file
  TString outputFileName = resultsPath+"root/"+label+"_1and2sigmasBands.root";
  TFile*  outputFile     = new TFile(outputFileName,"RECREATE");

  TGraph* gsv0sigma  = new TGraph(nmass-shift,massval+shift,sv0sigma+shift);
  TGraph* gsv1sigmaL = new TGraph(nmass-shift,massval+shift,sv1sigmaL+shift);
  TGraph* gsv1sigmaR = new TGraph(nmass-shift,massval+shift,sv1sigmaR+shift);
  TGraph* gsv2sigmaL = new TGraph(nmass-shift,massval+shift,sv2sigmaL+shift);
  TGraph* gsv2sigmaR = new TGraph(nmass-shift,massval+shift,sv2sigmaR+shift);
  gsv0sigma->SetName("gsv0sigma");
  gsv1sigmaL->SetName("gsv1sigmaL");
  gsv1sigmaR->SetName("gsv1sigmaR");
  gsv2sigmaL->SetName("gsv2sigmaL");
  gsv2sigmaR->SetName("gsv2sigmaR");
  gsv0sigma->Write();
  gsv1sigmaL->Write();
  gsv1sigmaR->Write();
  gsv2sigmaL->Write();
  gsv2sigmaR->Write();

  cout << endl << "Bands saved in file: " << outputFileName << endl;
  delete outputFile;

  // PLOTS
  if(makePlots)
    {
      // style
      setDefaultStyle();
      gStyle->SetPadRightMargin(0.12);
      gStyle->SetPadLeftMargin(0.12);

      // canvas for distributions
      const Int_t nlines = 5;  // number of lines in parabolas canvas
      Int_t ncols = TMath::Ceil((nmass-shift)/Float_t(nlines)); // number of columns in parabolas canvas
      TCanvas* c = new TCanvas("canvas","limits distributions",ncols*250,nlines*250);
      c->Divide(ncols,nlines);

      // loop over masses and plot histos
      for(Int_t imass=shift;imass<nmass;imass++)
        {         
          c->cd(imass+1);
          hsvLimDist[imass]->DrawCopy("e");
          gPad->SetLogy();
          gPad->SetLogx();
        }
      c->Print(resultsPath+"root/"+label+"_distributionLimits.root");
      c->Print(resultsPath+"pdf/"+label+"_distributionLimits.pdf");

      // bands
      Int_t npnts = gsv0sigma->GetN();
      TGraph* band1sigma = new TGraph(2*npnts+1);
      TGraph* band2sigma = new TGraph(2*npnts+1);

      for(Int_t imass=0;imass<npnts;imass++)
        {
          band1sigma->SetPoint(imass,      gsv1sigmaL->GetX()[imass],gsv1sigmaL->GetY()[imass]);
          band1sigma->SetPoint(npnts+imass,gsv1sigmaR->GetX()[npnts-imass-1],gsv1sigmaR->GetY()[npnts-imass-1]);
          band2sigma->SetPoint(imass,      gsv2sigmaL->GetX()[imass],gsv2sigmaL->GetY()[imass]);
          band2sigma->SetPoint(npnts+imass,gsv2sigmaR->GetX()[npnts-imass-1],gsv2sigmaR->GetY()[npnts-imass-1]);
        }
      band1sigma->SetPoint(2*npnts,gsv1sigmaL->GetX()[0],gsv1sigmaL->GetY()[0]);
      band2sigma->SetPoint(2*npnts,gsv2sigmaL->GetX()[0],gsv2sigmaL->GetY()[0]);

      // result
      TFile* f = TFile::Open(dataPath+label+"_Data_limits.root", "READ");
      TCanvas* lim = (TCanvas*)f->Get("limcanvas");
      TGraph* grsvlim = (TGraph*)lim->GetPrimitive(isDecay? "grtaulim": "grsvlim");
      grsvlim->SetName(isDecay? "grtaulim": "grsvlim");
      grsvlim->SetLineColor(kBlack);
      grsvlim->SetLineWidth(2);

      band1sigma->SetFillColorAlpha(3,0.35);
      band2sigma->SetFillColorAlpha(5,0.35);
      gsv0sigma->SetLineColor(4);
      gsv0sigma->SetLineStyle(2);
      gsv0sigma->SetLineWidth(1);

      // canvas
      TCanvas* limcanvas = new TCanvas("limcanvas","1 and 2 sigma bands",1000,800);
      limcanvas->cd();

      TH1I *dummylim = new TH1I("dummylim","",1,massval[0],massval[nmass-1]);
      dummylim->SetStats(0);
      dummylim->SetMinimum((isDecay? 1e23 : 1e-28));
      dummylim->SetMaximum((isDecay? 1e27 : 1e-20));
      dummylim->SetXTitle("m_{DM} [GeV]");
      dummylim->SetYTitle(Form("95%% %s [%s]",(isDecay? "#tau_{DM}^{LL}" : "<#sigma v>^{UL}"),(isDecay? "s" : "cm^{3}/s")));
      dummylim->DrawCopy();
      gPad->SetLogx();
      gPad->SetLogy();
      gPad->SetGrid();

      band2sigma->Draw("f");
      band1sigma->Draw("f");
      gsv0sigma->Draw("c");
      grsvlim->Draw("c");

      // thermal relic density taken from Steigman G., Dasgupta B, and Beacom J. F., 
      // Precise relic WIMP abundance and its impact onsearches for dark matter annihilation, 
      // Phys.Rev. D86(2012) 023506, [arXiv:1204.3622]
      const Int_t nThermalRelicSv = 22;
      Double_t thermalRelicMass[nThermalRelicSv] = {1.00e-01, 1.78e-01, 3.16e-01, 5.62e-01, 1.00e+00, 1.78e+00, 3.16e+00, 5.62e+00, 1.00e+01, 1.78e+01, 3.16e+01, 5.62e+01, 1.00e+02, 1.78e+02, 3.16e+02, 5.62e+02, 1.00e+03, 1.78e+03,3.16e+03, 5.62e+03, 1.00e+04, 1.00e+05};
      Double_t thermalRelicSigmav[nThermalRelicSv] = {4.8e-26, 4.9e-26, 5.1e-26, 5.0e-26, 4.7e-26, 4.5e-26, 3.9e-26, 2.8e-26, 2.5e-26, 2.3e-26, 2.2e-26, 2.2e-26, 2.2e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26,2.4e-26, 2.4e-26};
      TGraph *relicDensity = new TGraph(nThermalRelicSv, thermalRelicMass, thermalRelicSigmav);
      relicDensity->SetTitle("Thermal relic cross section");
      relicDensity->SetLineColor(kRed);
      relicDensity->SetLineWidth(3);
      relicDensity->SetLineStyle(5);
      relicDensity->SetMarkerColor(kRed);
      if(!isDecay) relicDensity->Draw("same");

      // legend
      TLegend* limleg;
      limleg = new TLegend(0.35, 0.14, 0.85, 0.34);
      limleg->SetTextSize(0.032);
      limleg->SetMargin(0.40);
      limleg->SetBorderSize(0);

      if(grsvlim)    limleg->AddEntry(grsvlim,"Limit","l");
      if(gsv0sigma)  limleg->AddEntry(gsv0sigma,"H_{0} median","l");
      if(band1sigma) limleg->AddEntry(band1sigma,"H_{0} 68% containment","f");
      if(band2sigma) limleg->AddEntry(band2sigma,"H_{0} 95% containment","f");
      if(!isDecay)   limleg->AddEntry(relicDensity,"Thermal relic cross section","l");
      limleg->Draw();

      TLatex* txchannel = new TLatex(0.75,0.4,strchannel);
      txchannel->SetTextSize(0.05);
      txchannel->SetNDC();
      txchannel->Draw();

      limcanvas->Print(resultsPath+"root/"+label+"_bands.root");
      limcanvas->Print(resultsPath+"pdf/"+label+"_bands.pdf");
    }
}

Int_t GetNSkippedMasses(Int_t nm,const Double_t* vm,Double_t minm)
{
  Int_t im=0;
  for(im=0;im<nm;im++)
    if(vm[im]>minm)
      break;
  return im;
}

void setDefaultStyle()
{
  // general settings
  //gStyle->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTextSize(1);
  gStyle->SetTitleSize(0.04, "xy");
  gStyle->SetTitleOffset(1.3, "x");
  gStyle->SetTitleOffset(1.6,"y");
  gStyle->SetTitleFillColor(0);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetOptStat(111110);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  //gStyle->SetColorModelPS(1);
  gStyle->SetPalette(1,0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetGridWidth(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetTickLength(0.04, "xy");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetLegendBorderSize(4);
  gStyle->SetGridColor(13);
}
