// test input unbinned fits files and create an IactEventListIrf from them
void testInputFITS(Bool_t checkInputs=kTRUE){
    TString dataFileName = "../data/fits/obs_id_23523_unbinned_spectrum.fits";

    // assign a TFITSHDU to each header data unit in the file
    TFITSHDU* hduHeader = new TFITSHDU(dataFileName + "[0]");
    TFITSHDU* hduOn = new TFITSHDU(dataFileName + "[1]");
    TFITSHDU* hduOff = new TFITSHDU(dataFileName + "[2]");
    TFITSHDU* hduAeff = new TFITSHDU(dataFileName + "[3]");
    TFITSHDU* hduEdisp = new TFITSHDU(dataFileName + "[4]");

    // initialise the empty IactEventListIRF
    IactEventListIrf* iactDataSet = new IactEventListIrf();

    // fetch effective time and Off to On exposure ratio (tau)
    Double_t obsTime = hduHeader->GetKeywordValue("TEFF").Atof();
    Double_t numberOnRegions = hduHeader->GetKeywordValue("ACC").Atof();
    Double_t numberOffRegions = hduHeader->GetKeywordValue("ACC_OFF").Atof();
    Double_t tau = numberOffRegions / numberOnRegions;
    iactDataSet->SetObsTime(obsTime);
    iactDataSet->SetTau(tau);

    // fetch the list of On events (energies in TeV)
    TVectorD* onEnergies = hduOn->GetTabRealVectorColumn("ENERGY");
    TVectorD* onTimes = hduOn->GetTabRealVectorColumn("TIME");
    for(Int_t i=0; i<onEnergies->GetNrows(); ++i){
    	// fOnSample->Fill(E, pointRA, pointDEC, dRA, dDEC, t, had);
    	// cout << "E / TeV : " << onEnergies[0][i] << " -- time / MJD : " << onTimes[0][i] << endl;
        iactDataSet->FillOnEvent(onEnergies[0][i] * 1e3, 0., 0., 0., 0., onTimes[0][i], 0.);
    }

    // fetch the list of Off events (energies in TeV)
    TVectorD* offEnergies = hduOff->GetTabRealVectorColumn("ENERGY");
    TVectorD* offTimes = hduOff->GetTabRealVectorColumn("TIME");
    for(Int_t i=0; i<offEnergies->GetNrows(); ++i){
    	// fOffSample->Fill(E,pointRA,pointDEC,dRA,dDEC,t,had);
    	// cout << "E / TeV : " << offEnergies[0][i] << " -- time / MJD : " << offTimes[0][i] << endl;
        iactDataSet->FillOffEvent(offEnergies[0][i] * 1e3, 0., 0., 0., 0., offTimes[0][i], 0.);
    }

    // fetch the effective area
    TVectorD* lowEtrueEdges = hduAeff->GetTabRealVectorColumn("ENERG_LO");
    TVectorD* highEtrueEdges = hduAeff->GetTabRealVectorColumn("ENERG_HI");
    TVectorD* aeffValues = hduAeff->GetTabRealVectorColumn("AEFF");
    Int_t nBinsAeff = aeffValues->GetNrows();
    // fill an array with the log10(E / GeV) of the bin edges (gLike input)
    Double_t log10EtrueEdges[nBinsAeff+1];
    for (Int_t i=0; i<nBinsAeff; ++i){
        log10EtrueEdges[i] = TMath::Log10(1e3 * lowEtrueEdges[0][i]);
    } 
    // the last bin edge is the last element of the highEtrueEdges
    log10EtrueEdges[nBinsAeff] = TMath::Log10(1e3 * highEtrueEdges[0][nBinsAeff-1]);
    TH1F* hAeff = new TH1F();
    hAeff->SetBins(nBinsAeff, log10EtrueEdges);
    for (Int_t i=0; i<nBinsAeff; ++i){
    	hAeff->SetBinContent(i, aeffValues[0][i]);
    }
    // set the effective area int the IactEventListIrf
    iactDataSet->SetHAeff(hAeff);
    delete hAeff;

    // fetch the migration matrix, set the bias and resolution TGraphs
    TVectorD* eTrue = hduEdisp->GetTabRealVectorColumn("E_TRUE");
    TVectorD* bias = hduEdisp->GetTabRealVectorColumn("BIAS");
    TVectorD* resolution = hduEdisp->GetTabRealVectorColumn("RES");
    TGraph* biasTGraph = new TGraph();
    TGraph* resolutionTGraph = new TGraph();
    Int_t nBinsEdisp = bias->GetNrows();
    // fill the TGraph remembering gLike works with log10(E / GeV)
    for (Int_t i=0; i<nBinsEdisp; ++i){
        biasTGraph->SetPoint(i, TMath::Log10(1e3 * eTrue[0][i]), bias[0][i]);
        resolutionTGraph->SetPoint(i, TMath::Log10(1e3 * eTrue[0][i]), resolution[0][i]);
    }
    iactDataSet->SetGEResoAndBias(resolutionTGraph, biasTGraph);
    delete biasTGraph, resolutionTGraph;

    if (checkInputs){
        cout << "checking that all the events and IRFs have been properly set ..." << endl;
        cout << "-> checking ON events:" << endl;
        iactDataSet->GetOnSample()->Print();
        cout << "-> checking OFF events:" << endl;
        iactDataSet->GetOffSample()->Print();
        cout << "-> checking effective area:" << endl;
        TCanvas* c1 = new TCanvas("c1", "aeff", 800, 600);
        c1->SetLogy();
        TH1F* hAeff = iactDataSet->GetHAeff();
        hAeff->GetXaxis()->SetTitle("log10(E  / GeV)");
        hAeff->GetYaxis()->SetTitle("A_{eff} / cm^{2}");
        hAeff->SetLineWidth(2);
        hAeff->Draw();
        cout << "-> checking resolution and bias:" << endl;
        TCanvas* c2 = new TCanvas("c2", "edisp", 800, 600);
        TGraph* biasTGraph = iactDataSet->GetGEbias();
        TGraph* resolutionTGraph = iactDataSet->GetGEreso();
        biasTGraph->SetLineColor(2);
        biasTGraph->SetLineWidth(2);
        resolutionTGraph->SetLineColor(4);
        resolutionTGraph->SetLineWidth(2);
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(biasTGraph, "bias");
        legend->AddEntry(resolutionTGraph, "resolution");
        biasTGraph->GetXaxis()->SetTitle("log10(E  / GeV)");
        biasTGraph->Draw();
        resolutionTGraph->Draw("same");
        legend->Draw("same");
    }
}
