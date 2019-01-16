///////////////////////////////////////////////////////////////////////////
//
// rootlogon.C
// ===========
//
// Loading libgLike.so when using root in this directory.
//
///////////////////////////////////////////////////////////////////////////

bool load(TString &dir, const TString &libfilename)
{
    cout << "\033[33m\033[1m" << "Loading " << dir << libfilename << " \033[0m" << flush;

    if (dir.IsNull()) dir = "../lib/";

    if (gSystem->Load(dir+libfilename)!=0) {
      cout << "\033[33m\033[1m" << "Error" << endl;
      return false;
    }
    else {
      cout << "\033[33m\033[1m" << "done." << endl;
      return true;
    }
}

void rootlogon()
{
    cout << endl;

    const TString libfilename = "libgLike.so"; // default
    const Bool_t  fileexist   = !gSystem->AccessPathName(libfilename.Data(), kFileExists);
    TString       gLikeDir    = fileexist ? "" : gSystem->Getenv("GLIKESYS");
    
    TString libdir = gLikeDir + "/lib/";
    TString outdir = gLikeDir + "/out/";
    gSystem->AddDynamicPath(outdir);

    if (!libdir.IsNull()) {
     cout << "\033[34m\033[1m" << "Searching gLike in " << libdir << "\033[0m" << endl ;
      if (!libdir.EndsWith("/")) libdir += "/";
    }

    if (!load(libdir,libfilename)) return;

    gInterpreter->AddIncludePath(gLikeDir+"/include");
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    cout << "\033[32m" << "Welcome to the gLike Root environment." << "\033[0m" << endl;
    cout << endl;
}
