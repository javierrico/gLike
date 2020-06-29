///////////////////////////////////////////////////////////////////////////
//
// rootlogon.C
// ===========
//
// Loading libgLike.so (.dylib) when using root in this directory.
//
///////////////////////////////////////////////////////////////////////////

void rootlogon()
{
    TString gLikeBuildDir = gSystem->ExpandPathName("$GLIKESYS");
    cout << "\033[34m\033[1m" << "Searching gLike in " << gLikeBuildDir << "/lib" << "\033[0m" << endl ;
    TString libFileNameUnix = "libgLike.so";
    TString libFileNameMacOs = "libgLike.dylib";
    // search if unix library exists, empty string if not
    TString libUnixPath =  gSystem->FindFile(gLikeBuildDir + "/lib/", libFileNameUnix);
    // search if MacOs library exists, empty string if not
    TString libMacOsPath =  gSystem->FindFile(gLikeBuildDir + "/lib/", libFileNameMacOs);
    
    if (!libUnixPath.IsNull()) {
        cout << "\033[33m\033[1m" << "Loading " << libUnixPath << " \033[0m" << flush;
        gSystem->Load(libUnixPath);
        cout << "\033[33m\033[1m" << "done." << endl;   
    } else if(!libMacOsPath.IsNull()) {
        cout << "\033[33m\033[1m" << "Loading " << libMacOsPath << " \033[0m" << flush;
        gSystem->Load(libMacOsPath);
        cout << "\033[33m\033[1m" << "done." << endl;  
    } else {
        cout << "\033[33m\033[1m" << "Error: gLike library not found" << endl;
        return;
    }
    
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    cout << "\033[32m" << "Welcome to the gLike Root environment." << "\033[0m" << endl;
    cout << endl;
}
