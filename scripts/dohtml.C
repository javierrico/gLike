#include "THtml.h"
///////////////////////////////////////////////////////////////////////////
//
// dohtml.C
// ========
//
// This is a service macro used to create the html-documentation from source code (THtml)
//
// Add here all directories in which files are stored from which the documentation
// should be created an add all macros which should be converted to HTML.
//
///////////////////////////////////////////////////////////////////////////

void dohtml()
{
    TString gLikeBuildDir = gSystem->ExpandPathName("$GLIKESYS");
    //Don't forget that the shared object must have been loaded
    TString libFileNameUnix = "libgLike.so";
    TString libFileNameMacOs = "libgLike.dylib";
    // search if unix library exists, empty string if it does not
    TString libUnixPath =  gSystem->FindFile(gLikeBuildDir + "/lib/", libFileNameUnix);
    // search if MacOs library exists, empty string if it does not
    TString libMacOsPath =  gSystem->FindFile(gLikeBuildDir + "/lib/", libFileNameMacOs);
    
    if (!libUnixPath.IsNull()) {
        gSystem->Load(libUnixPath);
    } else if(!libMacOsPath.IsNull()) {
        gSystem->Load(libMacOsPath);
    } else {
        return;
    }
    //Do not print 'Info' messages from the root system such like TCanvas::Print
    gErrorIgnoreLevel=kWarning;

    //Create the html document class
    THtml html;

    html.SetSourceDir(".:include:");
    html.SetOutputDir("htmldoc");
    html.SetProductName("gLike");
    html.SetAuthorTag("J. Rico & J. Aleksic");
    html.SetCopyrightTag("See Lkl.h for terms of usage");

    html.MakeAll(kTRUE);
}
