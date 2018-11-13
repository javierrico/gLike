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
    //Don't forget that the shared object must have been loaded
    gSystem->Load("lib/libgLike.so");

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
