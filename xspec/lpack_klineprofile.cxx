//
// Code auto-generated by initpackage (XSPEC12 local model package 
// code generator).  Do not edit
// Package: klineprofile Created :
// Initializer: klineprofile.cxx

#include <XSUser/UserInterface/xstcl.h>

#include  "klineprofileFunctionMap.h"
#include  <XSFunctions/Utilities/XSModelFunction.h>

extern "C" int Klineprofile_Init(Tcl_Interp* tclInterp);
extern "C" int Klineprofile_SafeInit(Tcl_Interp* tclInterp);

int Klineprofile_Init(Tcl_Interp* tclInterp)
{
        return Klineprofile_SafeInit(tclInterp);
}

int Klineprofile_SafeInit(Tcl_Interp* tclInterp)
{

        char PACKAGE[] = "klineprofile";
        char VERSION[] = "1.0";
        Tcl_PkgProvide(tclInterp, PACKAGE, VERSION);
        createklineprofileFunctionMap();
        XSModelFunction::updateComponentList
              ("/home/heasoft/datasets/lmodel.dat");
        return TCL_OK;

}
