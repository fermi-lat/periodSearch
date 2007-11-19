/** \file PeriodicityTestApp.cxx
    \brief Implmentation of PeriodicityTestApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PeriodicityTestApp.h"

#include "st_app/AppParGroup.h"

static const std::string s_cvs_id = "$Name:  $";

PeriodicityTestApp::PeriodicityTestApp(): m_os("PeriodicityTestApp", "", 2) {
  setName("gtptest");
  setVersion(s_cvs_id);
}

PeriodicityTestApp::~PeriodicityTestApp() throw() {}

void PeriodicityTestApp::run() {
  m_os.setMethod("run()");
  st_app::AppParGroup & pars(getParGroup("gtptest"));

  // Prompt and save.
  pars.Prompt();
  pars.Save();
}
