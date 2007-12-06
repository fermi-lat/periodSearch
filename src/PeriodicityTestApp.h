/** \file PeriodicityTestApp.h
    \brief Declaration of PeriodicityTestApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodicityTestApp_h
#define periodSearch_PeriodicityTestApp_h

#include "st_app/StApp.h"

#include "st_stream/StreamFormatter.h"

namespace st_app {
  class AppParGroup;
}

class PeriodicityTestArray;

class PeriodicityTestApp : public st_app::StApp {
  public:
    PeriodicityTestApp();
    virtual ~PeriodicityTestApp() throw();
    virtual void run();

  private:
    st_stream::StreamFormatter m_os;
};

#endif
