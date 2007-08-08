/** \file PowerSpectrumApp.h
    \brief Declaration of PowerSpectrumApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PowerSpectrumApp_h
#define periodSearch_PowerSpectrumApp_h

#include "pulsarDb/PulsarToolApp.h"

#include "st_stream/StreamFormatter.h"

namespace st_app {
  class AppParGroup;
}

class PowerSpectrumApp : public pulsarDb::PulsarToolApp {
  public:
    PowerSpectrumApp();
    virtual ~PowerSpectrumApp() throw();
    virtual void run();

    virtual void prompt(st_app::AppParGroup & pars);

  private:
    st_stream::StreamFormatter m_os;
};

#endif
