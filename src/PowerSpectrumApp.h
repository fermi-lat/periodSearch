/** \file PowerSpectrumApp.h
    \brief Declaration of PowerSpectrumApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PowerSpectrumApp_h
#define periodSearch_PowerSpectrumApp_h

#include "pulsarDb/PulsarToolApp.h"

#include "st_app/AppParGroup.h"

#include "st_stream/StreamFormatter.h"

class PowerSpectrumApp : public pulsarDb::PulsarToolApp {
  public:
    PowerSpectrumApp();
    virtual ~PowerSpectrumApp() throw();
    virtual void run();

  private:
    st_stream::StreamFormatter m_os;

    virtual void prompt(st_app::AppParGroup & pars);
};

#endif
