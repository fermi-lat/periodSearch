/** \file PowerSpectrumApp.h
    \brief Declaration of PowerSpectrumApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PowerSpectrumApp_h
#define periodSearch_PowerSpectrumApp_h

#include "st_app/StApp.h"

#include "st_stream/StreamFormatter.h"

#include <string>

namespace st_app {
  class AppParGroup;
}

class FourierAnalysis;

class PowerSpectrumApp : public st_app::StApp {
  public:
    PowerSpectrumApp();
    virtual ~PowerSpectrumApp() throw();
    virtual void run();

    virtual void prompt(st_app::AppParGroup & pars);

    const std::string & getDataDir();

  private:
    st_stream::StreamFormatter m_os;
    std::string m_data_dir;
    FourierAnalysis * m_test;
};

#endif
