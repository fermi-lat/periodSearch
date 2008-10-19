/** \file PeriodSearchApp.h
    \brief Declaration of PeriodSearchApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodSearchApp_h
#define periodSearch_PeriodSearchApp_h

#include "pulsarDb/PulsarToolApp.h"

#include "st_app/AppParGroup.h"

#include "st_stream/StreamFormatter.h"

class PeriodSearchApp : public pulsarDb::PulsarToolApp {
  public:
    PeriodSearchApp();
    virtual ~PeriodSearchApp() throw();
    virtual void runApp();

  private:
    st_stream::StreamFormatter m_os;

    virtual void prompt(st_app::AppParGroup & pars);
};

#endif
