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

/** \class PeriodSearchApp
    \brief Main application class for a periodicity search by a statistical test.
*/
class PeriodSearchApp : public pulsarDb::PulsarToolApp {
  public:
    /// \brief Construct a PeriodSearchApp object.
    PeriodSearchApp();

    /// \brief Destruct this PeriodSearchApp object.
    virtual ~PeriodSearchApp() throw();

    /// \brief Run the application.
    virtual void runApp();

  private:
    st_stream::StreamFormatter m_os;

    /** \brief Ask a user to input parameter values.
        \param pars Parameter set to be asked.
    */
    virtual void prompt(st_app::AppParGroup & pars);
};

#endif
