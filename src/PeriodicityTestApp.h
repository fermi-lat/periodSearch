/** \file PeriodicityTestApp.h
    \brief Declaration of PeriodicityTestApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodicityTestApp_h
#define periodSearch_PeriodicityTestApp_h

#include "st_app/StApp.h"

#include "st_stream/StreamFormatter.h"

/** \class PeriodicityTestApp
    \brief Main application class for various types of statistical test for periodicity.
*/
class PeriodicityTestApp : public st_app::StApp {
  public:
    /// \brief Construct a PeriodicityTestApp object.
    PeriodicityTestApp();

    /// \brief Destruct this PeriodicityTestApp object.
    virtual ~PeriodicityTestApp() throw();

    /// \brief Run the application.
    virtual void run();

  private:
    st_stream::StreamFormatter m_os;
};

#endif
