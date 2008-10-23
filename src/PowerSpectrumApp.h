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

/** \class PowerSpectrumApp
    \brief Main application class for a pulsation search by a power spectrum.
*/
class PowerSpectrumApp : public pulsarDb::PulsarToolApp {
  public:
    /// \brief Construct a PowerSpectrumApp object.
    PowerSpectrumApp();

    /// \brief Destruct this PowerSpectrumApp object.
    virtual ~PowerSpectrumApp() throw();

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
