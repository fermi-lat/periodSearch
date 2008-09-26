/** \file gtpspec.cxx
    \brief Period search tool that uses FFT technique to compute power spectrum density.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PowerSpectrumApp.h"

#include "st_app/StAppFactory.h"

#include "timeSystem/EventTimeHandler.h"
#include "timeSystem/GlastTimeHandler.h"

// List supported event file format(s).
timeSystem::EventTimeHandlerFactory<timeSystem::GlastScTimeHandler> glast_handler;

st_app::StAppFactory<PowerSpectrumApp> g_factory("gtpspec");
