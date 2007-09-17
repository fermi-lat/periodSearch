/** \file gtpspec.cxx
    \brief Period search tool that uses FFT technique to compute power spectrum density.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PowerSpectrumApp.h"

#include "st_app/StAppFactory.h"

#include "timeSystem/EventTimeHandler.h"
#include "timeSystem/GlastTimeHandler.h"

// List supported mission(s).
timeSystem::EventTimeHandlerFactory<timeSystem::GlastTimeHandler> glast_handler;

st_app::StAppFactory<PowerSpectrumApp> g_factory("gtpspec");
