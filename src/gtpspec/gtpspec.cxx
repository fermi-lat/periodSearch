/** \file gtpspec.cxx
    \brief Period search tool that uses FFT technique to compute power spectrum density.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PowerSpectrumApp.h"
#include "st_app/StAppFactory.h"

st_app::StAppFactory<PowerSpectrumApp> g_factory("gtpspec");
