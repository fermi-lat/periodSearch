/** \file gtptest.cxx
    \brief Periodicity test tool that computes statistical significance of periodicity in pulse phases assigned to events.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PeriodicityTestApp.h"

#include "st_app/StAppFactory.h"

st_app::StAppFactory<PeriodicityTestApp> g_factory("gtptest");
