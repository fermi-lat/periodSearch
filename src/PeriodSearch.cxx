/** \file PeriodSearch.cxx
    \brief Implementation of PeriodSearch class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include "periodSearch/PeriodSearch.h"

#include <cmath>

namespace periodSearch {

  const double PeriodSearch::s_2pi = 2. * 4. * atan(1.0);

  PeriodSearch::PeriodSearch(size_type num_bins): m_freq(num_bins), m_spec(num_bins) {}

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearch & test) {
    return test.write(os);
  }

}
