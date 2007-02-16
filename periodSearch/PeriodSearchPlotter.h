/** \file PeriodSearchPlotter.h
    \brief Declaration of PeriodSearchPlotter class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodSearchPlotter_h
#define periodSearch_PeriodSearchPlotter_h

#include <string>

namespace periodSearch {

  class PeriodSearch;

  /** \class PeriodSearchPlotter
      \brief Base class for various statistical tests used to determine frequency of pulsation
             when an approximate frequency is known.
  */
  class PeriodSearchPlotter {
    public:
      /** \brief Create a plotter object for plotting period search objects.
      */
      PeriodSearchPlotter();

      virtual ~PeriodSearchPlotter();

      /** \brief Display plot of statistic as a function of frequency for entire frequency range.
          \param title The title to display on the plot. (Purely cosmetic.)
          \param freq_unit The units to display on the x axis. (Purely cosmetic.)
      */
      virtual void plot(const PeriodSearch & search, const std::string & title, const std::string & freq_unit) const;

      /** \brief Display plot of statistics as a function of frequency over the given range.
          \param title The title to display on the plot. (Purely cosmetic.)
          \param freq_unit The units to display on the x axis. (Purely cosmetic.)
          \param min_freq The minimum frequency in the range (if negative, do not constrain minimum frequency.)
          \param max_freq The maximum frequency in the range (if negative, do not constrain maximum frequency.)
      */
      virtual void plotRange(const PeriodSearch & search, const std::string & title, const std::string & freq_unit,
        double min_freq = -1., double max_freq = -1.) const;
  };

}

#endif
