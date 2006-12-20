/** \file PeriodSearch.h
    \brief Declaration of PeriodSearch class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodSearch_h
#define periodSearch_PeriodSearch_h

#include <complex>
#include <string>
#include <utility>
#include <vector>

#include "st_stream/Stream.h"

namespace periodSearch {
  /** \class PeriodSearch
      \brief Base class for various statistical tests used to determine frequency of pulsation
             when an approximate frequency is known.
  */
  class PeriodSearch {
    public:
      typedef std::vector<double> cont_type;
      typedef cont_type::size_type size_type;

      PeriodSearch(size_type num_bins);

      virtual ~PeriodSearch() {}

      /** \brief Fill given time into histograms.
          \param evt_time The time of the event.
      */
      virtual void fill(double evt_time) = 0;

      /** \brief Use the trials as currently filled by data to compute statistics for this test. Details
                 depend on the specific test being performed in the subclass.
      */
      virtual const std::vector<double> & computeStats() = 0;

      /** \brief Display plot of statistics.
          \param title The title to display on the plot. (Purely cosmetic.)
          \param freq_unit The units to display on the x axis. (Purely cosmetic.)
      */
      virtual void plotStats(const std::string & title, const std::string & freq_unit) const = 0;

      /** \brief Find the frequency for which the statistic is maximized. Return the frequency and the value of
                 the statistic, as a pair.
      */
      virtual std::pair<double, double> findMax() const = 0;

      /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
          \param stat The value of the statistic.
      */
      virtual std::pair<double, double> chanceProb(double stat) const = 0;

      /** \brief Write the test values to the given stream.
          \param os The stream.
      */
      virtual st_stream::OStream & write(st_stream::OStream & os) const = 0;

    protected:
      // For convenience, compute once the value of 2 * pi.
      static const double s_2pi;

      // The frequencies forming the domain of the search/test.
      cont_type m_freq;

      // The statistical measure of the validity of each trial.
      cont_type m_spec;
  };

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearch & test);

}

#endif
