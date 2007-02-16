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

      /** \brief Find the frequency for which the statistic is maximized in a given frequency range. Return the
                 frequency and the value of the statistic, as a pair.
          \param min_freq The minimum frequency in the range.
          \param max_freq The maximum frequency in the range.
      */
      virtual std::pair<double, double> findMax(double min_freq = -1., double max_freq = -1.) const;

      /** \brief Return the number of independent trials for this search method.
      */
      virtual size_type numIndepTrials() const = 0;

      /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
          \param stat The value of the statistic.
      */
      virtual std::pair<double, double> chanceProbOneTrial(double stat) const = 0;

      /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
          \param stat The value of the statistic.
      */
      virtual std::pair<double, double> chanceProb(double stat) const;

      /** \brief Write statistical data as a function of frequency to the given stream.
          \param os The stream.
      */
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

      //* \brief Get frequency data associated with this search object.
      const cont_type getFreq() const;
      
      //* \brief Get spectral data associated with this search object.
      const cont_type getSpec() const;
      
      /* \brief Compute the probability that an event occurs at least once in N statistically
                independent trials, given that the probability of the event occurring in a single trial is p.
         \param prob_one_trial The probability p of the event occuring in one trial.
         \param num_indep_trial The number N of statistically independent trials.
      */
      static double chanceProbMultiTrial(double prob_one_trial, size_type num_indep_trial);

      /** \brief Write data over a specified frequency range as a function of frequency to the given stream.
          \param os The stream.
          \param min_freq The minimum frequency in the range.
          \param max_freq The maximum frequency in the range.
      */
      virtual st_stream::OStream & writeRange(st_stream::OStream & os, double min_freq = -1., double max_freq = -1.) const;

      /** \brief Given a frequency range, determine the indices of (inclusive) lower and upper bounds.
          \param min_freq The minimum frequency.
          \param max_freq The maximum frequency.
      */
      std::pair<size_type, size_type> getRangeIndex(double min_freq, double max_freq) const;

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
