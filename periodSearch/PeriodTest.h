/** \file PeriodTest.h
    \brief Declaration of PeriodTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodTest_h
#define periodSearch_PeriodTest_h

#include <complex>
#include <string>
#include <utility>
#include <vector>

#include "st_stream/Stream.h"

namespace periodSearch {
  /** \class PeriodTest
      \brief Base class for various statistical tests used to determine frequency of pulsation
             when an approximate frequency is known.
  */
  class PeriodTest {
    public:
      /// \brief Container type used to store histograms and/or Fourier transformed data.
      typedef std::vector<std::vector<std::complex<double> > > HistCont_t;

      virtual ~PeriodTest() {}

      /** \brief Fill given time into histograms.
          \param evt_time The time of the event.
      */
      virtual void fill(double evt_time);

      /** \brief For a given phase, compute a single statistical trial.
          \param phase The phase.
          \param trial A single trial array, whose exact interpretation depends on the type of test
                       being performed in the subclass.
      */
      virtual void fillOneTrial(double phase, std::vector<std::complex<double> > & trial) const = 0;

      /** \brief Use the trials as currently filled by data to compute statistics for this test. Details
                 depend on the specific test being performed in the subclass.
      */
      virtual const std::vector<double> & computeStats() = 0;

      /** \brief Display plot of statistics.
          \param title The title to display on the plot. (Purely cosmetic.)
          \param freq_unit The units to display on the x axis. (Purely cosmetic.)
      */
      virtual void plotStats(const std::string & title, const std::string & freq_unit) const;

      /** \brief Find the frequency for which the statistic is maximized. Return the frequency and the value of
                 the statistic, as a pair.
      */
      virtual std::pair<double, double> findMax() const;

      /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
          \param stat The value of the statistic.
      */
      virtual std::pair<double, double> chanceProb(double stat) const = 0;

      /** \brief Write the test values to the given stream.
          \param os The stream.
      */
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

    protected:
      /** \brief Construct a test object using given trial information.
          \param center The central value to test.
          \param step The step size to use.
          \param num_trials The number of trials in the test scan range.
          \param epoch The global time offset defining the origin for purposes of this computation.
          \param num_bins The number of bins used in each trial. Depending on the specific test, this
                          may be an upper limit on the actual number of bins used.
          \param duration The total time duration (only used if chance probability will be computed).
      */
      PeriodTest(double center, double step, long num_trials, double epoch, int num_bins, double duration);

      // For convenience, compute once the value of 2 * pi.
      static const double s_2pi;

      // The container of statistical trials.
      HistCont_t m_trial_hist;

      // The frequencies forming the domain of the search/test.
      std::vector<double> m_freqs;

      // The statistical measure of the validity of each trial.
      std::vector<double> m_stats;

      // Center of scan.
      double m_center;

      // Temporal origin.
      double m_epoch;

      // Duration of data.
      double m_duration;

      // Number of bins used internally for histograms and/or Fourier transformed data.
      int m_num_bins;

      // Number of events, used to normalize the trials.
      int m_num_events;

      // Number of independent trials.
      int m_num_indep_trials;
  };

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodTest & test);

}

#endif
