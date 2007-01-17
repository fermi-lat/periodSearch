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

#include "periodSearch/PeriodSearch.h"

#include "st_stream/Stream.h"

namespace periodSearch {
  /** \class PeriodTest
      \brief Base class for various statistical tests used to determine frequency of pulsation
             when an approximate frequency is known.
  */
  // TODO Rename to PeriodicityTest?
  class PeriodTest : public PeriodSearch {
    public:
      /// \brief Container type used to store histograms and/or Fourier transformed data.
      typedef std::vector<std::vector<std::complex<double> > > HistCont_t;

      virtual ~PeriodTest() {}

      /** \brief Fill given time into histograms.
          \param evt_time The time of the event.
      */
      virtual void fill(double evt_time);

      /** \brief Return the number of independent trials for this search method.
      */
      virtual size_type numIndepTrials() const;

      /** \brief For a given phase, compute a single statistical trial.
          \param phase The phase.
          \param trial A single trial array, whose exact interpretation depends on the type of test
                       being performed in the subclass.
      */
      virtual void fillOneTrial(double phase, std::vector<std::complex<double> > & trial) const = 0;

    protected:
      /** \brief Construct a test object using given trial information.
          \param center The central value to test.
          \param step The step size to use.
          \param num_trials The number of trials in the test scan range.
          \param epoch The global time offset defining the origin for purposes of this computation.
          \param array_size The size of array used for each trial. Depending on the specific test, this
                 may be an upper limit on the actual number of elements used.
          \param duration The total time duration (only used if chance probability will be computed).
      */
      PeriodTest(double center, double step, size_type num_trials, double epoch, size_type array_size, double duration);

      // The container of statistical trials.
      HistCont_t m_trial_hist;

      // Center of scan.
      double m_center;

      // Temporal origin.
      double m_epoch;

      // Duration of data.
      double m_duration;

      // Number of events, used to normalize the trials.
      int m_num_events;

      // Number of independent trials.
      int m_num_indep_trials;
  };

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodTest & test);

}

#endif
