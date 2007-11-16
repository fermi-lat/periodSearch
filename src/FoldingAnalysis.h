/** \file FoldingAnalysis.h
    \brief Declaration of FoldingAnalysis class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_FoldingAnalysis_h
#define periodSearch_FoldingAnalysis_h

#include <string>
#include <utility>
#include <vector>

#include "periodSearch/PeriodSearch.h"

namespace periodSearch {
  class PeriodicityTest;

  /** \class FoldingAnalysis
      \brief Class for periodicity search based on various statistical tests used to determine frequency of pulsation
             when an approximate frequency is known, using the epoch-folding technique.
  */
  class FoldingAnalysis : public PeriodSearch {
    public:
      /// \brief Container type used to store periodicity tests for trial frequencies.
      typedef std::vector<PeriodicityTest *> TestCont_t;

      /** \brief Construct a periodicity search object using given search criteria and a statistical test object.
          \param center The central value to test.
          \param step The step size to use.
          \param num_trials The number of trials in the test scan range.
          \param epoch The global time offset defining the origin for purposes of this computation.
          \param duration The total time duration (only used if chance probability will be computed).
          \param test The statistical test object to be used for this periodicity search.
      */
      FoldingAnalysis(double center, double step, size_type num_trials, double epoch, double duration,
                      const PeriodicityTest & test);

      virtual ~FoldingAnalysis() {}

      /** \brief Fill given time into periodicity test objects.
          \param evt_time The time of the event.
      */
      virtual void fill(double evt_time);

      /** \brief Compute test statistics for currently filled event times at all trial frequencies.
                 Details depend on the a specific test object given upon construction of this object.
      */
      // TODO: Need to return an array?  Remove the return value if not.
      virtual const std::vector<double> & computeStats();

      /** \brief Return the number of independent trials for this search method.
      */
      virtual size_type numIndepTrials(double min_freq = -1., double max_freq = -1.) const;

      /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
          \param stat The value of the statistic.
      */
      virtual std::pair<double, double> chanceProbOneTrial(double stat) const; 

      /** \brief Return a description of this search.
      */
      virtual std::string getDescription() const;

    protected:
      // The container of statistical trials.
      TestCont_t m_test_cont;

      // Sampling frequency (frequency step).
      double m_step;

      // Temporal origin.
      double m_epoch;

      // Fourier resolution of search.
      double m_fourier_res;

      // Number of events, used to normalize the trials.
      int m_num_events;
  };
}

#endif
