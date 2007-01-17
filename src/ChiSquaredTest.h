/** \file ChiSquaredTest.h
    \brief Declaration of ChiSquaredTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_ChiSquaredTest_h
#define periodSearch_ChiSquaredTest_h

#include <complex>
#include <utility>
#include <vector>

#include "periodSearch/PeriodTest.h"

/** \class ChiSquaredTest
    \brief PeriodTest subclass which uses a Chi^2 statistic for its search/test.
*/
class ChiSquaredTest : public periodSearch::PeriodTest {
  public:
    /// \brief Container type used to store trial histograms.
    typedef periodSearch::PeriodTest::HistCont_t HistCont_t;

    /** \brief Construct a test object using given trial information.
        \param center The central value to test.
        \param step The step size to use.
        \param num_trials The number of trials in the test scan range.
        \param epoch The global time offset defining the origin for purposes of this computation.
        \param num_phase_bins The number of bins used in each trial. Depending on the specific test, this
               may be an upper limit on the actual number of bins used.
        \param duration The total time duration (only used if chance probability will be computed).
    */
    ChiSquaredTest(double center, double step, size_type num_trials, double epoch, size_type num_phase_bins, double duration);

    /** \brief For a given phase, compute a single statistical trial.
        \param phase The phase.
        \param trial A single trial array, in this case a simple histogram. Only the real part of the
                     trial array is used.
    */
    virtual void fillOneTrial(double phase, std::vector<std::complex<double> > & trial) const;

    /** \brief Compute standard Chi^2 statistics for each trial histogram. The Chi^2 value is
               not reduced.
    */
    virtual const std::vector<double> & computeStats();

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProbOneTrial(double stat) const; 

  private:
    size_type m_num_phase_bins;
};

#endif
