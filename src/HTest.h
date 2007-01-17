/** \file HTest.h
    \brief Declaration of HTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_HTest_h
#define periodSearch_HTest_h

#include <utility>

#include "periodSearch/PeriodTest.h"

/** \class HTest
    \brief PeriodTest subclass which uses H statistic for its search/test.
*/
class HTest : public periodSearch::PeriodTest {
  public:
    /// \brief Container type used to store trials.
    typedef periodSearch::PeriodTest::HistCont_t HistCont_t;

    /** \brief Construct a test object using given trial information.
        \param center The central value to test.
        \param step The step size to use.
        \param num_trials The number of trials in the test scan range.
        \param epoch The global time offset defining the origin for purposes of this computation.
        \param max_harmonics The maximum number of harmonics used.
        \param duration The total time duration (only used if chance probability will be computed).
    */
    HTest(double center, double step, size_type num_trials, double epoch, size_type max_harmonics, double duration);

    /** \brief For a given phase, compute a single statistical trial.
        \param phase The phase.
        \param trial A single trial array, in this case a Fourier transform.
    */
    virtual void fillOneTrial(double phase, std::vector<std::complex<double> > & trial) const;

    /** \brief Compute H statistic (Fourier power) for each trial histogram.
    */
    virtual const std::vector<double> & computeStats();

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProbOneTrial(double stat) const; 

  private:
    size_type m_max_harm;
};

#endif
