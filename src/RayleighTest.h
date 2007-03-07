/** \file RayleighTest.h
    \brief Declaration of RayleighTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_RayleighTest_h
#define periodSearch_RayleighTest_h

#include "Z2nTest.h"

#include <sstream>

/** \class RayleighTest
    \brief PeriodTest subclass which uses Rayleigh statistic (=== Z2n with n == 1) for its search/test.
*/
class RayleighTest : public Z2nTest {
  public:
    /** \brief Construct a test object using given trial information.
        \param center The central value to test.
        \param step The step size to use.
        \param num_trials The number of trials in the test scan range.
        \param epoch The global time offset defining the origin for purposes of this computation.
        \param duration The total time duration (only used if chance probability will be computed).
    */
    RayleighTest(double center, double step, long num_trials, double epoch, double duration):
      Z2nTest(center, step, num_trials, epoch, 1, duration) {}

    /** \brief Return a description of this search.
    */
    virtual std::string getDescription() const {
      std::ostringstream os;
      os << PeriodTest::getDescription() << "\n" <<
        "Type of test: Rayleigh Test (Z2n with n = 1 harmonic)\n" << 
        "Probability distribution: Chi-squared, " << 2 << " degrees of freedom";
      return os.str();
    }
};

#endif
