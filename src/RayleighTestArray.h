/** \file RayleighTestArray.h
    \brief Declaration of RayleighTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_RayleighTestArray_h
#define periodSearch_RayleighTestArray_h

#include <string>
#include <sstream>

#include "Z2nTestArray.h"

/** \class RayleighTestArray
    \brief PeriodTest subclass which uses Rayleigh statistic (which is equal to Z2n statistic with n == 1) for its search/test.
*/
class RayleighTestArray : public Z2nTestArray {
  public:
    RayleighTestArray(size_type array_size): Z2nTestArray(array_size, 1) {}

    /** \brief Return a description of this test array.
    */
    virtual std::string getDescription() const {
      std::ostringstream os;
      os << "Type of test: " << getTestName() << " (Z2n with n = 1 harmonic)" << std::endl
         << "Probability distribution: Chi-squared, 2 degrees of freedom";
      return os.str();
    }

    /** \brief Return the name of this periodicity test.
    */
    virtual std::string getTestName() const { return "Rayleigh Test"; }
};

#endif
