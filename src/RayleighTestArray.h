/** \file RayleighTestArray.h
    \brief Declaration of RayleighTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_RayleighTestArray_h
#define periodSearch_RayleighTestArray_h

#include "Z2nTestArray.h"

#include <sstream>

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
      os << "Type of test: Rayleigh Test (Z2n with n = 1 harmonic)\n" << 
        "Probability distribution: Chi-squared, " << 2 << " degrees of freedom";
      return os.str();
    }
};

#endif
