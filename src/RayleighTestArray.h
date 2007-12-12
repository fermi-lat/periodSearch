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
    RayleighTestArray(size_type array_size): Z2nTestArray(array_size, 1) {
      // Set description of this statistical test.
      std::ostringstream os_dist;
      os_dist << "Chi-squared, 2 degrees of freedom";
      setDescription("Rayleigh Test", "", os_dist.str());
    }
};

#endif
