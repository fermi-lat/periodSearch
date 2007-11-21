/** \file PeriodicityTestArray.h
    \brief Declaration of PeriodicityTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodicityTestArray_h
#define periodSearch_PeriodicityTestArray_h

#include <complex>
#include <string>
#include <utility>
#include <vector>

/** \class PeriodicityTestArray
    \brief Base class to represent an array of various statistical tests used to evaluate significance of pulsation.
*/
class PeriodicityTestArray {
  public:
    // TODO: Is there any other way to do this?
    typedef int size_type;

    virtual ~PeriodicityTestArray() {}

    /** \brief Fill a given pulse phase into a given element of the array of periodicity test.
        \param phase The pulse phase to fill.
        \param array_index The index of the element of the periodicity test array, to which a given phase is filled.
    */
    virtual void fill(double phase, size_type array_index = 0) = 0;

    /** \brief Compute a test statistic for pulse phases currently filled in this object. Details
               depend on the specific test being performed in the subclass.
        \param array_index The index of the element of the periodicity test array, of which a test statistic is to be computed.
    */
    virtual double testStat(size_type array_index = 0) const = 0;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProb(double stat) const = 0;

    /** \brief Return a description of this periodicity test array.
    */
    virtual std::string getDescription() const = 0;

    /** \brief Return the size of this periodicity test array.
    */
    virtual size_type size() const = 0;

  protected:
    /** \brief Construct a periodicity test object.
    */
    PeriodicityTestArray() {};
};

#endif
