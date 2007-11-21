/** \file ChiSquaredTestArray.h
    \brief Declaration of ChiSquaredTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_ChiSquaredTestArray_h
#define periodSearch_ChiSquaredTestArray_h

#include <complex>
#include <utility>
#include <vector>

#include "PeriodicityTestArray.h"

/** \class ChiSquaredTestArray
    \brief PeriodicityTestArray subclass which uses a Chi^2 statistic for its test.
*/
class ChiSquaredTestArray : public PeriodicityTestArray {
  public:
    typedef std::vector<long> data_type;
    typedef std::vector<data_type> cont_type;

    ChiSquaredTestArray(size_type array_size, data_type::size_type num_phase_bins);

    virtual ~ChiSquaredTestArray() {}

    /** \brief Fill a given pulse phase into this periodicity test object.
        \param phase The pulse phase.
    */
    virtual void fill(double phase, size_type array_index = 0);

    /** \brief Compute a test statistic for pulse phases currently filled in this object. Details
               depend on the specific test being performed in the subclass.
    */
    virtual double testStat(size_type array_index = 0) const;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProb(double stat) const;

    /** \brief Return a description of this search.
    */
    virtual std::string getDescription() const;

    /** \brief Return the size of this periodicity test array.
    */
    virtual size_type size() const = 0;

  private:
    // The number of phase bins.
    size_type m_num_phase_bins;

    // The container for a folded light curve.
    cont_type m_curve_cont;

    // The number of events filled.
    long m_num_events;
};

#endif
