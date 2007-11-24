/** \file ChiSquaredTestArray.h
    \brief Declaration of ChiSquaredTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_ChiSquaredTestArray_h
#define periodSearch_ChiSquaredTestArray_h

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

    /** \brief Compute an S-value of this chi-squared test for pulse phases currently filled in this object.
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
    virtual size_type size() const;

    /** \brief Display a folded light curve used in this chi-squared test.
        \param title The plot title.
        \param array_index The index of the element of the periodicity test array, of which a plot is to be created.
    */
    virtual void plot(const std::string & title, size_type array_index = 0) const;

  private:
    // The number of phase bins.
    size_type m_num_phase_bins;

    // The container for a folded light curve.
    cont_type m_curve_cont;

    // The number of events filled for each element of this test array.
    std::vector<long> m_num_events;
};

#endif
