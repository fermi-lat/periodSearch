/** \file ChiSquaredTestArray.h
    \brief Declaration of ChiSquaredTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_ChiSquaredTestArray_h
#define periodSearch_ChiSquaredTestArray_h

#include <string>
#include <utility>
#include <vector>

#include "PeriodicityTestArray.h"
#include "StatisticViewer.h"

/** \class ChiSquaredTestArray
    \brief PeriodicityTestArray subclass which uses a Chi^2 statistic for its test.
*/
class ChiSquaredTestArray : public PeriodicityTestArray {
  public:
    typedef std::vector<long> data_type;
    typedef std::vector<data_type> cont_type;

    /** \brief Construct an array object of the chi-squared test.
        \param array_size The size of this test array.
        \param num_phase_bins The number of phase bins for the chi-squared test.
    */
    ChiSquaredTestArray(size_type array_size, data_type::size_type num_phase_bins);

    virtual ~ChiSquaredTestArray() {}

    /** \brief Fill a given pulse phase into this periodicity test object.
        \param array_index The index of the element of the periodicity test array, to which a given phase is filled.
        \param phase The pulse phase to fill.
    */
    virtual void fill(size_type array_index, double phase);

    /** \brief Compute an S-value of this chi-squared test for pulse phases currently filled in this object.
        \param array_index The index of the element of the periodicity test array, of which an S-value is to be computed.
    */
    virtual double testStat(size_type array_index) const;

    /** \brief Compute an S-value of this chi-squared test for pulse phases currently filled in this object, and set
               the folded light curve to the internal statistic viewer.
        \param array_index The index of the element of the periodicity test array, of which an S-value is to be computed.
    */
    virtual void updateViewer(size_type array_index);

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProb(double stat) const;

    /** \brief Return the size of this periodicity test array.
    */
    virtual size_type size() const;

  private:
    // The number of phase bins.
    size_type m_num_phase_bins;

    // The container for a folded light curve.
    cont_type m_curve_cont;

    // The number of events filled for each element of this test array.
    std::vector<long> m_num_events;
};

#endif
