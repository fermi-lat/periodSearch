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

    /** \brief Construct an array object of the chi-squared test.
        \param array_size The size of this test array.
        \param num_phase_bins The number of phase bins for the chi-squared test.
    */
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

    /** \brief Return a description of this test array.
    */
    virtual std::string getDescription() const;

    /** \brief Return the size of this periodicity test array.
    */
    virtual size_type size() const;

    /** \brief Return a pair of data arrays that represents a folded light curve for theis chi-squared test.
               The first array of the pair is an array of pulse phase values, and the second an array of the number
               of photons in each phase bin.
        \param array_index The index of the element of the periodicity test array, of which a data array is to be created.
    */
    virtual std::pair<std::vector<double>, std::vector<double> > getPlotData(size_type array_index = 0) const;

    /** \brief Return a pair of axis labels, each of which can be used as an X- and Y-axis label, respectively.
    */
    virtual std::pair<std::string, std::string> getPlotLabel() const;

    /** \brief Return a plot title that can be used with a return value of getPlotData method.
    */
    virtual std::string getPlotTitle() const;

  private:
    // The number of phase bins.
    size_type m_num_phase_bins;

    // The container for a folded light curve.
    cont_type m_curve_cont;

    // The number of events filled for each element of this test array.
    std::vector<long> m_num_events;
};

#endif
