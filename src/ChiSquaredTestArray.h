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

    /** \brief Fill a pair of given arrays with data that represents a folded light curve for theis chi-squared test.
        \param array_index The index of the element of the periodicity test array, of which a data array is to be created.
        \param phase The output array for pulse phases, which may be used an X-axis of a folded light curve.
        \param count The output array for the number of photons in a phase bins, which may be used an Y-axis of a folded light curve.
    */
    virtual void getPlotData(size_type array_index, std::vector<double> & phase, std::vector<double> & count) const;

    /** \brief Assign axis labels to a pair of given strings, each of which can be used as an X- and Y-axis label, respectively.
        \param x_data The output string that contains the label for X-axis of a plot to display.
        \param y_data The output string that contains the label for Y-axis of a plot to display.
    */
    virtual void getPlotLabel(std::string & x_label, std::string & y_label) const;

    /** \brief Return a plot title that can be used with a return value of getPlotData method.
    */
    virtual std::string getPlotTitle() const;

    /** \brief Return the name of this periodicity test.
    */
    virtual std::string getTestName() const;

    /** \brief Create a statistic viewer for an object of this class.
        \param array_index The index of the element of the periodicity test array, for which a viewer is to be created.
    */
    virtual StatisticViewer getViewer(size_type array_index) const;

  private:
    // The number of phase bins.
    size_type m_num_phase_bins;

    // The container for a folded light curve.
    cont_type m_curve_cont;

    // The number of events filled for each element of this test array.
    std::vector<long> m_num_events;
};

#endif
