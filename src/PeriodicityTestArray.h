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

#include "StatisticViewer.h"

/** \class PeriodicityTestArray
    \brief Base class to represent an array of various statistical tests used to evaluate significance of pulsation.
*/
class PeriodicityTestArray {
  public:
    // TODO: Is there any other way to do this?
    typedef int size_type;

    virtual ~PeriodicityTestArray() {}

    /** \brief Fill a given pulse phase into a given element of the array of periodicity test.
        \param array_index The index of the element of the periodicity test array, to which a given phase is filled.
        \param phase The pulse phase to fill.
    */
    virtual void fill(size_type array_index, double phase) = 0;

    /** \brief Compute a test statistic for pulse phases currently filled in this object. Details
               depend on the specific test being performed in the subclass.
        \param array_index The index of the element of the periodicity test array, of which a test statistic is to be computed.
    */
    virtual double testStat(size_type array_index) const = 0;
    // TODO: Rename the above to computeStat for consistency.

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProb(double stat) const = 0;
    // TODO: Rename the above to computeChanceProb for consistency.

    /** \brief Return a description of this periodicity test array.
    */
    virtual std::string getDescription() const = 0;

    /** \brief Return the size of this periodicity test array.
    */
    virtual size_type size() const = 0;

    /** \brief Fill a pair of given arrays with data that represents internal state of this periodicity test, such as a folded
               light curve for the chi-squared test, such that they can be used as X- and Y-axis of a plot to display. The given
               arrays may be resized if necessary. Details depend on the specific test being performed in the subclass.
        \param array_index The index of the element of the periodicity test array, of which a data array is to be created.
        \param x_data The output array that contains X-values for a plot to display.
        \param y_data The output array that contains Y-values for a plot to display.
    */
    virtual void getPlotData(size_type array_index, std::vector<double> & x_data, std::vector<double> & y_data) const = 0;

    /** \brief Assign axis labels to a pair of given strings, each of which can be used as an X- and Y-axis label, respectively.
        \param x_data The output string that contains the label for X-axis of a plot to display.
        \param y_data The output string that contains the label for Y-axis of a plot to display.
    */
    virtual void getPlotLabel(std::string & x_label, std::string & y_label) const = 0;

    /** \brief Return a plot title that can be used with a return value of getPlotData method.
    */
    virtual std::string getPlotTitle() const = 0;

    /** \brief Return the name of periodicity test being performed in the subclass.
    */
    virtual std::string getTestName() const = 0;

    /** \brief Get a reference to an internal statistic viewer for an object of this class.
        \param array_index The index of the element of the periodicity test array, for which a viewer is to be configured.
    */
    virtual StatisticViewer & getViewer(size_type array_index) = 0;

  protected:
    /** \brief Construct a periodicity test array object.
    */
    PeriodicityTestArray(StatisticViewer::index_type num_axis, StatisticViewer::data_type::size_type num_element):
      m_viewer(num_axis, num_element) {};

    StatisticViewer m_viewer;
};

#endif
