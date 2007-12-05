/** \file Z2nTestArray.h
    \brief Declaration of Z2nTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_Z2nTestArray_h
#define periodSearch_Z2nTestArray_h

#include <utility>
#include <vector>

#include "PeriodicityTestArray.h"
#include "StatisticViewer.h"

/** \class Z2nTestArray
    \brief PeriodicityTestArray subclass which uses a Z2n statistic for its test.
*/
class Z2nTestArray : public PeriodicityTestArray {
  public:
    typedef std::vector<double> data_type;
    typedef std::vector<data_type> cont_type;

    /** \brief Construct an array object of the Z2n test.
        \param array_size The size of this test array.
        \param num_phase_bins The number of harmonics to sum up for the Z2n test.
    */
    Z2nTestArray(size_type array_size, data_type::size_type num_harmonics);

    virtual ~Z2nTestArray() {}

    /** \brief Fill a given pulse phase into this periodicity test object.
        \param array_index The index of the element of the periodicity test array, to which a given phase is filled.
        \param phase The pulse phase to fill.
    */
    virtual void fill(size_type array_index, double phase);

    /** \brief Compute a Z2n value of this Z2n test for pulse phases currently filled in this object.
        \param array_index The index of the element of the periodicity test array, of which a Z2n value is to be computed.
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

    /** \brief Fill a pair of given arrays with data that contains the Fourier power (i.e., squared sum of sine and cosine component)
               for each harmonic number.
        \param array_index The index of the element of the periodicity test array, of which a data array is to be created.
        \param harmonic The output array for harmonic numbers, which may be used an X-axis of a plot.
        \param power The output array for the Fourier power for each harmonic number, which may be used an Y-axis of a plot.
    */
    virtual void getPlotData(size_type array_index, std::vector<double> & harmonic, std::vector<double> & power) const;

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

    /** \brief Get a reference to an internal statistic viewer for an object of this class.
        \param array_index The index of the element of the periodicity test array, for which a viewer is to be configured.
    */
    virtual StatisticViewer & getViewer(size_type array_index);

  protected:
    /** \brief Compute the Fourier power (i.e., squared sum of sine and cosine component) for each harmonic numbers.
        \param array_index The index of the element of the periodicity test array, for which the Fourier power is computed.
        \param power The container of the output Fourier powers.
    */
    void computePower(size_type array_index, data_type & power) const;

  private:
    // The number of harmonics to sum up.
    size_type m_num_harm;

    // The container for sine and cosine components.
    cont_type m_sine_cont;
    cont_type m_cosine_cont;

    // The number of events filled for each element of this test array.
    std::vector<long> m_num_events;
};

#endif
