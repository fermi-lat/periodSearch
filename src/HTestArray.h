/** \file HTestArray.h
    \brief Declaration of HTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_HTestArray_h
#define periodSearch_HTestArray_h

#include <utility>
#include <vector>

#include "Z2nTestArray.h"

/** \class HTestArray
    \brief PeriodicityTestArray subclass which uses a H statistic for its test.
*/
class HTestArray : public Z2nTestArray {
  public:
    typedef std::vector<double> data_type;
    typedef std::vector<data_type> cont_type;

    /** \brief Construct an array object of H test.
        \param array_size The size of this test array.
        \param num_phase_bins The maximum number of harmonics for the H test.
    */
    HTestArray(size_type array_size, data_type::size_type max_harmonics);

    virtual ~HTestArray() {}

    /** \brief Compute an H-value of this H test for pulse phases currently filled in this object.
        \param array_index The index of the element of the periodicity test array, of which an H-value is to be computed.
    */
    virtual double testStat(size_type array_index) const;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProb(double stat) const;

    /** \brief Return a description of this test array.
    */
    virtual std::string getDescription() const;

    /** \brief Fill a pair of given arrays with data that contains the candidate H values for each harmonic number,
               that are the Z2n values subtracted by 4*(m-1), where m is the harmonic number.
        \param array_index The index of the element of the periodicity test array, of which a data array is to be created.
        \param harmonic The output array for harmonic numbers, which may be used an X-axis of a plot.
        \param power The output array for the candidate H values for each harmonic number, which may be used an Y-axis of a plot.
    */
    virtual void getPlotData(size_type array_index, std::vector<double> & harmonic, std::vector<double> & H_value) const;

    /** \brief Assign axis labels to a pair of given strings, each of which can be used as an X- and Y-axis label, respectively.
        \param x_data The output string that contains the label for X-axis of a plot to display.
        \param y_data The output string that contains the label for Y-axis of a plot to display.
    */
    virtual void getPlotLabel(std::string & x_label, std::string & y_label) const;

    /** \brief Return a plot title that can be used with a return value of getPlotData method.
    */
    virtual std::string getPlotTitle() const;

  private:
    /** \brief Compute candidates for H-value for each harmonic numbers, from the given Fourier powers (i.e., squared sum
               of sine and cosine component).
        \param power The container of the input Fourier powers.
        \param H_candidate the container of the candidates of H-value (output).
    */
    void computeCandidate(data_type & power, data_type & H_candidate) const;

    // The maximum number of harmonics (to be used only in getDescription method).
    size_type m_max_harm;
};

#endif
