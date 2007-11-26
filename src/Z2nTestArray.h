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
        \param phase The pulse phase.
    */
    virtual void fill(double phase, size_type array_index = 0);

    /** \brief Compute a Z2n value of this Z2n test for pulse phases currently filled in this object.
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

    /** \brief Return a pair of data arrays that contains the Fourier power (i.e., squared sum of sine and cosine component)
               for each harmonic number. The first array of the pair is an array of harmonic numbers, and the second an array
               of the Fourier powers for each harmonic number.
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
    // The number of harmonics to sum up.
    size_type m_num_harm;

    // The container for sine and cosine components.
    cont_type m_sine_cont;
    cont_type m_cosine_cont;

    // The number of events filled for each element of this test array.
    std::vector<long> m_num_events;
};

#endif
