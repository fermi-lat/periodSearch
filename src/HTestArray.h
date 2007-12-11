/** \file HTestArray.h
    \brief Declaration of HTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_HTestArray_h
#define periodSearch_HTestArray_h

#include <string>
#include <utility>
#include <vector>

#include "Z2nTestArray.h"

class StatisticViewer;

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

    /** \brief Compute an H-value of this H test for pulse phases currently filled in this object, and set
               candidate H values to the internal statistic viewer.
        \param array_index The index of the element of the periodicity test array, of which an H-value is to be computed.
    */
    virtual double testStat(size_type array_index) const;

    /** \brief Compute an H-value of this H test for pulse phases currently filled in this object.
        \param array_index The index of the element of the periodicity test array, of which an H-value is to be computed.
    */
    virtual void updateViewer(size_type array_index);

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProb(double stat) const;

    /** \brief Return a description of this test array.
    */
    virtual std::string getDescription() const;

    /** \brief Return the name of this periodicity test.
    */
    virtual std::string getTestName() const;

  private:
    /** \brief Compute candidates for H-value for each harmonic numbers, from the given Fourier powers (i.e., squared sum
               of sine and cosine component).
        \param power The container of the input Fourier powers.
        \param H_candidate the container of the candidates of H-value (output).
    */
    void computeCandidate(data_type & power, data_type & H_candidate) const;
    void computeCandidate();

    // The maximum number of harmonics (to be used only in getDescription method).
    size_type m_max_harm;
};

#endif
