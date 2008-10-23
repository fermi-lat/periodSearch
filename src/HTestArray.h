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

    /** \brief Construct an array object for the H test.
        \param array_size The size of this test array.
        \param num_phase_bins The maximum number of harmonics for the H test.
    */
    HTestArray(size_type array_size, data_type::size_type max_harmonics);

    /// \brief Destruct this array object for the H test.
    virtual ~HTestArray() {}

    /** \brief Compute an H-value of this H test for pulse phases currently filled in this object, and set
               candidate H values to the internal statistic viewer.
        \param array_index The index of the element of the periodicity test array, of which an H-value is to be computed.
    */
    virtual double computeStat(size_type array_index) const;

    /** \brief Compute an H-value of this H test for pulse phases currently filled in this object.
        \param array_index The index of the element of the periodicity test array, of which an H-value is to be computed.
    */
    virtual void updateViewer(size_type array_index);

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> computeChanceProb(double stat) const;

  private:
    /** \brief Compute candidates for H-value for each harmonic numbers, from the given Fourier powers (i.e., squared sum
               of sine and cosine component), and return the H-value.
        \param array_index The index of the element of the periodicity test array, of which an H-value is to be computed.
        \param candidate the container of the candidates of H-value (output).
    */
    double computeH(size_type array_index, data_type & candidate) const;
};

#endif
