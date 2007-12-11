/** \file PeriodicityTestArray.h
    \brief Declaration of PeriodicityTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodicityTestArray_h
#define periodSearch_PeriodicityTestArray_h

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

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

    /** \brief Compute a test statistic for pulse phases currently filled in this object, and set viewable data in
               the internal statistic viewer. Details depend on the specific test being performed in the subclass.
        \param array_index The index of the element of the periodicity test array, of which a test statistic is to be computed.
    */
    virtual void updateViewer(size_type array_index) = 0;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProb(double stat) const = 0;
    // TODO: Rename the above to computeChanceProb for consistency.

    /** \brief Return a description of this periodicity test array.
    */
    virtual std::string getDescription() const = 0;

    /** \brief Return a summary of this periodicity test array.
        \param array_index The index of the element of the periodicity test array, of which a summary is to be returned.
    */
    virtual std::string getSummary(size_type array_index) {
      double test_stat = testStat(array_index);
      std::pair<double, double> chance_prob = chanceProb(test_stat);

      std::ostringstream os;
      os.precision(std::numeric_limits<double>::digits10);
      os << getDescription() << std::endl;
      os << "Test Statistic: " << test_stat << std::endl;
      os << "Chance Probability Range: " << "(" << chance_prob.first << ", " << chance_prob.second << ")";

      return os.str();
    }

    /** \brief Return the size of this periodicity test array.
    */
    virtual size_type size() const = 0;

    /** \brief Return the name of periodicity test being performed in the subclass.
    */
    virtual std::string getTestName() const = 0;

    /** \brief Get a reference to an internal statistic viewer for an object of this class.
    */
    virtual StatisticViewer & getViewer() { return m_viewer; };

  protected:
    /** \brief Construct a periodicity test array object.
    */
    PeriodicityTestArray(StatisticViewer::index_type num_axis, StatisticViewer::data_type::size_type num_element):
      m_viewer(num_axis, num_element) {};

    StatisticViewer m_viewer;
};

#endif
