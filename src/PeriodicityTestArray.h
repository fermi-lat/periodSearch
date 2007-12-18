/** \file PeriodicityTestArray.h
    \brief Declaration of PeriodicityTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodicityTestArray_h
#define periodSearch_PeriodicityTestArray_h

#include <cstddef>
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
    typedef std::size_t size_type;

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
    virtual double computeStat(size_type array_index) const = 0;

    /** \brief Compute a test statistic for pulse phases currently filled in this object, and set viewable data in
               the internal statistic viewer. Details depend on the specific test being performed in the subclass.
        \param array_index The index of the element of the periodicity test array, of which a test statistic is to be computed.
    */
    virtual void updateViewer(size_type array_index) = 0;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the test statistic.
    */
    virtual std::pair<double, double> computeChanceProb(double stat) const = 0;

    /** \brief Return the size of this periodicity test array.
    */
    virtual size_type size() const = 0;

    /** \brief Return the name of periodicity test being performed in the subclass.
    */
    std::string getTestName() const { return m_name; };

    /** \brief Return a description of this periodicity test array.
    */
    std::string getDescription() const { return m_description; };

    /** \brief Get a reference to an internal statistic viewer for an object of this class.
    */
    StatisticViewer & getViewer() { return m_viewer; };

  protected:
    /** \brief Construct a periodicity test array object.
        \param num_axis The number of axis for the internal data storage.
        \param num_element The number of element for the internal data storage.
    */
    PeriodicityTestArray(StatisticViewer::index_type num_axis, StatisticViewer::data_type::size_type num_element):
      m_name(), m_description(), m_viewer(num_axis, num_element) {};

    /** \brief Set the name and the description of this periodicity test.
        \param test_name The name of this periodicity test.
        \param test_cond The description of test condition, such as "10 phase bins" for the chi-squared test.
        \param prob_dist The description of the probability distribution that the test statistics follow.
    */
    void setDescription(const std::string & test_name, const std::string & test_cond, const std::string & prob_dist) {
      // Set the test name.
      m_name = test_name;

      // Set the test description.
      std::ostringstream os;
      os << "Type of test: " << m_name;
      if (!test_cond.empty()) os << ", " << test_cond;
      os << std::endl;
      os << "Probability distribution: " << prob_dist;
      m_description = os.str();
    };

    /** \brief Return a summary of this periodicity test array.
        \param stat The value of the test statistic.
    */
    std::string createSummary(double stat) {
      std::pair<double, double> chance_prob = computeChanceProb(stat);

      std::ostringstream os;
      os.precision(std::numeric_limits<double>::digits10);
      os << getDescription() << std::endl;
      os << "Test Statistic: " << stat << std::endl;
      os << "Chance Probability Range: " << "(" << chance_prob.first << ", " << chance_prob.second << ")";

      return os.str();
    }

  private:
    std::string m_name;
    std::string m_description;
    StatisticViewer m_viewer;
};

#endif
