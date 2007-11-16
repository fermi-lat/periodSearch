/** \file PeriodicityTest.h
    \brief Declaration of PeriodicityTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodicityTest_h
#define periodSearch_PeriodicityTest_h

#include <complex>
#include <string>
#include <utility>
#include <vector>

namespace periodSearch {
  /** \class PeriodicityTest
      \brief Base class for various statistical tests used to evaluate significance of pulsation.
  */
  class PeriodicityTest {
    public:
      /// \brief Container type used to store test-specific histograms, such as a folded light curve.
      typedef std::vector<std::complex<double> > DataCont_t;

      virtual ~PeriodicityTest() {}

      /** \brief Create a clone object of this object.
      */
      virtual PeriodicityTest * clone() const = 0;

      /** \brief Fill a given pulse phase into this periodicity test object.
          \param phase The pulse phase.
      */
      virtual void fill(double phase) const = 0;

      /** \brief Compute a test statistic for pulse phases currently filled in this object. Details
                 depend on the specific test being performed in the subclass.
      */
      virtual const double & testStat() = 0;

      /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
          \param stat The value of the statistic.
      */
      virtual std::pair<double, double> chanceProb(double stat) const = 0;

      /** \brief Return a description of this search.
      */
      virtual std::string getDescription() const = 0;

    protected:
      /** \brief Construct a periodicity test object.
      */
      PeriodicityTest() {};

      // The container of test-specific data for this statistical test.
      DataCont_t m_data_cont;
  };
}

#endif
