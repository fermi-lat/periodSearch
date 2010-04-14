/** \file FoldingAnalysis.h
    \brief Declaration of FoldingAnalysis class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_FoldingAnalysis_h
#define periodSearch_FoldingAnalysis_h

#include <string>
#include <utility>
#include <vector>

#include "PeriodSearch.h"

class PeriodicityTestArray;

/** \class FoldingAnalysis
    \brief Class for periodicity search based on various statistical tests used to determine frequency of pulsation
           when an approximate frequency is known, using the epoch-folding technique. 

           Note on the units of time and frequency: This class does not assume a specific time unit. Instead, all
           times (both time series data and time durations) given to an object of this class must be in the same time
           units throughout the lifetime of the object, and all frequencies given to it must be in the units of the
           inverse of the time unit. All results from the object will be in the same units of time or frequency.
*/
class FoldingAnalysis : public PeriodSearch {
  public:
    /** \brief Construct a FoldingAnalysis object, which performs periodicity search using given search criteria
               and a statistical test object.
        \param center The central frequency to test.
        \param step The step size of frequency sampling.
        \param test_array The statistical test array object to be used for this periodicity search.
        \param epoch The global time offset defining the origin for purposes of this computation.
        \param duration The total time duration (only used if chance probability will be computed).
        \param freq_unit Character string representing the unit of frequency, or the inverse of the unit of times given
               to this period search.
    */
    FoldingAnalysis(PeriodicityTestArray * test_array, double center, double step, double epoch, double duration,
      const std::string & freq_unit);

    /// \brief Destruct this FoldingAnalysis object.
    virtual ~FoldingAnalysis() {}

    /** \brief Fill given time into periodicity test objects.
        \param evt_time The time of the event.
    */
    virtual void fill(double evt_time);

    /** \brief Compute test statistics for currently filled event times at all trial frequencies.
               Details depend on the a specific test object given upon construction of this object.
    */
    virtual void computeStat();

    /** \brief Compute the number of independent trials for this search method.
        \param min_freq The minimum frequency.
        \param max_freq The maximum frequency.
    */
    virtual size_type computeNumIndepTrials(double min_freq = -1., double max_freq = -1.) const;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> computeChanceProbOneTrial(double stat) const;

  protected:
    // The container of statistical trials.
    PeriodicityTestArray * m_test_array;

    // Sampling frequency (frequency step).
    double m_step;

    // Temporal origin.
    double m_epoch;

    // Fourier resolution of search.
    double m_fourier_res;
};

#endif
