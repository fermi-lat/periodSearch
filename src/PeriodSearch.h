/** \file PeriodSearch.h
    \brief Declaration of PeriodSearch class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodSearch_h
#define periodSearch_PeriodSearch_h

#include <string>
#include <utility>
#include <vector>

#include "st_stream/Stream.h"

#include "StatisticViewer.h"

/** \class PeriodSearch
    \brief Base class for various statistical tests used to determine frequency of pulsation
           when an approximate frequency is known.
*/
class PeriodSearch {
  public:
    typedef std::vector<double> cont_type;
    typedef cont_type::size_type size_type;

    /** \brief Construct a PeriodSearch object.
        \param num_bins The number of frequency bins.
        \param freq_unit The unit of frequency, or the inverse of the unit of times given to this period search.
    */
    PeriodSearch(size_type num_bins, const std::string & freq_unit);

    /// \brief Destruct this PeriodSearch object.
    virtual ~PeriodSearch() {}

    /** \brief Fill given time into histograms.
        \param evt_time The time of the event.
    */
    virtual void fill(double evt_time) = 0;

    /** \brief Use the trials as currently filled by data to compute statistics for this test. Details
               depend on the specific test being performed in the subclass.
    */
    virtual void computeStat() = 0;

    /** \brief Perform a period search and set the result to the internal statistic viewer.
        \param min_freq The minimum frequency in the range.
        \param max_freq The maximum frequency in the range.
    */
    virtual void updateViewer(double min_freq = -1., double max_freq = -1.);

    /** \brief Find the frequency for which the statistic is maximized in a given frequency range. Return the
               frequency and the value of the statistic, as a pair.
        \param min_freq The minimum frequency in the range.
        \param max_freq The maximum frequency in the range.
    */
    virtual std::pair<double, double> findMax(double min_freq = -1., double max_freq = -1.) const;

    /** \brief Compute the number of independent trials for this search method.
        \param min_freq The minimum frequency.
        \param max_freq The maximum frequency.
    */
    virtual size_type computeNumIndepTrials(double min_freq = -1., double max_freq = -1.) const = 0;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> computeChanceProbOneTrial(double stat) const = 0;

    /** \brief Output a description of this search.
        \param os The stream.
    */
    virtual st_stream::OStream & write(st_stream::OStream & os) const;

    /* \brief Compute the probability that an event occurs at least once in N statistically
              independent trials, given that the probability of the event occurring in a single trial is p.
       \param prob_one_trial The probability p of the event occuring in one trial.
       \param num_indep_trial The number N of statistically independent trials.
    */
    static double computeChanceProbMultiTrial(double prob_one_trial, size_type num_indep_trial);

    /** \brief Given a frequency range, determine the indices of (inclusive) lower and upper bounds.
        \param min_freq The minimum frequency.
        \param max_freq The maximum frequency.
    */
    std::pair<size_type, size_type> getRangeIndex(double min_freq, double max_freq) const;

    /** \brief Return a description of this search.
    */
    std::string getDescription() const;

    /** \brief Get a reference to an internal statistic viewer for an object of this class.
    */
    StatisticViewer & getViewer();
    const StatisticViewer & getViewer() const;

  protected:
    /** \brief Set the name and the description of this periodicity test.
        \param search_type The type of this period search, such as "Folding Analysis".
        \param fourier_res The Fourier resolution of this period search.
        \param sampling_step The sampling step of this period search.
        \param search_info The detailed information about this period search.
    */
    void setDescription(const std::string & search_type, double fourier_res, double sampling_step, const std::string & search_info);

  private:
    // The unit of frequency.
    std::string m_freq_unit;

    // The description of this period search.
    std::string m_description;

    // The statistical measure of the validity of each trial.
    StatisticViewer m_viewer;
};

st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearch & test);

#endif
