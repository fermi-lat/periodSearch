/** \file PeriodSearch.h
    \brief Declaration of PeriodSearch class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodSearch_h
#define periodSearch_PeriodSearch_h

#include <complex>
#include <string>
#include <utility>
#include <vector>

#include "st_stream/Stream.h"

#include "StatisticViewer.h"

namespace periodSearch {

  /** \class PeriodSearchResult
      \brief Class encapsulating result of a PeriodSearch.
  */
  class PeriodSearchResult {
    public:
      typedef std::vector<double> cont_type;
      typedef cont_type::size_type size_type;

      /** \brief Create a search result.
          \param description Description of the probability distribution of the statistic.
          \param min_freq The minimum frequency in the range.
          \param max_freq The maximum frequency in the range.
          \param num_freq_bin The number of trial frequencies used.
          \param num_indep_trial The number of independent trial frequencies used.
          \param max_stat The frequency where the maximum statistic occurs and the value of the statistic at that frequency.
          \param chance_prob Range of chance probabilities.
      */
      PeriodSearchResult(const std::string & description, double min_freq, double max_freq, size_type num_freq_bin,
        size_type num_indep_trial, const std::pair<double, double> & max_stat, const std::pair<double, double> & chance_prob);

      /** \brief Write this search result to the given stream.
          \param os The stream.
      */
      template <typename StreamType>
      StreamType & write(StreamType & os) const;

    private:
      std::string m_description;
      double m_min_freq;
      double m_max_freq;
      size_type m_num_freq_bin;
      size_type m_num_indep_trial;
      std::pair<double, double> m_max_stat;
      std::pair<double, double> m_chance_prob;
  };

  template <typename StreamType>
  inline StreamType & PeriodSearchResult::write(StreamType & os) const {
    std::streamsize orig_precision = os.precision();
    os.precision(std::numeric_limits<double>::digits10);
    os << m_description << "\n"
       << "Search Range (Hz): [" << m_min_freq << ", " << m_max_freq << "]\n"
       << "Number of Trial Frequencies: " << m_num_freq_bin << "\n"
       << "Number of Independent Trials: " << m_num_indep_trial << "\n"
       << "Maximum Statistic: " << m_max_stat.second << " at " << m_max_stat.first << " Hz\n"
       << "Chance Probability Range: " << "(" << m_chance_prob.first << ", " << m_chance_prob.second << ")";
    os.precision(orig_precision);
    return os;
  }

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearchResult & result);

  /** \class PeriodSearch
      \brief Base class for various statistical tests used to determine frequency of pulsation
             when an approximate frequency is known.
  */
  class PeriodSearch {
    public:
      typedef std::vector<double> cont_type;
      typedef cont_type::size_type size_type;

      PeriodSearch(size_type num_bins);

      virtual ~PeriodSearch() {}

      /** \brief Fill given time into histograms.
          \param evt_time The time of the event.
      */
      virtual void fill(double evt_time) = 0;

      /** \brief Use the trials as currently filled by data to compute statistics for this test. Details
                 depend on the specific test being performed in the subclass.
      */
      virtual const std::vector<double> & computeStats() = 0;

      /** \brief Perform a period search and return the result of the test. */
      virtual PeriodSearchResult search(double min_freq = -1., double max_freq = -1.) const;

      /** \brief Find the frequency for which the statistic is maximized in a given frequency range. Return the
                 frequency and the value of the statistic, as a pair.
          \param min_freq The minimum frequency in the range.
          \param max_freq The maximum frequency in the range.
      */
      virtual std::pair<double, double> findMax(double min_freq = -1., double max_freq = -1.) const;

      /** \brief Return the number of independent trials for this search method.
      */
      virtual size_type numIndepTrials(double min_freq = -1., double max_freq = -1.) const = 0;
      // TODO: Rename the above to getNumIndepTrials for consistency.

      /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
          \param stat The value of the statistic.
      */
      virtual std::pair<double, double> chanceProbOneTrial(double stat) const = 0;
      // TODO: Rename the above to computeChanceProb for consistency.

      /** \brief Output a description of this search.
          \param os The stream.
      */
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

      /** \brief Return a description of this search.
      */
      virtual std::string getDescription() const = 0;

      //* \brief Get frequency data associated with this search object.
      const cont_type getFreq() const;
      
      //* \brief Get spectral data associated with this search object.
      const cont_type getSpec() const;
      
      /* \brief Compute the probability that an event occurs at least once in N statistically
                independent trials, given that the probability of the event occurring in a single trial is p.
         \param prob_one_trial The probability p of the event occuring in one trial.
         \param num_indep_trial The number N of statistically independent trials.
      */
      static double chanceProbMultiTrial(double prob_one_trial, size_type num_indep_trial);
      // TODO: Rename the above to computeChanceProb for consistency.

      /** \brief Given a frequency range, determine the indices of (inclusive) lower and upper bounds.
          \param min_freq The minimum frequency.
          \param max_freq The maximum frequency.
      */
      std::pair<size_type, size_type> getRangeIndex(double min_freq, double max_freq) const;

      /** \brief Get a reference to an internal statistic viewer for an object of this class.
          \param min_freq The minimum frequency to view.
          \param max_freq The maximum frequency to view.
      */
      virtual StatisticViewer & getViewer(double min_freq = -1., double max_freq = -1.);

    protected:

      // For convenience, compute once the value of 2 * pi.
      static const double s_2pi;

      // The frequencies forming the domain of the search/test.
      // TODO: Remove the below after ChiSquaredTest.{h,cxx} are cvs-removed.
      cont_type m_freq;

      // The statistical measure of the validity of each trial.
      // TODO: Remove the below after ChiSquaredTest.{h,cxx} are cvs-removed.
      cont_type m_spec;

      // The statistical measure of the validity of each trial.
      StatisticViewer m_viewer;
  };

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearch & test);

}

#endif
