/** \file PeriodSearch.cxx
    \brief Implementation of PeriodSearch class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include "periodSearch/PeriodSearch.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace periodSearch {

  const double PeriodSearch::s_2pi = 2. * 4. * atan(1.0);

  PeriodSearch::PeriodSearch(size_type num_bins): m_freq(num_bins), m_spec(num_bins) {}

  std::pair<double, double> PeriodSearch::findMax(double min_freq, double max_freq) const {
    bool found_max = false;
    size_type max_idx = 0;
    double max = 0.;

    // Impose range limits.
    std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
    size_type begin_index = indices.first;
    size_type end_index = indices.second;

    for (size_type ii = begin_index; ii < end_index; ++ii) {
      // If the value is larger than the current maximum, replace it.
      if (m_spec[ii] > max) {
        max = m_spec[ii];
        max_idx = ii;
        found_max = true;
      }
    }

    // Make sure a valid maximum was found.
    if (!found_max) {
      std::ostringstream os;
      os << "PeriodSearch::findMax cannot find any trial frequency in range [" << min_freq << ", " << max_freq << "]";
      throw std::runtime_error(os.str());
    }
    return std::pair<double, double>(m_freq[max_idx], max);
  }

  std::pair<double, double> PeriodSearch::chanceProb(double stat) const {
    // Compute probability for one trial.
    std::pair<double, double> chance_prob = chanceProbOneTrial(stat);

    // Get "N".
    size_type num_indep_trials = numIndepTrials();

    // Multiply by N (small probability approximation).
    chance_prob.first *= num_indep_trials;
    chance_prob.second *= num_indep_trials;

    // Limit to legally bounded range for probability.
    chance_prob.first = chance_prob.first > 1. ? 1. : chance_prob.first;
    chance_prob.second = chance_prob.second > 1. ? 1. : chance_prob.second;

    return chance_prob;
  }

  st_stream::OStream & PeriodSearch::write(st_stream::OStream & os) const {
    return writeRange(os);
  }

  const PeriodSearch::cont_type PeriodSearch::getFreq() const { return m_freq; }
      
  const PeriodSearch::cont_type PeriodSearch::getSpec() const { return m_spec; }

  double PeriodSearch::chanceProbMultiTrial(double prob_one_trial, size_type num_indep_trial) {
    static double epsilon = std::numeric_limits<double>::epsilon();
    double chance_prob = 0.;

    // Handle special cases.
    if (0 == num_indep_trial) chance_prob = 0.;
    else if (1 == num_indep_trial) chance_prob = prob_one_trial;
    else if (1. == prob_one_trial) chance_prob = 1.;
    else if (0. == prob_one_trial) chance_prob = 0.;
    else {
      // Compute x == -N * ln(1 - p) where N == num_indep_trial and p == prob_one_trial.
      double xx = 0.;
      if (.1 <= prob_one_trial) {
        // For large p use the formula directly.
        xx = -(signed long)(num_indep_trial) * std::log(1. - prob_one_trial);
      } else {
        // For small values of p, use iterative power series expansion.
        // Helper: p^n.
        double p_to_the_n = prob_one_trial;
        // qn == p^n/n -> q1 = p / 1.
        double q_sub_n = p_to_the_n;
        // Sn == sum qn -> Initial sum is 0.
        double S_sub_n = 0.;
        // Iterate over nn starting with 1.
        size_type nn = 1;
        do {
          // Compute Sn.
          S_sub_n += q_sub_n;
          // Prepare the next qn term and helper p_to_the_n:
          // n becomes "n+1"
          ++nn;
          // Compute p^"n+1"
          p_to_the_n *= prob_one_trial;
          // Compute q_sub_"n+1"
          q_sub_n = p_to_the_n / nn;
        } while (epsilon < q_sub_n / S_sub_n); // Note that this compares q_sub_n for the current nn to S_sub_n for the previous nn.

        // x == N * Sn for last n.
        xx = num_indep_trial * S_sub_n;
      }

      // Compute chance_prob == Q == 1. - exp(-x).
      if (.1 <= xx) {
        // For large x use the formula directly.
        chance_prob = 1. - std::exp(-xx);
      } else {
        // For small values of x, use iterative power series expansion.
        // Helper: leading term is x^(2n+1) / (2n + 1)!. For n = 0 this is x^1/1 = x.
        double leading_term = xx;
        // qn == leading term * (1 - xx/(2n + 2)). For n = 0 this is leading term * (1 - x/2)
        double q_sub_n = leading_term * (1. - xx / 2.);
        // Sn == sum qn -> Initial sum is 0.
        double S_sub_n = 0.;
        // Iterate over nn starting with 1.
        size_type nn = 0;
        do {
          // Compute Sn.
          S_sub_n += q_sub_n;
          // Prepare the next helper leading_term and qn.
          // n becomes "n+1"
          ++nn;
          // Compute next leading term.
          leading_term *= xx * xx / (2 * nn * (2 * nn + 1));
          // Compute q_sub_"n+1"
          q_sub_n = leading_term * (1. - xx / (2 * nn + 2));
        } while (epsilon < q_sub_n / S_sub_n); // Note that this compares q_sub_n for the current nn to S_sub_n for the previous nn.

        // Q == Sn for last n.
        chance_prob = S_sub_n;
      }
    }

    return chance_prob;
  }

  std::pair<PeriodSearch::size_type, PeriodSearch::size_type> PeriodSearch::getRangeIndex(double min_freq,
    double max_freq) const {
    size_type begin_index = 0;
    size_type end_index = m_freq.size();

    if (0. <= min_freq) {
      // Find first element whose frequency is not less than the minimum frequency.
      for (begin_index = 0; begin_index < m_freq.size() && m_freq[begin_index] < min_freq; ++begin_index);
    }

    if (0. <= max_freq) {
      // Find last element whose frequency is not greater than the maximum frequency.
      for (end_index = m_freq.size(); end_index > begin_index && m_freq[end_index - 1] > max_freq; --end_index);
    }

    return std::make_pair(begin_index, end_index);
  }

  st_stream::OStream & PeriodSearch::writeRange(st_stream::OStream & os, double min_freq, double max_freq) const {
    using namespace std;

    // Get info about the maximum.
    std::pair<double, double> max = findMax(min_freq, max_freq);

    // Chance probability.
    size_type num_indep_trials = numIndepTrials();
    std::pair<double, double> chance_prob = chanceProb(max.second);

    // Save current precision.
    int save_precision = os.precision();

    os.precision(15);

    // Write out the results.
    os << "Maximum at: " << max.first << std::endl << "Statistic: " << max.second << std::endl;
    os << "Chance probability range: (" << chance_prob.first << ", " << chance_prob.second << ")" << std::endl;
    os << "Number of statistically independent trials: " << num_indep_trials << std::endl;
    os << "Frequency\tStatistic";

    // Impose range limits.
    std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
    size_type begin_index = indices.first;
    size_type end_index = indices.second;

    // Write out the statistics.
    for (size_type ii = begin_index; ii < end_index; ++ii) os << std::endl << m_freq[ii] << "\t" << m_spec[ii];

    // Restore original precision.
    os.precision(save_precision);

    return os;
  }

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearch & test) {
    return test.write(os);
  }

}
