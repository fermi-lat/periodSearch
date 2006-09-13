/** \file ChiSquaredTest.cxx
    \brief Implementation of ChiSquaredTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#include <stdexcept>

#include "ChiSquaredProb.h"
#include "ChiSquaredTest.h"

ChiSquaredTest::ChiSquaredTest(double center, double step, long num_trials, double epoch, int num_bins, double duration):
  periodSearch::PeriodTest(center, step, num_trials, epoch, num_bins, duration) {
}

void ChiSquaredTest::fillOneTrial(double phase, std::vector<std::complex<double> > & trial) const {
  // Bin phase; it runs from [0, 1), so multiply by the number of bins to determine
  // bin for this phase.
  int bin_id = int(phase * m_num_bins);

  // Increment the count in that bin.
  trial[bin_id] += 1.;
}

const std::vector<double> & ChiSquaredTest::computeStats() {
  m_stats.assign(m_stats.size(), 0.);
  if (0 == m_num_events) return m_stats;

  // Pseudocode taken from work by M. Hirayama, which may be seen at:
  // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/pc_chi2test.txt
  //     avg = N_event / N_bin;
  double avg = double(m_num_events) / m_num_bins;

  //     for (i=0; i<N_trial; i++) {
  //       S_value[i] = 0.0;
  //       for (j=0; j<N_bin; j++) {
  //         value = [get the number of entries in the j-th bin of hist[i]];
  //         S_value[i] += (value - avg)*(value - avg)/avg;
  //       }
  //     }

  // Iterate over the number of trials.
  int num_trials = m_trial_hist.size();
  for (int ii = 0; ii < num_trials; ++ii) {
    // Reset statistics for this trial each time this is called.
    m_stats[ii] = 0.;

    // Iterate over bins in each trial.
    for (int jj = 0; jj < m_num_bins; ++jj) {
      // Compute deviation. Imaginary part of trial container is ignored.
      double dev = (m_trial_hist[ii][jj] - avg).real();

      // Sum squares of deviation divided by the average value.
      m_stats[ii] += dev * dev / avg;
    }

  }

  return m_stats;

}

std::pair<double,double> ChiSquaredTest::chanceProb(double stat) const {
  //    /* Leahy et al. 1983, ApJ 266, 160 */
  //    chance_prob = chi2prob(S_value[imax], N_bin-1) * N_Fourier;
  //    [where function chi2prob(chisq, dof) returns the chi-squared
  //     distribution for "dof" degrees of freedom, integrated from "chisq"
  //     to infinity];
  periodSearch::ChiSquaredProb prob(m_num_bins - 1);
  return prob(stat);
}
