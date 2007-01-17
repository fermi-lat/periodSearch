/** \file ChiSquaredTest.cxx
    \brief Implementation of ChiSquaredTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#include <stdexcept>

#include "ChiSquaredProb.h"
#include "ChiSquaredTest.h"

ChiSquaredTest::ChiSquaredTest(double center, double step, size_type num_trials, double epoch, size_type num_phase_bins,
  double duration): periodSearch::PeriodTest(center, step, num_trials, epoch, num_phase_bins, duration),
  m_num_phase_bins(num_phase_bins) {
}

void ChiSquaredTest::fillOneTrial(double phase, std::vector<std::complex<double> > & trial) const {
  // Bin phase; it runs from [0, 1), so multiply by the number of bins to determine
  // bin for this phase.
  size_type bin_id = size_type(phase * m_num_phase_bins);

  // Increment the count in that bin.
  trial[bin_id] += 1.;
}

const std::vector<double> & ChiSquaredTest::computeStats() {
  m_spec.assign(m_spec.size(), 0.);
  if (0 == m_num_events) return m_spec;

  // Pseudocode taken from work by M. Hirayama, which may be seen at:
  // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/pc_chi2test.txt
  //     avg = N_event / N_bin;
  double avg = double(m_num_events) / m_num_phase_bins;

  //     for (i=0; i<N_trial; i++) {
  //       S_value[i] = 0.0;
  //       for (j=0; j<N_bin; j++) {
  //         value = [get the number of entries in the j-th bin of hist[i]];
  //         S_value[i] += (value - avg)*(value - avg)/avg;
  //       }
  //     }

  // Iterate over the number of trials.
  size_type num_trials = m_trial_hist.size();
  for (size_type ii = 0; ii < num_trials; ++ii) {
    // Reset statistics for this trial each time this is called.
    m_spec[ii] = 0.;

    // Iterate over bins in each trial.
    for (size_type jj = 0; jj < m_num_phase_bins; ++jj) {
      // Compute deviation. Imaginary part of trial container is ignored.
      double dev = (m_trial_hist[ii][jj] - avg).real();

      // Sum squares of deviation divided by the average value.
      m_spec[ii] += dev * dev / avg;
    }

  }

  return m_spec;

}

std::pair<double,double> ChiSquaredTest::chanceProbOneTrial(double stat) const {
  //    /* Leahy et al. 1983, ApJ 266, 160 */
  //    chance_prob = chi2prob(S_value[imax], N_bin-1) * N_Fourier;
  //    [where function chi2prob(chisq, dof) returns the chi-squared
  //     distribution for "dof" degrees of freedom, integrated from "chisq"
  //     to infinity];
  periodSearch::ChiSquaredProb prob(m_num_phase_bins - 1);
  return prob(stat);
}
