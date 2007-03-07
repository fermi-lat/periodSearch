/** \file Z2nTest.cxx
    \brief Implementation of Z2nTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <sstream>
#include <stdexcept>

#include "ChiSquaredProb.h"
#include "Z2nTest.h"

Z2nTest::Z2nTest(double center, double step, size_type num_trials, double epoch, size_type num_harmonics, double duration):
  periodSearch::PeriodTest(center, step, num_trials, epoch, num_harmonics, duration), m_num_harm(num_harmonics) {
}

void Z2nTest::fillOneTrial(double phase, std::vector<std::complex<double> > & trial) const {
  // For each phase, the complex Fourier component is computed for each trial harmonic.
  for (size_type jj = 0; jj < m_num_harm; ++jj)
    trial[jj] += exp(std::complex<double>(0., s_2pi * (jj + 1) * phase));
}

const std::vector<double> & Z2nTest::computeStats() {
  m_spec.assign(m_spec.size(), 0.);
  if (0 == m_num_events) return m_spec;

  // Pseudocode taken from work by M. Hirayama, which may be seen at:
  // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/pc_z2ntest.txt
  //     norm = 2.0 / N_event;
  double fourier_norm = 2. / m_num_events;

  //     for (i=0; i<N_trial; i++) {
  //       for (j=0; j<N_harm; j++) {
  //         powspec[i][j] = norm * (sine_comp[i][j])^2 + cosine_comp[i][j])^2);
  //       }
  //     }
  // The above block of pseudo-code is built into the complex class.

  //     for (i=0; i<N_trial; i++) {
  //       test_stat[i] = 0.0;
  //       for (j=0; j<N_harm; j++) {
  //         test_stat[i] += powspec[i][j];
  //       }
  //     }

  // Iterate over the number of trials.
  size_type num_trials = m_trial_hist.size();
  for (size_type ii = 0; ii < num_trials; ++ii) {
    // Reset statistics for this trial each time this is called.
    m_spec[ii] = 0.;

    // Iterate over bins in each trial.
    for (size_type jj = 0; jj < m_num_harm; ++jj) {
      // Compute the power spectrum: just the sum of squared norms of all bins. (Note: norm of complex === norm^2)
      m_spec[ii] += norm(m_trial_hist[ii][jj]);
    }
    m_spec[ii] *= fourier_norm;
  }

  return m_spec;
}

std::pair<double, double> Z2nTest::chanceProbOneTrial(double stat) const {
  //      /* Leahy et al. 1983, ApJ 266, 160 */
  //      chance_prob = chi2prob(test_stat[imax], 2*N_harm) * N_Fourier;
  //      [where function chi2prob(chi2, dof) returns the chi-squared
  //       distribution for "dof" degrees of freedom, integrated from "chi2"
  //       to infinity];
  periodSearch::ChiSquaredProb prob(2 * m_num_harm);
  return prob(stat);
}

std::string Z2nTest::getDescription() const {
  std::ostringstream os;
  os << PeriodTest::getDescription() << "\n" <<
    "Type of test: Z2n Test, " << m_num_harm << " harmonics\n" <<
    "Probability distribution: Chi-squared, " << 2 * m_num_harm << " degrees of freedom";
  return os.str();
}
