/** \file HTest.cxx
    \brief Implementation of HTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <sstream>
#include <stdexcept>

#include "HTest.h"

HTest::HTest(double center, double step, size_type num_trials, double epoch, size_type max_harmonics, double duration):
  periodSearch::PeriodTest(center, step, num_trials, epoch, max_harmonics, duration), m_max_harm(max_harmonics) {
}

void HTest::fillOneTrial(double phase, std::vector<std::complex<double> > & trial) const {
  // For each phase, the complex Fourier component is computed for each trial harmonic.
  for (size_type jj = 0; jj < m_max_harm; ++jj)
    trial[jj] += exp(std::complex<double>(0., s_2pi * (jj + 1) * phase));
}

const std::vector<double> & HTest::computeStats() {
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
  //       z2_value = 0.0;
  //       test_stat[i] = 0.0;
  //       for (j=0; j<N_harm; j++) {
  //         z2_value += powspec[i][j];
  //         H_value = z2_value - 4.0*(double)j;
  //         if (H_value > test_stat[i]) {
  //           test_stat[i] = H_value;
  //         }
  //       }
  //     }

  // Iterate over the number of trials.
  size_type num_trials = m_trial_hist.size();
  for (size_type ii = 0; ii < num_trials; ++ii) {
    // Reset statistics for this trial each time this is called.
    m_spec[ii] = 0.;

    double z2_value = 0.;
    // Iterate over bins in each trial.
    for (size_type jj = 0; jj < m_max_harm; ++jj) {
      // Compute coefficient of power spectrum for each harmonic.
      z2_value += fourier_norm * norm(m_trial_hist[ii][jj]);
      double H_value = z2_value - 4. * jj;
      if (H_value > m_spec[ii]) {
        // Keep only the harmonic with the highest H_value.
        m_spec[ii] = H_value;
      }
    }
  }

  return m_spec;
}

std::pair<double, double> HTest::chanceProbOneTrial(double stat) const {
  /* De Jager et al. 1989, A&A 221, 180 */
  double lower_limit;
  double upper_limit;
  bool chance_prob_exact = true;

  if (stat <= 23.0) {
     double a = 0.9999755;
     double b = 0.39802;
     upper_limit = a * exp(-b * stat);
  } else if (stat < 50.0) {
     double c = 1.210597;
     double d = 0.45901;
     double e = 0.0022900;
     upper_limit = c * exp(-d * stat + e * stat * stat);
  } else {
     upper_limit = 4.0e-8; /* or less */
     chance_prob_exact = false;
  }

  upper_limit = upper_limit < 1. ? upper_limit : 1.;

  lower_limit = chance_prob_exact ? upper_limit : 0.;

  return std::make_pair(lower_limit, upper_limit);
}

std::string HTest::getDescription() const {
  std::ostringstream os;
  os << PeriodTest::getDescription() << "\n" <<
    "Type of test: H Test, " << m_max_harm << " maximum harmonics\n" <<
    "Probability distribution: H Test-specific";
  return os.str();
}
