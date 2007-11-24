/** \file HTestArray.cxx
    \brief Implementation of HTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <sstream>
#include <stdexcept>

#include "HTestArray.h"

HTestArray::HTestArray(size_type array_size, data_type::size_type max_harmonics):
  m_max_harm(max_harmonics), m_sine_cont(array_size, data_type(max_harmonics, 0.)),
  m_cosine_cont(array_size, data_type(max_harmonics, 0.)), m_num_events(array_size, 0) {}

void HTestArray::fill(double phase, size_type array_index) {
  // Define two pi (for convenience and clarity).
  static const double s_2pi = 2. * 4. * std::atan(1.0);

  // Get the storage for sine and consine component.
  data_type & sine_array = m_sine_cont.at(array_index);
  data_type & cosine_array = m_cosine_cont.at(array_index);

  // For each phase, the complex Fourier component is computed for each trial harmonic.
  for (size_type jj = 0; jj < m_max_harm; ++jj) {
    double phase_angle = s_2pi * (jj + 1) * phase;
    sine_array[jj] += std::sin(phase_angle);
    cosine_array[jj] += std::cos(phase_angle);
  }

  // Increment the number of events filled.
  ++(m_num_events.at(array_index));
}

double HTestArray::testStat(size_type array_index) const {
  // Get the storage for sine and consine component.
  const data_type & sine_array = m_sine_cont.at(array_index);
  const data_type & cosine_array = m_cosine_cont.at(array_index);

  // Compute normalization.
  double fourier_norm = 2. / m_num_events.at(array_index);

  // Compute H value.
  double highest_H = 0.;
  double z2_value = 0.;
  // Iterate over bins in each trial.
  for (size_type jj = 0; jj < m_max_harm; ++jj) {
    // Compute coefficient of power spectrum for each harmonic.
    z2_value += fourier_norm * (sine_array[jj] * sine_array[jj] + cosine_array[jj] * cosine_array[jj]);
    double H_value = z2_value - 4. * jj;
    if (H_value > highest_H) {
      // Keep only the harmonic with the highest H_value.
      highest_H = H_value;
    }
  }

  return highest_H;
}

std::pair<double, double> HTestArray::chanceProb(double stat) const {
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

std::string HTestArray::getDescription() const {
  std::ostringstream os;
  os << "Type of test: H Test, " << m_max_harm << " maximum harmonics\n" <<
    "Probability distribution: H Test-specific";
  return os.str();
}

HTestArray::size_type HTestArray::size() const {
  return m_sine_cont.size();
}

void HTestArray::plot(const std::string & /* title */, size_type /* array_index */) const {
  // TODO: Implement this method, plotting sine and cosine component against the harmonic number.
  throw std::runtime_error("HTestArray::plot is not implemented yet.");
}
