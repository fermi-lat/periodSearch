/** \file Z2nTestArray.cxx
    \brief Implementation of Z2nTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <sstream>
#include <stdexcept>

#include "ChiSquaredProb.h"
#include "Z2nTestArray.h"

Z2nTestArray::Z2nTestArray(size_type array_size, data_type::size_type num_harmonics):
  m_num_harm(num_harmonics), m_sine_cont(array_size, data_type(num_harmonics, 0.)),
  m_cosine_cont(array_size, data_type(num_harmonics, 0.)), m_num_events(array_size, 0) {}

void Z2nTestArray::fill(double phase, size_type array_index) {
  // Define two pi (for convenience and clarity).
  static const double s_2pi = 2. * 4. * std::atan(1.0);

  // Get the storage for sine and consine component.
  data_type & sine_array = m_sine_cont.at(array_index);
  data_type & cosine_array = m_cosine_cont.at(array_index);

  // For each phase, the complex Fourier component is computed for each trial harmonic.
  for (size_type jj = 0; jj < m_num_harm; ++jj) {
    double phase_angle = s_2pi * (jj + 1) * phase;
    sine_array[jj] += std::sin(phase_angle);
    cosine_array[jj] += std::cos(phase_angle);
  }

  // Increment the number of events filled.
  ++(m_num_events.at(array_index));
}

double Z2nTestArray::testStat(size_type array_index) const {
  // Get the storage for sine and consine component.
  const data_type & sine_array = m_sine_cont.at(array_index);
  const data_type & cosine_array = m_cosine_cont.at(array_index);

  // Compute normalization.
  double fourier_norm = 2. / m_num_events.at(array_index);

  // Compute the Fourier powers and sum them up over harmonics.
  double summed_power = 0.;
  for (size_type jj = 0; jj < m_num_harm; ++jj) {
    summed_power += sine_array[jj] * sine_array[jj];
    summed_power += cosine_array[jj] * cosine_array[jj];
  }

  // Return normalized summed power.
  return summed_power * fourier_norm;
}

std::pair<double, double> Z2nTestArray::chanceProb(double stat) const {
  //      /* Leahy et al. 1983, ApJ 266, 160 */
  //      chance_prob = chi2prob(test_stat[imax], 2*N_harm) * N_Fourier;
  //      [where function chi2prob(chi2, dof) returns the chi-squared
  //       distribution for "dof" degrees of freedom, integrated from "chi2"
  //       to infinity];
  periodSearch::ChiSquaredProb prob(2 * m_num_harm);
  return prob(stat);
}

std::string Z2nTestArray::getDescription() const {
  std::ostringstream os;
  os << "Type of test: Z2n Test, " << m_num_harm << " harmonics\n" <<
    "Probability distribution: Chi-squared, " << 2 * m_num_harm << " degrees of freedom";
  return os.str();
}

Z2nTestArray::size_type Z2nTestArray::size() const {
  return m_sine_cont.size();
}

void Z2nTestArray::plot(const std::string & /* title */, size_type /* array_index */) const {
  // TODO: Implement this method, plotting sine and cosine component against the harmonic number.
  throw std::runtime_error("Z2nTestArray::plot is not implemented yet.");
}
