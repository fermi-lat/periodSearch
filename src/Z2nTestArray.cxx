/** \file Z2nTestArray.cxx
    \brief Implementation of Z2nTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#include "Z2nTestArray.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "ChiSquaredProb.h"
#include "StatisticViewer.h"

Z2nTestArray::Z2nTestArray(size_type array_size, data_type::size_type num_harmonics):
  PeriodicityTestArray(2, num_harmonics), m_num_harm(num_harmonics), m_sine_cont(array_size, data_type(num_harmonics, 0.)),
  m_cosine_cont(array_size, data_type(num_harmonics, 0.)), m_num_events(array_size, 0) {
  // Set description of this statistical test.
  std::ostringstream os_cond;
  os_cond << m_num_harm << " harmonics";
  std::ostringstream os_dist;
  os_dist << "Chi-squared, " << 2 * m_num_harm << " degrees of freedom";
  setDescription("Z2n Test", os_cond.str(), os_dist.str());

  // Set harmonic numbers, starting with one (1).
  double harmonic_number = 1.;
  StatisticViewer::data_type & harmonic = m_viewer.getData(0);
  for (std::vector<double>::iterator itor = harmonic.begin(); itor != harmonic.end(); ++itor) *itor = harmonic_number++;

  // Set default labels to the viewer.
  m_viewer.setLabel(0, "HARMONIC_NUMBER");
  m_viewer.setLabel(1, "POWER");

  // Set default title to the viewer.
  m_viewer.setTitle("Fourier Powers");
}

void Z2nTestArray::fill(size_type array_index, double phase) {
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

void Z2nTestArray::computePower(size_type array_index, data_type & power) const {
  // Initialize the container of the Fourier powers to return.
  power.resize(m_num_harm);
  power.assign(m_num_harm, 0.);

  // Get the storage for sine and consine component.
  const data_type & sine_array = m_sine_cont.at(array_index);
  const data_type & cosine_array = m_cosine_cont.at(array_index);

  // Compute normalization.
  double fourier_norm = 2. / m_num_events.at(array_index);

  // Compute the Fourier powers.
  for (size_type jj = 0; jj < m_num_harm; ++jj) {
    power[jj] = fourier_norm * (sine_array[jj] * sine_array[jj] + cosine_array[jj] * cosine_array[jj]);
  }
}

double Z2nTestArray::computeZ2n(size_type array_index, data_type & power) const {
  // Compute the Fourier powers.
  computePower(array_index, power);

  // Sum up the Fourier powers over harmonics.
  double summed_power = 0.;
  for (data_type::const_iterator itor = power.begin(); itor != power.end(); ++itor) summed_power += *itor;

  // Return the summed power.
  return summed_power;
}

double Z2nTestArray::computeStat(size_type array_index) const {
  // Compute and return Z2n-value, without changing the viewer contents.
  data_type power;
  return computeZ2n(array_index, power);
}

void Z2nTestArray::updateViewer(size_type array_index) {
  // Compute Z2n-value, while keeping the Fourier powers in the viewer.
  StatisticViewer::data_type & power = m_viewer.getData(1);
  double Z2n_value = computeZ2n(array_index, power);

  // Set caption to the viewer.
  m_viewer.setCaption(createSummary(Z2n_value));
}

std::pair<double, double> Z2nTestArray::computeChanceProb(double stat) const {
  //      /* Leahy et al. 1983, ApJ 266, 160 */
  //      chance_prob = chi2prob(test_stat[imax], 2*N_harm) * N_Fourier;
  //      [where function chi2prob(chi2, dof) returns the chi-squared
  //       distribution for "dof" degrees of freedom, integrated from "chi2"
  //       to infinity];
  periodSearch::ChiSquaredProb prob(2 * m_num_harm);
  return prob(stat);
}

Z2nTestArray::size_type Z2nTestArray::size() const {
  return m_sine_cont.size();
}
