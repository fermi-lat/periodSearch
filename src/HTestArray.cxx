/** \file HTestArray.cxx
    \brief Implementation of HTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <sstream>
#include <stdexcept>

#include "HTestArray.h"

HTestArray::HTestArray(size_type array_size, data_type::size_type max_harmonics): Z2nTestArray(array_size, max_harmonics),
  m_max_harm(max_harmonics) {}

void HTestArray::computeCandidate(data_type & power, data_type & H_candidate) const {
  // Initialize output array.
  data_type::size_type array_size = power.size();
  H_candidate.resize(array_size);
  H_candidate.assign(array_size, 0.);

  // Compute candidates for the H value.
  double z2_value = 0.;
  for (data_type::size_type jj = 0; jj < array_size; ++jj) {
    z2_value += power[jj];
    H_candidate[jj] = z2_value - 4. * jj;
  }
}

double HTestArray::testStat(size_type array_index) const {
  // Compute the Fourier powers.
  data_type power;
  computePower(array_index, power);

  // Compute H values.
  data_type H_candidate;
  computeCandidate(power, H_candidate);

  // Keep only the harmonic with the highest H-value.
  double highest_H = 0.;
  for (data_type::const_iterator itor = H_candidate.begin(); itor != H_candidate.end(); ++itor) {
    if (*itor > highest_H) highest_H = *itor;
  }

  // Return the maximum value of H-value candidates.
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
  os << "Type of test: " << getTestName() << ", " << m_max_harm << " maximum harmonics\n" <<
    "Probability distribution: H Test-specific";
  return os.str();
}

std::string HTestArray::getTestName() const {
  return "H Test";
}

StatisticViewer & HTestArray::getViewer(size_type array_index) {
  // Let Z2nTestArray compute the Fourier powers.
  StatisticViewer & viewer = Z2nTestArray::getViewer(array_index);

  // Copy the computed Fourier powers to a local variable for safety.
  StatisticViewer::data_type power = viewer.getData(1);

  // Convert the Fourier powers to H-value candidates.
  StatisticViewer::data_type & candidate = viewer.getData(1);
  computeCandidate(power, candidate);

  // Overwrite Y-label in the viewer.
  viewer.setLabel(1, "CANDIDATE_VALUE");

  // Overwrite title in the viewer.
  viewer.setTitle("Candidates for H value");

  // Return the viewer.
  return viewer;
}
