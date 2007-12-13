/** \file HTestArray.cxx
    \brief Implementation of HTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#include "HTestArray.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

HTestArray::HTestArray(size_type array_size, data_type::size_type max_harmonics): Z2nTestArray(array_size, max_harmonics) {
  // Set description of this statistical test.
  std::ostringstream os_cond;
  os_cond << max_harmonics << " maximum harmonics";
  std::ostringstream os_dist;
  os_dist << "H Test-specific";
  setDescription("H Test", os_cond.str(), os_dist.str());

  // Overwrite Y-label in the viewer.
  StatisticViewer & viewer = getViewer();
  viewer.setLabel(1, "CANDIDATE_VALUE");

  // Overwrite title in the viewer.
  viewer.setTitle("Candidates for H value");
}

double HTestArray::computeH(size_type array_index, data_type & candidate) const {
  // Compute the Fourier powers.
  data_type power;
  computePower(array_index, power);

  // Initialize output array.
  data_type::size_type array_size = power.size();
  candidate.resize(array_size);
  candidate.assign(array_size, 0.);

  // Compute candidates for the H value.
  double z2_value = 0.;
  for (data_type::size_type jj = 0; jj < array_size; ++jj) {
    z2_value += power[jj];
    candidate[jj] = z2_value - 4. * jj;
  }

  // Find the largest candidate.
  double H_value = 0.;
  for (data_type::const_iterator itor = candidate.begin(); itor != candidate.end(); ++itor) {
    if (*itor > H_value) H_value = *itor;
  }

  // Return the maximum value of H-value candidates.
  return H_value;
}

double HTestArray::computeStat(size_type array_index) const {
  // Compute and return the H value, without changing the viewer contents.
  data_type candidate;
  return computeH(array_index, candidate);
}

void HTestArray::updateViewer(size_type array_index) {
  // Compute H value, while keeping candidate H values in the viewer.
  StatisticViewer & viewer = getViewer();
  StatisticViewer::data_type & candidate = viewer.getData(1);
  double H_value = computeH(array_index, candidate);

  // Set caption to the viewer.
  viewer.setCaption(createSummary(H_value));
}

std::pair<double, double> HTestArray::computeChanceProb(double stat) const {
  /* De Jager et al. 1989, A&A 221, 180 */
  double lower_limit;
  double upper_limit;
  bool chance_prob_exact = true;

  if (stat <= 23.0) {
     double a = 0.9999755;
     double b = 0.39802;
     upper_limit = a * std::exp(-b * stat);
  } else if (stat < 50.0) {
     double c = 1.210597;
     double d = 0.45901;
     double e = 0.0022900;
     upper_limit = c * std::exp(-d * stat + e * stat * stat);
  } else {
     upper_limit = 4.0e-8; /* or less */
     chance_prob_exact = false;
  }

  upper_limit = upper_limit < 1. ? upper_limit : 1.;

  lower_limit = chance_prob_exact ? upper_limit : 0.;

  return std::make_pair(lower_limit, upper_limit);
}
