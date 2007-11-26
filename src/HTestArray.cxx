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

double HTestArray::testStat(size_type array_index) const {
  // Compute the Fourier powers.
  data_type power;
  computePower(array_index, power);

  // Compute H value.
  double highest_H = 0.;
  double z2_value = 0.;
  // Iterate over bins in each trial.
  for (size_type jj = 0; jj < size_type(power.size()); ++jj) {
    // Compute coefficient of power spectrum for each harmonic.
    z2_value += power[jj];
    double H_value = z2_value - 4. * jj;

    // Keep only the harmonic with the highest H_value.
    if (H_value > highest_H) highest_H = H_value;
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

std::pair<std::vector<double>, std::vector<double> > HTestArray::getPlotData(size_type array_index) const {
  // TODO: Implement this method, plotting sine and cosine component against the harmonic number.
  throw std::runtime_error("HTestArray::getPlotData is not implemented yet.");
}

std::pair<std::string, std::string> HTestArray::getPlotLabel() const {
  return std::make_pair(std::string("Harmonic Number"), std::string("H"));
}

std::string HTestArray::getPlotTitle() const {
  return "Candidate H values";
}
