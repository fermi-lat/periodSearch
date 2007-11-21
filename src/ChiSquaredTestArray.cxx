/** \file ChiSquaredTestArray.cxx
    \brief Implementation of ChiSquaredTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#include <sstream>
#include <stdexcept>

#include "ChiSquaredProb.h"
#include "ChiSquaredTestArray.h"

ChiSquaredTestArray::ChiSquaredTestArray(size_type array_size, data_type::size_type num_phase_bins):
  m_num_phase_bins(num_phase_bins), m_curve_cont(array_size, data_type(num_phase_bins, 0)), m_num_events(0) {}

void ChiSquaredTestArray::fill(double phase, size_type array_index) {
  // Bin phase; it runs from [0, 1), so multiply by the number of bins to determine
  // bin for this phase.
  size_type bin_id = size_type(phase * m_num_phase_bins);

  // Round up the bin ID, just in case a give phase is out of range.
  bin_id %= m_num_phase_bins;

  // Increment the count in that bin.
  ++(m_curve_cont.at(array_index)[bin_id]);

  // Increment the number of events filled.
  ++m_num_events;
}

double ChiSquaredTestArray::testStat(size_type array_index) const {
  // Compute average count rate.
  double avg = double(m_num_events) / m_num_phase_bins;

  // Compute S-value.
  double S_value = 0.;
  const data_type & light_curve = m_curve_cont.at(array_index);
  for (data_type::const_iterator itor = light_curve.begin(); itor != light_curve.end(); ++itor) {
    // Compute deviation. Imaginary part of trial container is ignored.
    double dev = *itor - avg;

    // Sum squares of deviation divided by the average value.
    S_value += dev * dev / avg;
  }

  // Return S-value.
  return S_value;
}

std::pair<double, double> ChiSquaredTestArray::chanceProb(double stat) const {
  //    /* Leahy et al. 1983, ApJ 266, 160 */
  //    chance_prob = chi2prob(S_value[imax], N_bin-1) * N_Fourier;
  //    [where function chi2prob(chisq, dof) returns the chi-squared
  //     distribution for "dof" degrees of freedom, integrated from "chisq"
  //     to infinity];
  periodSearch::ChiSquaredProb prob(m_num_phase_bins - 1);
  return prob(stat);
}

std::string ChiSquaredTestArray::getDescription() const {
  std::ostringstream os;
  os << "Type of test: CHI2 Test, " << m_num_phase_bins << " phase bins\n"
     << "Probability distribution: Chi-squared, " << m_num_phase_bins - 1 << " degrees of freedom";
  return os.str();
}

ChiSquaredTestArray::size_type ChiSquaredTestArray::size() const {
  return m_curve_cont.size();
}