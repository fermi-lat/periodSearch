/** \file ChiSquaredTestArray.cxx
    \brief Implementation of ChiSquaredTestArray class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#include "ChiSquaredTestArray.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "ChiSquaredProb.h"
#include "StatisticViewer.h"

ChiSquaredTestArray::ChiSquaredTestArray(size_type array_size, data_type::size_type num_phase_bins):
  PeriodicityTestArray(2, num_phase_bins), m_num_phase_bins(num_phase_bins), m_curve_cont(array_size, data_type(num_phase_bins, 0)),
  m_num_events(array_size, 0) {
  // Check argument values.
  if (0 >= num_phase_bins) throw std::logic_error("ChiSquaredTestArray constructor was passed a non-positive number of phase bins");

  // Set description of this statistical test.
  std::ostringstream os_cond;
  os_cond << m_num_phase_bins << " phase bins";
  std::ostringstream os_dist;
  os_dist << "Chi-squared, " << m_num_phase_bins - 1 << " degrees of freedom";
  setDescription("Chi-squared Test", os_cond.str(), os_dist.str());

  // Set pulse phase values to the viewer.
  StatisticViewer & viewer = getViewer();
  StatisticViewer::data_type & phase = viewer.getData(0);
  for (StatisticViewer::data_type::size_type ii=0; ii < StatisticViewer::data_type::size_type(m_num_phase_bins); ++ii) {
    phase[ii] = (ii + 0.5) / m_num_phase_bins;
  }

  // Set default labels to the viewer.
  viewer.setLabel(0, "PULSE_PHASE");
  viewer.setLabel(1, "COUNTS");

  // Set default title to the viewer.
  viewer.setTitle("Folded Light Curve");
}

void ChiSquaredTestArray::fill(size_type array_index, double phase) {
  // Bin phase; it runs from [0, 1), so multiply by the number of bins to determine
  // bin for this phase.
  size_type bin_id = size_type(phase * m_num_phase_bins);

  // Round up the bin ID, just in case a give phase is out of range.
  bin_id %= m_num_phase_bins;

  // Increment the count in that bin.
  ++(m_curve_cont.at(array_index)[bin_id]);

  // Increment the number of events filled.
  ++(m_num_events.at(array_index));
}

double ChiSquaredTestArray::computeStat(size_type array_index) const {
  // Check the number of events.
  long num_events = m_num_events.at(array_index);
  if (num_events == 0) {
    std::ostringstream os;
    os << "ChiSquaredTestArray::computeStat was called with no events filled for test #" << array_index << std::endl;
    throw std::logic_error(os.str());
  }

  // Compute average count rate.
  double avg = double(num_events) / m_num_phase_bins;

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

void ChiSquaredTestArray::updateViewer(size_type array_index) {
  // Copy a light curve to the viewer.
  const data_type & light_curve = m_curve_cont.at(array_index);
  StatisticViewer & viewer = getViewer();
  StatisticViewer::data_type & viewer_curve = viewer.getData(1);
  for (StatisticViewer::data_type::size_type ii=0; ii < StatisticViewer::data_type::size_type(m_num_phase_bins); ++ii) {
    viewer_curve[ii] = light_curve[ii];
  }

  // Compute the S-value.
  double S_value = computeStat(array_index);

  // Set caption to the viewer.
  viewer.setCaption(createSummary(S_value));
}

std::pair<double, double> ChiSquaredTestArray::computeChanceProb(double stat) const {
  //    /* Leahy et al. 1983, ApJ 266, 160 */
  //    chance_prob = chi2prob(S_value[imax], N_bin-1) * N_Fourier;
  //    [where function chi2prob(chisq, dof) returns the chi-squared
  //     distribution for "dof" degrees of freedom, integrated from "chisq"
  //     to infinity];
  periodSearch::ChiSquaredProb prob(m_num_phase_bins - 1);
  return prob(stat);
}

ChiSquaredTestArray::size_type ChiSquaredTestArray::size() const {
  return m_curve_cont.size();
}
