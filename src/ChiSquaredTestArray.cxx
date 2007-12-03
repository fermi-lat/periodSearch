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
  m_num_phase_bins(num_phase_bins), m_curve_cont(array_size, data_type(num_phase_bins, 0)), m_num_events(array_size, 0),
  m_X_data(num_phase_bins, 0.), m_Y_data(num_phase_bins, 0.) {}

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

double ChiSquaredTestArray::testStat(size_type array_index) const {
  // Compute average count rate.
  double avg = double(m_num_events.at(array_index)) / m_num_phase_bins;

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
  os << "Type of test: " << getTestName() << ", " << m_num_phase_bins << " phase bins\n"
     << "Probability distribution: Chi-squared, " << m_num_phase_bins - 1 << " degrees of freedom";
  return os.str();
}

ChiSquaredTestArray::size_type ChiSquaredTestArray::size() const {
  return m_curve_cont.size();
}

void ChiSquaredTestArray::getPlotData(size_type array_index, std::vector<double> & phase, std::vector<double> & count) const {
  // Initialize the output arrays.
  phase.resize(m_num_phase_bins);
  phase.assign(m_num_phase_bins, 0.);
  count.resize(m_num_phase_bins);
  count.assign(m_num_phase_bins, 0.);

  // Create a light curve to plot.
  const data_type & curve = m_curve_cont.at(array_index);
  for (std::vector<double>::size_type ii=0; ii < std::vector<double>::size_type(m_num_phase_bins); ++ii) {
    phase[ii] = (ii + 0.5) / m_num_phase_bins;
    count[ii] = curve[ii];
  }
}

void ChiSquaredTestArray::getPlotLabel(std::string & x_label, std::string & y_label) const {
  x_label = "Pulse Phase";
  y_label = "Counts";
}

std::string ChiSquaredTestArray::getPlotTitle() const {
  return "Folded Light Curve";
}

std::string ChiSquaredTestArray::getTestName() const {
  return "Chi-squared Test";
}

StatisticViewer ChiSquaredTestArray::getViewer() const {
  // Create a viewer object to return.
  StatisticViewer viewer(2, m_num_phase_bins);

  // Copy data to the viewer.
  viewer.setData(0, m_X_data.begin(), true);
  viewer.setData(1, m_Y_data.begin(), true);

  // Set label to the viewer.
  viewer.setLabel(0, "Pulse Phase");
  viewer.setLabel(1, "Counts");

  // Set title and caption to the viewer.
  viewer.setTitle("Folded Light Curve");
  viewer.setCaption(getDescription());

  // Return the viewer.
  return viewer;
}

void ChiSquaredTestArray::computeViewerData(size_type array_index) {
  // Create a light curve to plot.
  const data_type & curve = m_curve_cont.at(array_index);
  for (std::vector<double>::size_type ii=0; ii < StatisticViewer::data_type::size_type(m_num_phase_bins); ++ii) {
    m_X_data[ii] = (ii + 0.5) / m_num_phase_bins;
    m_Y_data[ii] = curve[ii];
  }
}
