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

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "ChiSquaredProb.h"

ChiSquaredTestArray::ChiSquaredTestArray(size_type array_size, data_type::size_type num_phase_bins):
  m_num_phase_bins(num_phase_bins), m_curve_cont(array_size, data_type(num_phase_bins, 0)), m_num_events(array_size, 0) {}

void ChiSquaredTestArray::fill(double phase, size_type array_index) {
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
  os << "Type of test: CHI2 Test, " << m_num_phase_bins << " phase bins\n"
     << "Probability distribution: Chi-squared, " << m_num_phase_bins - 1 << " degrees of freedom";
  return os.str();
}

ChiSquaredTestArray::size_type ChiSquaredTestArray::size() const {
  return m_curve_cont.size();
}

void ChiSquaredTestArray::plot(const std::string & title, size_type array_index) const {
  using namespace st_graph;

  // Create a light curve to plot.
  const data_type & curve = m_curve_cont.at(array_index);
  typedef std::vector<double> hist_type;
  hist_type phase_value(m_num_phase_bins, 0.);
  hist_type light_curve(m_num_phase_bins, 0.);
  for (hist_type::size_type ii=0; ii < hist_type::size_type(m_num_phase_bins); ++ii) {
    phase_value[ii] = (ii + 0.5) / m_num_phase_bins;
    light_curve[ii] = curve[ii];
  }

  try {
    // Get graphics engine to set up graph.
    Engine & engine(Engine::instance());

    // Typedef for readability.
    typedef st_graph::ValueSequence<std::vector<double>::const_iterator> ValueSeq_t;

    // TODO Add text output (statistics, etc.)to a text box on the plot, and/or in a GUI output window.
    std::auto_ptr<IPlot> plot(engine.createPlot(title, 800, 600, "hist",
      ValueSeq_t(phase_value.begin(), phase_value.end()),
      ValueSeq_t(light_curve.begin(), light_curve.end())));

    // Set axes titles.
    std::vector<Axis> & axes(plot->getAxes());
    axes[0].setTitle("Pulse Phase");
    axes[1].setTitle("Counts");

    // Display plot.
    engine.run();

  } catch (const std::exception & x) {
    std::cerr << x.what() << std::endl;
    std::cerr << "Warning: ChiSquaredTestArray::plot could not display plot." << std::endl;
  }
}
