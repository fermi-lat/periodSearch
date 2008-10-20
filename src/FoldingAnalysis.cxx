/** \file FoldingAnalysis.cxx
    \brief Implementation of FoldingAnalysis class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "FoldingAnalysis.h"
#include "PeriodicityTestArray.h"
#include "StatisticViewer.h"

FoldingAnalysis::FoldingAnalysis(PeriodicityTestArray * test_array, double center, double step, double epoch, double duration,
  const std::string & freq_unit): PeriodSearch(test_array->size(), freq_unit), m_test_array(test_array), m_step(step),
  m_epoch(epoch), m_fourier_res(0.) {
  // Get the array size.
  size_type num_trials = test_array->size();

  // Check central frequency of the scan.
  if (0. >= center) {
    std::ostringstream os;
    os << "Non-positive value is given for the central frequency of a periodicity search: " << center;
    throw std::runtime_error(os.str());
  }
  if (0. >= m_step) {
    std::ostringstream os;
    os << "Non-positive value is given for the frequency step in a periodicity search: " << m_step;
    throw std::runtime_error(os.str());
  }
  if (0  >= num_trials) {
    std::ostringstream os;
    os << "Non-positive number is given for the number of trial frequencies in a periodicity search: " << num_trials;
    throw std::runtime_error(os.str());
  }
  if (0. >= duration) {
    std::ostringstream os;
    os << "Non-positive value is given for the time span of the observation for a periodicity search: " << duration;
    throw std::runtime_error(os.str());
  }

  // Create vector containing the trial frequencies.
  size_type ii_cent = num_trials / 2;
  double min = center - ii_cent * m_step;

  // Check whether the step was too big, leading to a negative frequency.
  if (0. >= min) {
    std::ostringstream os;
    os << "Trial frequencies for a periodicity search include a non-positive value such as " << min;
    throw std::runtime_error(os.str());
  }

  // Step from minimum frequency on up, populating internal arrays.
  StatisticViewer & viewer = getViewer();
  for (size_type ii = 0; ii < num_trials; ++ii) {
    // Populating frequency array.
    StatisticViewer::data_type & freq = viewer.getData(0);
    freq[ii] = min + ii * m_step;
  }

  // Compute Fourier resolution.
  m_fourier_res = 1. / duration;

  // Add/modify plot title.
  viewer.setTitle("Folding Analysis: " + m_test_array->getTestName());

  // Set description of this period search.
  setDescription("Folding Analysis", m_fourier_res, m_step, m_test_array->getDescription());
}

void FoldingAnalysis::fill(double evt_time) {
  // Pseudocode taken from work by M. Hirayama, which may be seen at:
  // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/pc_chi2test.txt
  //     dt = evtime - epoch;
  double dt = evt_time - m_epoch;

  //     for (i=0; i<N_trial; i++) {
  //       phase = dt / period_array[i];
  //       phase -= floor(phase);
  //
  //       [fill one entry into hist[i] at "phase"];
  //     }

  // Iterate over the number of trial frequencies.
  StatisticViewer & viewer = getViewer();
  size_type num_trials = m_test_array->size();
  for (size_type ii = 0; ii < num_trials; ++ii) {
    // For each frequency, compute the phase.
    const StatisticViewer::data_type & freq = viewer.getData(0);
    double phase = dt * freq[ii];
    phase -= floor(phase);

    // Use this phase information to fill in the corresponding trial.
    m_test_array->fill(ii, phase);
  }
}

void FoldingAnalysis::computeStat() {
  // Prepare a returning array.
  StatisticViewer & viewer = getViewer();
  StatisticViewer::data_type & spec = viewer.getData(1);
  spec.assign(spec.size(), 0.);

  // Iterate over the number of trials.
  size_type num_trials = m_test_array->size();
  for (size_type ii = 0; ii < num_trials; ++ii) {
    spec[ii] = m_test_array->computeStat(ii);
  }
}

PeriodSearch::size_type FoldingAnalysis::computeNumIndepTrials(double min_freq, double max_freq) const {
  std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);

  // Reset min/max frequency if either bound was not explicitly specified (negative).
  const StatisticViewer & viewer = getViewer();
  const StatisticViewer::data_type & freq = viewer.getData(0);
  if (0. > min_freq) min_freq = freq[indices.first];
  if (0. > max_freq && indices.second > 0) max_freq = freq[indices.second - 1];

  //    N_Fourier = (stop - start) / Fourier_step
  size_type n_fourier = size_type(ceil(fabs((max_freq - min_freq)) / m_fourier_res));

  // Compute also the number of bins.
  size_type num_bins = indices.second >= indices.first ? indices.second - indices.first : 0;

  // Whichever is smaller is the number of independent trials.
  return (n_fourier < num_bins) ? n_fourier : num_bins;
}

std::pair<double, double> FoldingAnalysis::computeChanceProbOneTrial(double stat) const {
  // Delegate the computation of a chance probablity.
  return m_test_array->computeChanceProb(stat);
}
