/** \file FoldingAnalysis.cxx
    \brief Implementation of FoldingAnalysis class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <cmath>
#include <iostream>
//#include <memory>
#include <sstream>
#include <stdexcept>

//#include "st_graph/Axis.h"
//#include "st_graph/Engine.h"
//#include "st_graph/IPlot.h"
//#include "st_graph/Sequence.h"

#include "FoldingAnalysis.h"
#include "PeriodicityTest.h"
#include "PeriodicityTestArray.h"

namespace periodSearch {

  FoldingAnalysis::FoldingAnalysis(PeriodicityTestArray * test_array, double center, double step, double epoch, double duration):
    PeriodSearch(test_array->size()), m_test_array(test_array), m_step(step), m_epoch(epoch), m_fourier_res(0.) {
    // Get the array size.
    size_type num_trials = test_array->size();

    // Make certain there is no error in the input.
    if (0. >= center) throw std::logic_error("FoldingAnalysis constructor was passed a non-positive center");
    if (0. >= m_step) throw std::logic_error("FoldingAnalysis constructor was passed a non-positive step");
    if (0  >= num_trials) throw std::logic_error("FoldingAnalysis constructor was passed a non-positive num_trials");
    if (0. >= duration) throw std::logic_error("FoldingAnalysis constructor was passed a non-positive duration");

    // Create vector containing the trial frequencies.
    size_type ii_cent = num_trials / 2;
    double min = center - ii_cent * m_step;

    // Check whether the step was too big, leading to a negative frequency.
    if (0. >= min) throw std::logic_error("FoldingAnalysis constructor computed a non-positive trial frequency");

    // Step from minimum frequency on up, populating internal arrays.
    for (size_type ii = 0; ii < num_trials; ++ii) {
      // Populating frequency array.
      m_freq[ii] = min + ii * m_step;
    }

    // Compute Fourier resolution.
    m_fourier_res = 1. / duration;
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
    size_type num_trials = m_test_array->size();
    for (size_type ii = 0; ii < num_trials; ++ii) {
      // For each frequency, compute the phase.
      double phase = dt * m_freq[ii];
      phase -= floor(phase);

      // Use this phase information to fill in the corresponding trial.
      m_test_array->fill(phase, ii);
    }
  }

  const std::vector<double> & FoldingAnalysis::computeStats() {
    // Prepare a returning array.
    m_spec.assign(m_spec.size(), 0.);

    // Iterate over the number of trials.
    size_type num_trials = m_test_array->size();
    for (size_type ii = 0; ii < num_trials; ++ii) {
      m_spec[ii] = m_test_array->testStat(ii);
    }

    // Return the result.
    return m_spec;
  }

  PeriodSearch::size_type FoldingAnalysis::numIndepTrials(double min_freq, double max_freq) const {
    std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);

    // Reset min/max frequency if either bound was not explicitly specified (negative).
    if (0. > min_freq) min_freq = m_freq[indices.first];
    if (0. > max_freq && indices.second > 0) max_freq = m_freq[indices.second - 1];

    //    N_Fourier = (stop - start) / Fourier_step
    size_type n_fourier = size_type(ceil(fabs((max_freq - min_freq)) / m_fourier_res));

    // Compute also the number of bins.
    size_type num_bins = indices.second >= indices.first ? indices.second - indices.first : 0;

    // Whichever is smaller is the number of independent trials.
    return (n_fourier < num_bins) ? n_fourier : num_bins;
  }

  std::pair<double, double> FoldingAnalysis::chanceProbOneTrial(double stat) const {
    // Delegate the computation of a chance probablity.
    return m_test_array->chanceProb(stat);
  }

  std::string FoldingAnalysis::getDescription() const {
    // Write out common parameters.
    std::ostringstream os;
    os << "Search Type: Folding Analysis\n"
       << "Fourier Resolution: " << m_fourier_res << " Hz\n"
       << "Sampling Frequency: " << m_step << " Hz";

    // Write out the test-specific output.
    os << "\n" << m_test_array->getDescription();

    // Return the string.
    return os.str();
  }
}
