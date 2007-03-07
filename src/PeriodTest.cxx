/** \file PeriodTest.cxx
    \brief Implementation of PeriodTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "periodSearch/PeriodTest.h"

namespace periodSearch {

  PeriodTest::PeriodTest(double center, double step, size_type num_trials, double epoch, size_type array_size, double duration):
    PeriodSearch(num_trials), m_trial_hist(num_trials, std::vector<std::complex<double> >(array_size, 0.)),
    m_step(step), m_epoch(epoch), m_fourier_res(0.), m_num_events(0) {
    // Make certain there is no error in the input.
    if (0. >= center) throw std::logic_error("PeriodTest constructor was passed a non-positive center");
    if (0. >= m_step) throw std::logic_error("PeriodTest constructor was passed a non-positive step");
    if (0 == num_trials) throw std::logic_error("PeriodTest constructor was passed a non-positive num_trials");
    if (0. >= duration) throw std::logic_error("PeriodTest constructor was passed a non-positive duration");
    if (0 == array_size) throw std::logic_error("PeriodTest constructor was passed a non-positive array_size");

    // Create vector containing the trial frequencies.
    size_type ii_cent = num_trials / 2;
    double min = center - ii_cent * m_step;

    // Check whether the step was too big, leading to a negative frequency.
    if (0. >= min) throw std::logic_error("PeriodTest constructor computed a non-positive trial frequency");

    // Step from minimum frequency on up, populating frequency array.
    for (size_type ii = 0; ii < num_trials; ++ii) {
      m_freq[ii] = min + ii * m_step;
    }

    // Compute Fourier resolution.
    m_fourier_res = 1. / duration;
  }

  void PeriodTest::fill(double evt_time) {
    // Pseudocode taken from work by M. Hirayama, which may be seen at:
    // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/pc_chi2test.txt
    //     dt = evtime - epoch;
    double dt = evt_time - m_epoch;

    //     N_event++;
    ++m_num_events;

    //     for (i=0; i<N_trial; i++) {
    //       phase = dt / period_array[i];
    //       phase -= floor(phase);
    //
    //       [fill one entry into hist[i] at "phase"];
    //     }

    // Iterate over the number of trial frequencies.
    size_type num_trials = m_trial_hist.size();
    for (size_type ii = 0; ii < num_trials; ++ii) {
      // For each frequency, compute the phase.
      double phase = dt * m_freq[ii];
      phase -= floor(phase);

      // Use this phase information to fill in the corresponding trial.
      fillOneTrial(phase, m_trial_hist[ii]);
    }
  }

  PeriodSearch::size_type PeriodTest::numIndepTrials(double min_freq, double max_freq) const {
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

  std::string PeriodTest::getDescription() const {
    std::ostringstream os;
    os << "Search Type: Folding Analysis\n"
       << "Fourier Resolution: " << m_fourier_res << " Hz\n"
       << "Sampling Frequency: " << m_step << " Hz";
    return os.str();
  }

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodTest & test) {
    return test.write(os);
  }
}
