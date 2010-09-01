/** \file FourierAnalysis.cxx
    \brief Implementation of FourierAnalysis class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "ChiSquaredProb.h"
#include "FourierAnalysis.h"
#include "StatisticViewer.h"

#include "fftw/fftw3.h"

FourierAnalysis::FourierAnalysis(double t_start, double t_stop, double width, PeriodSearch::size_type num_bins,
  const std::string & freq_unit, int /* num_events */):
  PeriodSearch(num_bins / 2 + 1, freq_unit), m_index(), m_t_start(t_start), m_t_stop(t_stop), m_width(width),
  m_fourier_res(0.), m_num_segments(0) , m_num_bins(num_bins) {
  // Check time interval.
  if (t_start > t_stop) {
    std::ostringstream os;
    os << "Bad time interval for Discrete Fast Fourier Transform (DFFT): [" << t_start << ", " << t_stop <<
        "] (start time is later than stop time)";
    throw std::runtime_error(os.str());
  }

  // Check width of time bin.
  if (0. >= width) {
    std::ostringstream os;
    os << "Non-positive value is given for the width of time bin: " << width;
    throw std::runtime_error(os.str());
  }

  // Check the number of time bins to be FFT'ed at a time.
  if (0  >= num_bins) {
    std::ostringstream os;
    os << "Non-positive number is given for the number of time bins to be Fourier-transfered at a time: " << num_bins;
    throw std::runtime_error(os.str());
  }
  // m_index.reserve(num_events);

  // Set up frequency array.
  StatisticViewer & viewer = getViewer();
  StatisticViewer::data_type & freq = viewer.getData(0);
  m_fourier_res = 1. / (m_width * m_num_bins);
  for (std::size_t ii = 0; ii < freq.size(); ++ii) {
    freq[ii] = ii * m_fourier_res;
  }

  // Add/modify plot title.
  viewer.setTitle("Fourier Analysis: Power Spectrum");

  // Set description of this period search.
  // Note: Description changes with events being filled, as the number of segments grows with events.
  //       Give generic information here, and update with more information at FFT computation.
  std::ostringstream os;
  os << "Data Binning: " << m_num_bins << " time bins to be Fourier-transfered at a time";
  setDescription("Fourier Analysis", m_fourier_res, m_fourier_res, os.str());
}

void FourierAnalysis::fill(double evt_time) {
  if (m_t_start <= evt_time && evt_time <= m_t_stop) {
    // Compute index as if we had one huge array with all the segments.
    size_type global_idx = size_type(std::floor((evt_time - m_t_start) / m_width) + .5);
    // Determine segment for this event.
    size_type segment_idx = global_idx / m_num_bins;
    // Determine bin within the segment for this event.
    size_type bin_idx = global_idx % m_num_bins;

    // Keep track of the last segment seen.
    m_num_segments = std::max(m_num_segments, segment_idx + 1);
    // Track the segment and bin index in the m_index member.
    m_index.insert(std::make_pair(segment_idx, bin_idx));
  }
}

void FourierAnalysis::computeStat() {
  // Set description of this period search.
  // Note: This is the only place to correctly display the number of segments because it grows with events.
  std::ostringstream os;
  os << "Data Binning: " << m_num_segments << " segments with " << m_num_bins << " time bins in each segment" << std::endl
     << "Probability Distribution: Chi Squared with " << 2 * m_num_segments << " degrees of freedom";
  setDescription("Fourier Analysis", m_fourier_res, m_fourier_res, os.str());

  // Prepare statistic viewer.
  StatisticViewer & viewer = getViewer();
  const StatisticViewer::data_type & freq = viewer.getData(0);
  StatisticViewer::data_type & spec = viewer.getData(1);

  // Prepare variables for FFTW.
  double * in = 0;
  fftw_complex * out = 0;
  fftw_plan p = 0;
  std::size_t num_cpx_elements = freq.size();
  std::size_t num_dbl_elements = 2 * num_cpx_elements;

  // Allocate array for fftw input/output.
  in = (double *) fftw_malloc(sizeof(double) * num_dbl_elements);
  if (0 == in) throw std::runtime_error("Could not allocate a memory space for Discrete Fast Fourier Transform (DFFT)");

  // Perform transform in place.
  out = (fftw_complex *) in;

  // Set up to perform fftw-style transform.
  // p = fftw_plan_dft_1d(m_num_bins, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  p = fftw_plan_dft_r2c_1d(m_num_bins, in, out, FFTW_ESTIMATE);

  // Iterate over segments.
  for (std::size_t seg_idx = 0; seg_idx < m_num_segments; ++seg_idx) {
    // Initialize array to be all zero.
    std::memset(in, '\0', sizeof(double) * num_dbl_elements);

    // Populate the real portion of the input array.
    double num_events = 0.;
    for (index_map_type::iterator itor = m_index.lower_bound(seg_idx); itor != m_index.upper_bound(seg_idx); ++itor) {
      // Increment the number of counts observed in this bin.
      ++in[itor->second];
      // Increment the total number of events observed in this segment.
      ++num_events;
    }

    if (num_events == 0) {
      // Pack zeros into the array holding the power density.
      for (std::size_t ii = 0; ii < num_cpx_elements; ++ii) spec[ii] = 0.;

    } else {
      // Shift origin by the average number of events per bin.
      double events_per_bin = num_events / m_num_bins;
      for (std::size_t ii = 0; ii < m_num_bins; ++ii) {
        //in[ii][0] -= events_per_bin;
        in[ii] -= events_per_bin;
      }

      // Do the transformation.
      fftw_execute(p);

      // Pack results into the array holding the power density.
      for (std::size_t ii = 0; ii < num_cpx_elements; ++ii) {
        const double & real = out[ii][0];
        const double & imag = out[ii][1];
        spec[ii] += (real * real + imag * imag) * 2. / num_events;
      }
    }
  }

  fftw_destroy_plan(p);
  fftw_free(in);
}

PeriodSearch::size_type FourierAnalysis::computeNumIndepTrials(double min_freq, double max_freq) const {
  std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
  size_type num_indep_trials = indices.second >= indices.first ? indices.second - indices.first : 0;
  return num_indep_trials;
}

std::pair<double,double> FourierAnalysis::computeChanceProbOneTrial(double stat) const {
  ChiSquaredProb prob(2 * m_num_segments);
  return prob(stat);
}
