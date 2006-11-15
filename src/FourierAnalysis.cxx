/** \file FourierAnalysis.cxx
    \brief Implementation of FourierAnalysis class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "FourierAnalysis.h"

#include "fftw/fftw3.h"

namespace periodSearch {

  FourierAnalysis::FourierAnalysis(double t_start, double t_stop, double width, size_type num_bins, int /* num_events */):
    m_freq(num_bins / 2 + 1), m_spec(num_bins / 2 + 1), m_index(), m_t_start(t_start), m_t_stop(t_stop),
    m_width(width), m_num_segments(0) , m_num_bins(num_bins) {
    if (t_start > t_stop) throw std::runtime_error("FourierAnalysis: start time is > stop time");
    // m_index.reserve(num_events);

    // Set up frequency array.
    double freq_step = 1. / (m_width * m_num_bins);
    for (size_t ii = 0; ii < m_freq.size(); ++ii) {
      m_freq[ii] = ii * freq_step;
    }
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

  const std::vector<double> & FourierAnalysis::computeStats() {
    double * in = 0;
    fftw_complex * out = 0;
    fftw_plan p = 0;
    size_t num_cpx_elements = m_freq.size();
    size_t num_dbl_elements = 2 * num_cpx_elements;

    // Allocate array for fftw input/output.
    in = (double *) fftw_malloc(sizeof(double) * num_dbl_elements);
    if (0 == in) throw std::runtime_error("FourierAnalysis could not allocate array");

    // Perform transform in place.
    out = (fftw_complex *) in;

    // Set up to perform fftw-style transform.
    // p = fftw_plan_dft_1d(m_num_bins, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    p = fftw_plan_dft_r2c_1d(m_num_bins, in, out, FFTW_ESTIMATE);

    // Iterate over segments.
    for (size_t seg_idx = 0; seg_idx < m_num_segments; ++seg_idx) {
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

      // Shift origin by the average number of events per bin.
      double events_per_bin = num_events / m_num_bins;
      for (size_t ii = 0; ii < m_num_bins; ++ii) {
        //in[ii][0] -= events_per_bin;
        in[ii] -= events_per_bin;
      }

      // Do the transformation.
      fftw_execute(p);

      // Pack results into the array holding the power density.
      for (size_t ii = 0; ii < num_cpx_elements; ++ii) {
        const double & real = out[ii][0];
        const double & imag = out[ii][1];
        m_spec[ii] += real * real + imag * imag * 2. / num_events;
      }
    }

    fftw_destroy_plan(p);
    fftw_free(in);

    return m_spec;
  }

  void FourierAnalysis::plot(const std::string & title, const std::string & freq_unit, double min_freq, double max_freq) const {
    using namespace st_graph;

    try {
      // Display value of maximum frequency/statistic in title.
      std::ostringstream os;
      std::pair<double, double> max = findMaxRange(min_freq, max_freq);
      std::pair<double, double> chance_prob = chanceProb(max.second);

      os << title << ", max at: " << max.first << ", stat: " << max.second;

      // Massage display: if difference between min and max is small enough just use max.
      os.setf(std::ios::scientific);
      os.precision(2); // 3 digits -> < 1. e -4. limit in next line.
      if ((chance_prob.second - chance_prob.first) / chance_prob.second < 1.e-4)
        os << ", chance prob: " << chance_prob.second;
      else
        os << ", chance prob < " << chance_prob.second;

      // Get graphics engine to set up graph.
      Engine & engine(Engine::instance());

      // Typedef for readability.
      typedef st_graph::ValueSequence<std::vector<double>::const_iterator> ValueSeq_t;

      // Impose range limits.
      std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
      size_type begin_index = indices.first;
      size_type end_index = indices.second;

      // Create plot, using m_freq as x, and m_spec as y.
      std::auto_ptr<IPlot> plot(engine.createPlot(os.str(), 800, 600, "hist",
        ValueSeq_t(m_freq.begin() + begin_index, m_freq.begin() + end_index),
        ValueSeq_t(m_spec.begin() + begin_index, m_spec.begin() + end_index)));

      // Set axes titles.
      std::vector<Axis> & axes(plot->getAxes());
      axes[0].setTitle("Frequency " + freq_unit);
      axes[1].setTitle("Test Statistic");

      // Display plot.
      engine.run();

    } catch (const std::exception & x) {
      std::cerr << x.what() << std::endl;
      std::cerr << "Warning: FourierAnalysis::plotStats could not display plot." << std::endl;
    }
  }

  void FourierAnalysis::plotStats(const std::string & title, const std::string & freq_unit) const {
    plot(title, freq_unit);
  }

  std::pair<double, double> FourierAnalysis::findMax() const {
    return findMaxRange();
  }

  std::pair<double, double> FourierAnalysis::findMaxRange(double min_freq, double max_freq) const {
    long max_idx = -1;
    double max = 0.;

    // Impose range limits.
    std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
    size_type begin_index = indices.first;
    size_type end_index = indices.second;

    for (unsigned long ii = begin_index; ii < end_index; ++ii) {
      // If the value is larger than the current maximum, replace it.
      if (m_spec[ii] > max) {
        max = m_spec[ii];
        max_idx = ii;
      }
    }

    // Make sure a valid maximum was found.
    if (0 > max_idx) {
      std::ostringstream os;
      os << "FourierAnalysis::findMaxRange cannot find any trial frequency in range [" << min_freq << ", " << max_freq << "]";
      throw std::runtime_error(os.str());
    }
    return std::pair<double, double>(m_freq[max_idx], max);
  }

  std::pair<double,double> FourierAnalysis::chanceProb(double /* stat */) const {
    // TODO Use stat to compute chanceProb.
    return std::pair<double, double>(.001, .001);
  }

  st_stream::OStream & FourierAnalysis::write(st_stream::OStream & os) const {
    return writeRange(os);
  }

  st_stream::OStream & FourierAnalysis::writeRange(st_stream::OStream & os, double min_freq, double max_freq) const {
    using namespace std;

    // Get info about the maximum.
    std::pair<double, double> max = findMaxRange(min_freq, max_freq);

    // Chance probability.
    std::pair<double, double> chance_prob = chanceProb(max.second);

    // Save current precision.
    int save_precision = os.precision();

    os.precision(15);

    // Write out the results.
    os << "Maximum at: " << max.first << std::endl << "Statistic: " << max.second << std::endl;
    os << "Chance probability range: (" << chance_prob.first << ", " << chance_prob.second << ")" << std::endl;
    os << "Frequency\tStatistic";

    // Impose range limits.
    std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
    size_type begin_index = indices.first;
    size_type end_index = indices.second;

    // Write out the statistics.
    for (size_type ii = begin_index; ii < end_index; ++ii) os << std::endl << m_freq[ii] << "\t" << m_spec[ii];

    // Restore original precision.
    os.precision(save_precision);

    return os;
  }

  std::pair<FourierAnalysis::size_type, FourierAnalysis::size_type> FourierAnalysis::getRangeIndex(double min_freq,
    double max_freq) const {
    size_type begin_index = 0;
    size_type end_index = m_freq.size();

    if (0. <= min_freq) {
      // Find first element whose frequency is not less than the minimum frequency.
      for (begin_index = 0; begin_index < m_freq.size() && m_freq[begin_index] < min_freq; ++begin_index);
    }

    if (0. <= max_freq) {
      // Find last element whose frequency is not greater than the maximum frequency.
      for (end_index = m_freq.size(); end_index > begin_index && m_freq[end_index - 1] > max_freq; --end_index);
    }

    return std::make_pair(begin_index, end_index);
  }

}
