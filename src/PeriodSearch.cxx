/** \file PeriodSearch.cxx
    \brief Implementation of PeriodSearch class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include "periodSearch/PeriodSearch.h"

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace periodSearch {

  const double PeriodSearch::s_2pi = 2. * 4. * atan(1.0);

  PeriodSearch::PeriodSearch(size_type num_bins): m_freq(num_bins), m_spec(num_bins) {}

  void PeriodSearch::plot(const std::string & title, const std::string & freq_unit, double min_freq, double max_freq) const {
    using namespace st_graph;

    try {
      // Display value of maximum frequency/statistic in title.
      std::ostringstream os;
      std::pair<double, double> max = findMax(min_freq, max_freq);
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
      std::cerr << "Warning: PeriodSearch::plot could not display plot." << std::endl;
    }
  }

  std::pair<double, double> PeriodSearch::findMax(double min_freq, double max_freq) const {
    bool found_max = false;
    size_type max_idx = 0;
    double max = 0.;

    // Impose range limits.
    std::pair<size_type, size_type> indices = getRangeIndex(min_freq, max_freq);
    size_type begin_index = indices.first;
    size_type end_index = indices.second;

    for (size_type ii = begin_index; ii < end_index; ++ii) {
      // If the value is larger than the current maximum, replace it.
      if (m_spec[ii] > max) {
        max = m_spec[ii];
        max_idx = ii;
        found_max = true;
      }
    }

    // Make sure a valid maximum was found.
    if (!found_max) {
      std::ostringstream os;
      os << "PeriodSearch::findMax cannot find any trial frequency in range [" << min_freq << ", " << max_freq << "]";
      throw std::runtime_error(os.str());
    }
    return std::pair<double, double>(m_freq[max_idx], max);
  }

  st_stream::OStream & PeriodSearch::write(st_stream::OStream & os, double min_freq, double max_freq) const {
    using namespace std;

    // Get info about the maximum.
    std::pair<double, double> max = findMax(min_freq, max_freq);

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

  std::pair<PeriodSearch::size_type, PeriodSearch::size_type> PeriodSearch::getRangeIndex(double min_freq,
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

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearch & test) {
    return test.write(os);
  }

}
