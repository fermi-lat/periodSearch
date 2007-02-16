/** \file PeriodSearchPlotter.cxx
    \brief Implementation of PeriodSearchPlotter class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include "periodSearch/PeriodSearch.h"
#include "periodSearch/PeriodSearchPlotter.h"

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace periodSearch {

  PeriodSearchPlotter::PeriodSearchPlotter() {}

  PeriodSearchPlotter::~PeriodSearchPlotter() {}

  void PeriodSearchPlotter::plot(const PeriodSearch & search, const std::string & title, const std::string & freq_unit) const {
    plotRange(search, title, freq_unit);
  }

  void PeriodSearchPlotter::plotRange(const PeriodSearch & search, const std::string & title, const std::string & freq_unit,
    double min_freq, double max_freq) const {
    using namespace st_graph;
    typedef PeriodSearch::cont_type cont_type;
    typedef PeriodSearch::size_type size_type;

    try {
      // Display value of maximum frequency/statistic in title.
      std::ostringstream os;
      std::pair<double, double> max = search.findMax(min_freq, max_freq);
      std::pair<double, double> chance_prob = search.chanceProb(max.second);

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
      std::pair<size_type, size_type> indices = search.getRangeIndex(min_freq, max_freq);
      size_type begin_index = indices.first;
      size_type end_index = indices.second;

      // Create plot, using frequency as x, and spectrum/statistic as y.
      const cont_type & freq(search.getFreq());
      const cont_type & spec(search.getSpec());
      std::auto_ptr<IPlot> plot(engine.createPlot(os.str(), 800, 600, "hist",
        ValueSeq_t(freq.begin() + begin_index, freq.begin() + end_index),
        ValueSeq_t(spec.begin() + begin_index, spec.begin() + end_index)));

      // Set axes titles.
      std::vector<Axis> & axes(plot->getAxes());
      axes[0].setTitle("Frequency " + freq_unit);
      axes[1].setTitle("Test Statistic");

      // Display plot.
      engine.run();

    } catch (const std::exception & x) {
      std::cerr << x.what() << std::endl;
      std::cerr << "Warning: PeriodSearchPlotter::plot could not display plot." << std::endl;
    }
  }

}
