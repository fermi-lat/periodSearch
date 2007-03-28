/** \file PeriodSearchViewer.cxx
    \brief Implementation of PeriodSearchViewer class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include "periodSearch/PeriodSearch.h"
#include "periodSearch/PeriodSearchViewer.h"

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "st_stream/Stream.h"

#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace periodSearch {

  PeriodSearchViewer::PeriodSearchViewer(const PeriodSearch & search, double min_freq, double max_freq):
    m_search(&search), m_min_freq(min_freq), m_max_freq(max_freq) {}

  PeriodSearchViewer::~PeriodSearchViewer() {}

  void PeriodSearchViewer::plot(const std::string & title, const std::string & freq_unit) const {
    using namespace st_graph;
    typedef PeriodSearch::cont_type cont_type;
    typedef PeriodSearch::size_type size_type;

    try {
      // Get graphics engine to set up graph.
      Engine & engine(Engine::instance());

      // Typedef for readability.
      typedef st_graph::ValueSequence<std::vector<double>::const_iterator> ValueSeq_t;

      // Impose range limits.
      // TODO Move interpretation of min/max frequency into constructor and store the selected index locating
      // the min and max to avoid recomputing min/max index each time this is called.
      std::pair<size_type, size_type> indices = m_search->getRangeIndex(m_min_freq, m_max_freq);
      size_type begin_index = indices.first;
      size_type end_index = indices.second;

      // Create plot, using frequency as x, and spectrum/statistic as y.
      const cont_type & freq(m_search->getFreq());
      const cont_type & spec(m_search->getSpec());
      // TODO Add output from PeriodSearchResult::write(...) to a text box on the plot, and/or in a
      // GUI output window.
      std::auto_ptr<IPlot> plot(engine.createPlot(title, 800, 600, "hist",
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
      std::cerr << "Warning: PeriodSearchViewer::plot could not display plot." << std::endl;
    }
  }

  st_stream::OStream & PeriodSearchViewer::writeSummary(st_stream::OStream & os) const {
    os << m_search->search(m_min_freq, m_max_freq);
    return os;
  }

  st_stream::OStream & PeriodSearchViewer::writeData(st_stream::OStream & os) const {
    using namespace std;

    // Save current precision.
    int save_precision = os.precision();

    os.precision(std::numeric_limits<double>::digits10);

    // Write out the data.
    os << "Frequency\tStatistic";

    // Impose range limits.
    // TODO Move interpretation of min/max frequency into constructor and store the selected index locating
    // the min and max to avoid recomputing min/max index each time this is called.
    std::pair<PeriodSearch::size_type, PeriodSearch::size_type> indices = m_search->getRangeIndex(m_min_freq, m_max_freq);
    PeriodSearch::size_type begin_index = indices.first;
    PeriodSearch::size_type end_index = indices.second;

    const std::vector<double> & freq(m_search->getFreq());
    const std::vector<double> & spec(m_search->getSpec());

    // Write out the statistics.
    for (PeriodSearch::size_type ii = begin_index; ii < end_index; ++ii) os << std::endl << freq[ii] << "\t" << spec[ii];

    // Restore original precision.
    os.precision(save_precision);

    return os;
  }

}
