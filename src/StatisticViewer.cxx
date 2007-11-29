/** \file StatisticViewer.cxx
    \brief Implementation of StatisticViewer class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#include "StatisticViewer.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"

#include "tip/Header.h"
#include "tip/Table.h"

StatisticViewer::StatisticViewer(index_type num_axis, data_type::size_type num_element): m_begin_cont(num_axis),
  m_num_element(num_element), m_label_cont(num_axis), m_unit_cont(num_axis), m_title(), m_caption() {}

void StatisticViewer::setData(index_type axis_index, const data_type::const_iterator & begin) {
  m_begin_cont.at(axis_index) = begin;
}

void StatisticViewer::setLabel(index_type axis_index, const std::string & label) {
  m_label_cont.at(axis_index) = label;
}

void StatisticViewer::setUnit(index_type axis_index, const std::string & unit) {
  m_unit_cont.at(axis_index) = unit;
}

void StatisticViewer::setTitle(const std::string & title) {
  m_title = title;
}

void StatisticViewer::setCaption(const std::string & caption) {
  m_caption = caption;
}

void StatisticViewer::plot(index_type x_axis_index, index_type y_axis_index) const {
  using namespace st_graph;

  // Get data to plot here, in order to check the indicies before creating a plot.
  data_type::const_iterator x_begin = m_begin_cont.at(x_axis_index);
  data_type::const_iterator y_begin = m_begin_cont.at(y_axis_index);

  // Create axis labels here, in order to check the indicies before creating a plot.
  std::string x_label = m_label_cont.at(x_axis_index) + " " + m_unit_cont.at(x_axis_index);
  std::string y_label = m_label_cont.at(y_axis_index) + " " + m_unit_cont.at(y_axis_index);

  try {
    // Get graphics engine to set up graph.
    Engine & engine(Engine::instance());

    // Typedef for readability.
    typedef st_graph::ValueSequence<std::vector<double>::const_iterator> ValueSeq_t;

    // Create plot, using frequency as x, and spectrum/statistic as y.
    // TODO: Add m_caption in a text box on the plot, and/or in a GUI output window.
    std::auto_ptr<IPlot> plot(engine.createPlot(m_title, 800, 600, "hist",
      ValueSeq_t(x_begin, x_begin + m_num_element), ValueSeq_t(y_begin, y_begin + m_num_element)));

    // Set axes titles.
    std::vector<Axis> & axes(plot->getAxes());
    axes[0].setTitle(x_label);
    axes[1].setTitle(y_label);

    // Display plot.
    engine.run();

  } catch (const std::exception & x) {
    std::cerr << x.what() << std::endl;
    std::cerr << "Warning: StatisticViewer::plot could not display a plot." << std::endl;
  }
}

st_stream::StreamFormatter & StatisticViewer::write(st_stream::StreamFormatter & os) const {
  // Select chatness levels for caption and data.
  st_stream::OStream & os_caption = os.info(eIncludeCaption);
  st_stream::OStream & os_data = os.info(eIncludeData);

  // Write out the caption.
  os_caption << m_caption << std::endl;

  // Get the number of axes.
  index_type num_axis = m_begin_cont.size();

  // Write out axis labels and units.
  for (index_type ii = 0; ii < num_axis; ++ii) {
    if (ii != 0) os_data << "\t";
    os_data << m_label_cont[ii] << "(" << m_unit_cont[ii] << ")";
  }
  os_data << std::endl;

  // Save current precision, and set desired precision in this method.
  int save_precision = os_data.precision();
  os_data.precision(std::numeric_limits<double>::digits10);

  // Write out the statistics.
  for (data_type::size_type diff = 0; diff < m_num_element; ++diff) {
    for (index_type ii = 0; ii < num_axis; ++ii) {
      if (ii != 0) os_data << "\t";
      os_data << *(m_begin_cont[ii] + diff);
    }
    os_data << std::endl;
  }

  // Restore original precision.
  os_data.precision(save_precision);

  // Return the stream.
  return os;
}

tip::Table & StatisticViewer::write(tip::Table & table) const {
  // Write description of this search into the header.
  std::stringstream ss;
  ss << m_caption;
  while (ss.good()) {
    const unsigned int buf_size = 1024;
    char buf[buf_size];
    ss.getline(buf, buf_size);
    table.getHeader().addComment(buf);
  }

  // Resize the table to accomodate all the data.
  table.setNumRecords(m_num_element);

  // Start at the beginning of the table.
  tip::Table::Iterator itor = table.begin();

  // Get the number of axes.
  index_type num_axis = m_begin_cont.size();

  // Write out the statistics.
  for (data_type::size_type diff = 0; diff < m_num_element; ++diff, ++itor) {
    for (index_type ii = 0; ii < num_axis; ++ii) {
      std::string label = m_label_cont[ii];
      double value = *(m_begin_cont[ii] + diff);
      (*itor)[label].set(value);
    }
  }

  return table;
}
