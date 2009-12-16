/** \file StatisticViewer.cxx
    \brief Implementation of StatisticViewer class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#include "StatisticViewer.h"

#include <iostream>
#include <limits>
#include <memory>
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
#include "tip/TipException.h"

StatisticViewer::StatisticViewer(index_type num_axis, data_type::size_type num_element): m_num_axis(num_axis),
  m_num_element(num_element), m_data_cont(num_axis, data_type(num_element, 0.)), m_begin_index(0), m_end_index(num_element),
  m_label_cont(num_axis), m_unit_cont(num_axis), m_title(), m_caption() {}

const StatisticViewer::data_type & StatisticViewer::getData(index_type axis_index) const {
  return m_data_cont.at(axis_index);
}

StatisticViewer::data_type & StatisticViewer::getData(index_type axis_index) {
  return m_data_cont.at(axis_index);
}

void StatisticViewer::selectData(data_type::size_type begin_index, data_type::size_type end_index) {
  // Check whether the desired range is found or not.
  if (begin_index >= end_index) {
    std::ostringstream os;
    os << "No bins contained in the given index range [" << begin_index << ", " << end_index << "]";
    throw std::runtime_error(os.str());

  } else {
    // Set the indices to the data member.
    m_begin_index = begin_index;
    m_end_index = end_index;
  }
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
  const data_type & x_data = m_data_cont.at(x_axis_index);
  const data_type & y_data = m_data_cont.at(y_axis_index);

  // Get units here, in order to check the indicies before creating a plot.
  const std::string & x_unit = m_unit_cont.at(x_axis_index);
  const std::string & y_unit = m_unit_cont.at(y_axis_index);

  // Create axis labels here, in order to check the indicies before creating a plot.
  std::string x_label = m_label_cont.at(x_axis_index);
  if (!x_unit.empty()) x_label += " (" + x_unit + ")";
  std::string y_label = m_label_cont.at(y_axis_index);
  if (!y_unit.empty()) y_label += " (" + y_unit + ")";

  try {
    // Get graphics engine to set up graph.
    Engine & engine(Engine::instance());

    // Typedef for readability.
    typedef st_graph::ValueSequence<std::vector<double>::const_iterator> ValueSeq_t;

    // Create plot, using frequency as x, and spectrum/statistic as y.
    // TODO: Add m_caption in a text box on the plot, and/or in a GUI output window.
    std::auto_ptr<IPlot> plot(engine.createPlot(m_title, 800, 600, "hist",
      ValueSeq_t(x_data.begin() + m_begin_index, x_data.begin() + m_end_index),
      ValueSeq_t(y_data.begin() + m_begin_index, y_data.begin() + m_end_index)));

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

st_stream::StreamFormatter & StatisticViewer::write(st_stream::StreamFormatter & os, int chat_caption, int chat_data) const {
  // Write out the caption.
  os.info(chat_caption) << m_caption << std::endl;

  // Write out axis labels and units.
  for (index_type axis_index = 0; axis_index < m_num_axis; ++axis_index) {
    if (axis_index != 0) os.info(chat_data) << "\t";
    os.info(chat_data) << m_label_cont[axis_index];
    const std::string & unit_string = m_unit_cont[axis_index];
    if (!unit_string.empty()) os.info(chat_data) << "(" << unit_string << ")";
  }
  os.info(chat_data) << std::endl;

  // Save current precision, and set desired precision in this method.
  int save_precision = os.info(chat_data).precision();
  os.info(chat_data).precision(std::numeric_limits<double>::digits10);

  // Write out the statistics.
  for (data_type::size_type elem_index = m_begin_index; elem_index < m_end_index; ++elem_index) {
    for (index_type axis_index = 0; axis_index < m_num_axis; ++axis_index) {
      if (axis_index != 0) os.info(chat_data) << "\t";
      os.info(chat_data) << m_data_cont[axis_index][elem_index];
    }
    os.info(chat_data) << std::endl;
  }

  // Restore original precision.
  os.info(chat_data).precision(save_precision);

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

  // Append FITS columns if missing.
  for (index_type axis_index = 0; axis_index < m_num_axis; ++axis_index) {
    const std::string & field_name = m_label_cont[axis_index];
    try {
      table.getFieldIndex(field_name);
    } catch (const tip::TipException &) {
      table.appendField(field_name, "1D");
    }

    // Set the field unit to each column, if unit is not blank.
    const std::string & field_unit = m_unit_cont[axis_index];
    if (!field_unit.empty()) {
      tip::IColumn * column = table.getColumn(table.getFieldIndex(field_name));
      column->getColumnKeyword("TUNIT").set(field_unit);
    }
  }

  // Resize the table to accomodate all the data.
  table.setNumRecords(m_end_index - m_begin_index);

  // Start at the beginning of the table.
  tip::Table::Iterator itor = table.begin();

  // Write out the statistics.
  for (data_type::size_type elem_index = m_begin_index; elem_index < m_end_index; ++elem_index, ++itor) {
    for (index_type axis_index = 0; axis_index < m_num_axis; ++axis_index) {
      std::string label = m_label_cont[axis_index];
      double value = m_data_cont[axis_index][elem_index];
      (*itor)[label].set(value);
    }
  }

  return table;
}
