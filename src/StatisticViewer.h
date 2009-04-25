/** \file StatisticViewer.h
    \brief Declaration of StatisticViewer class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_StatisticViewer_h
#define periodSearch_StatisticViewer_h

#include <string>
#include <vector>

namespace st_stream {
  class StreamFormatter;
}

namespace tip {
  class Table;
}

/** \class StatisticViewer
    \brief Handler of a graphical plot, a text output, and a FITS file output of numeric data, designed to display
           statistical data obtained from pulsation searches and periodicity tests.
*/
class StatisticViewer {
  public:
    typedef unsigned long index_type;
    typedef std::vector<double> data_type;

    /** \brief Create a viewer object for plotting and writing designated data for human viewing.
        \param num_axis The number of axes of data to be viewed. For one-dimensional histogram,
               for example, give 2 for num_axis (1 for X-axis and 1 for Y-axis).
        \param num_element The number of data elements to be viewed (common to all axes).
    */
    StatisticViewer(index_type num_axis, data_type::size_type num_element);

    /** \brief Get a reference to an internal data storage for a given axis.
        \param axis_index The index of axis, for which a reference to the data storage is to be returned.
    */
    const data_type & getData(index_type axis_index) const;
    data_type & getData(index_type axis_index);

    /** \brief Mark a part of the internal data storage as "selected" for viewing.
        \param begin_index The index for the internal storage which points to the first element of the selected range.
        \param end_index The index for the internal storage which points to one past the last element of the selected range.
    */
    void selectData(data_type::size_type begin_index, data_type::size_type end_index);

    /** \brief Set an axis label.
        \param axis_index The index of axis, for which an axis label is to be set.
        \param label The axis label to set.
    */
    void setLabel(index_type axis_index, const std::string & label);

    /** \brief Set a unit of an axis.
        \param axis_index The index of axis, for which a unit is to be set.
        \param label The unit for the axis to set.
    */
    void setUnit(index_type axis_index, const std::string & unit);

    /** \brief Set a title of a plot, a text output, and a FITS output.
        \param title The title to set.
    */
    void setTitle(const std::string & title);

    /** \brief Set a caption of a plot, a text output, and a FITS output.
        \param caption The caption to set.
    */
    void setCaption(const std::string & caption);

    /** \brief Display a graphical plot of designated data.
        \param x_index The axis index of the data to be used for X-axis of the plot.
        \param y_index The axis index of the data to be used for Y-axis of the plot.
    */
    void plot(index_type x_axis_index = 0, index_type y_axis_index = 1) const;

    /** \brief Write a text output to an output stream. The output will be controlled by chatness levels given.
        \param os StreamFormatter object, to which the text output is to be forwarded.
        \param chat_caption Minimum chatness level to write the caption.
        \param chat_data Minimum chatness level to write the data contents.
    */
    st_stream::StreamFormatter & write(st_stream::StreamFormatter & os, int chat_caption = 2, int chat_data = 3) const;

    /** \brief Write a text data (title, caption, etc.) into a FITS header and numerical data into a FITS table.
        \param table FITS table, which the data are wirtten into.
    */
    tip::Table & write(tip::Table & table) const;

  private:
    data_type::size_type m_num_axis;
    data_type::size_type m_num_element;
    std::vector<data_type> m_data_cont;
    data_type::size_type m_begin_index;
    data_type::size_type m_end_index;
    std::vector<std::string> m_label_cont;
    std::vector<std::string> m_unit_cont;
    std::string m_title;
    std::string m_caption;
};

#endif
