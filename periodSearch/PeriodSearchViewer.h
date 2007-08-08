/** \file PeriodSearchViewer.h
    \brief Declaration of PeriodSearchViewer class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_PeriodSearchViewer_h
#define periodSearch_PeriodSearchViewer_h

#include <string>

namespace st_stream {
  class OStream;
}

namespace tip {
  class Header;
  class Table;
}

namespace periodSearch {

  class PeriodSearch;

  /** \class PeriodSearchViewer
      \brief Base class for various statistical tests used to determine frequency of pulsation
             when an approximate frequency is known.
  */
  class PeriodSearchViewer {
    public:
      /** \brief Create a viewer object for viewing period search objects (plotting and/or writing output).
          \param search The search being viewed.
          \param min_freq The minimum frequency in the range (if negative, do not constrain minimum frequency.)
          \param max_freq The maximum frequency in the range (if negative, do not constrain maximum frequency.)
      */
      PeriodSearchViewer(const PeriodSearch & search, double min_freq = -1., double max_freq = -1.);

      virtual ~PeriodSearchViewer();

      /** \brief Display plot of statistic as a function of frequency for entire frequency range.
          \param title The title to display on the plot. (Purely cosmetic.)
          \param freq_unit The units to display on the x axis. (Purely cosmetic.)
      */
      virtual void plot(const std::string & title, const std::string & freq_unit) const;

      /** \brief Write description of this search to the given stream.
          \param os The stream to which to send output.
      */
      virtual st_stream::OStream & writeSummary(st_stream::OStream & os) const;

      /** \brief Write description of this search to the given tip header.
          \param header The tip header to which to write the summary.
      */
      virtual tip::Header & writeSummary(tip::Header & header) const;

      /** \brief Write data used to perform this search to the given stream.
          \param os The stream to which to send output.
      */
      virtual st_stream::OStream & writeData(st_stream::OStream & os) const;

      /** \brief Write data used to perform this search to the given tip table.
          \param table The tip table to which to send output.
      */
      virtual tip::Table & writeData(tip::Table & table) const;

    private:
      const PeriodSearch * m_search;
      double m_min_freq;
      double m_max_freq;
  };

}

#endif
