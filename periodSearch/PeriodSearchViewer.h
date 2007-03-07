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

      /** \brief Write description of this test to the given stream.
          \param os The stream to which to send output.
      */
      virtual st_stream::OStream & write(st_stream::OStream & os) const;

      // TODO Add fits file output method.
      //  o Summary goes to COMMENT header keywords (human readable).
      //  o Data goes to FITS columns (FREQUENCY and POWER, to match Xronos? How do we handle errors?).

    private:
      const PeriodSearch * m_search;
      double m_min_freq;
      double m_max_freq;
  };

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodSearchViewer & viewer);
}

#endif
