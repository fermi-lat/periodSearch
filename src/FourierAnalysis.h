/** \file FourierAnalysis.h
    \brief Declaration of FourierAnalysis class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef periodSearch_FourierAnalysis_h
#define periodSearch_FourierAnalysis_h

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "periodSearch/PeriodSearch.h"

// TODO: Put this in a namespace.

/** \class FourierAnalysis
*/
class FourierAnalysis : public periodSearch::PeriodSearch {
  public:
    /// \brief Container type used to store histograms.
//    typedef periodSearch::PeriodTest::HistCont_t HistCont_t;

    /** \brief Construct a FourierAnalysis.
        \param t_start Time lower boundary.
        \param t_stop Time upper boundary.
        \param width Width of one time bin in one data subset to be transformed.
        \param num_bins The number of bins used in each FFT. Depending on the specific test, this
                        may be an upper limit on the actual number of bins used.
        \param num_events Hint giving the anticipated number of events to be filled.
    */
    FourierAnalysis(double t_start, double t_stop, double width, size_type num_bins, int num_events = 0);

    /** \brief Fill given time into histograms.
        \param evt_time The time of the event.
    */
    virtual void fill(double evt_time);

    /** \brief
    */
    virtual const std::vector<double> & computeStats();

    /** \brief Display plot of statistics.
        \param title The title to display on the plot. (Purely cosmetic.)
        \param freq_unit The units to display on the x axis. (Purely cosmetic.)
        \param min_freq The minimum frequency in the range.
        \param max_freq The maximum frequency in the range.
    */
    virtual void plot(const std::string & title, const std::string & freq_unit, double min_freq = -1., double max_freq = -1.) const;

    /** \brief Display plot of statistics.
        \param title The title to display on the plot. (Purely cosmetic.)
        \param freq_unit The units to display on the x axis. (Purely cosmetic.)
    */
    virtual void plotStats(const std::string & title, const std::string & freq_unit) const;

    /** \brief Find the frequency for which the statistic is maximized. Return the frequency and the value of
               the statistic, as a pair.
    */
    virtual std::pair<double, double> findMax() const;

    /** \brief Find the frequency for which the statistic is maximized in a given frequency range. Return the
               frequency and the value of the statistic, as a pair.
        \param min_freq The minimum frequency in the range.
        \param max_freq The maximum frequency in the range.
    */
    virtual std::pair<double, double> findMaxRange(double min_freq = -1., double max_freq = -1.) const;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProb(double stat) const; 

    /** \brief Write data as a function of frequency to the given stream.
        \param os The stream.
    */
    virtual st_stream::OStream & write(st_stream::OStream & os) const;

    /** \brief Write data over a specified frequency range as a function of frequency to the given stream.
        \param os The stream.
        \param min_freq The minimum frequency in the range.
        \param max_freq The maximum frequency in the range.
    */
    virtual st_stream::OStream & writeRange(st_stream::OStream & os, double min_freq = -1., double max_freq = -1.) const;

  private:
    /** \brief Given a frequency range, determine the indices of (inclusive) lower and upper bounds.
        \param min_freq The minimum frequency.
        \param max_freq The maximum frequency.
    */
    std::pair<size_type, size_type> getRangeIndex(double min_freq, double max_freq) const;

    typedef std::multimap<size_type, size_type> index_map_type;
    index_map_type m_index;
    double m_t_start;
    double m_t_stop;
    double m_width;
    size_type m_num_segments;
    size_type m_num_bins;
};

#endif
