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

/** \class FourierAnalysis
*/
class FourierAnalysis : public periodSearch::PeriodSearch {
  public:
    /** \brief Construct a FourierAnalysis.
        \param t_start Time lower boundary.
        \param t_stop Time upper boundary.
        \param width Width of one time bin in one data subset to be transformed.
        \param num_bins The number of bins used in each FFT. Depending on the specific test, this
                        may be an upper limit on the actual number of bins used.
        \param num_events Hint giving the anticipated number of events to be filled.
    */
    // TODO: put minimum frequency in constructor low_freq_cutoff.
    FourierAnalysis(double t_start, double t_stop, double width, size_type num_bins, int num_events = 0);

    /** \brief Fill given time into histograms.
        \param evt_time The time of the event.
    */
    virtual void fill(double evt_time);

    /** \brief
    */
    virtual const std::vector<double> & computeStats();

    /** \brief Return the number of independent trials for this search method.
    */
    virtual size_type numIndepTrials() const;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> chanceProbOneTrial(double stat) const; 

    /** \brief Write data over a specified frequency range as a function of frequency to the given stream.
        \param os The stream.
        \param min_freq The minimum frequency in the range.
        \param max_freq The maximum frequency in the range.
    */
    virtual st_stream::OStream & writeRange(st_stream::OStream & os, double min_freq = -1., double max_freq = -1.) const;

  private:
    typedef std::multimap<size_type, size_type> index_map_type;
    index_map_type m_index;
    double m_t_start;
    double m_t_stop;
    double m_width;
    size_type m_num_segments;
    size_type m_num_bins;
};

#endif
