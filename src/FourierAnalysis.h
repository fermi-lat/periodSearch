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

#include "PeriodSearch.h"

/** \class FourierAnalysis
    \brief Class for pulsation search by computing a power spectrum of a given time series data using Fast Fourier
           Transform (FFT) technique.

           Note on the units of time and frequency: This class does not assume a specific time unit. Instead, all
           times (both time series data and time durations) given to an object of this class must be in the same time
           units throughout the lifetime of the object, and all frequencies given to it must be in the units of the
           inverse of the time unit. All results from the object will be in the same units of time or frequency.
*/
class FourierAnalysis : public PeriodSearch {
  public:
    /** \brief Construct a FourierAnalysis object, which computes a power spectrum of a given time series data using FFT.
        \param t_start Time lower boundary.
        \param t_stop Time upper boundary.
        \param width Width of one time bin in one data subset to be transformed.
        \param num_bins The number of bins used in each FFT. Depending on the specific test, this
                        may be an upper limit on the actual number of bins used.
        \param num_events Hint giving the anticipated number of events to be filled.
        \param freq_unit Character string representing the unit of frequency, or the inverse of the unit of times given
               to this period search.
    */
    FourierAnalysis(double t_start, double t_stop, double width, size_type num_bins, const std::string & freq_unit, int num_events = 0);

    /// \brief Destruct this FourierAnalysis object.
    virtual ~FourierAnalysis() {}

    /** \brief Fill given time into the internal storage of this object.
        \param evt_time The time of the event.
    */
    virtual void fill(double evt_time);

    /** \brief Compute Fourier powers at trial frequencies with Discrete Fast Fourier Transform algorighm.
    */
    virtual void computeStat();

    /** \brief Compute the number of independent trials for this search method.
        \param min_freq The minimum frequency.
        \param max_freq The maximum frequency.
    */
    virtual size_type computeNumIndepTrials(double min_freq = -1., double max_freq = -1.) const;

    /** \brief Compute the chance probability for the given parameters. Return pair with lower, upper limit.
        \param stat The value of the statistic.
    */
    virtual std::pair<double, double> computeChanceProbOneTrial(double stat) const; 

  private:
    typedef std::multimap<size_type, size_type> index_map_type;
    index_map_type m_index;
    double m_t_start;
    double m_t_stop;
    double m_width;
    double m_fourier_res;
    size_type m_num_segments;
    size_type m_num_bins;
};

#endif
