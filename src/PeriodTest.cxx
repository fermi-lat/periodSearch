/** \file PeriodTest.cxx
    \brief Implementation of PeriodTest class.
    \authors Masaharu Hirayama, GSSC
             James Peachey, HEASARC/GSSC
*/

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "periodSearch/PeriodTest.h"

namespace periodSearch {

  const double PeriodTest::s_2pi = 2. * 4. * atan(1.0);

  PeriodTest::PeriodTest(double center, double step, long num_trials, double epoch, int num_bins, double duration):
    m_trial_hist(num_trials, std::vector<std::complex<double> >(num_bins, 0.)), m_freqs(num_trials), m_stats(num_trials),
    m_center(center), m_epoch(epoch), m_duration(duration), m_num_bins(num_bins), m_num_events(0), m_num_indep_trials(0) {
    // Make certain there is no error in the input.
    if (0. >= m_duration) throw std::logic_error("PeriodTest constructor was passed a non-positive duration");
    if (0 >= m_num_bins) throw std::logic_error("PeriodTest constructor was passed a non-positive number of bins");

    // Create vector containing the trial frequencies.
    long ii_cent = num_trials / 2;
    for (long ii = 0; ii < num_trials; ++ii) {
      m_freqs[ii] = center + (ii - ii_cent) * step;
      if (0. >= m_freqs[ii]) throw std::logic_error("PeriodTest constructor computed a non-positive trial frequency");
    }

    // Compute number of independent trials.
    double fourier_res = 1. / m_duration;

    //    N_Fourier = (stop - start) / Fourier_step
    unsigned int n_fourier = int(ceil(fabs((m_freqs.back() - m_freqs.front()) / fourier_res)));

    m_num_indep_trials = (n_fourier < m_trial_hist.size()) ? n_fourier : m_trial_hist.size();
  }

  void PeriodTest::fill(double evt_time) {
    // Pseudocode taken from work by M. Hirayama, which may be seen at:
    // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/pc_chi2test.txt
    //     dt = evtime - epoch;
    double dt = evt_time - m_epoch;

    //     N_event++;
    ++m_num_events;

    //     for (i=0; i<N_trial; i++) {
    //       phase = dt / period_array[i];
    //       phase -= floor(phase);
    //
    //       [fill one entry into hist[i] at "phase"];
    //     }

    // Iterate over the number of trial frequencies.
    int num_trials = m_trial_hist.size();
    for (int ii = 0; ii < num_trials; ++ii) {
      // For each frequency, compute the phase.
      double phase = dt * m_freqs[ii];
      phase -= floor(phase);

      // Use this phase information to fill in the corresponding trial.
      fillOneTrial(phase, m_trial_hist[ii]);
    }
  }

  void PeriodTest::plotStats(const std::string & title, const std::string & freq_unit) const {
    using namespace st_graph;

    try {
      // Display value of maximum frequency/statistic in title.
      std::ostringstream os;
      std::pair<double, double> max = findMax();
      std::pair<double, double> chance_prob = chanceProb(max.second);

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

      // Create plot, using m_freqs as x, and m_stats as y.
      std::auto_ptr<IPlot> plot(engine.createPlot(os.str(), 800, 600, "hist", ValueSeq_t(m_freqs.begin(), m_freqs.end()),
        ValueSeq_t(m_stats.begin(), m_stats.end())));

      // Set axes titles.
      std::vector<Axis> & axes(plot->getAxes());
      axes[0].setTitle("Frequency " + freq_unit);
      axes[1].setTitle("Test Statistic");

      // Display plot.
      engine.run();

    } catch (const std::exception & x) {
      std::cerr << x.what() << std::endl;
      std::cerr << "Warning: PeriodTest::plotStats could not display plot." << std::endl;
      return;
    }
  }

  std::pair<double, double> PeriodTest::findMax() const {
    unsigned long max_idx = 0;
    double max = 0.;
    for (unsigned long ii = 0; ii < m_stats.size(); ++ii) {
      if (m_stats[ii] > max) {
        max = m_stats[ii];
        max_idx = ii;
      }
    }
    return std::pair<double, double>(m_freqs[max_idx], max);
  }

  st_stream::OStream & PeriodTest::write(st_stream::OStream & os) const {
    using namespace std;

    // Get info about the maximum.
    std::pair<double, double> max = findMax();

    // Chance probability.
    std::pair<double, double> chance_prob = chanceProb(max.second);

    // Save current precision.
    int save_precision = os.precision();

    os.precision(15);

    // Write out the results.
    os << "Maximum at: " << max.first << std::endl << "Statistic: " << max.second << std::endl;
    os << "Chance probability range: (" << chance_prob.first << ", " << chance_prob.second << ")" << std::endl;
    os << "Frequency\tStatistic";

    // Write out the statistics.
    for (unsigned long ii = 0; ii < m_stats.size(); ++ii) os << std::endl << m_freqs[ii] << "\t" << m_stats[ii];

    // Restore original precision.
    os.precision(save_precision);

    return os;
  }

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodTest & test) {
    return test.write(os);
  }

}
