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

  PeriodTest::PeriodTest(double center, double step, size_type num_trials, double epoch, size_type array_size, double duration):
    PeriodSearch(num_trials), m_trial_hist(num_trials, std::vector<std::complex<double> >(array_size, 0.)),
    m_center(center), m_epoch(epoch), m_duration(duration), m_num_events(0), m_num_indep_trials(0) {
    // Make certain there is no error in the input.
    if (0. >= center) throw std::logic_error("PeriodTest constructor was passed a non-positive center");
    if (0. >= step) throw std::logic_error("PeriodTest constructor was passed a non-positive step");
    if (0 == num_trials) throw std::logic_error("PeriodTest constructor was passed a non-positive num_trials");
    if (0. >= m_duration) throw std::logic_error("PeriodTest constructor was passed a non-positive duration");
    if (0 == array_size) throw std::logic_error("PeriodTest constructor was passed a non-positive array_size");

    // Create vector containing the trial frequencies.
    size_type ii_cent = num_trials / 2;
    double min = center - ii_cent * step;

    // Check whether the step was too big, leading to a negative frequency.
    if (0. >= min) throw std::logic_error("PeriodTest constructor computed a non-positive trial frequency");

    // Step from minimum frequency on up, populating frequency array.
    for (size_type ii = 0; ii < num_trials; ++ii) {
      m_freq[ii] = min + ii * step;
    }

    // Compute number of independent trials.
    double fourier_res = 1. / m_duration;

    //    N_Fourier = (stop - start) / Fourier_step
    size_type n_fourier = size_type(ceil(fabs((m_freq.back() - m_freq.front()) / fourier_res)));

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
    size_type num_trials = m_trial_hist.size();
    for (size_type ii = 0; ii < num_trials; ++ii) {
      // For each frequency, compute the phase.
      double phase = dt * m_freq[ii];
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

      // Create plot, using m_freq as x, and m_spec as y.
      std::auto_ptr<IPlot> plot(engine.createPlot(os.str(), 800, 600, "hist", ValueSeq_t(m_freq.begin(), m_freq.end()),
        ValueSeq_t(m_spec.begin(), m_spec.end())));

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
    size_type max_idx = 0;
    double max = 0.;
    for (size_type ii = 0; ii < m_spec.size(); ++ii) {
      if (m_spec[ii] > max) {
        max = m_spec[ii];
        max_idx = ii;
      }
    }
    return std::pair<double, double>(m_freq[max_idx], max);
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
    for (size_type ii = 0; ii < m_spec.size(); ++ii) os << std::endl << m_freq[ii] << "\t" << m_spec[ii];

    // Restore original precision.
    os.precision(save_precision);

    return os;
  }

  st_stream::OStream & operator <<(st_stream::OStream & os, const PeriodTest & test) {
    return test.write(os);
  }

}
