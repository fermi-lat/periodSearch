/** \file test_psearch.cxx
    \brief Period search tool.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/commonUtilities.h"

#include "pulsarDb/PdotCanceler.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/PulsarTestApp.h"
#include "timeSystem/TimeInterval.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "ChiSquaredTestArray.h"
#include "FoldingAnalysis.h"
#include "FourierAnalysis.h"
#include "HTestArray.h"
#include "PeriodicityTestArray.h"
#include "PeriodSearch.h"
#include "RayleighTestArray.h"
#include "Z2nTestArray.h"

static const std::string s_cvs_id = "$Name:  $";

class TestPeriodSearchApp : public timeSystem::PulsarTestApp {
  public:
    TestPeriodSearchApp();
    virtual ~TestPeriodSearchApp() throw() {}
    virtual void run();

    void testPeriodSearch();

    void testChanceProb();

    void testChiSquaredTestArray();

    void testZ2nTestArray();

    void testHTestArray();

    void testRayleighTestArray();

  private:
    void testAllStats(const std::string & prefix, const std::vector<double> & events, double t_start, double t_stop,
      double center, double step, long num_trials, double epoch, int num_bins,
      double fourier_width, int fourier_num_bins, double fourier_min_freq, double fourier_max_freq, bool plot);

    void testOneSearch(const std::vector<double> & events, PeriodSearch & search, const std::string & plot_title,
      const std::string & out_file, bool plot, double min_freq = -1., double max_freq = -1.);
};

TestPeriodSearchApp::TestPeriodSearchApp(): PulsarTestApp("periodSearch") {
  setName("test_stpsearch");
  setVersion(s_cvs_id);
}

void TestPeriodSearchApp::run() {
  // Test PeriodSearch subclasses.
  testPeriodSearch();

  // Test computations of chance probability.
  testChanceProb();

  // Test PeriodicityTestArary subclasses.
  testChiSquaredTestArray();
  testZ2nTestArray();
  testHTestArray();
  testRayleighTestArray();

  // Check test status.
  reportStatus();
}

void TestPeriodSearchApp::testPeriodSearch() {
  setMethod("testPeriodSearch");

  // Trick up some fake events.
  int num_evts = 1000;

  std::vector<double> fake_evts(num_evts);
  setPrecision(16);

  for (int ii = 0; ii < num_evts; ++ii)
    fake_evts[ii] = ii + 1;

  // Set up search parameters.
  double central = 1.;
  double step = .5e-3;
  int num_pds = 61;
  double epoch = .2;
  int num_bins = 10;
  double duration = 1000.;

  bool plot = getParGroup()["plot"];

  // First do simple test with this highly artificial data.
  // Note for Fourier test: width of .1 s -> Nyquist = 1/.2s = 5 Hz.
  testAllStats("artificial", fake_evts, 0., duration, central, step, num_pds, epoch, num_bins, .1, 10000, .9, 1.1, plot);

  // Data taken from M. Hirayama's work with modified ASCA data.
  // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/existing.html#tryout003
  num_pds = 101;

  central = 1. / 50.41843041e-3;
  step = .168e-7 * central * central;

  epoch = 212380785.922;

  fake_evts.clear();

  // Read some real data.
  std::string event_file = facilities::commonUtilities::joinPath(getDataPath(), "step-01.fits");
  std::auto_ptr<const tip::Table> evt_table(tip::IFileSvc::instance().readTable(event_file, "EVENTS"));

  // Get gti_table.
  std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(event_file, "GTI"));

  // Make the array big enough to hold these events.
  fake_evts.resize(evt_table->getNumRecords());

  // Correct event times for changed MJDREF.
  const double event_time_offset = (54101 - 51910) * 86400.;
  std::vector<double>::iterator event_itor = fake_evts.begin();
  for (tip::Table::ConstIterator itor = evt_table->begin(); itor != evt_table->end(); ++itor, ++event_itor) {
    *event_itor = (*itor)["TIME"].get() + event_time_offset;
  }

  double tstart = 0.;
  double tstop = 0.;

  const tip::Header & header(evt_table->getHeader());

  // If possible, get tstart and tstop from first and last interval in GTI extension.
  tip::Table::ConstIterator gti_itor = gti_table->begin();
  if (gti_itor != gti_table->end()) {
    // TSTART is the start of the first interval.
    tstart = (*gti_itor)["START"].get();
    // TSTOP is from the stop of the last interval.
    gti_itor = gti_table->end();
    --gti_itor;
    tstop = (*gti_itor)["STOP"].get();
  } else {
    header["TSTART"].get(tstart);
    header["TSTOP"].get(tstop);
  }

  // Correct tstart and tstop for changed MJDREF.
  tstart += event_time_offset;
  tstop += event_time_offset;

  // Repeat simple test with this somewhat less artificial data.
  // Note for Fourier test: width of .01 s -> Nyquist = 1/.02s = 50 Hz.
  testAllStats("psrb0540", fake_evts, tstart, tstop, central, step, num_pds, epoch, num_bins, .01, 1000000, 19.82, 19.85, plot);

  // Now test pdot correction.
  timeSystem::AbsoluteTime glast_origin("TDB", 51910, 0.);
  timeSystem::AbsoluteTime abs_epoch = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, epoch));
  double pdot = 4.7967744e-13;
  std::vector<double> fdot_ratio(1, -pdot * central);
  pulsarDb::PdotCanceler canceler("TDB", abs_epoch, fdot_ratio);

  // Correct the data.
  timeSystem::AbsoluteTime evt_time("TDB", 0, 0.);
  for (std::vector<double>::iterator itor = fake_evts.begin(); itor != fake_evts.end(); ++itor) {
    evt_time = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, *itor));
    canceler.cancelPdot(evt_time);
    *itor = (evt_time - glast_origin).computeDuration("TDB", "Sec");
  }

  // Cancel pdot in tstart, tstop and epoch to be consistent.
  evt_time = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, tstart));
  canceler.cancelPdot(evt_time);
  tstart = (evt_time - glast_origin).computeDuration("TDB", "Sec");

  evt_time = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, tstop));
  canceler.cancelPdot(evt_time);
  tstop = (evt_time - glast_origin).computeDuration("TDB", "Sec");

  evt_time = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, epoch));
  canceler.cancelPdot(evt_time);
  epoch = (evt_time - glast_origin).computeDuration("TDB", "Sec");

  // Repeat test with the pdot corrected data.
  testAllStats("psrb0540-pdot", fake_evts, tstart, tstop, central, step, num_pds, epoch, num_bins, .01, 1000000,
    19.82, 19.85, plot);
}

void TestPeriodSearchApp::testChanceProb() {
  setMethod("testChanceProb");

  // Vector to hold array of number of statistically independent trials to test chanceProb.
  std::vector<PeriodSearch::size_type>::size_type trial_size = 11;
  std::vector<PeriodSearch::size_type> num_indep_trial(trial_size, 0);
  num_indep_trial[1] = 1;
  num_indep_trial[2] = 2;
  num_indep_trial[3] = 10;
  for (std::vector<PeriodSearch::size_type>::size_type idx = 4; idx != trial_size; ++idx) {
    num_indep_trial[idx] = 10 * num_indep_trial[idx - 1];
  }

  // Vector to hold array of single-trial probabilities used to test chanceProb calculation.
  std::vector<double>::size_type prob_size = 201;
  std::vector<double> prob_one_trial(prob_size, 0.);
  for (std::vector<double>::size_type idx = 1; idx != prob_size; ++idx) {
    prob_one_trial[idx] = std::pow(.9, double(prob_size - (idx + 1)));
  }

  // Populate array with approximate answers using a standard math library call. Note that this is
  // inaccurate especially for probabilities near 0, and/or for large numbers of trials.
  std::vector<std::vector<double> > approx_chance_prob(trial_size, std::vector<double>(prob_size, 0.));
  for (std::vector<PeriodSearch::size_type>::size_type ii = 0; ii != trial_size; ++ii) {
    for (std::vector<double>::size_type jj = 0; jj != prob_size; ++jj) {
      approx_chance_prob[ii][jj] = 1. - std::pow(1. - prob_one_trial[jj], double(num_indep_trial[ii]));
    }
  }

  // Require the agreement between the approximate simple formula and the form used in the PeriodSearch class
  // to be to about 6.5 digits. Note that this limit cannot be refined because the approximate values are
  // not sufficiently accurate.
  double epsilon = 1.e-7;

  for (std::vector<PeriodSearch::size_type>::size_type ii = 0; ii != trial_size; ++ii) {
    for (std::vector<double>::size_type jj = 0; jj != prob_size; ++jj) {
      double chance_prob = PeriodSearch::computeChanceProbMultiTrial(prob_one_trial[jj], num_indep_trial[ii]);
      if (0. > chance_prob) {
        err() << "ERROR: computeChanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] <<
          ") unexpectedly returned " << chance_prob << ", which is < 0." << std::endl;
      } else if (1. < chance_prob) {
        err() << "ERROR: computeChanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] <<
          ") unexpectedly returned " << chance_prob << ", which is > 1." << std::endl;
      } else if ((0. == approx_chance_prob[ii][jj] && 0. != chance_prob) ||
        (0. != approx_chance_prob[ii][jj] && epsilon < std::fabs(chance_prob / approx_chance_prob[ii][jj] - 1.))) {
        err() << "ERROR: computeChanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] << ") returned " <<
          chance_prob << ", not " << approx_chance_prob[ii][jj] << ", as expected." << std::endl;
      }
    }
  }
}

void TestPeriodSearchApp::testChiSquaredTestArray() {
  setMethod("testChiSquaredTestArray");

  // Prepare a ChiSquaredTestArray object.
  const int num_trial = 3;
  const int num_bin = 4;
  ChiSquaredTestArray chi2_test(num_trial, num_bin);
  const StatisticViewer & chi2_viewer(chi2_test.getViewer());

  // Fill events.
  const int num_events = 10;
  for (int ii = 0; ii < num_events; ++ii) {
    double phase = (ii + 0.1) / num_events;
    chi2_test.fill(0, phase);
    chi2_test.fill(1, phase); chi2_test.fill(1, phase);
    chi2_test.fill(2, phase); chi2_test.fill(2, phase); chi2_test.fill(2, phase);
  }

  // Check size.
  int test_size = chi2_test.size();
  if (test_size != num_trial) {
    err() << "Size of the ChiSquaredTestArray was reported as " << test_size << ", not " << num_trial << "." << std::endl;
  }

  // Set the comparison precision.
  const double epsilon = 1.e-12;

  // Check the folded light curve.
  const StatisticViewer::data_type & result_counts = chi2_viewer.getData(1);
  const double expected_counts[] = {3., 2., 3., 2.};
  for (int trial = 0; trial < num_trial; ++trial) {
    chi2_test.updateViewer(trial);
    for (int ii = 0; ii < num_bin; ++ii) {
      double result = result_counts[ii];
      double expected = expected_counts[ii] * (trial + 1);
      if (std::fabs(result/expected - 1.) > epsilon) {
        err() << "Photon count in phase bin " << ii << " of trial " << trial << " was " << result
          << ", not " << expected << "." << std::endl;
      }
    }
  }

  // Check S-value.
  const double expected_values[] = {0.4, 0.8, 1.2};
  for (int trial = 0; trial < num_trial; ++trial) {
    double result = chi2_test.computeStat(trial);
    double expected = expected_values[trial];
    if (std::fabs(result/expected - 1.) > epsilon) {
      err() << "S-value for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }
}

void TestPeriodSearchApp::testZ2nTestArray() {
  setMethod("testZ2nTestArray");

  // Prepare a Z2nTestArray object.
  const int num_trial = 3;
  const int num_harm = 4;
  Z2nTestArray z2n_test(num_trial, num_harm);
  const StatisticViewer & z2n_viewer(z2n_test.getViewer());

  // Fill events.
  const int num_events = 12;
  for (int ii = 0; ii < num_events; ++ii) {
    z2n_test.fill(0, 0.123);
    double phase = (ii % 3 + 1) / 8.; // 45, 90, and 135 degrees.
    z2n_test.fill(1, phase);
    phase = ii / 12.; // one whole round.
    z2n_test.fill(2, phase);
  }

  // Check size.
  int test_size = z2n_test.size();
  if (test_size != num_trial) {
    err() << "Size of the Z2nTestArray was reported as " << test_size << ", not " << num_trial << "." << std::endl;
  }

  // Set the comparison precision.
  const double epsilon = 1.e-12;

  // Check the power spectrum.
  const StatisticViewer::data_type & result_powers = z2n_viewer.getData(1);
  const double norm = 2. / num_events;
  double power10 = 4. * (std::sqrt(2.) + 1.);
  power10 = norm * power10 * power10;
  double power11 = norm * 4. * 4.;
  double power12 = 4. * (std::sqrt(2.) - 1.);
  power12 = norm * power12 * power12;
  double power13 = norm * 4. * 4.;
  const double expected_powers[][num_harm] = { {24., 24., 24., 24.}, {power10, power11, power12, power13}, {0., 0., 0., 0.} };
  for (int trial = 0; trial < num_trial; ++trial) {
    z2n_test.updateViewer(trial);
    for (int ii = 0; ii < num_harm; ++ii) {
      double result = result_powers[ii];
      double expected = expected_powers[trial][ii];
      if ((expected == 0. && std::fabs(result) > epsilon)
          || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
        err() << "Power for harmonic nubmer " << ii << " of trial " << trial << " was " << result
          << ", not " << expected << "." << std::endl;
      }
    }
  }

  // Check Z2n-value.
  double expected_values[] = {96., power10 + power11 + power12 + power13, 0.};
  for (int trial = 0; trial < num_trial; ++trial) {
    double result = z2n_test.computeStat(trial);
    double expected = expected_values[trial];
    if ((expected == 0. && std::fabs(result) > epsilon)
        || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
      err() << "Z2n-value for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }
}

void TestPeriodSearchApp::testHTestArray() {
  setMethod("testHTestArray");

  // Prepare a HTestArray object.
  const int num_trial = 3;
  const int num_harm = 4;
  HTestArray h_test(num_trial, num_harm);
  const StatisticViewer & h_viewer(h_test.getViewer());

  // Fill events.
  const int num_events = 12;
  for (int ii = 0; ii < num_events; ++ii) {
    h_test.fill(0, 0.123);
    double phase = (ii % 3 + 1) / 8.; // 45, 90, and 135 degrees.
    h_test.fill(1, phase);
    phase = ii / 12.; // one whole round.
    h_test.fill(2, phase);
  }

  // Check size.
  int test_size = h_test.size();
  if (test_size != num_trial) {
    err() << "Size of the HTestArray was reported as " << test_size << ", not " << num_trial << "." << std::endl;
  }

  // Set the comparison precision.
  const double epsilon = 1.e-12;

  // Check the power spectrum.
  const StatisticViewer::data_type & result_powers = h_viewer.getData(1);
  const double norm = 2. / num_events;
  double power10 = 4. * (std::sqrt(2.) + 1.);
  power10 = norm * power10 * power10;
  double power11 = norm * 4. * 4.;
  double power12 = 4. * (std::sqrt(2.) - 1.);
  power12 = norm * power12 * power12;
  double power13 = norm * 4. * 4.;
  const double expected_powers[][num_harm] = {
    {24., 44., 64., 84.},
    {power10, power10 + power11 - 4., power10 + power11 + power12 - 8., power10 + power11 + power12 + power13 - 12.},
    {0., -4., -8., -12.}
  };
  for (int trial = 0; trial < num_trial; ++trial) {
    h_test.updateViewer(trial);
    for (int ii = 0; ii < num_harm; ++ii) {
      double result = result_powers[ii];
      double expected = expected_powers[trial][ii];
      if ((expected == 0. && std::fabs(result) > epsilon)
          || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
        err() << "Candidate H value for harmonic nubmer " << ii << " of trial " << trial << " was " << result
          << ", not " << expected << "." << std::endl;
      }
    }
  }

  // Check H-value.
  double expected_values[] = {84., power10, 0.};
  for (int trial = 0; trial < num_trial; ++trial) {
    double result = h_test.computeStat(trial);
    double expected = expected_values[trial];
    if ((expected == 0. && std::fabs(result) > epsilon)
        || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
      err() << "H-value for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }
}

void TestPeriodSearchApp::testRayleighTestArray() {
  setMethod("testRayleighTestArray");

  // Prepare a RayleighTestArray object.
  const int num_trial = 3;
  RayleighTestArray rayleigh_test(num_trial);
  const StatisticViewer & rayleigh_viewer(rayleigh_test.getViewer());

  // Fill events.
  const int num_events = 12;
  for (int ii = 0; ii < num_events; ++ii) {
    rayleigh_test.fill(0, 0.123);
    double phase = (ii % 3 + 1) / 8.; // 45, 90, and 135 degrees.
    rayleigh_test.fill(1, phase);
    phase = ii / 12.; // one whole round.
    rayleigh_test.fill(2, phase);
  }

  // Check size.
  int test_size = rayleigh_test.size();
  if (test_size != num_trial) {
    err() << "Size of the RayleighTestArray was reported as " << test_size << ", not " << num_trial << "." << std::endl;
  }

  // Set the comparison precision.
  const double epsilon = 1.e-12;

  // Check the power spectrum.
  const StatisticViewer::data_type & result_powers = rayleigh_viewer.getData(1);
  const double norm = 2. / num_events;
  double power10 = 4. * (std::sqrt(2.) + 1.);
  power10 = norm * power10 * power10;
  const double expected_powers[] = {24., power10, 0.};
  for (int trial = 0; trial < num_trial; ++trial) {
    rayleigh_test.updateViewer(trial);
    double result = result_powers[0];
    double expected = expected_powers[trial];
    if ((expected == 0. && std::fabs(result) > epsilon)
        || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
      err() << "Viewer content for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }

  // Check Rayleigh-value.
  double expected_values[] = {24., power10, 0.};
  for (int trial = 0; trial < num_trial; ++trial) {
    double result = rayleigh_test.computeStat(trial);
    double expected = expected_values[trial];
    if ((expected == 0. && std::fabs(result) > epsilon)
        || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
      err() << "Rayleigh statistic for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }
}

void TestPeriodSearchApp::testAllStats(const std::string & prefix, const std::vector<double> & events, double t_start, double t_stop,
  double center, double step, long num_trials, double epoch, int num_bins,
  double fourier_width, int fourier_num_bins, double fourier_min_freq, double fourier_max_freq, bool plot) {
  setMethod("testAllStats");

  double duration = t_stop - t_start;

  // Test ChiSquared case.
  ChiSquaredTestArray chi2_test(num_trials, num_bins);
  FoldingAnalysis chi2_search(&chi2_test, center, step, epoch, duration, "Hz");
  testOneSearch(events, chi2_search, "Folding Analysis: Chi Squared Statistic", prefix + "-chi-sq.fits", plot);

  // Test Z2n case.
  Z2nTestArray z2n_test(num_trials, num_bins);
  FoldingAnalysis z2n_search(&z2n_test, center, step, epoch, duration, "Hz");

  testOneSearch(events, z2n_search, "Folding Analysis: Z2n Statistic", prefix + "-z2n.fits", plot);

  // Test Rayleigh case.
  RayleighTestArray rayleigh_test(num_trials);
  FoldingAnalysis rayleigh_search(&rayleigh_test, center, step, epoch, duration, "Hz");

  testOneSearch(events, rayleigh_search, "Folding Analysis: Rayleigh Statistic", prefix + "-rayleigh.fits", plot);

  // Test H case.
  HTestArray h_test(num_trials, num_bins);
  FoldingAnalysis h_search(&h_test, center, step, epoch, duration, "Hz");

  testOneSearch(events, h_search, "Folding Analysis: H Statistic", prefix + "-h.fits", plot);

  // Create analysis object.
  FourierAnalysis fourier_search(t_start, t_stop, fourier_width, fourier_num_bins, "Hz", events.size());

  testOneSearch(events, fourier_search, "Fourier Analysis: Power Spectrum", prefix + "-fourier.fits", plot,
    fourier_min_freq, fourier_max_freq);
}

void TestPeriodSearchApp::testOneSearch(const std::vector<double> & events, PeriodSearch & search, const std::string & plot_title,
  const std::string & out_file, bool plot, double min_freq, double max_freq) {
  // Get the viewer.
  StatisticViewer & viewer = search.getViewer();

  // Fill the data into the search object.
  for (std::vector<double>::const_iterator itor = events.begin(); itor != events.end(); ++itor) {
    search.fill(*itor);
  }

  // Perform the search operation.
  search.computeStat();
  search.updateViewer(min_freq, max_freq);

  // Find the template file.
  std::string template_file = facilities::commonUtilities::joinPath(getDataPath(), "period-search-out.tpl");

  // Create output file.
  tip::IFileSvc::instance().createFile(out_file, template_file, true);

  // Open output file.
  tip::Table * out_table(tip::IFileSvc::instance().editTable(out_file, "POWER_SPECTRUM"));

  // Write the summary to the output header, and the data to the output table.
  viewer.write(*out_table);

  // Close output file. This is necessary before checking output file against reference file.
  delete out_table;

  // Check the result against its reference file in data/outref/ directory.
  checkOutputFits(out_file, true);

  // Plot if requested.
  viewer.setTitle(plot_title);
  viewer.setLabel(0, "Hz");
  if (plot) viewer.plot();
}

st_app::StAppFactory<TestPeriodSearchApp> g_factory("test_periodSearch");
