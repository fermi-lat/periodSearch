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

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Env.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/TimeRep.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "periodSearch/PeriodTest.h"
#include "periodSearch/PeriodSearchViewer.h"
#include "ChiSquaredTest.h"
#include "FourierAnalysis.h"
#include "HTest.h"
#include "RayleighTest.h"
#include "Z2nTest.h"

static const std::string s_cvs_id = "$Name:  $";

using namespace periodSearch;

class PSearchTestApp : public st_app::StApp {
  public:

    PSearchTestApp(): m_os("PSearchTestApp", "", 2), m_failed(false) {
      setName("test_stpsearch");
      setVersion(s_cvs_id);
    }

    virtual ~PSearchTestApp() throw() {}
    virtual void run();

    void testAllStats(const std::string & prefix, const std::vector<double> & events, double t_start, double t_stop,
      double center, double step, long num_trials, double epoch, int num_bins,
      double fourier_width, int fourier_num_bins, double fourier_min_freq, double fourier_max_freq, bool plot);

    void testFindMax(const periodSearch::PeriodTest & test);

    void testChooseEph(const std::string & ev_file, const std::string & eph_file, const std::string & pulsar_name, double epoch);

    void testChanceProb();

    std::string findFile(const std::string & file_name);

  private:
    void testOneSearch(const std::vector<double> & events, PeriodSearch & search,
      const std::string & text_title, const std::string & plot_title, const std::string & out_file,
      bool plot, double min_freq = -1., double max_freq = -1.);

    st_stream::StreamFormatter m_os;
    bool m_failed;
};

void PSearchTestApp::run() {

  // Trick up some fake events.
  int num_evts = 1000;

  std::vector<double> fake_evts(num_evts);
  m_os.out().precision(16);
  m_os.err().precision(16);

  for (int ii = 0; ii < num_evts; ++ii)
    fake_evts[ii] = ii + 1;

  // Set up search parameters.
  double central = 1.;
  double step = .5e-3;
  int num_pds = 61;
  double epoch = .2;
  int num_bins = 10;
  double duration = 1000.;

  // Test process of picking the ephemeris.
  testChooseEph(findFile("ft1tiny.fits"), findFile("groD4-dc2v4.fits"), "crab", 212380785.922);



  if (m_failed) throw std::runtime_error("UNIT TEST FAILED");

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
  std::auto_ptr<const tip::Table> evt_table(tip::IFileSvc::instance().readTable(findFile("step-01.fits"), "EVENTS"));

  // Get gti_table.
  std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(findFile("step-01.fits"), "GTI"));

  // Make the array big enough to hold these events.
  fake_evts.resize(evt_table->getNumRecords());

  // Correct event times for changed MJDREF.
  timeSystem::MetRep orig_glast_time("TDB", 54101, 0., 0.);
  timeSystem::MetRep current_glast_time("TDB", 51910, 0., 0.);
  std::vector<double>::iterator event_itor = fake_evts.begin();
  for (tip::Table::ConstIterator itor = evt_table->begin(); itor != evt_table->end(); ++itor, ++event_itor) {
    orig_glast_time.setValue((*itor)["TIME"].get());
    current_glast_time = timeSystem::AbsoluteTime(orig_glast_time);
    *event_itor = current_glast_time.getValue();
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
  orig_glast_time.setValue(tstart);
  current_glast_time = timeSystem::AbsoluteTime(orig_glast_time);
  tstart = current_glast_time.getValue();
  orig_glast_time.setValue(tstop);
  current_glast_time = timeSystem::AbsoluteTime(orig_glast_time);
  tstop = current_glast_time.getValue();

  // Repeat simple test with this somewhat less artificial data.
  // Note for Fourier test: width of .01 s -> Nyquist = 1/.02s = 50 Hz.
  testAllStats("psrb0540", fake_evts, tstart, tstop, central, step, num_pds, epoch, num_bins, .01, 1000000, 19.82, 19.85, plot);

  // Now test pdot correction.
  double phi0 = 0.; // Ignored for these purposes anyway.
  double pdot = 4.7967744e-13;
  double p2dot = 0.; // Not available.

  using namespace pulsarDb;
  using namespace timeSystem;

  MetRep glast_tdb("TDB", 51910, 0., 0.);
  glast_tdb.setValue(epoch);
  AbsoluteTime abs_epoch(glast_tdb);

  PeriodEph eph("TDB", abs_epoch, abs_epoch, abs_epoch, phi0, 1. / central, pdot, p2dot);
  TimingModel timing_model;

  // Correct the data.
  AbsoluteTime evt_time(glast_tdb);
  for (std::vector<double>::iterator itor = fake_evts.begin(); itor != fake_evts.end(); ++itor) {
    glast_tdb.setValue(*itor);
    evt_time = glast_tdb;
    timing_model.cancelPdot(eph, evt_time);
    glast_tdb = evt_time;
    *itor = glast_tdb.getValue();
  }

  // Cancel pdot in tstart, tstop and epoch to be consistent.
  glast_tdb.setValue(tstart);
  evt_time = glast_tdb;
  timing_model.cancelPdot(eph, evt_time);
  glast_tdb = evt_time;
  tstart = glast_tdb.getValue();

  glast_tdb.setValue(tstop);
  evt_time = glast_tdb;
  timing_model.cancelPdot(eph, evt_time);
  glast_tdb = evt_time;
  tstop = glast_tdb.getValue();

  glast_tdb.setValue(epoch);
  evt_time = glast_tdb;
  timing_model.cancelPdot(eph, evt_time);
  glast_tdb = evt_time;
  epoch = glast_tdb.getValue();

  // Repeat test with the pdot corrected data.
  testAllStats("psrb0540-pdot", fake_evts, tstart, tstop, central, step, num_pds, epoch, num_bins, .01, 1000000,
    19.82, 19.85, plot);

  // Test process of picking the ephemeris.
  testChooseEph(findFile("ft1tiny.fits"), findFile("groD4-dc2v4.fits"), "crab", epoch);

  // Test computations of chance probability.
  testChanceProb();
}

void PSearchTestApp::testAllStats(const std::string & prefix, const std::vector<double> & events, double t_start, double t_stop,
  double center, double step, long num_trials, double epoch, int num_bins,
  double fourier_width, int fourier_num_bins, double fourier_min_freq, double fourier_max_freq, bool plot) {
  using namespace periodSearch;

  m_os.setMethod("testAllStats");

  double duration = t_stop - t_start;

  // Test ChiSquared case.
  ChiSquaredTest test(center, step, num_trials, epoch, num_bins, duration);

  testOneSearch(events, test, "Chi Squared Statistic", "Folding Analysis: Chi Squared Statistic", prefix + "-chi-sq.fits",
    plot);

  // Test Z2n case.
  Z2nTest test_z2n(center, step, num_trials, epoch, num_bins, duration);

  testOneSearch(events, test_z2n, "Z2n Statistic", "Folding Analysis: Z2n Statistic", prefix + "-z2n.fits",
    plot);

  // Test Rayleigh case.
  RayleighTest test_rayleigh(center, step, num_trials, epoch, duration);

  testOneSearch(events, test_rayleigh, "Rayleigh Statistic", "Folding Analysis: Rayleigh Statistic", prefix + "-rayleigh.fits",
    plot);

  // Test H case.
  HTest test_h(center, step, num_trials, epoch, num_bins, duration);

  testOneSearch(events, test_h, "H Statistic", "Folding Analysis: H Statistic", prefix + "-h.fits",
    plot);

  // Create analysis object.
  FourierAnalysis test_fourier(t_start, t_stop, fourier_width, fourier_num_bins, events.size());

  testOneSearch(events, test_fourier, "Fourier Power", "Fourier Analysis: Power Spectrum", prefix + "-fourier.fits",
    plot, fourier_min_freq, fourier_max_freq);
}

void PSearchTestApp::testChooseEph(const std::string & ev_file, const std::string & eph_file, const std::string & pulsar_name,
  double epoch) {
  using namespace pulsarDb;
  using namespace timeSystem;
  using namespace tip;

  m_os.setMethod("testChooseEph");

  // Open event file.
  std::auto_ptr<const Table> ev_table(IFileSvc::instance().readTable(ev_file, "EVENTS"));

  // Need some keywords.
  const Header & header(ev_table->getHeader());
  double mjdref = 0.L;
  header["MJDREF"].get(mjdref);

  // Create a timing model object from which to compute the frequency.
  TimingModel model;
  SloppyEphChooser chooser;

  // Create a computer.
  EphComputer computer(model, chooser);

  // Get database access.
  PulsarDb db(eph_file);

  // Limit database to this pulsar only.
  db.filterName(pulsar_name);

  // Load ephemerides into computer.
  computer.load(db);

  MetRep glast_tdb("TDB", 51910, 0., epoch);
  FrequencyEph freq = computer.calcPulsarEph(AbsoluteTime(glast_tdb));

  const double epsilon = 1.e-8;

  double correct_f0 = 29.93633350069171;
  if (fabs(correct_f0/freq.f0() - 1.) > epsilon) {
    m_failed = true;
    m_os.err() << "f0 was computed to be " << freq.f0() << ", not " << correct_f0 << std::endl;
  }

  double correct_f1 = -3.772519499263467e-10;
  if (fabs(correct_f1/freq.f1() - 1.) > epsilon) {
    m_failed = true;
    m_os.err() << "f1 was computed to be " << freq.f1() << ", not " << correct_f1 << std::endl;
  }

  // Select the best ephemeris for this time.
  glast_tdb.setValue(epoch);
  const PulsarEph & chosen_eph(chooser.choose(computer.getPulsarEphCont(), AbsoluteTime(glast_tdb)));

  double correct_f2 = chosen_eph.f2();
  if (fabs(correct_f2/freq.f2() - 1.) > epsilon) {
    m_failed = true;
    m_os.err() << "ERROR: in testChooseEph, f2 was computed to be " << freq.f2() << ", not " << correct_f2 << std::endl;
  }
}

void PSearchTestApp::testChanceProb() {
  using namespace periodSearch;

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
      double chance_prob = PeriodSearch::chanceProbMultiTrial(prob_one_trial[jj], num_indep_trial[ii]);
      if (0. > chance_prob) {
        m_failed = true;
        m_os.err() << "ERROR: chanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] <<
          ") unexpectedly returned " << chance_prob << ", which is < 0." << std::endl;
      } else if (1. < chance_prob) {
        m_failed = true;
        m_os.err() << "ERROR: chanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] <<
          ") unexpectedly returned " << chance_prob << ", which is > 1." << std::endl;
      } else if ((0. == approx_chance_prob[ii][jj] && 0. != chance_prob) ||
        (0. != approx_chance_prob[ii][jj] && epsilon < std::fabs(chance_prob / approx_chance_prob[ii][jj] - 1.))) {
        m_failed = true;
        m_os.err() << "ERROR: chanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] << ") returned " <<
          chance_prob << ", not " << approx_chance_prob[ii][jj] << ", as expected." << std::endl;
      }
    }
  }
}


std::string PSearchTestApp::findFile(const std::string & file_name) {
  using namespace st_facilities;
  return Env::appendFileName(Env::getDataDir("periodSearch"), file_name);
}

void PSearchTestApp::testOneSearch(const std::vector<double> & events, PeriodSearch & search,
  const std::string & text_title, const std::string & plot_title, const std::string & out_file,
  bool plot, double min_freq, double max_freq) {

  PeriodSearchViewer viewer(search, min_freq, max_freq);

  // Fill the data into the search object.
  for (std::vector<double>::const_iterator itor = events.begin(); itor != events.end(); ++itor) {
    search.fill(*itor);
  }

  // Perform the search operation.
  search.computeStats();

  // Find the template file.
  std::string template_file = findFile("period-search-out.tpl");

  // Create output file.
  tip::IFileSvc::instance().createFile(out_file, template_file, true);

  // Open output file and get reference to its header.
  std::auto_ptr<tip::Table> out_table(tip::IFileSvc::instance().editTable(out_file, "POWER_SPECTRUM"));
  tip::Header & out_header(out_table->getHeader());

  // Write the summary to the output header, and the data to the output table.
  viewer.writeSummary(out_header);
  viewer.writeData(*out_table);

  enum ChatLevel {
    eIncludeSummary= 2,
    eAllDetails = 3
  };

  // Write the stats to the screen.
  m_os.info(eIncludeSummary) << text_title << std::endl;
  viewer.writeSummary(m_os.info(eIncludeSummary)) << std::endl;

  // Write details of test result if chatter is high enough.
  viewer.writeData(m_os.info(eAllDetails)) << std::endl;

  if (plot) viewer.plot(plot_title, "(Hz)");
}

st_app::StAppFactory<PSearchTestApp> g_factory("test_periodSearch");
