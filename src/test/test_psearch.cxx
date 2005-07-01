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
#include "pulsarDb/GlastTime.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Env.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "periodSearch/PeriodTest.h"
#include "ChiSquaredTest.h"
#include "HTest.h"
#include "RayleighTest.h"
#include "Z2nTest.h"

static const double s_sec_per_day = 86400.;
static const std::string s_cvs_id = "$Name:  $";

class PSearchTestApp : public st_app::StApp {
  public:
    PSearchTestApp(): m_os("PSearchTestApp", "", 2), m_data_dir(), m_failed(false) {
      setName("test_stpsearch");
      setVersion(s_cvs_id);
    }

    virtual ~PSearchTestApp() throw() {}
    virtual void run();

    void testAllStats(double center, double step, long num_trials, double epoch, int num_bins, const std::vector<double> & events,
      double duration, const std::string & unit);

    void testFindMax(const periodSearch::PeriodTest & test);

    void testChooseEph(const std::string & ev_file, const std::string & eph_file, const std::string pulsar_name, double epoch);

    const std::string & getDataDir();

    std::string findFile(const std::string & file_name);

    std::string makeTitle(const periodSearch::PeriodTest & test, const std::string & init_title);

  private:
    st_stream::StreamFormatter m_os;
    std::string m_data_dir;
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
  std::string unit = "(/s)";

  // Test process of picking the ephemeris.
  testChooseEph(findFile("ft1tiny.fits"), findFile("groD4-dc2v4.fits"), "crab", 23078385.922);

  if (m_failed) throw std::runtime_error("UNIT TEST FAILED");

  // First do simple test with this highly artificial data.
  testAllStats(central, step, num_pds, epoch, num_bins, fake_evts, duration, unit);

  // Data taken from M. Hirayama's work with modified ASCA data.
  // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/existing.html#tryout003
  num_pds = 101;

  central = 1. / 50.41843041e-3;
  step = .168e-7 * central * central;

  epoch = 23078385.922;

  fake_evts.clear();

  // Read some real data.
  std::auto_ptr<const tip::Table> evt_table(tip::IFileSvc::instance().readTable(findFile("step-01.fits"), "EVENTS"));

  // Read telapse for duration from GTI.
  std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(findFile("step-01.fits"), "GTI"));
  duration = 0.;
  gti_table->getHeader()["TELAPSE"].get(duration);

  double valid_since;
  double valid_until;

  gti_table->getHeader()["TSTART"].get(valid_since);
  gti_table->getHeader()["TSTOP"].get(valid_until);

  // Make the array big enough to hold these events.
  fake_evts.resize(evt_table->getNumRecords());

  std::vector<double>::iterator event_itor = fake_evts.begin();
  for (tip::Table::ConstIterator itor = evt_table->begin(); itor != evt_table->end(); ++itor, ++event_itor) {
    *event_itor = (*itor)["TIME"].get();
  }

  // Repeat simple test with this somewhat less artificial data.
  testAllStats(central, step, num_pds, epoch, num_bins, fake_evts, duration, unit);

  // Now test pdot correction.
  double phi0 = 0.; // Ignored for these purposes anyway.
  double pdot = 4.7967744e-13;
  double p2dot = 0.; // Not available.

  using namespace pulsarDb;

  // The following declarator looks like a function prototype.
  // PeriodEph eph(GlastTdbTime(valid_since), GlastTdbTime(valid_until), GlastTdbTime(epoch), phi0, 1. / central, pdot, p2dot);
  // Resolve the misunderstanding by using a temporary variable for the first argument.
  GlastTdbTime vs(valid_since);
  PeriodEph eph(vs, GlastTdbTime(valid_until), GlastTdbTime(epoch), phi0, 1. / central, pdot, p2dot);
  TimingModel timing_model;

  // Correct the data.
  for (std::vector<double>::iterator itor = fake_evts.begin(); itor != fake_evts.end(); ++itor) {
    GlastTdbTime tdb(*itor);
    timing_model.correctPdot(eph, tdb);
    *itor = tdb.elapsed();
  }

  // Repeat test with the pdot corrected data.
  testAllStats(central, step, num_pds, epoch, num_bins, fake_evts, duration, unit);

  // Test process of picking the ephemeris.
  testChooseEph(findFile("ft1tiny.fits"), findFile("groD4-dc2v4.fits"), "crab", epoch);
}

void PSearchTestApp::testAllStats(double center, double step, long num_trials, double epoch, int num_bins,
  const std::vector<double> & events, double duration, const std::string & unit) {
  m_os.setMethod("testAllStats");

  // Test ChiSquared case.
  ChiSquaredTest test(center, step, num_trials, epoch, num_bins, duration);

  // Iterate over the fake events.
  for (std::vector<double>::const_iterator itor = events.begin(); itor != events.end(); ++itor) {
    test.fill(*itor);
  }

  test.computeStats();

  m_os.out() << "Chi Squared Statistic" << std::endl;
  m_os.out() << test << std::endl;
  test.plotStats("Chi Squared Statistic", unit);

  // Test Z2n case.
  Z2nTest test_z2n(center, step, num_trials, epoch, num_bins, duration);

  // Iterate over the fake events.
  for (std::vector<double>::const_iterator itor = events.begin(); itor != events.end(); ++itor) {
    test_z2n.fill(*itor);
  }

  test_z2n.computeStats();

  m_os.out() << "Z2n Statistic" << std::endl;
  m_os.out() << test_z2n << std::endl;
  test_z2n.plotStats("Z2n Statistic", unit);

  // Test Rayleigh case.
  RayleighTest test_rayleigh(center, step, num_trials, epoch, duration);

  // Iterate over the fake events.
  for (std::vector<double>::const_iterator itor = events.begin(); itor != events.end(); ++itor) {
    test_rayleigh.fill(*itor);
  }

  test_rayleigh.computeStats();

  m_os.out() << "Rayleigh Statistic" << std::endl;
  m_os.out() << test_rayleigh << std::endl;
  test_rayleigh.plotStats("Rayleigh Statistic", unit);

  // Test H case.
  HTest test_h(center, step, num_trials, epoch, num_bins, duration);

  // Iterate over the fake events.
  for (std::vector<double>::const_iterator itor = events.begin(); itor != events.end(); ++itor) {
    test_h.fill(*itor);
  }

  test_h.computeStats();

  m_os.out() << "H Statistic" << std::endl;
  m_os.out() << test_h << std::endl;
  test_h.plotStats("H Statistic", unit);
}

void PSearchTestApp::testChooseEph(const std::string & ev_file, const std::string & eph_file, const std::string pulsar_name,
  double epoch) {
  using namespace pulsarDb;
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

  FrequencyEph freq = computer.calcPulsarEph(GlastTdbTime(epoch));

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
  const PulsarEph & chosen_eph(chooser.choose(computer.getPulsarEphCont(), GlastTdbTime(epoch)));

  double correct_f2 = chosen_eph.f2();
  if (fabs(correct_f2/freq.f2() - 1.) > epsilon) {
    m_failed = true;
    m_os.err() << "ERROR: in testChooseEph, f2 was computed to be " << freq.f2() << ", not " << correct_f2 << std::endl;
  }
}

const std::string & PSearchTestApp::getDataDir() {
  m_data_dir = st_facilities::Env::getDataDir("periodSearch");
  return m_data_dir;
}

std::string PSearchTestApp::findFile(const std::string & file_name) {
  return st_facilities::Env::appendFileName(getDataDir(), file_name);
}

std::string PSearchTestApp::makeTitle(const periodSearch::PeriodTest & test, const std::string & init_title) {
  std::ostringstream os;
  std::pair<double, double> max = test.findMax();
  std::pair<double, double> chance_prob = test.chanceProb(max.second);

  os << init_title << ", Max at: " << max.first << ", stat: " << max.second;

  // Massage display: if difference between min and max is small enough just use max.
  os.setf(std::ios::scientific);
  os.precision(2); // 2 digits -> < 1. e -4. limit in next line.
  if ((chance_prob.second - chance_prob.first) / chance_prob.second < 1.e-4) 
    os << ", chance prob: " << chance_prob.second;
  else
    os << ", chance prob < " << chance_prob.second;

  return os.str();
}

st_app::StAppFactory<PSearchTestApp> g_factory("gtpsearch");
