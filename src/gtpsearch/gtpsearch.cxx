/** \file gtpsearch.cxx
    \brief Period search tool.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"

#include "pulsePhase/TimingModel.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Env.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"

#include "periodSearch/PeriodTest.h"
#include "ChiSquaredTest.h"
#include "HTest.h"
#include "RayleighTest.h"
#include "Z2nTest.h"

using namespace periodSearch;

class PSearchApp : public st_app::StApp {
  public:
    PSearchApp();
    virtual ~PSearchApp() throw();
    virtual void run();

    virtual void prompt(st_app::AppParGroup & pars);

    pulsePhase::TimingModel::FrequencyCoeff estimateFrequency(const std::string & psrdb_file, const std::string & psr_name,
      double epoch, double mjdref);

    const std::string & getDataDir();

  private:
    st_stream::StreamFormatter m_os;
    std::string m_data_dir;
    PeriodTest * m_test;
};

PSearchApp::PSearchApp(): m_os("PSearchApp", "", 2), m_data_dir(), m_test(0) {
}

PSearchApp::~PSearchApp() throw() {
  delete m_test;
}

void PSearchApp::run() {
  m_os.setMethod("run()");
  st_app::AppParGroup & pars(getParGroup("gtpsearch"));

  // Prompt for all parameters and save them.
  prompt(pars);

  // Get parameters.
  std::string event_file = pars["evfile"];
  std::string event_extension = pars["evtable"];
  double scan_step = pars["scanstep"];
  long num_trials = pars["numtrials"];
  double epoch = pars["epoch"];
  long num_bins = pars["numbins"];
  std::string time_col = pars["timecol"];
  bool plot = pars["plot"];
  std::string title = pars["title"];
  bool correct_pdot = pars["correctpdot"];

  using namespace pulsePhase;

  // Store frequency info in a timing model
  std::auto_ptr<TimingModel> timing_model(0);

  // Ignored but needed for timing model.
  double phi0 = 0.;

  // Open the test file.
  std::auto_ptr<const tip::Table> event_table(tip::IFileSvc::instance().readTable(event_file, event_extension));

  // Read GTI.
  std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(event_file, "GTI"));

  // Get necessary keywords.
  double mjdref = 0.;
  double duration = 0.;
  const tip::Header & header(gti_table->getHeader());
  header["MJDREF"].get(mjdref);
  header["TELAPSE"].get(duration);

  // Compute frequency step from scan step and the Fourier resolution == 1. / duration.
  if (0. >= duration) throw std::runtime_error("TELAPSE for data is not positive!");
  double f_step = scan_step / duration;

  // Handle either period or frequency-style input.
  std::string eph_style = pars["ephstyle"];
  if (eph_style == "FREQ") {
    double f0 = pars["f0"];
    double f1 = pars["f1"];
    double f2 = pars["f2"];
    timing_model.reset(new TimingModel(epoch, phi0, TimingModel::FrequencyCoeff(f0, f1, f2)));
  } else if (eph_style == "PER") {
    double p0 = pars["p0"];
    double p1 = pars["p1"];
    double p2 = pars["p2"];
    timing_model.reset(new TimingModel(epoch, phi0, TimingModel::PeriodCoeff(p0, p1, p2)));
  } else if (eph_style == "FILE") {
    std::string psrdb_file = pars["psrdbfile"];
    std::string psr_name = pars["psrname"];
    timing_model.reset(new TimingModel(epoch, phi0, estimateFrequency(psrdb_file, psr_name, epoch, mjdref)));
  }

  // Choose which kind of test to create.
  std::string algorithm = pars["algorithm"];

  // Create the proper test.
  if (algorithm == "Chi2") 
    m_test = new ChiSquaredTest(timing_model->getF0(), f_step, num_trials, epoch, num_bins, duration);
  else if (algorithm == "H") 
    m_test = new HTest(timing_model->getF0(), f_step, num_trials, epoch, num_bins, duration);
  else if (algorithm == "Z2n") 
    m_test = new Z2nTest(timing_model->getF0(), f_step, num_trials, epoch, num_bins, duration);
  else throw std::runtime_error("PSearchApp: invalid test algorithm " + algorithm);

  // Iterate over table, filling the search/test object with temporal data.
  for (tip::Table::ConstIterator itor = event_table->begin(); itor != event_table->end(); ++itor) {
    // Get value from the table.
    double evt_time = (*itor)[time_col].get();

    // Perform pdot correction if so desired.
    if (correct_pdot) evt_time = timing_model->calcPdotCorr(evt_time);

    // Fill into the test.
    m_test->fill(evt_time);
  }

  // Compute the statistics.
  m_test->computeStats();

  // Use default title if user did not specify one.
  if (title == "DEFAULT") title = algorithm + " Test";

  // Write the stats to the screen.
  m_os.out() << title << std::endl;
  m_os.out() << *m_test << std::endl;

  // TODO: When tip supports getting units from a column, replace the following:
  std::string unit = "(/s)";
  // with:
  // std::string unit = "(/" + event_table->getColumn(time_col)->getUnit() + ")";
  // Display a plot, if desired.
  if (plot) m_test->plotStats(title, unit);
}

void PSearchApp::prompt(st_app::AppParGroup & pars) {
  // Prompt for most parameters automatically.
  pars.Prompt("algorithm");
  pars.Prompt("evfile");
  pars.Prompt("evtable");
  pars.Prompt("ephstyle");

  std::string eph_style = pars["ephstyle"];
  if (eph_style == "FREQ")
    pars.Prompt("f0");
  else if (eph_style == "PER")
    pars.Prompt("p0");
  else if (eph_style == "FILE") {
    pars.Prompt("psrdbfile");
    pars.Prompt("psrname");
  } else
    throw std::runtime_error("Unknown ephemeris style " + eph_style);

  pars.Prompt("scanstep");
  pars.Prompt("correctpdot");

  // Only prompt for f1 & f2 / p1 & p2 if pdot correction is selected.
  if (true == bool(pars["correctpdot"])) {
    if (eph_style == "FREQ") {
      pars.Prompt("f1");
      pars.Prompt("f2");
    } else if (eph_style == "PER") {
      pars.Prompt("p1");
      pars.Prompt("p2");
    }
    // if eph_style == "FILE", coeffs will be determined.
  }

  pars.Prompt("numtrials");
  pars.Prompt("epoch");
  pars.Prompt("numbins");
  pars.Prompt("timecol");
  pars.Prompt("plot");
  pars.Prompt("title");

  // Save current values of the parameters.
  pars.Save();
}

pulsePhase::TimingModel::FrequencyCoeff PSearchApp::estimateFrequency(const std::string & psrdb_file, const std::string & psr_name,
  double epoch, double mjdref) {
  using namespace pulsePhase;
  using namespace pulsarDb;
  using namespace tip;
  static const double s_sec_per_day = 86400.;

  // Get database access.
  PulsarDb db(psrdb_file);

  // Limit database to this pulsar only.
  db.filterName(psr_name);

  // Select the best ephemeris for this time.
  PulsarEph chosen_eph(db.chooseEph(mjdref + (long double)(epoch) / s_sec_per_day, true));

  // Create a timing model object from which to compute the frequency.
  TimingModel model((chosen_eph.m_t0 - mjdref) * s_sec_per_day, 0., chosen_eph.m_f0, chosen_eph.m_f1, chosen_eph.m_f2);

  // Get estimated frequency coefficients.
  return model.calcFreq(epoch);
}

const std::string & PSearchApp::getDataDir() {
  m_data_dir = st_facilities::Env::getDataDir("gtpsearch");
  return m_data_dir;
}

st_app::StAppFactory<PSearchApp> g_factory("gtpsearch");
