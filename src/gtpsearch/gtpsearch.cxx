/** \file gtpsearch.cxx
    \brief Period search tool.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "pulsarDb/AbsoluteTime.h"
#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/GlastTime.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

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

static const std::string s_cvs_id = "$Name:  $";

class PSearchApp : public st_app::StApp {
  public:
    PSearchApp();
    virtual ~PSearchApp() throw();
    virtual void run();

    virtual void prompt(st_app::AppParGroup & pars);

    pulsarDb::FrequencyEph estimateFrequency(const std::string & psrdb_file, const std::string & psr_name,
      double epoch, double mjdref);

    const std::string & getDataDir();

  private:
    st_stream::StreamFormatter m_os;
    std::string m_data_dir;
    PeriodTest * m_test;
};

PSearchApp::PSearchApp(): m_os("PSearchApp", "", 2), m_data_dir(), m_test(0) {
  st_app::AppParGroup & pars(getParGroup("gtpsearch"));

  setName("gtpsearch");
  setVersion(s_cvs_id);

  pars.setSwitch("ephstyle");
  pars.setSwitch("correctpdot");
  pars.setCase("ephstyle", "FREQ", "f0");
  pars.setCase("ephstyle", "FREQ", "correctpdot");
  pars.setCase("ephstyle", "PER", "correctpdot");
  pars.setCase("correctpdot", "true", "f1");
  pars.setCase("ephstyle", "FREQ", "f1");
  pars.setCase("ephstyle", "FREQ", "f2");
  pars.setCase("ephstyle", "PER", "p0");
  pars.setCase("ephstyle", "PER", "p1");
  pars.setCase("correctpdot", "true", "p1");
  pars.setCase("ephstyle", "PER", "p2");
  pars.setCase("ephstyle", "DB", "psrname");
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
  std::string demod_binary_string = pars["demodbin"];

  using namespace pulsarDb;

  // Ignored but needed for timing model.
  double phi0 = 0.;

  // Find the pulsar database.
  std::string psrdb_file = pars["psrdbfile"];
  std::string psrdb_file_uc = psrdb_file;
  for (std::string::iterator itor = psrdb_file_uc.begin(); itor != psrdb_file_uc.end(); ++itor) *itor = toupper(*itor);
  if ("DEFAULT" == psrdb_file_uc) {
    using namespace st_facilities;
    psrdb_file = Env::appendFileName(Env::getDataDir("pulsePhase"), "master_pulsardb.fits");
  }
  std::string psr_name = pars["psrname"];
  std::string demod_bin_string = pars["demodbin"];
  for (std::string::iterator itor = demod_bin_string.begin(); itor != demod_bin_string.end(); ++itor) *itor = toupper(*itor);
  
  // Open the test file.
  std::auto_ptr<const tip::Table> event_table(tip::IFileSvc::instance().readTable(event_file, event_extension));

  // Read GTI.
  std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(event_file, "GTI"));

  // Get necessary keywords.
  double mjdref = 0.;
  double duration = 0.;
  double valid_since = 0.;
  double valid_until = 0.;
  std::string telescope;
  std::string time_sys;
  const tip::Header & header(gti_table->getHeader());
  header["MJDREF"].get(mjdref);
  header["TELAPSE"].get(duration);
  header["TSTART"].get(valid_since);
  header["TSTOP"].get(valid_until);
  header["TELESCOP"].get(telescope);
  header["TIMESYS"].get(time_sys);

  if (telescope != "GLAST") throw std::runtime_error("Only GLAST supported for now");

  std::auto_ptr<AbsoluteTime> abs_epoch(0);
  std::auto_ptr<AbsoluteTime> abs_valid_since(0);
  std::auto_ptr<AbsoluteTime> abs_valid_until(0);
  if (time_sys == "TDB") {
    abs_epoch.reset(new GlastTdbTime(epoch));
    abs_valid_since.reset(new GlastTdbTime(valid_since));
    abs_valid_until.reset(new GlastTdbTime(valid_until));
  } else if (time_sys == "TT") {
    abs_epoch.reset(new GlastTtTime(epoch));
    abs_valid_since.reset(new GlastTtTime(valid_since));
    abs_valid_until.reset(new GlastTtTime(valid_until));
  } else {
    throw std::runtime_error("Only TDB or TT time systems are supported for now");
  }

  // Compute frequency step from scan step and the Fourier resolution == 1. / duration.
  if (0. >= duration) throw std::runtime_error("TELAPSE for data is not positive!");
  double f_step = scan_step / duration;

  // A TimingModel will be needed for several steps below.
  TimingModel model;
  SloppyEphChooser chooser;
  EphComputer computer(model, chooser);

  // Open the database.
  PulsarDb database(psrdb_file);

  // Select only ephemerides for this pulsar.
  database.filterName(psr_name);

  // Load the selected ephemerides.
  computer.load(database);

  // Determine whether to perform binary demodulation.
  bool demod_bin = false;
  if (demod_bin_string != "NO") {
    // User selected not "no", so attempt to perform demodulation
    if (!computer.getOrbitalEphCont().empty()) {
      demod_bin = true;
    } else if (demod_bin_string == "YES") {
      throw std::runtime_error("Binary demodulation was required by user, but no orbital ephemeris was found");
    }
  }

  // Handle either period or frequency-style input.
  std::string eph_style = pars["ephstyle"];
  if (eph_style == "FREQ") {
    double f0 = pars["f0"];
    double f1 = pars["f1"];
    double f2 = pars["f2"];

    if (0. >= f0) throw std::runtime_error("Frequency must be positive.");

    // Override any ephemerides which may have been found in the database with the ephemeris the user provided.
    PulsarEphCont & ephemerides(computer.getPulsarEphCont());
    ephemerides.clear();
    ephemerides.push_back(FrequencyEph(*abs_valid_since, *abs_valid_until, *abs_epoch, phi0, f0, f1, f2).clone());
  } else if (eph_style == "PER") {
    double p0 = pars["p0"];
    double p1 = pars["p1"];
    double p2 = pars["p2"];

    if (0. >= p0) throw std::runtime_error("Period must be positive.");

    // Override any ephemerides which may have been found in the database with the ephemeris the user provided.
    PulsarEphCont & ephemerides(computer.getPulsarEphCont());
    ephemerides.clear();
    ephemerides.push_back(PeriodEph(*abs_valid_since, *abs_valid_until, *abs_epoch, phi0, p0, p1, p2).clone());
  } else if (eph_style == "DB") {
    // No action needed.
  }

  // Choose which kind of test to create.
  std::string algorithm = pars["algorithm"];

  for (std::string::iterator itor = algorithm.begin(); itor != algorithm.end(); ++itor) *itor = toupper(*itor);

  // Compute an ephemeris at abs_epoch to use for the test.
  std::auto_ptr<PulsarEph> eph(computer.calcPulsarEph(*abs_epoch).clone());

  // Create the proper test.
  if (algorithm == "CHI2") 
    m_test = new ChiSquaredTest(eph->f0(), f_step, num_trials, epoch, num_bins, duration);
  else if (algorithm == "H") 
    m_test = new HTest(eph->f0(), f_step, num_trials, epoch, num_bins, duration);
  else if (algorithm == "Z2N") 
    m_test = new Z2nTest(eph->f0(), f_step, num_trials, epoch, num_bins, duration);
  else throw std::runtime_error("PSearchApp: invalid test algorithm " + algorithm);

  // Iterate over table, filling the search/test object with temporal data.
  for (tip::Table::ConstIterator itor = event_table->begin(); itor != event_table->end(); ++itor) {
    // Get value from the table.
    double evt_time = (*itor)[time_col].get();

    if (time_sys == "TDB") {
      GlastTdbTime tdb(evt_time);
      // Perform binary correction if so desired.
      if (demod_bin) computer.demodulateBinary(tdb);

      // Perform pdot correction if so desired.
      // For efficiency use the TimingModel directly here, instead of using the EphComputer.
      if (correct_pdot) computer.cancelPdot(tdb);
      evt_time = tdb.elapsed();
    } else {
      GlastTtTime tt(evt_time);
      // Perform binary correction if so desired.
      if (demod_bin) computer.demodulateBinary(tt);

      // Perform pdot correction if so desired.
      // For efficiency use the TimingModel directly here, instead of using the EphComputer.
      if (correct_pdot) computer.cancelPdot(tt);
      evt_time = tt.elapsed();
    }

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
  pars.Prompt("psrdbfile");
  pars.Prompt("psrname");
  pars.Prompt("ephstyle");
  pars.Prompt("demodbin");

  std::string eph_style = pars["ephstyle"];
  if (eph_style == "FREQ") {
    pars.Prompt("epoch");
    pars.Prompt("f0");
  } else if (eph_style == "PER") {
    pars.Prompt("epoch");
    pars.Prompt("p0");
  } else if (eph_style == "DB") {
    // No action needed.
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
    // if eph_style == "DB", coeffs will be determined.
  }

  pars.Prompt("numtrials");
  pars.Prompt("numbins");
  pars.Prompt("timecol");
  pars.Prompt("plot");
  pars.Prompt("title");

  // Save current values of the parameters.
  pars.Save();
}

const std::string & PSearchApp::getDataDir() {
  m_data_dir = st_facilities::Env::getDataDir("periodSearch");
  return m_data_dir;
}

st_app::StAppFactory<PSearchApp> g_factory("gtpsearch");
