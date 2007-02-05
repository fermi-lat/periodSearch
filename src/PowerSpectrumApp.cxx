/** \file PowerSpectrumApp.cxx
    \brief Implmentation of PowerSpectrumApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PowerSpectrumApp.h"

#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "st_app/AppParGroup.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/GlastMetRep.h"
#include "timeSystem/TimeRep.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Env.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"
#include "timeSystem/TimeSystem.h"

#include "periodSearch/PeriodSearch.h"
#include "FourierAnalysis.h"

using namespace periodSearch;
using namespace timeSystem;

static const std::string s_cvs_id = "$Name:  $";

PowerSpectrumApp::PowerSpectrumApp(): m_os("PowerSpectrumApp", "", 2), m_data_dir(), m_test(0) {
  st_app::AppParGroup & pars(getParGroup("gtpspec"));
  pars.setSwitch("ephstyle");
  pars.setSwitch("cancelpdot");
  pars.setCase("ephstyle", "FREQ", "f0");
  pars.setCase("ephstyle", "DB", "cancelpdot");
  pars.setCase("ephstyle", "FREQ", "cancelpdot");
  pars.setCase("ephstyle", "PER", "cancelpdot");
  pars.setCase("cancelpdot", "true", "f1");
  pars.setCase("ephstyle", "FREQ", "f1");
  pars.setCase("ephstyle", "FREQ", "f2");
  pars.setCase("ephstyle", "FREQ", "ephepoch");
  pars.setCase("ephstyle", "FREQ", "timeformat");
  pars.setCase("ephstyle", "FREQ", "timesys");
  pars.setCase("ephstyle", "PER", "p0");
  pars.setCase("ephstyle", "PER", "p1");
  pars.setCase("ephstyle", "PER", "ephepoch");
  pars.setCase("ephstyle", "PER", "timeformat");
  pars.setCase("ephstyle", "PER", "timesys");
  pars.setCase("cancelpdot", "true", "p1");
  pars.setCase("ephstyle", "PER", "p2");
  pars.setCase("ephstyle", "DB", "psrname");
}

PowerSpectrumApp::~PowerSpectrumApp() throw() {
  delete m_test;
}

void PowerSpectrumApp::run() {
  m_os.setMethod("run()");
  st_app::AppParGroup & pars(getParGroup("gtpspec"));

  // Prompt for all parameters and save them.
  prompt(pars);

  // Get parameters.
  std::string event_file = pars["evfile"];
  std::string event_extension = pars["evtable"];
  std::string epoch = pars["ephepoch"];
  std::string epoch_time_format = pars["timeformat"];
  long num_bins = pars["numbins"];
  double low_f_cut = pars["lowfcut"];
  std::string time_field = pars["timefield"];
  bool plot = pars["plot"];
  std::string title = pars["title"];
  bool cancel_pdot = pars["cancelpdot"];
  std::string demod_binary_string = pars["demodbin"];
  std::string eph_style = pars["ephstyle"];

  // Open the event file.
  std::auto_ptr<const tip::Table> event_table(tip::IFileSvc::instance().readTable(event_file, event_extension));

  // Read GTI.
  std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(event_file, "GTI"));

  // Get necessary keywords.
  double tstart = 0.;
  double tstop = 0.;
  std::string telescope;
  std::string event_time_sys;

  // Find all other keywords from events extension.
  const tip::Header & header(event_table->getHeader());
  header["TSTART"].get(tstart);
  header["TSTOP"].get(tstop);
  header["TELESCOP"].get(telescope);
  header["TIMESYS"].get(event_time_sys);

  // Get the mjdref from the header, which is not as simple as just reading a single keyword.
  MjdRefDatabase mjd_ref_db;
  IntFracPair mjd_ref(mjd_ref_db(header));

  // Make names of time system and mission case insensitive.
  for (std::string::iterator itor = telescope.begin(); itor != telescope.end(); ++itor) *itor = std::toupper(*itor);
  for (std::string::iterator itor = event_time_sys.begin(); itor != event_time_sys.end(); ++itor) *itor = std::toupper(*itor);

  if (telescope != "GLAST") throw std::runtime_error("Only GLAST supported for now");

  // Handle leap seconds.
  std::string leap_sec_file = pars["leapsecfile"];
  timeSystem::TimeSystem::setDefaultLeapSecFileName(leap_sec_file);

  // Make time formats etc. case insensitive.
  for (std::string::iterator itor = epoch_time_format.begin(); itor != epoch_time_format.end(); ++itor) *itor = std::toupper(*itor);
  for (std::string::iterator itor = eph_style.begin(); itor != eph_style.end(); ++itor) *itor = std::toupper(*itor);

  // Determine the time system used for the ephemeris epoch.
  std::string epoch_time_sys;
  if (eph_style == "DB") epoch_time_sys = "TDB";
  else epoch_time_sys = pars["timesys"].Value();

  // Make time formats etc. case insensitive.
  for (std::string::iterator itor = epoch_time_sys.begin(); itor != epoch_time_sys.end(); ++itor) *itor = std::toupper(*itor);

  // Interpret FILE option: match epoch to event time system.
  if ("FILE" == epoch_time_sys) epoch_time_sys = event_time_sys;

  using namespace pulsarDb;

  // Ignored but needed for timing model.
  double phi0 = 0.;

  std::string psr_name = pars["psrname"];
  std::string demod_bin_string = pars["demodbin"];
  for (std::string::iterator itor = demod_bin_string.begin(); itor != demod_bin_string.end(); ++itor) *itor = std::toupper(*itor);

  std::string evt_time_sys;
  header["TIMESYS"].get(evt_time_sys);

  // Set up event time representation.
  MetRep evt_time_rep(evt_time_sys, mjd_ref, 0.);
  evt_time_rep.setValue(tstart);
  AbsoluteTime abs_tstart(evt_time_rep);
  evt_time_rep.setValue(tstop);
  AbsoluteTime abs_tstop(evt_time_rep);
  evt_time_rep.setValue(.5 * (tstart + tstop));
  AbsoluteTime abs_middle(evt_time_rep);

  // A TimingModel will be needed for several steps below.
  TimingModel model;
  SloppyEphChooser chooser;
  EphComputer computer(model, chooser);

  if (eph_style != "DB") {
    std::auto_ptr<TimeRep> epoch_rep(0);
    // Create representation for this time format and time system.
    if (epoch_time_format == "FILE") {
      epoch_rep.reset(new MetRep(epoch_time_sys, mjd_ref, 0.));
    } else if (epoch_time_format == "GLAST") {
      epoch_rep.reset(new GlastMetRep(epoch_time_sys, 0.));
    } else if (epoch_time_format == "MJD") {
      epoch_rep.reset(new MjdRep(epoch_time_sys, 0, 0.));
    } else {
      throw std::runtime_error("Time format \"" + epoch_time_format + "\" is not supported for manual ephemeris epoch");
    }

    // Assign the ephepoch supplied by the user to the representation.
    epoch_rep->assign(epoch);

    AbsoluteTime abs_epoch(*epoch_rep);

    // Handle either period or frequency-style input.
    if (eph_style == "FREQ") {
      double f0 = pars["f0"];
      double f1 = pars["f1"];
      double f2 = pars["f2"];

      if (0. >= f0) throw std::runtime_error("Frequency must be positive.");

      // Override any ephemerides which may have been found in the database with the ephemeris the user provided.
      PulsarEphCont & ephemerides(computer.getPulsarEphCont());
      ephemerides.push_back(FrequencyEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, phi0, f0, f1, f2).clone());
    } else if (eph_style == "PER") {
      double p0 = pars["p0"];
      double p1 = pars["p1"];
      double p2 = pars["p2"];

      if (0. >= p0) throw std::runtime_error("Period must be positive.");

      // Override any ephemerides which may have been found in the database with the ephemeris the user provided.
      PulsarEphCont & ephemerides(computer.getPulsarEphCont());
      ephemerides.push_back(PeriodEph(epoch_time_sys, abs_epoch, abs_epoch, abs_epoch, phi0, p0, p1, p2).clone());
    }
  }

  if (eph_style == "DB" || demod_bin_string != "NO") {
    // Find the pulsar database.
    std::string psrdb_file = pars["psrdbfile"];
    std::string psrdb_file_uc = psrdb_file;
    for (std::string::iterator itor = psrdb_file_uc.begin(); itor != psrdb_file_uc.end(); ++itor) *itor = std::toupper(*itor);
    if ("DEFAULT" == psrdb_file_uc) {
      using namespace st_facilities;
      psrdb_file = Env::appendFileName(Env::getDataDir("periodSearch"), "master_pulsardb.fits");
    }

    // Open the database.
    PulsarDb database(psrdb_file);

    // Select only ephemerides for this pulsar.
    database.filterName(psr_name);

    // Load the selected ephemerides.
    if (eph_style == "DB") computer.loadPulsarEph(database);
    computer.loadOrbitalEph(database);
  }

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

  double bin_width = pars["binwidth"];

  // Create the proper test.
  m_test = new FourierAnalysis(tstart, tstop, bin_width, num_bins);

  // Iterate over table, filling the search/test object with temporal data.
  for (tip::Table::ConstIterator itor = event_table->begin(); itor != event_table->end(); ++itor) {
    // Get value from the table.
    double evt_time = (*itor)[time_field].get();

    evt_time_rep.setValue(evt_time);
    AbsoluteTime abs_evt_time(evt_time_rep);
    // Perform binary correction if so desired.
    if (demod_bin) computer.demodulateBinary(abs_evt_time);

    // Perform pdot correction if so desired.
    // For efficiency use the TimingModel directly here, instead of using the EphComputer.
    if (cancel_pdot) computer.cancelPdot(abs_evt_time);
    evt_time_rep = abs_evt_time;
    evt_time = evt_time_rep.getValue();

    // Fill into the test.
    m_test->fill(evt_time);
  }

  // Compute the statistics.
  m_test->computeStats();

  enum ChatLevel {
    eDefault = 2,
    eIncludeSummary= 3,
    eAllDetails = 4,
  };

  // Use default title if user did not specify one.
  if (title == "DEFAULT") title = "FFT Test";

  // Write the stats to the screen.
  m_os.info(eIncludeSummary) << title << std::endl;
  m_test->writeRange(m_os.info(eAllDetails), low_f_cut);

  // TODO: When tip supports getting units from a column, replace the following:
  std::string unit = "(Hz)";
  // with:
  // std::string unit = "(/" + event_table->getColumn(time_field)->getUnit() + ")";
  // Display a plot, if desired.
  if (plot) m_test->plotRange(title, unit, low_f_cut);
}

void PowerSpectrumApp::prompt(st_app::AppParGroup & pars) {
  // Prompt for most parameters automatically.
  pars.Prompt("evfile");
  pars.Prompt("evtable");
  pars.Prompt("psrdbfile");
  pars.Prompt("psrname");
  pars.Prompt("binwidth");
  pars.Prompt("numbins");
  pars.Prompt("lowfcut");
  pars.Prompt("ephstyle");
  pars.Prompt("demodbin");

  std::string eph_style = pars["ephstyle"];
  if (eph_style == "FREQ") {
    pars.Prompt("ephepoch");
    pars.Prompt("timeformat");
    pars.Prompt("timesys");
    pars.Prompt("f0");
  } else if (eph_style == "PER") {
    pars.Prompt("ephepoch");
    pars.Prompt("timeformat");
    pars.Prompt("timesys");
    pars.Prompt("p0");
  } else if (eph_style == "DB") {
    // No action needed.
  } else
    throw std::runtime_error("Unknown ephemeris style " + eph_style);

  pars.Prompt("cancelpdot");

  // Only prompt for f1 & f2 / p1 & p2 if pdot correction is selected.
  if (true == bool(pars["cancelpdot"])) {
    if (eph_style == "FREQ") {
      pars.Prompt("f1");
      pars.Prompt("f2");
    } else if (eph_style == "PER") {
      pars.Prompt("p1");
      pars.Prompt("p2");
    }
    // if eph_style == "DB", coeffs will be determined.
  }

  pars.Prompt("timefield");
  pars.Prompt("plot");
  pars.Prompt("title");
  pars.Prompt("leapsecfile");

  // Save current values of the parameters.
  pars.Save();
}

const std::string & PowerSpectrumApp::getDataDir() {
  m_data_dir = st_facilities::Env::getDataDir("periodSearch");
  return m_data_dir;
}
