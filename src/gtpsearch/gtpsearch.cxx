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

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/PulsarDb.h"
#include "pulsarDb/PulsarEph.h"
#include "pulsarDb/TimingModel.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

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
#include "periodSearch/PeriodSearchViewer.h"
#include "ChiSquaredTest.h"
#include "HTest.h"
#include "RayleighTest.h"
#include "Z2nTest.h"

using namespace periodSearch;
using namespace timeSystem;

static const std::string s_cvs_id = "$Name:  $";

class PSearchApp : public st_app::StApp {
  public:
    PSearchApp();
    virtual ~PSearchApp() throw();
    virtual void run();

    virtual void prompt(st_app::AppParGroup & pars);

    const std::string & getDataDir();

  private:
    st_stream::StreamFormatter m_os;
    std::string m_data_dir;
    PeriodSearch * m_test;
};

PSearchApp::PSearchApp(): m_os("PSearchApp", "", 2), m_data_dir(), m_test(0) {
  st_app::AppParGroup & pars(getParGroup("gtpsearch"));

  setName("gtpsearch");
  setVersion(s_cvs_id);

  pars.setSwitch("ephstyle");
  pars.setSwitch("cancelpdot");
  pars.setSwitch("timeorigin");
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
  pars.setCase("timeorigin", "USER", "usertime");
  pars.setCase("timeorigin", "USER", "userformat");
  pars.setCase("timeorigin", "USER", "usersys");
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
  std::string epoch = pars["ephepoch"];
  std::string epoch_time_format = pars["timeformat"];
  long num_bins = pars["numbins"];
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
  double duration = 0.;
  double tstart = 0.;
  double tstop = 0.;
  std::string telescope;
  std::string event_time_sys;
  // Find TELAPSE from GTI extension.
  gti_table->getHeader()["TELAPSE"].get(duration);

  // Find all other keywords from events extension.
  const tip::Header & header(event_table->getHeader());
  header["TELESCOP"].get(telescope);
  header["TIMESYS"].get(event_time_sys);

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
  
  // Compute frequency step from scan step and the Fourier resolution == 1. / duration.
  if (0. >= duration) throw std::runtime_error("TELAPSE for data is not positive!");
  double f_step = scan_step / duration;

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

  // Set up event time representations.
  MetRep evt_time_rep(evt_time_sys, mjd_ref, 0.);

  // Set up start time of observation (data set). Apply necessary corrections.
  evt_time_rep.setValue(tstart);
  AbsoluteTime abs_tstart(evt_time_rep);
  if (demod_bin) computer.demodulateBinary(abs_tstart);
  if (cancel_pdot) computer.cancelPdot(abs_tstart);
  // Propagate corrected time back to the original tstart variable.
  evt_time_rep = abs_tstart;
  tstart = evt_time_rep.getValue();

  // Set up stop time of observation (data set). Apply necessary corrections.
  evt_time_rep.setValue(tstop);
  AbsoluteTime abs_tstop(evt_time_rep);
  if (demod_bin) computer.demodulateBinary(abs_tstop);
  if (cancel_pdot) computer.cancelPdot(abs_tstop);
  // Propagate corrected time back to the original tstop variable.
  evt_time_rep = abs_tstop;
  tstop = evt_time_rep.getValue();

  // Compute a time to resresent the center of the observation.
  // It is not necessary to perform binary demodulation or pdot cancellation again because
  // the tstart and tstop already had these corrections.
  evt_time_rep.setValue(.5 * (tstart + tstop));
  AbsoluteTime abs_middle(evt_time_rep);

  // Handle styles of origin input.
  std::string origin_style = pars["timeorigin"];
  for (std::string::iterator itor = origin_style.begin(); itor != origin_style.end(); ++itor) *itor = std::toupper(*itor);
  AbsoluteTime abs_origin("TDB", Duration(0, 0.), Duration(0, 0.));
  std::string origin_time_sys;
  if (origin_style == "START") {
    // Get time of origin and its time system from event file.
    abs_origin = abs_tstart;
    origin_time_sys = event_time_sys;
  } else if (origin_style == "STOP") {
    // Get time of origin and its time system from event file.
    abs_origin = abs_tstop;
    origin_time_sys = event_time_sys;
  } else if (origin_style == "MIDDLE") {
    // Get time of origin and its time system from event file.
    abs_origin = abs_middle;
    origin_time_sys = event_time_sys;
  } else if (origin_style == "USER") {
    // Get time of origin and its format and system from parameters.
    std::string origin_time = pars["usertime"];
    std::string origin_time_format = pars["userformat"];
    origin_time_sys = pars["usersys"].Value();

    // Make case insensitive.
    for (std::string::iterator itor = origin_time_format.begin(); itor != origin_time_format.end(); ++itor) *itor = std::toupper(*itor);
    for (std::string::iterator itor = origin_time_sys.begin(); itor != origin_time_sys.end(); ++itor) *itor = std::toupper(*itor);

    // Interpret FILE option: match origin to event time system.
    if ("FILE" == origin_time_sys) origin_time_sys = event_time_sys;

    // Set up the origin using the given time system.
    std::auto_ptr<TimeRep> origin_rep(0);
    // Create representation for this time format and time system.
    if (origin_time_format == "FILE") {
      origin_rep.reset(new MetRep(origin_time_sys, mjd_ref, 0.));
    } else if (origin_time_format == "GLAST") {
      origin_rep.reset(new GlastMetRep(origin_time_sys, 0.));
    } else if (origin_time_format == "MJD") {
      origin_rep.reset(new MjdRep(origin_time_sys, 0, 0.));
    } else {
      throw std::runtime_error("Time format \"" + origin_time_format + "\" is not supported for user time origin");
    }

    // Assign the user time supplied by the user to the representation.
    origin_rep->assign(origin_time);

    abs_origin = *origin_rep;
  } else {
    throw std::runtime_error("Unsupported origin style " + origin_style);
  }
  
  // Choose which kind of test to create.
  std::string algorithm = pars["algorithm"];

  for (std::string::iterator itor = algorithm.begin(); itor != algorithm.end(); ++itor) *itor = std::toupper(*itor);

  // Compute an ephemeris at abs_origin to use for the test.
  std::auto_ptr<PulsarEph> eph(computer.calcPulsarEph(abs_origin).clone());

  // Reset computer to contain only the corrected ephemeris which was just computed.
  PulsarEphCont & ephemerides(computer.getPulsarEphCont());
  ephemerides.clear();
  ephemerides.push_back(eph->clone());

  // Convert absolute origin to the time system demanded by event file.
  evt_time_rep = abs_origin;
  double origin = evt_time_rep.getValue();

  // Create the proper test.
  if (algorithm == "CHI2") 
    m_test = new ChiSquaredTest(eph->f0(), f_step, num_trials, origin, num_bins, duration);
  else if (algorithm == "H") 
    m_test = new HTest(eph->f0(), f_step, num_trials, origin, num_bins, duration);
  else if (algorithm == "Z2N") 
    m_test = new Z2nTest(eph->f0(), f_step, num_trials, origin, num_bins, duration);
  else throw std::runtime_error("PSearchApp: invalid test algorithm " + algorithm);

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
  if (title == "DEFAULT") title = "Folding Analysis: " + algorithm + " Test";

  // Create a viewer for plotting and writing output.
  periodSearch::PeriodSearchViewer viewer(*m_test);

  // Write the stats to the screen.
  m_os.info(eIncludeSummary) << title << std::endl;
  viewer.writeSummary(m_os.info(eIncludeSummary)) << std::endl;

  // Write details of test result if chatter is high enough.
  viewer.writeData(m_os.info(eAllDetails)) << std::endl;

  // TODO: When tip supports getting units from a column, replace the following:
  std::string unit = "(Hz)";
  // with:
  // std::string unit = "(/" + event_table->getColumn(time_field)->getUnit() + ")";
  // Display a plot, if desired.
  if (plot) viewer.plot(title, unit);
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

  pars.Prompt("scanstep");
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

  pars.Prompt("numtrials");
  pars.Prompt("numbins");

  pars.Prompt("timeorigin");
  std::string origin_style = pars["timeorigin"];
  for (std::string::iterator itor = origin_style.begin(); itor != origin_style.end(); ++itor) *itor = std::toupper(*itor);
  if (origin_style == "USER") {
    pars.Prompt("usertime");
    pars.Prompt("userformat");
    pars.Prompt("usersys");
  }

  pars.Prompt("timefield");
  pars.Prompt("plot");
  pars.Prompt("title");
  pars.Prompt("leapsecfile");

  // Save current values of the parameters.
  pars.Save();
}

const std::string & PSearchApp::getDataDir() {
  m_data_dir = st_facilities::Env::getDataDir("periodSearch");
  return m_data_dir;
}

st_app::StAppFactory<PSearchApp> g_factory("gtpsearch");
