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
    typedef std::vector<const tip::Table *> table_cont_type;

    PSearchApp();
    virtual ~PSearchApp() throw();
    virtual void run();

    virtual void prompt(st_app::AppParGroup & pars);

    virtual std::auto_ptr<TimeRep> createTimeRep(const std::string & time_format, const std::string & time_system,
      const std::string & time_value);

    virtual std::auto_ptr<TimeRep> createTimeRep(const std::string & time_format, const std::string & time_system,
      const std::string & time_value, const tip::Header & header);

    void openEventFile(const st_app::AppParGroup & pars, table_cont_type & event_table_cont, table_cont_type & gti_table_cont);

    void initEphComputer(const st_app::AppParGroup & pars, const tip::Header & header, pulsarDb::EphComputer & computer);

    void computeTimeBoundary(const table_cont_type & gti_table_cont, bool demod_bin, bool cancel_pdot,
      pulsarDb::EphComputer & computer, AbsoluteTime & abs_tstart, AbsoluteTime & abs_tstop);

    AbsoluteTime computeTimeOrigin(const st_app::AppParGroup & pars, const tip::Header & header,
      const AbsoluteTime & abs_tstart, const AbsoluteTime & abs_tstop, timeSystem::TimeRep & time_rep);

    void updateEphComputer(const AbsoluteTime & abs_time, pulsarDb::EphComputer & computer);

    void initTimeCorrection(const st_app::AppParGroup & pars, const pulsarDb::EphComputer & computer,
      bool & demod_bin, bool & cancel_pdot);

    AbsoluteTime applyTimeCorrection(double time_value, TimeRep & time_rep, const pulsarDb::EphComputer & computer,
      bool demod_bin, bool cancel_pdot);

    double computeTimeValue(const AbsoluteTime & abs_time, TimeRep & time_rep);

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
  std::string out_file = pars["outfile"];
  double scan_step = pars["scanstep"];
  long num_trials = pars["numtrials"];
  long num_bins = pars["numbins"];
  std::string time_field = pars["timefield"];
  bool plot = pars["plot"];
  std::string title = pars["title"];
  bool clobber = pars["clobber"];

  // Open the event file(s).
  table_cont_type event_table_cont;
  table_cont_type gti_table_cont;
  openEventFile(pars, event_table_cont, gti_table_cont);

  const tip::Table * reference_table = event_table_cont.at(0);

  // Identify mission and time system from events extension.
  const tip::Header & reference_header(reference_table->getHeader());

  // Handle leap seconds.
  std::string leap_sec_file = pars["leapsecfile"];
  timeSystem::TimeSystem::setDefaultLeapSecFileName(leap_sec_file);

  // Set up EphComputer for arrival time corrections.
  pulsarDb::TimingModel model;
  pulsarDb::SloppyEphChooser chooser;
  pulsarDb::EphComputer computer(model, chooser);
  initEphComputer(pars, reference_header, computer);

  // Use user input (parameters) together with computer to determine corrections to apply.
  bool demod_bin = false;
  bool cancel_pdot = false;
  initTimeCorrection(pars, computer, demod_bin, cancel_pdot);

// START HERE
// Note: "target time rep" is the common TimeRep used to compute the time series to be analyzed.
// o another helper to determine the target TimeRep: always using MetRep using TDB unless no corrections. If no
//     corrections, require all input files to have same time system and match target to input event time system.

  // Determine start/stop of the observation interval in AbsoluteTime.
  AbsoluteTime abs_tstart("TDB", Duration(0, 0.), Duration(0, 0.));
  AbsoluteTime abs_tstop("TDB", Duration(0, 0.), Duration(0, 0.));
  computeTimeBoundary(gti_table_cont, demod_bin, cancel_pdot, computer, abs_tstart, abs_tstop);

  // Set up target time representation, used to compute the time series to analyze.
  std::string target_time_system = "TDB";
  if (!demod_bin && !cancel_pdot) {
    // TODO When multiple files are used, check to ensure all have same time system.
    // For now one event file, automatically self-consistent:
    reference_header["TIMESYS"].get(target_time_system);
  }
  // TODO: compute mjd_ref from abs_tstart.
  Duration mjd_ref(51910, 64.184);
  IntFracPair mjd_ref_int_frac = mjd_ref.getValue(Day);
  std::auto_ptr<TimeRep> target_time_rep(new MetRep(target_time_system, mjd_ref_int_frac, 0.));

  // Compute time origin for periodicity search, both in AbsoluteTime and in double.
  AbsoluteTime abs_origin = computeTimeOrigin(pars, reference_header, abs_tstart, abs_tstop, *target_time_rep);
  double origin = computeTimeValue(abs_origin, *target_time_rep);

  // Compute spin ephemeris to be used in periodicity search and pdot cancellation, and replace PulsarEph in EphComputer with it.
  updateEphComputer(abs_origin, computer);

  // Get central frequency of periodicity search.
  double f_center = computer.choosePulsarEph(abs_origin).f0();

  // Compute frequency step from scan step and the Fourier resolution == 1. / duration,
  double duration = computeTimeValue(abs_tstop, *target_time_rep) - computeTimeValue(abs_tstart, *target_time_rep);
  double f_step = scan_step / duration;

  // Choose which kind of test to create.
  std::string algorithm = pars["algorithm"];
  for (std::string::iterator itor = algorithm.begin(); itor != algorithm.end(); ++itor) *itor = std::toupper(*itor);

  // Create the proper test.
  if (algorithm == "CHI2")
    m_test = new ChiSquaredTest(f_center, f_step, num_trials, origin, num_bins, duration);
  else if (algorithm == "H")
    m_test = new HTest(f_center, f_step, num_trials, origin, num_bins, duration);
  else if (algorithm == "Z2N")
    m_test = new Z2nTest(f_center, f_step, num_trials, origin, num_bins, duration);
  else throw std::runtime_error("PSearchApp: invalid test algorithm " + algorithm);

  for (table_cont_type::iterator table_itor = event_table_cont.begin(); table_itor != event_table_cont.end(); ++table_itor) {
    const tip::Table * event_table = *table_itor;
    const tip::Header & header(event_table->getHeader());
    std::auto_ptr<TimeRep> evt_time_rep(createTimeRep("FILE", "FILE", "0.", header));

    // Iterate over table, filling the search/test object with temporal data.
    for (tip::Table::ConstIterator itor = event_table->begin(); itor != event_table->end(); ++itor) {
      // Get value from the table.
      double evt_time = (*itor)[time_field].get();

      // Apply time corrections, convert back out to a double value.
      AbsoluteTime abs_evt_time(applyTimeCorrection(evt_time, *evt_time_rep, computer, demod_bin, cancel_pdot));
      double target_evt_time = computeTimeValue(abs_evt_time, *target_time_rep);

      // Fill into the test.
      m_test->fill(target_evt_time);
    }
  }

  // Compute the statistics.
  m_test->computeStats();

  enum ChatLevel {
    eIncludeSummary= 2,
    eAllDetails = 3
  };

  // Use default title if user did not specify one.
  if (title == "DEFAULT") title = "Folding Analysis: " + algorithm + " Test";

  // Create a viewer for plotting and writing output.
  periodSearch::PeriodSearchViewer viewer(*m_test);

  // Interpret output file parameter.
  std::string out_file_uc = out_file;
  for (std::string::iterator itor = out_file_uc.begin(); itor != out_file_uc.end(); ++itor) *itor = std::toupper(*itor);

  if ("NONE" != out_file_uc) {
    // Find the template file.
    using namespace st_facilities;
    std::string template_file = Env::appendFileName(Env::getDataDir("periodSearch"), "period-search-out.tpl");

    // Create output file.
    tip::IFileSvc::instance().createFile(out_file, template_file, clobber);

    // Open output file and get reference to its header.
    std::auto_ptr<tip::Table> out_table(tip::IFileSvc::instance().editTable(out_file, "POWER_SPECTRUM"));
    tip::Header & out_header(out_table->getHeader());

    // Write the summary to the output header, and the data to the output table.
    viewer.writeSummary(out_header);
    viewer.writeData(*out_table);
  }

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
  pars.Prompt("outfile");
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
  pars.Prompt("clobber");

  // Save current values of the parameters.
  pars.Save();
}

std::auto_ptr<TimeRep> PSearchApp::createTimeRep(const std::string & time_format, const std::string & time_system,
  const std::string & time_value) {
  std::auto_ptr<TimeRep> time_rep(0);

  // Make upper case copies of input for case insensitive comparisons.
  std::string time_format_uc(time_format);
  for (std::string::iterator itor = time_format_uc.begin(); itor != time_format_uc.end(); ++itor) *itor = std::toupper(*itor);

  std::string time_system_uc(time_system);
  for (std::string::iterator itor = time_system_uc.begin(); itor != time_system_uc.end(); ++itor) *itor = std::toupper(*itor);

  // Create representation for this time format and time system.
  if ("GLAST" == time_format_uc) {
    time_rep.reset(new GlastMetRep(time_system, 0.));
  } else if ("MJD" == time_format_uc) {
    time_rep.reset(new MjdRep(time_system, 0, 0.));
  } else {
    throw std::runtime_error("Time format \"" + time_format + "\" is not supported for ephemeris time");
  }

  // Assign the ephtime supplied by the user to the representation.
  time_rep->assign(time_value);

  return time_rep;
}

std::auto_ptr<TimeRep> PSearchApp::createTimeRep(const std::string & time_format, const std::string & time_system,
  const std::string & time_value, const tip::Header & header) {
  std::auto_ptr<TimeRep> time_rep(0);

  // Make upper case copies of input for case insensitive comparisons.
  std::string time_format_uc(time_format);
  for (std::string::iterator itor = time_format_uc.begin(); itor != time_format_uc.end(); ++itor) *itor = std::toupper(*itor);

  std::string time_system_uc(time_system);
  for (std::string::iterator itor = time_system_uc.begin(); itor != time_system_uc.end(); ++itor) *itor = std::toupper(*itor);

  // Make a local modifiable copy to hold the rationalized time system.
  std::string time_system_rat(time_system);

  // First check whether time system should be read from the tip::Header.
  if ("FILE" == time_system_uc) header["TIMESYS"].get(time_system_rat);

  // Create representation for this time format and time system.
  if ("FILE" == time_format_uc) {
    // Check TELESCOP keyword for supported missions.
    std::string telescope;
    header["TELESCOP"].get(telescope);
    for (std::string::iterator itor = telescope.begin(); itor != telescope.end(); ++itor) *itor = std::toupper(*itor);
    if (telescope != "GLAST") throw std::runtime_error("Only GLAST supported for now");

    // Get the mjdref from the header, which is not as simple as just reading a single keyword.
    MjdRefDatabase mjd_ref_db;
    IntFracPair mjd_ref(mjd_ref_db(header));
    time_rep.reset(new MetRep(time_system_rat, mjd_ref, 0.));
  } else {
    // Delegate to overload that does not use tip.
    return createTimeRep(time_format, time_system_rat, time_value);
  }

  // Assign the time supplied by the user to the representation.
  time_rep->assign(time_value);

  return time_rep;
}

void PSearchApp::openEventFile(const st_app::AppParGroup & pars, table_cont_type & event_table_cont,
  table_cont_type & gti_table_cont) {
  std::string event_file = pars["evfile"];
  std::string event_extension = pars["evtable"];

  // Clear out any gti tables already in gti_table_cont.
  for (table_cont_type::reverse_iterator itor = gti_table_cont.rbegin(); itor != gti_table_cont.rend(); ++itor) {
    delete *itor;
  }
  gti_table_cont.clear();

  // Clear out any event tables already in event_table_cont.
  for (table_cont_type::reverse_iterator itor = event_table_cont.rbegin(); itor != event_table_cont.rend(); ++itor) {
    delete *itor;
  }
  event_table_cont.clear();

  // Open the event table.
  const tip::Table * event_table(tip::IFileSvc::instance().readTable(event_file, event_extension));

  // Add the table to the container.
  event_table_cont.push_back(event_table);

  // Open the GTI table.
  const tip::Table * gti_table(tip::IFileSvc::instance().readTable(event_file, "GTI"));

  // Add the table to the container.
  gti_table_cont.push_back(gti_table);

}

void PSearchApp::initEphComputer(const st_app::AppParGroup & pars, const tip::Header & header, pulsarDb::EphComputer & computer) {
  using namespace periodSearch;
  using namespace pulsarDb;

  std::string eph_style = pars["ephstyle"];
  for (std::string::iterator itor = eph_style.begin(); itor != eph_style.end(); ++itor) *itor = std::toupper(*itor);

  // Determine the time system used for the ephemeris epoch.
  std::string epoch_time_sys;
  if (eph_style == "DB") epoch_time_sys = "TDB";
  else epoch_time_sys = pars["timesys"].Value();

  // Ignored but needed for timing model.
  double phi0 = 0.;

  std::string psr_name = pars["psrname"];

  if (eph_style != "DB") {
    std::string epoch_time_format = pars["timeformat"];
    std::string epoch = pars["ephepoch"];
    AbsoluteTime abs_epoch(*createTimeRep(epoch_time_format, epoch_time_sys, epoch, header));

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

  std::string demod_bin_string = pars["demodbin"];
  for (std::string::iterator itor = demod_bin_string.begin(); itor != demod_bin_string.end(); ++itor) *itor = std::toupper(*itor);
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
}

void PSearchApp::computeTimeBoundary(const PSearchApp::table_cont_type & gti_table_cont, bool demod_bin, bool cancel_pdot,
  pulsarDb::EphComputer & computer, AbsoluteTime & abs_tstart, AbsoluteTime & abs_tstop) {
  bool candidate_found = false;

  // First, look for first and last times in the GTI.
  for (table_cont_type::const_iterator itor = gti_table_cont.begin(); itor != gti_table_cont.end(); ++itor) {
    const tip::Table & gti_table = *(*itor);

    // If possible, get tstart and tstop from first and last interval in GTI extension.
    tip::Table::ConstIterator gti_itor = gti_table.begin();
    if (gti_itor != gti_table.end()) {
      // Get start of the first interval and stop of last interval in GTI.
      double gti_start = (*gti_itor)["START"].get();
      gti_itor = gti_table.end();
      --gti_itor;
      double gti_stop = (*gti_itor)["STOP"].get();

      // Set up event time representation using GTI header.
      std::auto_ptr<TimeRep> time_rep(createTimeRep("FILE", "FILE", "0.", gti_table.getHeader()));

      // Correct the time.
      AbsoluteTime abs_gti_start = applyTimeCorrection(gti_start, *time_rep, computer, demod_bin, cancel_pdot);
      AbsoluteTime abs_gti_stop = applyTimeCorrection(gti_stop, *time_rep, computer, demod_bin, cancel_pdot);

      if (candidate_found) {
        // See if current interval extends the observation.
        if (abs_gti_start < abs_tstart) abs_tstart = abs_gti_start;
        if (abs_gti_stop > abs_tstop) abs_tstop = abs_gti_stop;
      } else {
        // First candidate is the current interval.
        abs_tstart = abs_gti_start;
        abs_tstop = abs_gti_stop;
        candidate_found = true;
      }
    }
  }

  // In the unlikely event that there were no GTI files, no intervals in the GTI, and no event files, this is a
  // serious problem.
  if (!candidate_found) throw std::runtime_error("computeTimeBoundary: cannot determine start/stop of the observation interval");
}


AbsoluteTime PSearchApp::computeTimeOrigin(const st_app::AppParGroup & pars, const tip::Header & header,
  const AbsoluteTime & abs_tstart, const AbsoluteTime & abs_tstop, timeSystem::TimeRep & time_rep) {

  // Handle styles of origin input.
  std::string origin_style = pars["timeorigin"];
  for (std::string::iterator itor = origin_style.begin(); itor != origin_style.end(); ++itor) *itor = std::toupper(*itor);
  AbsoluteTime abs_origin("TDB", Duration(0, 0.), Duration(0, 0.));
  if (origin_style == "START") {
    // Get time of origin and its time system from event file.
    abs_origin = abs_tstart;
  } else if (origin_style == "STOP") {
    // Get time of origin and its time system from event file.
    abs_origin = abs_tstop;
  } else if (origin_style == "MIDDLE") {
    // Use the center of the observation as the time origin.
    double tstart = computeTimeValue(abs_tstart, time_rep);
    double tstop = computeTimeValue(abs_tstop, time_rep);
    time_rep.set("TIME", .5 * (tstart + tstop));
    abs_origin = time_rep;
  } else if (origin_style == "USER") {
    // Get time of origin and its format and system from parameters.
    std::string origin_time = pars["usertime"];
    std::string origin_time_format = pars["userformat"];
    std::string origin_time_sys = pars["usersys"].Value();

    abs_origin = *createTimeRep(origin_time_format, origin_time_sys, origin_time, header);

  } else {
    throw std::runtime_error("Unsupported origin style " + origin_style);
  }

  return abs_origin;
}

void PSearchApp::updateEphComputer(const AbsoluteTime & abs_origin, pulsarDb::EphComputer & computer) {
  using namespace pulsarDb;

  // Compute an ephemeris at abs_origin to use for the test.
  std::auto_ptr<PulsarEph> eph(computer.calcPulsarEph(abs_origin).clone());

  // Reset computer to contain only the corrected ephemeris which was just computed.
  PulsarEphCont & ephemerides(computer.getPulsarEphCont());
  ephemerides.clear();
  ephemerides.push_back(eph->clone());
}

void PSearchApp::initTimeCorrection(const st_app::AppParGroup & pars, const pulsarDb::EphComputer & computer,
  bool & demod_bin, bool & cancel_pdot) {
  std::string demod_bin_string = pars["demodbin"];
  for (std::string::iterator itor = demod_bin_string.begin(); itor != demod_bin_string.end(); ++itor) *itor = std::toupper(*itor);

  // Determine whether to perform binary demodulation.
  demod_bin = false;
  if (demod_bin_string != "NO") {
    // User selected not "no", so attempt to perform demodulation
    if (!computer.getOrbitalEphCont().empty()) {
      demod_bin = true;
    } else if (demod_bin_string == "YES") {
      throw std::runtime_error("Binary demodulation was required by user, but no orbital ephemeris was found");
    }
  }

  // Determine whether to cancel pdot.
  cancel_pdot = bool(pars["cancelpdot"]);

}

AbsoluteTime PSearchApp::applyTimeCorrection(double time_value, TimeRep & time_rep, const pulsarDb::EphComputer & computer,
  bool demod_bin, bool cancel_pdot) {
  // Assign the value to the time representation.
  time_rep.set("TIME", time_value);

  // Convert TimeRep into AbsoluteTime so that computer can perform the necessary corrections.
  AbsoluteTime abs_time(time_rep);

  // Apply selected corrections.
  if (demod_bin) computer.demodulateBinary(abs_time);
  if (cancel_pdot) computer.cancelPdot(abs_time);

  return abs_time;
}

double PSearchApp::computeTimeValue(const AbsoluteTime & abs_time, TimeRep & time_rep) {
  double time_value = 0.;

  // Assign the absolute time to the time representation.
  time_rep = abs_time;

  // Get value from the time representation.
  time_rep.get("TIME", time_value);

  return time_value;
}

const std::string & PSearchApp::getDataDir() {
  m_data_dir = st_facilities::Env::getDataDir("periodSearch");
  return m_data_dir;
}

st_app::StAppFactory<PSearchApp> g_factory("gtpsearch");
