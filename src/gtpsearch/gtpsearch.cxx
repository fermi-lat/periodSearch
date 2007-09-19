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
#include "pulsarDb/PulsarToolApp.h"
#include "pulsarDb/TimingModel.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/EventTimeHandler.h"
#include "timeSystem/GlastTimeHandler.h"
#include "timeSystem/TimeSystem.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Env.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"

#include "periodSearch/PeriodSearch.h"
#include "periodSearch/PeriodSearchViewer.h"
#include "ChiSquaredTest.h"
#include "HTest.h"
#include "RayleighTest.h"
#include "Z2nTest.h"

using namespace periodSearch;
using namespace timeSystem;

static const std::string s_cvs_id = "$Name:  $";

class PSearchApp : public pulsarDb::PulsarToolApp {
  public:
    PSearchApp();
    virtual ~PSearchApp() throw();
    virtual void run();

    virtual void prompt(st_app::AppParGroup & pars);

  private:
    st_stream::StreamFormatter m_os;
};

PSearchApp::PSearchApp(): m_os("PSearchApp", "", 2) {
  st_app::AppParGroup & pars(getParGroup("gtpsearch"));

  setName("gtpsearch");
  setVersion(s_cvs_id);

  pars.setSwitch("ephstyle");
  pars.setSwitch("cancelpdot");
  pars.setSwitch("timeorigin");
  pars.setCase("ephstyle", "FREQ", "ra");
  pars.setCase("ephstyle", "FREQ", "dec");
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
  pars.setCase("ephstyle", "PER", "ra");
  pars.setCase("ephstyle", "PER", "dec");
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

PSearchApp::~PSearchApp() throw() {}

void PSearchApp::run() {
  // Clean up from any previous runs.
  resetApp();

  m_os.setMethod("run()");
  st_app::AppParGroup & pars(getParGroup("gtpsearch"));

  // Prompt for all parameters and save them.
  prompt(pars);

  // Get parameters.
  std::string out_file = pars["outfile"];
  double scan_step = pars["scanstep"];
  long num_trials = pars["numtrials"];
  long num_bins = pars["numbins"];
  bool plot = pars["plot"];
  std::string title = pars["title"];
  bool clobber = pars["clobber"];

  // Open the event file(s).
  openEventFile(pars);

  // Handle leap seconds.
  std::string leap_sec_file = pars["leapsecfile"];
  timeSystem::TimeSystem::setDefaultLeapSecFileName(leap_sec_file);

  // Setup time correction mode.
  defineTimeCorrectionMode("NONE", SUPPRESSED, SUPPRESSED, SUPPRESSED);
  defineTimeCorrectionMode("AUTO", ALLOWED,    ALLOWED,    ALLOWED);
  defineTimeCorrectionMode("BARY", REQUIRED,   SUPPRESSED, SUPPRESSED);
  defineTimeCorrectionMode("BIN",  REQUIRED,   REQUIRED,   SUPPRESSED);
  defineTimeCorrectionMode("PDOT", REQUIRED,   SUPPRESSED, REQUIRED);
  defineTimeCorrectionMode("ALL",  REQUIRED,   REQUIRED,   REQUIRED);

  
  selectTimeCorrectionMode(pars);

  // Set up EphComputer for arrival time corrections.
  pulsarDb::TimingModel model;
  pulsarDb::SloppyEphChooser chooser;
  initEphComputer(pars, model, chooser);

  // Use user input (parameters) together with computer to determine corrections to apply.
  bool guess_pdot = true;
  initTimeCorrection(pars, guess_pdot);

  // Compute central frequency of periodicity search, which is an expected pulse frequency at the time origin for the search.
  double origin = 0.;
  AbsoluteTime abs_origin = computeAbsoluteTime(origin);
  double f_center = getEphComputer().calcPulsarEph(abs_origin).f0();

  // Compute frequency step from scan step and the Fourier resolution == 1. / duration,
  // using start/stop of the observation interval with all corrections applied.
  double duration = computeElapsedSecond(getStopTime()) - computeElapsedSecond(getStartTime());
  double f_step = scan_step / duration;

  // Choose which kind of test to create.
  std::string algorithm = pars["algorithm"];
  for (std::string::iterator itor = algorithm.begin(); itor != algorithm.end(); ++itor) *itor = std::toupper(*itor);

  // Create the proper test.
  std::auto_ptr<PeriodSearch> search(0);
  if (algorithm == "CHI2")
    search.reset(new ChiSquaredTest(f_center, f_step, num_trials, origin, num_bins, duration));
  else if (algorithm == "H")
    search.reset(new HTest(f_center, f_step, num_trials, origin, num_bins, duration));
  else if (algorithm == "Z2N")
    search.reset(new Z2nTest(f_center, f_step, num_trials, origin, num_bins, duration));
  else throw std::runtime_error("PSearchApp: invalid test algorithm " + algorithm);

  for (setFirstEvent(); !isEndOfEventList(); setNextEvent()) {
    // Get event time as AbsoluteTime.
    AbsoluteTime abs_evt_time(getEventTime());

    // Convert event time to target time representation.
    double target_evt_time = computeElapsedSecond(abs_evt_time);

    // Fill into the test.
    search->fill(target_evt_time);
  }

  // Compute the statistics.
  search->computeStats();

  enum ChatLevel {
    eIncludeSummary= 2,
    eAllDetails = 3
  };

  // Use default title if user did not specify one.
  if (title == "DEFAULT") title = "Folding Analysis: " + algorithm + " Test";

  // Create a viewer for plotting and writing output.
  periodSearch::PeriodSearchViewer viewer(*search);

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

  // Display a plot, if desired.
  if (plot) viewer.plot(title, "(Hz)");
}

void PSearchApp::prompt(st_app::AppParGroup & pars) {
  // Prompt for most parameters automatically.
  pars.Prompt("algorithm");
  pars.Prompt("evfile");
  pars.Prompt("scfile");
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
    pars.Prompt("ra");
    pars.Prompt("dec");
    pars.Prompt("f0");
  } else if (eph_style == "PER") {
    pars.Prompt("ephepoch");
    pars.Prompt("timeformat");
    pars.Prompt("timesys");
    pars.Prompt("ra");
    pars.Prompt("dec");
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

  pars.Prompt("evtable");
  pars.Prompt("sctable");
  pars.Prompt("timefield");
  pars.Prompt("plot");
  pars.Prompt("title");
  pars.Prompt("leapsecfile");
  pars.Prompt("clobber");

  // Save current values of the parameters.
  pars.Save();
}

// List supported mission(s).
timeSystem::EventTimeHandlerFactory<timeSystem::GlastTimeHandler> glast_handler;

st_app::StAppFactory<PSearchApp> g_factory("gtpsearch");
