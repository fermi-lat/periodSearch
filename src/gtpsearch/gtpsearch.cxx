/** \file gtpsearch.cxx
    \brief Period search tool.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include <cctype>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "facilities/commonUtilities.h"

#include "pulsarDb/EphChooser.h"
#include "pulsarDb/EphComputer.h"
#include "pulsarDb/EphStatus.h"
#include "pulsarDb/PulsarToolApp.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/EventTimeHandler.h"
#include "timeSystem/GlastTimeHandler.h"
#include "timeSystem/TimeSystem.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"

#include "ChiSquaredTestArray.h"
#include "FoldingAnalysis.h"
#include "HTestArray.h"
#include "PeriodSearch.h"
#include "RayleighTestArray.h"
#include "StatisticViewer.h"
#include "Z2nTestArray.h"

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
  setName("gtpsearch");
  setVersion(s_cvs_id);

  st_app::AppParGroup & pars(getParGroup());

  pars.setSwitch("algorithm");
  pars.setSwitch("ephstyle");
  pars.setSwitch("timeorigin");
  pars.setCase("algorithm", "CHI2", "numphase");
  pars.setCase("algorithm", "Z2N", "numharm");
  pars.setCase("algorithm", "H", "maxharm");
  pars.setCase("ephstyle", "FREQ", "ephepoch");
  pars.setCase("ephstyle", "FREQ", "timeformat");
  pars.setCase("ephstyle", "FREQ", "timesys");
  pars.setCase("ephstyle", "FREQ", "ra");
  pars.setCase("ephstyle", "FREQ", "dec");
  pars.setCase("ephstyle", "FREQ", "f0");
  pars.setCase("ephstyle", "FREQ", "f1");
  pars.setCase("ephstyle", "FREQ", "f2");
  pars.setCase("ephstyle", "PER", "ephepoch");
  pars.setCase("ephstyle", "PER", "timeformat");
  pars.setCase("ephstyle", "PER", "timesys");
  pars.setCase("ephstyle", "PER", "ra");
  pars.setCase("ephstyle", "PER", "dec");
  pars.setCase("ephstyle", "PER", "p0");
  pars.setCase("ephstyle", "PER", "p1");
  pars.setCase("ephstyle", "PER", "p2");
  pars.setCase("timeorigin", "USER", "usertime");
  pars.setCase("timeorigin", "USER", "userformat");
  pars.setCase("timeorigin", "USER", "usersys");
}

PSearchApp::~PSearchApp() throw() {}

void PSearchApp::run() {
  // Clean up from any previous runs.
  resetApp();

  m_os.setMethod("run()");
  st_app::AppParGroup & pars(getParGroup());

  // Prompt for all parameters and save them.
  prompt(pars);

  // Get parameters.
  std::string out_file = pars["outfile"];
  double scan_step = pars["scanstep"];
  long num_trials = pars["numtrials"];
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
  pulsarDb::SloppyEphChooser chooser;
  initEphComputer(pars, chooser);

  // Use user input (parameters) together with computer to determine corrections to apply.
  bool vary_ra_dec = true;
  bool guess_pdot = true;
  initTimeCorrection(pars, vary_ra_dec, guess_pdot);

  // Report ephemeris status.
  std::set<pulsarDb::EphStatusCodeType> code_to_report;
  code_to_report.insert(pulsarDb::Remarked);
  reportEphStatus(m_os.warn(), code_to_report);

  // Compute central frequency of periodicity search, which is an expected pulse frequency at the time origin for the search.
  double origin = 0.;
  timeSystem::AbsoluteTime abs_origin = computeAbsoluteTime(origin);
  const pulsarDb::PulsarEph & eph = getEphComputer().choosePulsarEph(abs_origin);
  double f_center = eph.calcFrequency(abs_origin);
  // Note: The time system in which the frequency (f_center) is measured is ignored here, because it is only a rough estimate
  // of the frequency at the time origin (abs_origin) for the purpose of finding a reasonable scan range.

  // Compute frequency step from scan step and the Fourier resolution == 1. / duration,
  // using start/stop of the observation interval with all corrections applied.
  double duration = computeElapsedSecond(getStopTime()) - computeElapsedSecond(getStartTime());
  double f_step = scan_step / duration;

  // Choose which kind of test to create.
  std::string algorithm = pars["algorithm"];
  for (std::string::iterator itor = algorithm.begin(); itor != algorithm.end(); ++itor) *itor = std::toupper(*itor);

  // Create the proper period search object..
  std::auto_ptr<PeriodSearch> search(0);
  std::auto_ptr<PeriodicityTestArray> test_array(0);
  if (algorithm == "CHI2") {
    long num_phase = pars["numphase"];
    test_array.reset(new ChiSquaredTestArray(num_trials, num_phase));
    search.reset(new FoldingAnalysis(test_array.get(), f_center, f_step, origin, duration, "Hz"));
  } else if (algorithm == "RAYLEIGH") {
    test_array.reset(new RayleighTestArray(num_trials));
    search.reset(new FoldingAnalysis(test_array.get(), f_center, f_step, origin, duration, "Hz"));
  } else if (algorithm == "Z2N") {
    long num_harm = pars["numharm"];
    test_array.reset(new Z2nTestArray(num_trials, num_harm));
    search.reset(new FoldingAnalysis(test_array.get(), f_center, f_step, origin, duration, "Hz"));
  } else if (algorithm == "H") {
    long max_harm = pars["maxharm"];
    test_array.reset(new HTestArray(num_trials, max_harm));
    search.reset(new FoldingAnalysis(test_array.get(), f_center, f_step, origin, duration, "Hz"));
  } else {
    throw std::runtime_error("PSearchApp: invalid test algorithm " + algorithm);
  }

  for (setFirstEvent(); !isEndOfEventList(); setNextEvent()) {
    // Get event time as AbsoluteTime.
    timeSystem::AbsoluteTime abs_evt_time(getEventTime());

    // Convert event time to target time representation.
    double target_evt_time = computeElapsedSecond(abs_evt_time);

    // Fill into the period search object.
    search->fill(target_evt_time);
  }

  // Compute the statistics.
  search->computeStat();

  // Update the statistic viewer in PeriodSearch object, and get a reference to it.
  search->updateViewer();
  StatisticViewer & viewer(search->getViewer());

  // Set a plot title: use default title if user did not specify one.
  std::string title_uc(title);
  for (std::string::iterator itor = title_uc.begin(); itor != title_uc.end(); ++itor) *itor = std::toupper(*itor);
  if (title_uc != "DEFAULT") viewer.setTitle(title);

  // Interpret output file parameter.
  std::string out_file_uc = out_file;
  for (std::string::iterator itor = out_file_uc.begin(); itor != out_file_uc.end(); ++itor) *itor = std::toupper(*itor);

  if ("NONE" != out_file_uc) {
    // Find the template file.
    using namespace facilities;
    std::string template_file = commonUtilities::joinPath(commonUtilities::getDataPath("periodSearch"), "period-search-out.tpl");

    // Create output file.
    tip::IFileSvc::instance().createFile(out_file, template_file, clobber);

    // Open output file.
    std::auto_ptr<tip::Table> out_table(tip::IFileSvc::instance().editTable(out_file, "POWER_SPECTRUM"));

    // Write the summary to the output header, and the data to the output table.
    viewer.write(*out_table);

    // Get a target name to put in the output file.
    std::string psr_name = pars["psrname"];
    for (std::string::iterator itor = psr_name.begin(); itor != psr_name.end(); ++itor) *itor = std::toupper(*itor);

    // Update header keywords.
    tip::Header & header(out_table->getHeader());
    tip::Header::KeyValCont_t keywords;
    keywords.push_back(tip::Header::KeyValPair_t("DATE", header.formatTime(time(0))));
    keywords.push_back(tip::Header::KeyValPair_t("CREATOR", getName() + " " + getVersion()));
    keywords.push_back(tip::Header::KeyValPair_t("OBJECT", psr_name));
    keywords.push_back(tip::Header::KeyValPair_t("DATASUM", "-1")); // Force update of DATASUM keyword.
    header.update(keywords);
  }

  // Write the stats to the screen.
  viewer.write(m_os);

  // Display a plot, if desired.
  if (plot) {
    viewer.setLabel(0, "Pulse Frequency");
    viewer.setLabel(1, "Test Statistic");
    viewer.plot();
  }
}

void PSearchApp::prompt(st_app::AppParGroup & pars) {
  // Prompt for most parameters automatically.
  pars.Prompt("evfile");
  pars.Prompt("evtable");
  pars.Prompt("timefield");
  pars.Prompt("scfile");
  pars.Prompt("sctable");
  pars.Prompt("psrdbfile");
  pars.Prompt("psrname");
  pars.Prompt("outfile");

  pars.Prompt("algorithm");
  std::string algorithm = pars["algorithm"];
  for (std::string::iterator itor = algorithm.begin(); itor != algorithm.end(); ++itor) *itor = std::toupper(*itor);
  if (algorithm == "CHI2") {
    pars.Prompt("numphase");
  } else if (algorithm == "Z2N") {
    pars.Prompt("numharm");
  } else if (algorithm == "H") {
    pars.Prompt("maxharm");
  }

  pars.Prompt("scanstep");
  pars.Prompt("numtrials");

  pars.Prompt("timeorigin");
  std::string origin_style = pars["timeorigin"];
  for (std::string::iterator itor = origin_style.begin(); itor != origin_style.end(); ++itor) *itor = std::toupper(*itor);
  if (origin_style == "USER") {
    pars.Prompt("usertime");
    pars.Prompt("userformat");
    pars.Prompt("usersys");
  }

  pars.Prompt("ephstyle");
  std::string eph_style = pars["ephstyle"];
  for (std::string::iterator itor = eph_style.begin(); itor != eph_style.end(); ++itor) *itor = std::toupper(*itor);
  if (eph_style == "FREQ") {
    pars.Prompt("ephepoch");
    pars.Prompt("timeformat");
    pars.Prompt("timesys");
    pars.Prompt("ra");
    pars.Prompt("dec");
    pars.Prompt("f0");
    pars.Prompt("f1");
    pars.Prompt("f2");
  } else if (eph_style == "PER") {
    pars.Prompt("ephepoch");
    pars.Prompt("timeformat");
    pars.Prompt("timesys");
    pars.Prompt("ra");
    pars.Prompt("dec");
    pars.Prompt("p0");
    pars.Prompt("p1");
    pars.Prompt("p2");
  } else if (eph_style == "DB") {
    // No action needed.
  } else
    throw std::runtime_error("Unknown ephemeris style " + eph_style);

  pars.Prompt("tcorrect");
  pars.Prompt("solareph");
  pars.Prompt("matchsolareph");
  pars.Prompt("angtol");
  pars.Prompt("plot");
  pars.Prompt("title");
  pars.Prompt("leapsecfile");
  pars.Prompt("reportephstatus");
  pars.Prompt("chatter");
  pars.Prompt("clobber");
  pars.Prompt("debug");
  pars.Prompt("gui");
  pars.Prompt("mode");

  // Save current values of the parameters.
  pars.Save();
}

// List supported event file format(s).
timeSystem::EventTimeHandlerFactory<timeSystem::GlastScTimeHandler> glast_handler;

st_app::StAppFactory<PSearchApp> g_factory("gtpsearch");
