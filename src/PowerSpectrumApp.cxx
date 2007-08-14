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
#include "periodSearch/PeriodSearchViewer.h"
#include "FourierAnalysis.h"

using namespace periodSearch;
using namespace timeSystem;

static const std::string s_cvs_id = "$Name:  $";

PowerSpectrumApp::PowerSpectrumApp(): m_os("PowerSpectrumApp", "", 2) {
  st_app::AppParGroup & pars(getParGroup("gtpspec"));
  pars.setSwitch("ephstyle");
  pars.setCase("ephstyle", "FREQ", "cancelpdot");
  pars.setCase("ephstyle", "PER", "cancelpdot");
  pars.setCase("ephstyle", "FREQ", "f1f0ratio");
  pars.setCase("ephstyle", "FREQ", "f2f0ratio");
  pars.setCase("ephstyle", "PER", "p1p0ratio");
  pars.setCase("ephstyle", "PER", "p2p0ratio");
}

PowerSpectrumApp::~PowerSpectrumApp() throw() {
}

void PowerSpectrumApp::run() {
  // Clean up from any previous runs.
  resetApp();

  m_os.setMethod("run()");
  st_app::AppParGroup & pars(getParGroup("gtpspec"));

  // Prompt for all parameters and save them.
  prompt(pars);

  // Get parameters.
  std::string out_file = pars["outfile"];
  long num_bins = pars["numbins"];
  double low_f_cut = pars["lowfcut"];
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
  initEphComputer(pars, model, chooser, "DB");

  // Use user input (parameters) together with computer to determine corrections to apply.
  bool guess_pdot = false;
  initTimeCorrection(pars, guess_pdot);

  // Compute start time of the data set.
  double tstart = computeElapsedSecond(getStartTime());
  double tstop = computeElapsedSecond(getStopTime());

  // Get binwidth parameter.
  double bin_width = pars["binwidth"];

  // Create the proper test.
  std::auto_ptr<FourierAnalysis> test(new FourierAnalysis(tstart, tstop, bin_width, num_bins));

  for (setFirstEvent(); !isEndOfEventList(); setNextEvent()) {
    // Get event time as AbsoluteTime.
    AbsoluteTime abs_evt_time(getEventTime());

    // Convert event time to target time representation.
    double target_evt_time = computeElapsedSecond(abs_evt_time);

    // Fill into the test.
    test->fill(target_evt_time);
  }

  // Compute the statistics.
  test->computeStats();

  enum ChatLevel {
    eIncludeSummary= 2,
    eAllDetails = 3
  };

  // Use default title if user did not specify one.
  if (title == "DEFAULT") title = "Fourier Analysis: Power Spectrum";

  // Create a viewer for plotting and writing output.
  periodSearch::PeriodSearchViewer viewer(*test, low_f_cut);

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

void PowerSpectrumApp::prompt(st_app::AppParGroup & pars) {
  // Prompt for most parameters automatically.
  pars.Prompt("evfile");
  pars.Prompt("outfile");
  pars.Prompt("evtable");
  pars.Prompt("binwidth");
  pars.Prompt("numbins");
  pars.Prompt("lowfcut");
  pars.Prompt("ephstyle");

  // Only prompt for f1 & f2 / p1 & p2 if pdot correction is selected.
  if (true == bool(pars["cancelpdot"])) {
    std::string eph_style = pars["ephstyle"];
    if (eph_style == "FREQ") {
      pars.Prompt("f1f0ratio");
      pars.Prompt("f2f0ratio");
    } else if (eph_style == "PER") {
      pars.Prompt("p1p0ratio");
      pars.Prompt("p2p0ratio");
    } else {
      throw std::runtime_error("Unknown ephemeris style " + eph_style);
    }
  }

  pars.Prompt("psrdbfile");
  pars.Prompt("psrname");
  pars.Prompt("cancelpdot");
  pars.Prompt("demodbin");

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
