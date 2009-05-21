/** \file test_psearch.cxx
    \brief Period search tool.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include <cmath>
#include <memory>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "PeriodicityTestApp.h"
#include "PeriodSearchApp.h"
#include "PowerSpectrumApp.h"

#include "pulsarDb/PdotCanceler.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_stream/Stream.h"
#include "st_stream/StreamFormatter.h"
#include "st_stream/st_stream.h"

#include "timeSystem/AbsoluteTime.h"
#include "timeSystem/Duration.h"
#include "timeSystem/ElapsedTime.h"
#include "timeSystem/EventTimeHandler.h"
#include "timeSystem/GlastTimeHandler.h"
#include "timeSystem/PulsarTestApp.h"
#include "timeSystem/TimeInterval.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "ChiSquaredTestArray.h"
#include "FoldingAnalysis.h"
#include "FourierAnalysis.h"
#include "HTestArray.h"
#include "PeriodicityTestApp.h"
#include "PeriodicityTestArray.h"
#include "PeriodSearch.h"
#include "PeriodSearchApp.h"
#include "PowerSpectrumApp.h"
#include "RayleighTestArray.h"
#include "Z2nTestArray.h"

static const std::string s_cvs_id = "$Name:  $";

/** \class PeriodSearchAppTester
    \brief Test PeriodSearchApp application (gtpsearch).
*/
class PeriodSearchAppTester: public timeSystem::PulsarApplicationTester {
  public:
  /** \brief Construct a PeriodSearchAppTester object.
      \param test_app Unit test appliction of pulsar tool package, under which this application tester is to run.
  */
  PeriodSearchAppTester(timeSystem::PulsarTestApp & test_app);

  /// \brief Destruct this PeriodSearchAppTester object.
  virtual ~PeriodSearchAppTester() throw() {}

  /// \brief Returns an application object to be tested.
  virtual st_app::StApp * createApplication() const;

  /** \brief Return a logical true if the given header keyword is determined correct, and a logical false otherwise.
      \param keyword_name Name of the header keyword to be verified.
      \param out_keyword Header keyword taken from the output file to be verified.
      \param ref_keyword Header keyword taken from the reference file which out_keyword is checked against.
      \param error_stream Output stream for this method to put an error messages when verification fails.
  */
  virtual bool verify(const std::string & keyword_name, const tip::KeyRecord & out_keyword,
    const tip::KeyRecord & ref_keyword, std::ostream & error_stream) const;

  /** \brief Return a logical true if the given table cell is considered correct, and a logical false otherwise.
      \param column_name Name of the FITS column that the given table cell belongs to.
      \param out_cell Table cell taken from the output file to be verified.
      \param ref_cell Table cell taken from the reference file which out_cell is checked against.
      \param error_stream Output stream for this method to put an error message when verification fails.
  */
  virtual bool verify(const std::string & column_name, const tip::TableCell & out_cell, const tip::TableCell & ref_cell,
    std::ostream & error_stream) const;

  /** \brief Return a logical true if the given character string is considered correct, and a logical false otherwise.
      \param out_string Character string taken from the output file to be verified.
      \param ref_string Character string taken from the reference file which out_string is checked against.
      \param error_stream Output stream for this method to put an error message when verification fails.
  */
  virtual bool verify(const std::string & out_string, const std::string & ref_string, std::ostream & error_stream) const;
};

PeriodSearchAppTester::PeriodSearchAppTester(timeSystem::PulsarTestApp & test_app): PulsarApplicationTester("gtpsearch", test_app) {}

st_app::StApp * PeriodSearchAppTester::createApplication() const {
  return new PeriodSearchApp();
}

bool PeriodSearchAppTester::verify(const std::string & keyword_name, const tip::KeyRecord & out_keyword,
  const tip::KeyRecord & ref_keyword, std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  // Compare values.
  if ("COMMENT" == keyword_name) {
    // Extract keyword values as character strings.
    std::string out_value = out_keyword.getComment();
    std::string ref_value = ref_keyword.getComment();

    // Compare exactly like text outputs.
    verified = verify(out_value, ref_value, error_stream);

  } else {
    // Extract keyword values as character strings.
    std::string out_value = out_keyword.getValue();
    std::string ref_value = ref_keyword.getValue();

    // Require an exact match.
    verified = (out_value == ref_value);
    if (!verified) error_stream << "Value \"" << out_value << "\" not identical to reference \"" << ref_value << "\".";
  }

  // Return the result.
  return verified;
}

bool PeriodSearchAppTester::verify(const std::string & column_name, const tip::TableCell & out_cell, const tip::TableCell & ref_cell,
  std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  // Extract cell values as floating-point numbers.
  double out_value;
  double ref_value;
  out_cell.get(out_value);
  ref_cell.get(ref_value);
  error_stream.precision(std::numeric_limits<double>::digits10);

  // Compare values.
  if ("STATISTIC" == column_name) {
    // Require them be close enough as floating-point numbers of type double whose value is of the order of unity.
    double abs_tol = 1.e-6;
    verified = (std::fabs(out_value - ref_value) <= abs_tol);
    if (!verified) {
      error_stream << "Test statistic " << out_value << " not equivalent to reference " << ref_value <<
        " with absolute tolerance of " << abs_tol << ".";
    }

  } else {
    // Compare as floating-point numbers.
    // NOTE: This requires out_value be exactly zero when ref_value is zero.
    double rel_tol = std::numeric_limits<double>::epsilon() * 1000.;
    verified = (std::fabs(out_value - ref_value) <= std::fabs(ref_value) * rel_tol);
    if (!verified) {
      error_stream << "Value " << out_value << " not equivalent to reference " << ref_value <<
        " with relative tolerance of " << rel_tol;
    }
  }

  // Return the result.
  return verified;
}

bool PeriodSearchAppTester::verify(const std::string & out_string, const std::string & ref_string, std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  if (ref_string.find("MJD") != std::string::npos) {
    // Require a match down to the 10th decimal point, which is approx. 10 microseconds.
    verified = equivalent(out_string, ref_string, 1.e-10, 0.);
    if (!verified) {
      error_stream << "MJD number not equivalent to reference with absolute tolerance of 1e-10 days (8.64 microseconds)." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else if (ref_string.find("Fourier Resolution") != std::string::npos || ref_string.find("Sampling Frequency") != std::string::npos ||
    ref_string.find("Search Range") != std::string::npos || ref_string.find("Maximum Statistic") != std::string::npos ||
    ref_string.find("Chance Probability Range") != std::string::npos) {
    // Compare as floating-point numbers.
    verified = equivalent(out_string, ref_string);
    if (!verified) {
      error_stream << "Floating-point number not equivalent to reference." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else {
    // Require an exact match.
    verified = (out_string == ref_string);
    if (!verified) {
      error_stream << "Line not identical to reference." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }
  }

  // Return the result.
  return verified;
}

/** \class PowerSpectrumAppTester
    \brief Test PowerSpectrumApp application (gtpspec).
*/
class PowerSpectrumAppTester: public timeSystem::PulsarApplicationTester {
  public:
  /** \brief Construct a PowerSpectrumAppTester object.
      \param test_app Unit test appliction of pulsar tool package, under which this application tester is to run.
  */
  PowerSpectrumAppTester(timeSystem::PulsarTestApp & test_app);

  /// \brief Destruct this PowerSpectrumAppTester object.
  virtual ~PowerSpectrumAppTester() throw() {}

  /// \brief Returns an application object to be tested.
  virtual st_app::StApp * createApplication() const;

  /** \brief Return a logical true if the given header keyword is determined correct, and a logical false otherwise.
      \param keyword_name Name of the header keyword to be verified.
      \param out_keyword Header keyword taken from the output file to be verified.
      \param ref_keyword Header keyword taken from the reference file which out_keyword is checked against.
      \param error_stream Output stream for this method to put an error messages when verification fails.
  */
  virtual bool verify(const std::string & keyword_name, const tip::KeyRecord & out_keyword,
    const tip::KeyRecord & ref_keyword, std::ostream & error_stream) const;

  /** \brief Return a logical true if the given table cell is considered correct, and a logical false otherwise.
      \param column_name Name of the FITS column that the given table cell belongs to.
      \param out_cell Table cell taken from the output file to be verified.
      \param ref_cell Table cell taken from the reference file which out_cell is checked against.
      \param error_stream Output stream for this method to put an error message when verification fails.
  */
  virtual bool verify(const std::string & column_name, const tip::TableCell & out_cell, const tip::TableCell & ref_cell,
    std::ostream & error_stream) const;

  /** \brief Return a logical true if the given character string is considered correct, and a logical false otherwise.
      \param out_string Character string taken from the output file to be verified.
      \param ref_string Character string taken from the reference file which out_string is checked against.
      \param error_stream Output stream for this method to put an error message when verification fails.
  */
  virtual bool verify(const std::string & out_string, const std::string & ref_string, std::ostream & error_stream) const;
};

PowerSpectrumAppTester::PowerSpectrumAppTester(timeSystem::PulsarTestApp & test_app): PulsarApplicationTester("gtpspec", test_app) {}

st_app::StApp * PowerSpectrumAppTester::createApplication() const {
  return new PowerSpectrumApp();
}

bool PowerSpectrumAppTester::verify(const std::string & keyword_name, const tip::KeyRecord & out_keyword,
  const tip::KeyRecord & ref_keyword, std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  // Extract keyword values as character strings.
  std::string out_value = out_keyword.getValue();
  std::string ref_value = ref_keyword.getValue();

  // Compare values.
  if ("COMMENT" == keyword_name) {
    // Compare exactly like text outputs.
    verified = verify(out_value, ref_value, error_stream);

  } else {
    // Require an exact match.
    verified = (out_value == ref_value);
    if (!verified) error_stream << "Value \"" << out_value << "\" not identical to reference \"" << ref_value << "\".";
  }

  // Return the result.
  return verified;
}

bool PowerSpectrumAppTester::verify(const std::string & column_name, const tip::TableCell & out_cell, const tip::TableCell & ref_cell,
  std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  // Extract cell values as floating-point numbers.
  double out_value;
  double ref_value;
  out_cell.get(out_value);
  ref_cell.get(ref_value);
  error_stream.precision(std::numeric_limits<double>::digits10);

  // Compare values.
  if ("STATISTIC" == column_name) {
    // Require them be close enough as floating-point numbers of type double whose value is of the order of unity.
    double abs_tol = 1.e-6;
    verified = (std::fabs(out_value - ref_value) <= abs_tol);
    if (!verified) {
      error_stream << "Test statistic " << out_value << " not equivalent to reference " << ref_value <<
        " with absolute tolerance of " << abs_tol << ".";
    }

  } else {
    // Compare as floating-point numbers.
    // NOTE: This requires out_value be exactly zero when ref_value is zero.
    double rel_tol = std::numeric_limits<double>::epsilon() * 1000.;
    verified = (std::fabs(out_value - ref_value) <= std::fabs(ref_value) * rel_tol);
    if (!verified) {
      error_stream << "Value " << out_value << " not equivalent to reference " << ref_value <<
        " with relative tolerance of " << rel_tol;
    }
  }

  // Return the result.
  return verified;
}

bool PowerSpectrumAppTester::verify(const std::string & out_string, const std::string & ref_string, std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  if (ref_string.find("MJD") != std::string::npos) {
    // Require a match down to the 10th decimal point, which is approx. 10 microseconds.
    verified = equivalent(out_string, ref_string, 1.e-10, 0.);
    if (!verified) {
      error_stream << "MJD number not equivalent to reference with absolute tolerance of 1e-10 days (8.64 microseconds)." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else if (ref_string.find("Fourier Resolution") != std::string::npos || ref_string.find("Sampling Frequency") != std::string::npos ||
    ref_string.find("Search Range") != std::string::npos || ref_string.find("Maximum Statistic") != std::string::npos ||
    ref_string.find("Chance Probability Range") != std::string::npos) {
    // Compare as floating-point numbers.
    verified = equivalent(out_string, ref_string);
    if (!verified) {
      error_stream << "Floating-point number not equivalent to reference." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else {
    // Require an exact match.
    verified = (out_string == ref_string);
    if (!verified) {
      error_stream << "Line not identical to reference." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }
  }

  // Return the result.
  return verified;
}

/** \class PeriodicityTestAppTester
    \brief Test PeriodicityTestApp application (gtptest).
*/
class PeriodicityTestAppTester: public timeSystem::PulsarApplicationTester {
  public:
  /** \brief Construct a PeriodicityTestAppTester object.
      \param test_app Unit test appliction of pulsar tool package, under which this application tester is to run.
  */
  PeriodicityTestAppTester(timeSystem::PulsarTestApp & test_app);

  /// \brief Destruct this PeriodicityTestAppTester object.
  virtual ~PeriodicityTestAppTester() throw() {}

  /// \brief Returns an application object to be tested.
  virtual st_app::StApp * createApplication() const;

  /** \brief Return a logical true if the given header keyword is determined correct, and a logical false otherwise.
      \param keyword_name Name of the header keyword to be verified.
      \param out_keyword Header keyword taken from the output file to be verified.
      \param ref_keyword Header keyword taken from the reference file which out_keyword is checked against.
      \param error_stream Output stream for this method to put an error messages when verification fails.
  */
  virtual bool verify(const std::string & keyword_name, const tip::KeyRecord & out_keyword,
    const tip::KeyRecord & ref_keyword, std::ostream & error_stream) const;

  /** \brief Return a logical true if the given table cell is considered correct, and a logical false otherwise.
      \param column_name Name of the FITS column that the given table cell belongs to.
      \param out_cell Table cell taken from the output file to be verified.
      \param ref_cell Table cell taken from the reference file which out_cell is checked against.
      \param error_stream Output stream for this method to put an error message when verification fails.
  */
  virtual bool verify(const std::string & column_name, const tip::TableCell & out_cell, const tip::TableCell & ref_cell,
    std::ostream & error_stream) const;

  /** \brief Return a logical true if the given character string is considered correct, and a logical false otherwise.
      \param out_string Character string taken from the output file to be verified.
      \param ref_string Character string taken from the reference file which out_string is checked against.
      \param error_stream Output stream for this method to put an error message when verification fails.
  */
  virtual bool verify(const std::string & out_string, const std::string & ref_string, std::ostream & error_stream) const;
};

PeriodicityTestAppTester::PeriodicityTestAppTester(timeSystem::PulsarTestApp & test_app): PulsarApplicationTester("gtptest", test_app) {}

st_app::StApp * PeriodicityTestAppTester::createApplication() const {
  return new PeriodicityTestApp();
}

bool PeriodicityTestAppTester::verify(const std::string & keyword_name, const tip::KeyRecord & out_keyword,
  const tip::KeyRecord & ref_keyword, std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  // Extract keyword values as character strings.
  std::string out_value = out_keyword.getValue();
  std::string ref_value = ref_keyword.getValue();

  // Compare values.
  if ("COMMENT" == keyword_name) {
    // Compare exactly like text outputs.
    verified = verify(out_value, ref_value, error_stream);

  } else {
    // Require an exact match.
    verified = (out_value == ref_value);
    if (!verified) error_stream << "Value \"" << out_value << "\" not identical to reference \"" << ref_value << "\".";
  }

  // Return the result.
  return verified;
}

bool PeriodicityTestAppTester::verify(const std::string & column_name, const tip::TableCell & out_cell, const tip::TableCell & ref_cell,
  std::ostream & error_stream) const {
  // Initialize return value.
  bool verified = false;

  // Extract cell values as floating-point numbers.
  double out_value;
  double ref_value;
  out_cell.get(out_value);
  ref_cell.get(ref_value);
  error_stream.precision(std::numeric_limits<double>::digits10);

  // Compare values.
  if ("POWER" == column_name || "CANDIDATE_VALUE" == column_name) {
    // Require them be close enough as floating-point numbers of type double whose value is of the order of unity.
    double abs_tol = 1.e-6;
    verified = (std::fabs(out_value - ref_value) <= abs_tol);
    if (!verified) {
      error_stream << "Test statistic " << out_value << " not equivalent to reference " << ref_value <<
        " with absolute tolerance of " << abs_tol << ".";
    }

  } else {
    // Compare as floating-point numbers.
    // NOTE: This requires out_value be exactly zero when ref_value is zero.
    double rel_tol = std::numeric_limits<double>::epsilon() * 1000.;
    verified = (std::fabs(out_value - ref_value) <= std::fabs(ref_value) * rel_tol);
    if (!verified) {
      error_stream << "Value " << out_value << " not equivalent to reference " << ref_value <<
        " with relative tolerance of " << rel_tol;
    }
  }

  // Return the result.
  return verified;
}

bool PeriodicityTestAppTester::verify(const std::string & out_string, const std::string & ref_string, std::ostream & error_stream)
  const {
  // Initialize return value.
  bool verified = false;

  if (ref_string.find("Test Statistic") != std::string::npos ||
    ref_string.find("Chance Probability Range") != std::string::npos) {
    // Compare as floating-point numbers.
    verified = equivalent(out_string, ref_string);
    if (!verified) {
      error_stream << "Floating-point number not equivalent to reference." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }

  } else {
    // Require an exact match.
    verified = (out_string == ref_string);
    if (!verified) {
      error_stream << "Line not identical to reference." <<
        std::endl << "[OUT] " << out_string << std::endl << "[REF] " << ref_string;
    }
  }

  // Return the result.
  return verified;
}

/** \class PeriodSearchTestApp
    \brief Test periodSearch package and applications in it.
*/
class PeriodSearchTestApp : public timeSystem::PulsarTestApp {
  public:
    /// \brief Construct a PeriodSearchTestApp object.
    PeriodSearchTestApp();

    /// \brief Destruct this PeriodSearchTestApp object.
    virtual ~PeriodSearchTestApp() throw() {}

    /// \brief Do all tests.
    virtual void runTest();

    /// \brief Test PeriodSearch class.
    void testPeriodSearch();

    /// \brief Test chance probability computations.
    void testChanceProb();

    /// \brief Test ChiSquaredTestArray class.
    void testChiSquaredTestArray();

    /// \brief Test Z2nTestArray class.
    void testZ2nTestArray();

    /// \brief Test HTestArray class.
    void testHTestArray();

    /// \brief Test RayleighTestArray class.
    void testRayleighTestArray();

    /// \brief Test PeriodSearchApp class.
    void testPeriodSearchApp();

    /// \brief Test PowerSpectrumApp class.
    void testPowerSpectrumApp();

    /// \brief Test PeriodicityTestApp class.
    void testPeriodicityTestApp();

  private:
    /** \brief Helper method for testPeriodSearch, to test all types of test statistics.
        \param prefix Character string to be used as a prefix of output file names.
        \param events List of event times to be tested periodicity.
        \param t_start Start time of the time series data.
        \param t_stop Stop time of the time series data.
        \param center Central frequency of a frequency scan in Hz.
        \param step Frequency step size in Hz.
        \param num_trials Number of trial frequencies.
        \param epoch Reference epoch of the time series data.
        \param num_bins The nmber of bins for the chi-squared test, the number of harmonics for the Z2n test,
                        and the maximumn number of harmonics for the H test.
        \param fourier_width Width of time bins for the Fourier transformation.
        \param fourier_num_bins Number of time bins for the Fourier transformation.
        \param fourier_min_freq Minimum frequency to be computed by the Fourier transformation.
        \param fourier_max_freq Maximum frequency to be computed by the Fourier transformation.
        \param plot Set to true if the result should be plotted on screen. Set to false otherwise.
    */
    void testAllStats(const std::string & prefix, const std::vector<double> & events, double t_start, double t_stop,
      double center, double step, long num_trials, double epoch, int num_bins,
      double fourier_width, int fourier_num_bins, double fourier_min_freq, double fourier_max_freq, bool plot);

    /** \brief Helper method for testAllStats, to test one type of test statistics.
        \param events List of event times to be tested periodicity.
        \param search PeriodSearch object to be tested.
        \param plot_title Title of a plot that displays a test result on screen.
        \param out_file Name of an output file that stores a test result.
        \param min_freq Minimum frequency to compute test statistics. Set to a negative value for not specifying the minimum.
        \param max_freq Maximum frequency to compute test statistics. Set to a negative value for not specifying the maximum.
    */
    void testOneSearch(const std::vector<double> & events, PeriodSearch & search, const std::string & plot_title,
      const std::string & out_file, bool plot, double min_freq = -1., double max_freq = -1.);
};

PeriodSearchTestApp::PeriodSearchTestApp(): PulsarTestApp("periodSearch") {
  setName("test_stpsearch");
  setVersion(s_cvs_id);
}

void PeriodSearchTestApp::runTest() {
  // Test PeriodSearch subclasses.
  testPeriodSearch();

  // Test computations of chance probability.
  testChanceProb();

  // Test PeriodicityTestArary subclasses.
  testChiSquaredTestArray();
  testZ2nTestArray();
  testHTestArray();
  testRayleighTestArray();

  // Test applications.
  testPeriodSearchApp();
  testPowerSpectrumApp();
  testPeriodicityTestApp();
}

void PeriodSearchTestApp::testPeriodSearch() {
  setMethod("testPeriodSearch");

  // Trick up some fake events.
  int num_evts = 1000;

  std::vector<double> fake_evts(num_evts);
  setPrecision(16);

  for (int ii = 0; ii < num_evts; ++ii)
    fake_evts[ii] = ii + 1;

  // Set up search parameters.
  double central = 1.;
  double step = .5e-3;
  int num_pds = 61;
  double epoch = .2;
  int num_bins = 10;
  double duration = 1000.;

  bool plot = getParGroup()["plot"];

  // First do simple test with this highly artificial data.
  // Note for Fourier test: width of .1 s -> Nyquist = 1/.2s = 5 Hz.
  testAllStats("artificial", fake_evts, 0., duration, central, step, num_pds, epoch, num_bins, .1, 10000, .9001, 1.1001, plot);

  // Data taken from M. Hirayama's work with modified ASCA data.
  // http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/existing.html#tryout003
  num_pds = 101;

  central = 1. / 50.41843041e-3;
  step = .168e-7 * central * central;

  epoch = 212380785.922;

  fake_evts.clear();

  // Read some real data.
  std::string event_file = prependDataPath("step-01.fits");
  std::auto_ptr<const tip::Table> evt_table(tip::IFileSvc::instance().readTable(event_file, "EVENTS"));

  // Get gti_table.
  std::auto_ptr<const tip::Table> gti_table(tip::IFileSvc::instance().readTable(event_file, "GTI"));

  // Make the array big enough to hold these events.
  fake_evts.resize(evt_table->getNumRecords());

  // Correct event times for changed MJDREF.
  const double event_time_offset = (54101 - 51910) * 86400.;
  std::vector<double>::iterator event_itor = fake_evts.begin();
  for (tip::Table::ConstIterator itor = evt_table->begin(); itor != evt_table->end(); ++itor, ++event_itor) {
    *event_itor = (*itor)["TIME"].get() + event_time_offset;
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
  tstart += event_time_offset;
  tstop += event_time_offset;

  // Repeat simple test with this somewhat less artificial data.
  // Note for Fourier test: width of .01 s -> Nyquist = 1/.02s = 50 Hz.
  testAllStats("psrb0540", fake_evts, tstart, tstop, central, step, num_pds, epoch, num_bins, .01, 1000000, 19.82001, 19.85001, plot);

  // Now test pdot correction.
  timeSystem::AbsoluteTime glast_origin("TDB", 51910, 0.);
  timeSystem::AbsoluteTime abs_epoch = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, epoch));
  double pdot = 4.7967744e-13;
  std::vector<double> fdot_ratio(1, -pdot * central);
  pulsarDb::PdotCanceler canceler("TDB", abs_epoch, fdot_ratio);

  // Correct the data.
  timeSystem::AbsoluteTime evt_time("TDB", 0, 0.);
  for (std::vector<double>::iterator itor = fake_evts.begin(); itor != fake_evts.end(); ++itor) {
    evt_time = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, *itor));
    canceler.cancelPdot(evt_time);
    *itor = (evt_time - glast_origin).computeDuration("TDB", "Sec");
  }

  // Cancel pdot in tstart, tstop and epoch to be consistent.
  evt_time = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, tstart));
  canceler.cancelPdot(evt_time);
  tstart = (evt_time - glast_origin).computeDuration("TDB", "Sec");

  evt_time = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, tstop));
  canceler.cancelPdot(evt_time);
  tstop = (evt_time - glast_origin).computeDuration("TDB", "Sec");

  evt_time = glast_origin + timeSystem::ElapsedTime("TDB", timeSystem::Duration(0, epoch));
  canceler.cancelPdot(evt_time);
  epoch = (evt_time - glast_origin).computeDuration("TDB", "Sec");

  // Repeat test with the pdot corrected data.
  testAllStats("psrb0540-pdot", fake_evts, tstart, tstop, central, step, num_pds, epoch, num_bins, .01, 1000000,
    19.82001, 19.85001, plot);
}

void PeriodSearchTestApp::testChanceProb() {
  setMethod("testChanceProb");

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
      double chance_prob = PeriodSearch::computeChanceProbMultiTrial(prob_one_trial[jj], num_indep_trial[ii]);
      if (0. > chance_prob) {
        err() << "ERROR: computeChanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] <<
          ") unexpectedly returned " << chance_prob << ", which is < 0." << std::endl;
      } else if (1. < chance_prob) {
        err() << "ERROR: computeChanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] <<
          ") unexpectedly returned " << chance_prob << ", which is > 1." << std::endl;
      } else if ((0. == approx_chance_prob[ii][jj] && 0. != chance_prob) ||
        (0. != approx_chance_prob[ii][jj] && epsilon < std::fabs(chance_prob / approx_chance_prob[ii][jj] - 1.))) {
        err() << "ERROR: computeChanceProbMultiTrial(" << prob_one_trial[jj] << ", " << num_indep_trial[ii] << ") returned " <<
          chance_prob << ", not " << approx_chance_prob[ii][jj] << ", as expected." << std::endl;
      }
    }
  }
}

void PeriodSearchTestApp::testChiSquaredTestArray() {
  setMethod("testChiSquaredTestArray");

  // Prepare a ChiSquaredTestArray object.
  const int num_trial = 3;
  const int num_bin = 4;
  ChiSquaredTestArray chi2_test(num_trial, num_bin);
  const StatisticViewer & chi2_viewer(chi2_test.getViewer());

  // Fill events.
  const int num_events = 10;
  for (int ii = 0; ii < num_events; ++ii) {
    double phase = (ii + 0.1) / num_events;
    chi2_test.fill(0, phase);
    chi2_test.fill(1, phase); chi2_test.fill(1, phase);
    chi2_test.fill(2, phase); chi2_test.fill(2, phase); chi2_test.fill(2, phase);
  }

  // Check size.
  int test_size = chi2_test.size();
  if (test_size != num_trial) {
    err() << "Size of the ChiSquaredTestArray was reported as " << test_size << ", not " << num_trial << "." << std::endl;
  }

  // Set the comparison precision.
  const double epsilon = 1.e-12;

  // Check the folded light curve.
  const StatisticViewer::data_type & result_counts = chi2_viewer.getData(1);
  const double expected_counts[] = {3., 2., 3., 2.};
  for (int trial = 0; trial < num_trial; ++trial) {
    chi2_test.updateViewer(trial);
    for (int ii = 0; ii < num_bin; ++ii) {
      double result = result_counts[ii];
      double expected = expected_counts[ii] * (trial + 1);
      if (std::fabs(result/expected - 1.) > epsilon) {
        err() << "Photon count in phase bin " << ii << " of trial " << trial << " was " << result
          << ", not " << expected << "." << std::endl;
      }
    }
  }

  // Check S-value.
  const double expected_values[] = {0.4, 0.8, 1.2};
  for (int trial = 0; trial < num_trial; ++trial) {
    double result = chi2_test.computeStat(trial);
    double expected = expected_values[trial];
    if (std::fabs(result/expected - 1.) > epsilon) {
      err() << "S-value for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }
}

void PeriodSearchTestApp::testZ2nTestArray() {
  setMethod("testZ2nTestArray");

  // Prepare a Z2nTestArray object.
  const int num_trial = 3;
  const int num_harm = 4;
  Z2nTestArray z2n_test(num_trial, num_harm);
  const StatisticViewer & z2n_viewer(z2n_test.getViewer());

  // Fill events.
  const int num_events = 12;
  for (int ii = 0; ii < num_events; ++ii) {
    z2n_test.fill(0, 0.123);
    double phase = (ii % 3 + 1) / 8.; // 45, 90, and 135 degrees.
    z2n_test.fill(1, phase);
    phase = ii / 12.; // one whole round.
    z2n_test.fill(2, phase);
  }

  // Check size.
  int test_size = z2n_test.size();
  if (test_size != num_trial) {
    err() << "Size of the Z2nTestArray was reported as " << test_size << ", not " << num_trial << "." << std::endl;
  }

  // Set the comparison precision.
  const double epsilon = 1.e-12;

  // Check the power spectrum.
  const StatisticViewer::data_type & result_powers = z2n_viewer.getData(1);
  const double norm = 2. / num_events;
  double power10 = 4. * (std::sqrt(2.) + 1.);
  power10 = norm * power10 * power10;
  double power11 = norm * 4. * 4.;
  double power12 = 4. * (std::sqrt(2.) - 1.);
  power12 = norm * power12 * power12;
  double power13 = norm * 4. * 4.;
  const double expected_powers[][num_harm] = { {24., 24., 24., 24.}, {power10, power11, power12, power13}, {0., 0., 0., 0.} };
  for (int trial = 0; trial < num_trial; ++trial) {
    z2n_test.updateViewer(trial);
    for (int ii = 0; ii < num_harm; ++ii) {
      double result = result_powers[ii];
      double expected = expected_powers[trial][ii];
      if ((expected == 0. && std::fabs(result) > epsilon)
          || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
        err() << "Power for harmonic nubmer " << ii << " of trial " << trial << " was " << result
          << ", not " << expected << "." << std::endl;
      }
    }
  }

  // Check Z2n-value.
  double expected_values[] = {96., power10 + power11 + power12 + power13, 0.};
  for (int trial = 0; trial < num_trial; ++trial) {
    double result = z2n_test.computeStat(trial);
    double expected = expected_values[trial];
    if ((expected == 0. && std::fabs(result) > epsilon)
        || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
      err() << "Z2n-value for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }
}

void PeriodSearchTestApp::testHTestArray() {
  setMethod("testHTestArray");

  // Prepare a HTestArray object.
  const int num_trial = 3;
  const int num_harm = 4;
  HTestArray h_test(num_trial, num_harm);
  const StatisticViewer & h_viewer(h_test.getViewer());

  // Fill events.
  const int num_events = 12;
  for (int ii = 0; ii < num_events; ++ii) {
    h_test.fill(0, 0.123);
    double phase = (ii % 3 + 1) / 8.; // 45, 90, and 135 degrees.
    h_test.fill(1, phase);
    phase = ii / 12.; // one whole round.
    h_test.fill(2, phase);
  }

  // Check size.
  int test_size = h_test.size();
  if (test_size != num_trial) {
    err() << "Size of the HTestArray was reported as " << test_size << ", not " << num_trial << "." << std::endl;
  }

  // Set the comparison precision.
  const double epsilon = 1.e-12;

  // Check the power spectrum.
  const StatisticViewer::data_type & result_powers = h_viewer.getData(1);
  const double norm = 2. / num_events;
  double power10 = 4. * (std::sqrt(2.) + 1.);
  power10 = norm * power10 * power10;
  double power11 = norm * 4. * 4.;
  double power12 = 4. * (std::sqrt(2.) - 1.);
  power12 = norm * power12 * power12;
  double power13 = norm * 4. * 4.;
  const double expected_powers[][num_harm] = {
    {24., 44., 64., 84.},
    {power10, power10 + power11 - 4., power10 + power11 + power12 - 8., power10 + power11 + power12 + power13 - 12.},
    {0., -4., -8., -12.}
  };
  for (int trial = 0; trial < num_trial; ++trial) {
    h_test.updateViewer(trial);
    for (int ii = 0; ii < num_harm; ++ii) {
      double result = result_powers[ii];
      double expected = expected_powers[trial][ii];
      if ((expected == 0. && std::fabs(result) > epsilon)
          || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
        err() << "Candidate H value for harmonic nubmer " << ii << " of trial " << trial << " was " << result
          << ", not " << expected << "." << std::endl;
      }
    }
  }

  // Check H-value.
  double expected_values[] = {84., power10, 0.};
  for (int trial = 0; trial < num_trial; ++trial) {
    double result = h_test.computeStat(trial);
    double expected = expected_values[trial];
    if ((expected == 0. && std::fabs(result) > epsilon)
        || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
      err() << "H-value for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }
}

void PeriodSearchTestApp::testRayleighTestArray() {
  setMethod("testRayleighTestArray");

  // Prepare a RayleighTestArray object.
  const int num_trial = 3;
  RayleighTestArray rayleigh_test(num_trial);
  const StatisticViewer & rayleigh_viewer(rayleigh_test.getViewer());

  // Fill events.
  const int num_events = 12;
  for (int ii = 0; ii < num_events; ++ii) {
    rayleigh_test.fill(0, 0.123);
    double phase = (ii % 3 + 1) / 8.; // 45, 90, and 135 degrees.
    rayleigh_test.fill(1, phase);
    phase = ii / 12.; // one whole round.
    rayleigh_test.fill(2, phase);
  }

  // Check size.
  int test_size = rayleigh_test.size();
  if (test_size != num_trial) {
    err() << "Size of the RayleighTestArray was reported as " << test_size << ", not " << num_trial << "." << std::endl;
  }

  // Set the comparison precision.
  const double epsilon = 1.e-12;

  // Check the power spectrum.
  const StatisticViewer::data_type & result_powers = rayleigh_viewer.getData(1);
  const double norm = 2. / num_events;
  double power10 = 4. * (std::sqrt(2.) + 1.);
  power10 = norm * power10 * power10;
  const double expected_powers[] = {24., power10, 0.};
  for (int trial = 0; trial < num_trial; ++trial) {
    rayleigh_test.updateViewer(trial);
    double result = result_powers[0];
    double expected = expected_powers[trial];
    if ((expected == 0. && std::fabs(result) > epsilon)
        || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
      err() << "Viewer content for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }

  // Check Rayleigh-value.
  double expected_values[] = {24., power10, 0.};
  for (int trial = 0; trial < num_trial; ++trial) {
    double result = rayleigh_test.computeStat(trial);
    double expected = expected_values[trial];
    if ((expected == 0. && std::fabs(result) > epsilon)
        || (expected != 0. && std::fabs(result/expected - 1.) > epsilon)) {
      err() << "Rayleigh statistic for trial " << trial << " was " << result << ", not " << expected << "." << std::endl;
    }
  }
}

void PeriodSearchTestApp::testPeriodSearchApp() {
  setMethod("testPeriodSearchApp");

  // Create an application tester object.
  PeriodSearchAppTester app_tester(*this);

  // List supported event file format(s).
  timeSystem::EventTimeHandlerFactory<timeSystem::GlastScTimeHandler> glast_sctime_handler;

  // Prepare variables to create application objects.
  std::list<std::string> test_name_cont;
  test_name_cont.push_back("par1");
  test_name_cont.push_back("par2");
  test_name_cont.push_back("par3");
  test_name_cont.push_back("par4");

  // Prepare files to be used in the tests.
  std::string ev_file = prependDataPath("my_pulsar_events_v3.fits");
  std::string sc_file = prependDataPath("my_pulsar_spacecraft_data_v3r1.fits");
  std::string master_pulsardb = prependDataPath("master_pulsardb_v2.fits");
  std::string ev_file_2gti = prependDataPath("my_pulsar_events_2gti.fits");
  std::string sc_file_bogus = prependDataPath("my_pulsar_spacecraft_data_bogus.fits");

  // Loop over parameter sets.
  for (std::list<std::string>::const_iterator test_itor = test_name_cont.begin(); test_itor != test_name_cont.end(); ++test_itor) {
    const std::string & test_name = *test_itor;
    std::string log_file(getMethod() + "_" + test_name + ".log");
    std::string log_file_ref(getMethod() + "_" + test_name + ".ref");
    std::string out_file(getMethod() + "_" + test_name + ".fits");
    std::string out_file_ref(prependOutrefPath(out_file));
    bool ignore_exception(false);

    // Set default parameters.
    st_app::AppParGroup pars(app_tester.getName());
    pars["evfile"] = "";
    pars["scfile"] = "";
    pars["psrdbfile"] = "";
    pars["psrname"] = "ANY";
    pars["outfile"] = "";
    pars["algorithm"] = "CHI2";
    pars["numphase"] = 10;
    pars["numharm"] = 10;
    pars["maxharm"] = 10;
    pars["scanstep"] = 0.5;
    pars["numtrials"] = 100;
    pars["timeorigin"] = "MIDDLE";
    pars["usertime"] = 0.;
    pars["userformat"] = "FILE";
    pars["usersys"] = "FILE";
    pars["ephstyle"] = "DB";
    pars["ephepoch"] = 0.;
    pars["timeformat"] = "FILE";
    pars["timesys"] = "FILE";
    pars["ra"] = 0.;
    pars["dec"] = 0.;
    pars["f0"] = 1.;
    pars["f1"] = 0.;
    pars["f2"] = 0.;
    pars["p0"] = 1.;
    pars["p1"] = 0.;
    pars["p2"] = 0.;
    pars["tcorrect"] = "AUTO";
    pars["solareph"] = "JPL DE405";
    pars["matchsolareph"] = "ALL";
    pars["angtol"] = 1.e-8;
    pars["evtable"] = "EVENTS";
    pars["timefield"] = "TIME";
    pars["sctable"] = "SC_DATA";
    pars["plot"] = "yes";
    pars["title"] = "DEFAULT";
    pars["leapsecfile"] = "DEFAULT";
    pars["reportephstatus"] = "yes";
    pars["chatter"] = 2;
    pars["clobber"] = "yes";
    pars["debug"] = "no";
    pars["gui"] = "no";
    pars["mode"] = "ql";

    // Set test-specific parameters.
    if ("par1" == test_name) {
      // Test standard computation with DB option.
      pars["algorithm"] = "Chi2";
      pars["evfile"] = ev_file;
      pars["evtable"] = "EVENTS";
      pars["scfile"] = sc_file;
      pars["outfile"] = out_file;
      pars["ephstyle"] = "DB";
      pars["psrdbfile"] = master_pulsardb;
      pars["psrname"] = "PSR B0540-69";
      pars["scanstep"] = .5;
      pars["tcorrect"] = "AutO";
      pars["matchsolareph"] = "nONe";
      pars["numtrials"] = 200;
      pars["timeorigin"] = "uSer";
      pars["usertime"] = 214380785.922;
      pars["usersys"] = "tDb";
      pars["userformat"] = "fIle";
      pars["numphase"] = 10;
      pars["timefield"] = "TIME";
      pars["plot"] = "nO";
      pars["title"] = "My statistical test";
      pars["gui"] = "No";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par2" == test_name) {
      // Test ephemeris status reporting.
      std::string summary_file("psrdb_summary.txt");
      remove(summary_file.c_str());
      std::ofstream ofs_summary(summary_file.c_str());
      ofs_summary << prependDataPath("psrdb_spin.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_remark.txt") << std::endl;
      ofs_summary.close();
      pars["algorithm"] = "Chi2";
      pars["evfile"] = ev_file_2gti;
      pars["evtable"] = "EVENTS";
      pars["scfile"] = sc_file_bogus;
      pars["outfile"] = out_file;
      pars["ephstyle"] = "DB";
      pars["psrdbfile"] = "@" + summary_file;
      pars["psrname"] = "PSR J0540-6919";
      pars["scanstep"] = .5;
      pars["tcorrect"] = "AutO";
      pars["matchsolareph"] = "nONe";
      pars["numtrials"] = 200;
      pars["timeorigin"] = "uSer";
      pars["usertime"] = 214380785.922;
      pars["usersys"] = "tDb";
      pars["userformat"] = "fIle";
      pars["numphase"] = 10;
      pars["timefield"] = "TIME";
      pars["plot"] = "nO";
      pars["title"] = "My statistical test";
      pars["gui"] = "No";

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      ofs << "gtpsearch: WARNING: The following pulsar ephemeris status are reported." << std::endl;
      ofs << "[1] Remarked \"Test remark entry No.2\" since 53990 MJD (TDB) until 54010 MJD (TDB)" << std::endl;
      ofs << "[2] Remarked \"Test remark entry No.3\" since 54025 MJD (TDB) until 54035 MJD (TDB)" << std::endl;
      ofs << "[3] Remarked \"Test remark entry No.4\" since 54050 MJD (TDB) until 54070 MJD (TDB)" << std::endl;
      ofs << "[4] Remarked \"Test remark entry No.6\" since 53990 MJD (TDB) until 54030 MJD (TDB)" << std::endl;
      ofs << "[5] Remarked \"Test remark entry No.7\" since 54030 MJD (TDB) until 54070 MJD (TDB)" << std::endl;
      ofs << "[6] Remarked \"Test remark entry No.8\" since 53990 MJD (TDB) until 54070 MJD (TDB)" << std::endl;
      std::logic_error error("Error while computing an S-value of chi-squared test: No events filled for test #0");
      app_tester.writeException(ofs, error);
      ofs << std::endl;
      ofs.close();

      out_file.erase();
      out_file_ref.erase();
      ignore_exception = true;

    } else if ("par3" == test_name) {
      // Test no reporting of ephemeris status with reportephstatus=no.
      std::string summary_file("psrdb_summary.txt");
      remove(summary_file.c_str());
      std::ofstream ofs_summary(summary_file.c_str());
      ofs_summary << prependDataPath("psrdb_spin.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_remark.txt") << std::endl;
      ofs_summary.close();
      pars["algorithm"] = "Chi2";
      pars["evfile"] = ev_file_2gti;
      pars["evtable"] = "EVENTS";
      pars["scfile"] = sc_file_bogus;
      pars["outfile"] = out_file;
      pars["ephstyle"] = "DB";
      pars["psrdbfile"] = "@" + summary_file;
      pars["psrname"] = "PSR J0540-6919";
      pars["scanstep"] = .5;
      pars["tcorrect"] = "AutO";
      pars["matchsolareph"] = "nONe";
      pars["numtrials"] = 200;
      pars["timeorigin"] = "uSer";
      pars["usertime"] = 214380785.922;
      pars["usersys"] = "tDb";
      pars["userformat"] = "fIle";
      pars["numphase"] = 10;
      pars["timefield"] = "TIME";
      pars["plot"] = "nO";
      pars["title"] = "My statistical test";
      pars["gui"] = "No";
      pars["reportephstatus"] = "no";

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      std::logic_error error("Error while computing an S-value of chi-squared test: No events filled for test #0");
      app_tester.writeException(ofs, error);
      ofs << std::endl;
      ofs.close();

      out_file.erase();
      out_file_ref.erase();
      ignore_exception = true;

    } else if ("par4" == test_name) {
      // Test reporting of database creation history.
      std::string summary_file("psrdb_summary.txt");
      remove(summary_file.c_str());
      std::ofstream ofs_summary(summary_file.c_str());
      ofs_summary << prependDataPath("psrdb_spin.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_remark.txt") << std::endl;
      ofs_summary.close();
      pars["algorithm"] = "Chi2";
      pars["evfile"] = ev_file_2gti;
      pars["evtable"] = "EVENTS";
      pars["scfile"] = sc_file_bogus;
      pars["outfile"] = out_file;
      pars["ephstyle"] = "DB";
      pars["psrdbfile"] = "@" + summary_file;
      pars["psrname"] = "PSR J0540-6919";
      pars["scanstep"] = .5;
      pars["tcorrect"] = "AutO";
      pars["matchsolareph"] = "nONe";
      pars["numtrials"] = 200;
      pars["timeorigin"] = "uSer";
      pars["usertime"] = 214380785.922;
      pars["usersys"] = "tDb";
      pars["userformat"] = "fIle";
      pars["numphase"] = 10;
      pars["timefield"] = "TIME";
      pars["plot"] = "nO";
      pars["title"] = "My statistical test";
      pars["gui"] = "No";
      pars["reportephstatus"] = "no";
      pars["chatter"] = 4;

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      ofs << "gtpsearch: INFO: ==========================" << std::endl;
      ofs << "gtpsearch: INFO: Pulsar ephemerides are loaded and/or filtered as follows:" << std::endl;
      ofs << "gtpsearch: INFO:    Load TEXTDB SPIN_PARAMETERS(FREQ) FILENAME='psrdb_spin.txt'" << std::endl;
      ofs << "gtpsearch: INFO:    Load TEXTDB REMARKS FILENAME='psrdb_remark.txt'" << std::endl;
      ofs << "gtpsearch: INFO:    Filter by pulsar name 'PSR J0540-6919'" << std::endl;
      ofs << "gtpsearch: INFO: ==========================" << std::endl;
      ofs << "gtpsearch: INFO: --------------------------" << std::endl;
      ofs << "gtpsearch: INFO: Arrival time corrections are applied as follows:" << std::endl;
      ofs << "gtpsearch: INFO:    Barycentric correction: Applied if necessary" << std::endl;
      ofs << "gtpsearch: INFO:    Binary demodulation: Not applied" << std::endl;
      ofs << "gtpsearch: INFO:    Pdot cancellation: Applied" << std::endl;
      ofs << "gtpsearch: INFO: Following time system(s) are listed for this task:" << std::endl;
      ofs << "gtpsearch: INFO:    Spin ephemeri(de)s are defined in: TDB(5)" << std::endl;
      ofs << "gtpsearch: INFO:    Orbital ephemeri(de)s are defined in: None" << std::endl;
      ofs << "gtpsearch: INFO:    Pdot cancellation will be performed in: TDB" << std::endl;
      ofs << "gtpsearch: INFO:    Time series analysis will be performed in: TDB" << std::endl;
      ofs << "gtpsearch: INFO: --------------------------" << std::endl;
      std::logic_error error("Error while computing an S-value of chi-squared test: No events filled for test #0");
      app_tester.writeException(ofs, error);
      ofs << std::endl;
      ofs.close();

      out_file.erase();
      out_file_ref.erase();
      ignore_exception = true;

    } else {
      // Skip this iteration.
      continue;
    }

    // Test the application.
    app_tester.test(pars, log_file, log_file_ref, out_file, out_file_ref, ignore_exception);
  }
}

void PeriodSearchTestApp::testPowerSpectrumApp() {
  setMethod("testPowerSpectrumApp");

  // Create an application tester object.
  PowerSpectrumAppTester app_tester(*this);

  // List supported event file format(s).
  timeSystem::EventTimeHandlerFactory<timeSystem::GlastScTimeHandler> glast_sctime_handler;

  // Prepare variables to create application objects.
  std::list<std::string> test_name_cont;
  test_name_cont.push_back("par1");
  test_name_cont.push_back("par2");
  test_name_cont.push_back("par3");
  test_name_cont.push_back("par4");

  // Prepare files to be used in the tests.
  std::string ev_file = prependDataPath("my_pulsar_events_v3.fits");
  std::string sc_file = prependDataPath("my_pulsar_spacecraft_data_v3r1.fits");
  std::string master_pulsardb = prependDataPath("master_pulsardb_v2.fits");
  std::string ev_file_2gti = prependDataPath("my_pulsar_events_2gti.fits");
  std::string sc_file_bogus = prependDataPath("my_pulsar_spacecraft_data_bogus.fits");

  // Loop over parameter sets.
  for (std::list<std::string>::const_iterator test_itor = test_name_cont.begin(); test_itor != test_name_cont.end(); ++test_itor) {
    const std::string & test_name = *test_itor;
    std::string log_file(getMethod() + "_" + test_name + ".log");
    std::string log_file_ref(getMethod() + "_" + test_name + ".ref");
    std::string out_file(getMethod() + "_" + test_name + ".fits");
    std::string out_file_ref(prependOutrefPath(out_file));
    bool ignore_exception(false);

    // Set default parameters.
    st_app::AppParGroup pars(app_tester.getName());
    pars["evfile"] = "";
    pars["scfile"] = "";
    pars["psrdbfile"] = "";
    pars["psrname"] = "ANY";
    pars["outfile"] = "";
    pars["binwidth"] = 1.e-2;
    pars["numbins"] = 1000000;
    pars["timeorigin"] = "MIDDLE";
    pars["usertime"] = 0.;
    pars["userformat"] = "FILE";
    pars["usersys"] = "FILE";
    pars["ra"] = 0.;
    pars["dec"] = 0.;
    pars["ephstyle"] = "FREQ";
    pars["f1f0ratio"] = 0.;
    pars["f2f0ratio"] = 0.;
    pars["p1p0ratio"] = 0.;
    pars["p2p0ratio"] = 0.;
    pars["tcorrect"] = "AUTO";
    pars["solareph"] = "JPL DE405";
    pars["matchsolareph"] = "ALL";
    pars["angtol"] = 1.e-8;
    pars["evtable"] = "EVENTS";
    pars["timefield"] = "TIME";
    pars["sctable"] = "SC_DATA";
    pars["lowfcut"] = .01;
    pars["plot"] = "yes";
    pars["title"] = "DEFAULT";
    pars["leapsecfile"] = "DEFAULT";
    pars["reportephstatus"] = "yes";
    pars["chatter"] = 2;
    pars["clobber"] = "yes";
    pars["debug"] = "no";
    pars["gui"] = "no";
    pars["mode"] = "ql";

    // Set test-specific parameters.
    if ("par1" == test_name) {
      // Test standard computation with DB option.
      pars["evfile"] = ev_file;
      pars["evtable"] = "EVENTS";
      pars["scfile"] = sc_file;
      pars["outfile"] = out_file;
      pars["ephstyle"] = "FREQ";
      pars["ra"] = 85.0482;
      pars["dec"] = -69.3319;
      pars["f1f0ratio"] = 0.0;
      pars["f2f0ratio"] = 0.0;
      pars["timeorigin"] = "mIDdLe";
      pars["psrdbfile"] = master_pulsardb;
      pars["psrname"] = "PSR B0540-69";
      pars["binwidth"] = 0.01;
      pars["numbins"] = 1000000;
      pars["tcorrect"] = "BaRY";
      pars["timefield"] = "TIME";
      pars["lowfcut"] = 0.01001;
      pars["plot"] = "No";
      pars["title"] = "My Fourier analysis";
      pars["gui"] = "No";
      log_file_ref = prependOutrefPath(log_file);

    } else if ("par2" == test_name) {
      // Test ephemeris status reporting.
      std::string summary_file("psrdb_summary.txt");
      remove(summary_file.c_str());
      std::ofstream ofs_summary(summary_file.c_str());
      ofs_summary << prependDataPath("psrdb_spin.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_binary.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_remark.txt") << std::endl;
      ofs_summary.close();
      pars["evfile"] = ev_file_2gti;
      pars["evtable"] = "EVENTS";
      pars["scfile"] = sc_file_bogus;
      pars["outfile"] = out_file;
      pars["ephstyle"] = "FREQ";
      pars["ra"] = 85.0482;
      pars["dec"] = -69.3319;
      pars["f1f0ratio"] = 0.0;
      pars["f2f0ratio"] = 0.0;
      pars["timeorigin"] = "mIDdLe";
      pars["psrdbfile"] = "@" + summary_file;
      pars["psrname"] = "PSR J0540-6919";
      pars["binwidth"] = 0.01;
      pars["numbins"] = 1000000;
      pars["tcorrect"] = "BaRY";
      pars["timefield"] = "TIME";
      pars["plot"] = "No";
      pars["title"] = "My Fourier analysis";
      pars["gui"] = "No";

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      ofs << "gtpspec: WARNING: The following pulsar ephemeris status are reported." << std::endl;
      ofs << "[1] Remarked \"Test remark entry No.2\" since 53990 MJD (TDB) until 54010 MJD (TDB)" << std::endl;
      ofs << "[2] Remarked \"Test remark entry No.3\" since 54025 MJD (TDB) until 54035 MJD (TDB)" << std::endl;
      ofs << "[3] Remarked \"Test remark entry No.4\" since 54050 MJD (TDB) until 54070 MJD (TDB)" << std::endl;
      ofs << "[4] Remarked \"Test remark entry No.6\" since 53990 MJD (TDB) until 54030 MJD (TDB)" << std::endl;
      ofs << "[5] Remarked \"Test remark entry No.7\" since 54030 MJD (TDB) until 54070 MJD (TDB)" << std::endl;
      ofs << "[6] Remarked \"Test remark entry No.8\" since 53990 MJD (TDB) until 54070 MJD (TDB)" << std::endl;
      std::runtime_error error("Could not find the maximum statistic at any trial frequency in range [0.01, -1]");
      app_tester.writeException(ofs, error);
      ofs.close();

      out_file.erase();
      out_file_ref.erase();
      ignore_exception = true;

    } else if ("par3" == test_name) {
      // Test no reporting of ephemeris status with reportephstatus=no.
      std::string summary_file("psrdb_summary.txt");
      remove(summary_file.c_str());
      std::ofstream ofs_summary(summary_file.c_str());
      ofs_summary << prependDataPath("psrdb_spin.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_binary.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_remark.txt") << std::endl;
      ofs_summary.close();
      pars["evfile"] = ev_file_2gti;
      pars["evtable"] = "EVENTS";
      pars["scfile"] = sc_file_bogus;
      pars["outfile"] = out_file;
      pars["ephstyle"] = "FREQ";
      pars["ra"] = 85.0482;
      pars["dec"] = -69.3319;
      pars["f1f0ratio"] = 0.0;
      pars["f2f0ratio"] = 0.0;
      pars["timeorigin"] = "mIDdLe";
      pars["psrdbfile"] = "@" + summary_file;
      pars["psrname"] = "PSR J0540-6919";
      pars["binwidth"] = 0.01;
      pars["numbins"] = 1000000;
      pars["tcorrect"] = "BaRY";
      pars["timefield"] = "TIME";
      pars["plot"] = "No";
      pars["title"] = "My Fourier analysis";
      pars["gui"] = "No";
      pars["reportephstatus"] = "no";

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      std::runtime_error error("Could not find the maximum statistic at any trial frequency in range [0.01, -1]");
      app_tester.writeException(ofs, error);
      ofs.close();

      out_file.erase();
      out_file_ref.erase();
      ignore_exception = true;

    } else if ("par4" == test_name) {
      // Test reporting of database creation history.
      std::string summary_file("psrdb_summary.txt");
      remove(summary_file.c_str());
      std::ofstream ofs_summary(summary_file.c_str());
      ofs_summary << prependDataPath("psrdb_spin.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_binary.txt") << std::endl;
      ofs_summary << prependDataPath("psrdb_remark.txt") << std::endl;
      ofs_summary.close();
      pars["evfile"] = ev_file_2gti;
      pars["evtable"] = "EVENTS";
      pars["scfile"] = sc_file_bogus;
      pars["outfile"] = out_file;
      pars["ephstyle"] = "FREQ";
      pars["ra"] = 85.0482;
      pars["dec"] = -69.3319;
      pars["f1f0ratio"] = 0.0;
      pars["f2f0ratio"] = 0.0;
      pars["timeorigin"] = "mIDdLe";
      pars["psrdbfile"] = "@" + summary_file;
      pars["psrname"] = "PSR J0540-6919";
      pars["binwidth"] = 0.01;
      pars["numbins"] = 1000000;
      pars["tcorrect"] = "AuTO";
      pars["timefield"] = "TIME";
      pars["plot"] = "No";
      pars["title"] = "My Fourier analysis";
      pars["gui"] = "No";
      pars["reportephstatus"] = "no";
      pars["chatter"] = 4;

      remove(log_file_ref.c_str());
      std::ofstream ofs(log_file_ref.c_str());
      ofs << "gtpspec: INFO: ==========================" << std::endl;
      ofs << "gtpspec: INFO: Pulsar ephemerides are loaded and/or filtered as follows:" << std::endl;
      ofs << "gtpspec: INFO:    Load TEXTDB SPIN_PARAMETERS(FREQ) FILENAME='psrdb_spin.txt'" << std::endl;
      ofs << "gtpspec: INFO:    Load TEXTDB ORBITAL_PARAMETERS(DD) FILENAME='psrdb_binary.txt'" << std::endl;
      ofs << "gtpspec: INFO:    Load TEXTDB REMARKS FILENAME='psrdb_remark.txt'" << std::endl;
      ofs << "gtpspec: INFO:    Filter by pulsar name 'PSR J0540-6919'" << std::endl;
      ofs << "gtpspec: INFO:    Filter by solar system ephemeris 'JPL DE405'" << std::endl;
      ofs << "gtpspec: INFO: ==========================" << std::endl;
      ofs << "gtpspec: INFO: --------------------------" << std::endl;
      ofs << "gtpspec: INFO: Arrival time corrections are applied as follows:" << std::endl;
      ofs << "gtpspec: INFO:    Barycentric correction: Applied if necessary" << std::endl;
      ofs << "gtpspec: INFO:    Binary demodulation: Not applied" << std::endl;
      ofs << "gtpspec: INFO:    Pdot cancellation: Applied" << std::endl;
      ofs << "gtpspec: INFO: Following time system(s) are listed for this task:" << std::endl;
      ofs << "gtpspec: INFO:    Spin ephemeri(de)s are defined in: None" << std::endl;
      ofs << "gtpspec: INFO:    Orbital ephemeri(de)s are defined in: None" << std::endl;
      ofs << "gtpspec: INFO:    Pdot cancellation will be performed in: TDB" << std::endl;
      ofs << "gtpspec: INFO:    Time series analysis will be performed in: TDB" << std::endl;
      ofs << "gtpspec: INFO: --------------------------" << std::endl;
      std::runtime_error error("Could not find the maximum statistic at any trial frequency in range [0.01, -1]");
      app_tester.writeException(ofs, error);
      ofs.close();

      out_file.erase();
      out_file_ref.erase();
      ignore_exception = true;

    } else {
      // Skip this iteration.
      continue;
    }

    // Test the application.
    app_tester.test(pars, log_file, log_file_ref, out_file, out_file_ref, ignore_exception);
  }
}

void PeriodSearchTestApp::testPeriodicityTestApp() {
  setMethod("testPeriodicityTestApp");

  // Create an application tester object.
  PeriodicityTestAppTester app_tester(*this);

  // List supported event file format(s).
  timeSystem::EventTimeHandlerFactory<timeSystem::GlastScTimeHandler> glast_sctime_handler;

  // Prepare variables to create application objects.
  std::list<std::string> test_name_cont;
  test_name_cont.push_back("par1");

  // Prepare files to be used in the tests.
  std::string ev_file = prependDataPath("my_pulsar_events_phase_v3r1.fits");

  // Loop over parameter sets.
  for (std::list<std::string>::const_iterator test_itor = test_name_cont.begin(); test_itor != test_name_cont.end(); ++test_itor) {
    const std::string & test_name = *test_itor;
    std::string log_file(getMethod() + "_" + test_name + ".log");
    std::string log_file_ref(getMethod() + "_" + test_name + ".ref");
    std::string out_file(getMethod() + "_" + test_name + ".fits");
    std::string out_file_ref(prependOutrefPath(out_file));
    bool ignore_exception(false);

    // Set default parameters.
    st_app::AppParGroup pars(app_tester.getName());
    pars["evfile"] = "";
    pars["outfile"] = "";
    pars["numphase"] = 10;
    pars["numharm"] = 10;
    pars["maxharm"] = 10;
    pars["evtable"] = "EVENTS";
    pars["pphasefield"] = "PULSE_PHASE";
    pars["plot"] = "yes";
    pars["title"] = "DEFAULT";
    pars["chatter"] = 2;
    pars["clobber"] = "yes";
    pars["debug"] = "no";
    pars["gui"] = "no";
    pars["mode"] = "ql";

    // Set test-specific parameters.
    if ("par1" == test_name) {
      // Test standard computation with DB option.
      pars["evfile"] = ev_file;
      pars["evtable"] = "EVENTS";
      pars["pphasefield"] = "PULSE_PHASE";
      pars["numphase"] = 10;
      pars["numharm"] = 3;
      pars["maxharm"] = 5;
      pars["outfile"] = out_file;
      pars["plot"] = "nO";
      pars["title"] = "All statistical test results";
      pars["gui"] = "No";
      log_file_ref = prependOutrefPath(log_file);

    } else {
      // Skip this iteration.
      continue;
    }

    // Test the application.
    app_tester.test(pars, log_file, log_file_ref, out_file, out_file_ref, ignore_exception);
  }
}

void PeriodSearchTestApp::testAllStats(const std::string & prefix, const std::vector<double> & events, double t_start, double t_stop,
  double center, double step, long num_trials, double epoch, int num_bins,
  double fourier_width, int fourier_num_bins, double fourier_min_freq, double fourier_max_freq, bool plot) {
  setMethod("testAllStats");

  double duration = t_stop - t_start;

  // Test ChiSquared case.
  ChiSquaredTestArray chi2_test(num_trials, num_bins);
  FoldingAnalysis chi2_search(&chi2_test, center, step, epoch, duration, "Hz");
  testOneSearch(events, chi2_search, "Folding Analysis: Chi Squared Statistic", prefix + "-chi-sq.fits", plot);

  // Test Z2n case.
  Z2nTestArray z2n_test(num_trials, num_bins);
  FoldingAnalysis z2n_search(&z2n_test, center, step, epoch, duration, "Hz");

  testOneSearch(events, z2n_search, "Folding Analysis: Z2n Statistic", prefix + "-z2n.fits", plot);

  // Test Rayleigh case.
  RayleighTestArray rayleigh_test(num_trials);
  FoldingAnalysis rayleigh_search(&rayleigh_test, center, step, epoch, duration, "Hz");

  testOneSearch(events, rayleigh_search, "Folding Analysis: Rayleigh Statistic", prefix + "-rayleigh.fits", plot);

  // Test H case.
  HTestArray h_test(num_trials, num_bins);
  FoldingAnalysis h_search(&h_test, center, step, epoch, duration, "Hz");

  testOneSearch(events, h_search, "Folding Analysis: H Statistic", prefix + "-h.fits", plot);

  // Create analysis object.
  FourierAnalysis fourier_search(t_start, t_stop, fourier_width, fourier_num_bins, "Hz", events.size());

  testOneSearch(events, fourier_search, "Fourier Analysis: Power Spectrum", prefix + "-fourier.fits", plot,
    fourier_min_freq, fourier_max_freq);
}

void PeriodSearchTestApp::testOneSearch(const std::vector<double> & events, PeriodSearch & search, const std::string & plot_title,
  const std::string & out_file, bool plot, double min_freq, double max_freq) {
  // Create an application tester object.
  PeriodSearchAppTester tester(*this);

  // Get the viewer.
  StatisticViewer & viewer = search.getViewer();

  // Fill the data into the search object.
  for (std::vector<double>::const_iterator itor = events.begin(); itor != events.end(); ++itor) {
    search.fill(*itor);
  }

  // Perform the search operation.
  search.computeStat();
  search.updateViewer(min_freq, max_freq);

  // Find the template file.
  std::string template_file = prependDataPath("period-search-out.tpl");

  // Create output file.
  tip::IFileSvc::instance().createFile(out_file, template_file, true);

  // Open output file.
  tip::Table * out_table(tip::IFileSvc::instance().editTable(out_file, "POWER_SPECTRUM"));

  // Write the summary to the output header, and the data to the output table.
  viewer.write(*out_table);

  // Close output file. This is necessary before checking output file against reference file.
  delete out_table;

  // Check the result against its reference file in data/outref/ directory.
  tester.checkOutputFits(out_file, prependOutrefPath(out_file));

  // Plot if requested.
  viewer.setTitle(plot_title);
  viewer.setLabel(0, "Hz");
  if (plot) viewer.plot();
}

st_app::StAppFactory<PeriodSearchTestApp> g_factory("test_periodSearch");
