/** \file PeriodicityTestApp.cxx
    \brief Implmentation of PeriodicityTestApp class.
    \author Masaharu Hirayama, GSSC
            James Peachey, HEASARC/GSSC
*/
#include "PeriodicityTestApp.h"

#include <limits>
#include <list>
#include <stdexcept>

#include "st_app/AppParGroup.h"

#include "st_facilities/FileSys.h"

#include "st_graph/Axis.h"
#include "st_graph/Engine.h"
#include "st_graph/IPlot.h"
#include "st_graph/Sequence.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "ChiSquaredTestArray.h"
#include "HTestArray.h"
#include "PeriodicityTestArray.h"
#include "RayleighTestArray.h"
#include "Z2nTestArray.h"

static const std::string s_cvs_id = "$Name:  $";

PeriodicityTestApp::PeriodicityTestApp(): m_os("PeriodicityTestApp", "", 2) {
  setName("gtptest");
  setVersion(s_cvs_id);
}

PeriodicityTestApp::~PeriodicityTestApp() throw() {}

void PeriodicityTestApp::run() {
  m_os.setMethod("run()");
  st_app::AppParGroup & pars(getParGroup("gtptest"));

  // Prompt and save.
  pars.Prompt();
  pars.Save();

  // Read parameters to open input event files.
  std::string event_file = pars["evfile"];
  std::string event_extension = pars["evtable"];
  std::string field_name = pars["pphasefield"];

  // Open the event table(s) for reading.
  typedef std::list<const tip::Table *> table_list_type;
  table_list_type table_list;
  st_facilities::FileSys::FileNameCont file_name_cont = st_facilities::FileSys::expandFileList(event_file);
  for (st_facilities::FileSys::FileNameCont::const_iterator itor = file_name_cont.begin(); itor != file_name_cont.end(); ++itor) {
    std::string file_name = *itor;
    const tip::Table * table = tip::IFileSvc::instance().readTable(file_name, event_extension);

    // Check whether pulse phase field exists.
    try {
      table->getFieldIndex(field_name);
    } catch (const tip::TipException &) {
      throw std::runtime_error("Input file \"" + file_name + "\" does not contain a pulse phase column named \"" + field_name + "\"");
    }

    // Add it to the table list.
    table_list.push_back(table);
  }

  // Create a list of periodicity tests.
  typedef std::list<PeriodicityTestArray *> test_list_type;
  test_list_type test_list;
  PeriodicityTestArray * test_array(0);

  // Add periodicity tests to the list.
  long num_phase = pars["numphase"];
  test_array = new ChiSquaredTestArray(1, num_phase);
  test_list.push_back(test_array);
  PeriodicityTestArray & test_to_plot = *test_array;

  test_array = new RayleighTestArray(1);
  test_list.push_back(test_array);

  long num_harm = pars["numharm"];
  test_array = new Z2nTestArray(1, num_harm);
  test_list.push_back(test_array);

  long max_harm = pars["maxharm"];
  test_array = new HTestArray(1, max_harm);
  test_list.push_back(test_array);

  // Loop over events in the event table(s).
  for (table_list_type::const_iterator table_itor = table_list.begin(); table_itor != table_list.end(); ++table_itor) {
    const tip::Table & table = **table_itor;
    for (tip::Table::ConstIterator event_itor = table.begin(); event_itor != table.end(); ++event_itor) {
      const tip::ConstTableRecord & record = *event_itor;

      // Read phase value from this event, as a signle variable of double type.
      double phase_value;
      record[field_name].get(phase_value);

      // Fill the phase value into the tests.
      for (test_list_type::iterator itor = test_list.begin(); itor != test_list.end(); ++itor) {
        PeriodicityTestArray & test_array = **itor;
        test_array.fill(0, phase_value);
      }
    }
  }

  // Get parameters for output.
  std::string out_file = pars["outfile"];
  bool plot = pars["plot"];
  std::string title = pars["title"];
  bool clobber = pars["clobber"];

  // Use default title if user did not specify one.
  std::string title_uc(title);
  for (std::string::iterator itor = title_uc.begin(); itor != title_uc.end(); ++itor) *itor = std::toupper(*itor);
  if (title_uc != "DEFAULT") test_to_plot.getViewer(0).setTitle(title);

  // Compute the statistical test results and write them to the screen.
  for (test_list_type::iterator itor = test_list.begin(); itor != test_list.end(); ++itor) {
    PeriodicityTestArray & test_array = **itor;
    test_array.getViewer(0).write(m_os);
  }

  // Display a plot, if desired.
  if (plot) test_to_plot.getViewer(0).plot();

  // Delete the event table(s).
  for (table_list_type::iterator itor = table_list.begin(); itor != table_list.end(); ++itor) delete *itor;

  // Delete the peridicity tests.
  for (test_list_type::iterator itor = test_list.begin(); itor != test_list.end(); ++itor) delete *itor;
}
