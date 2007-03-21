/**
    \mainpage periodSearch package

    \author  Masaharu Hirayama hirayama@jca.umbc.edu
    \author  James Peachey peachey@milkyway.gsfc.nasa.gov

    \section synopsis Synopsis
    This package contains a library and application gtpsearch which searches for pulsation frequencies
    near to a known, guessed or estimated reference frequency. It has
    several test statistics available, including Chi-squared (X2), Z-squared-n (Z2n),
    H, and Rayleigh (Z2n with n == 1). For details about how
    each of these is computed, see
    <a href="http://glast.gsfc.nasa.gov/ssc/dev/psr_tools/algorithm.html"> pulsar_algorithms </a>

    \section parameters Application Parameters

    \subsection key Key To Parameter Descriptions
\verbatim
Automatic parameters:
par_name [ = value ] type

Hidden parameters:
(par_name = value ) type

Where "par_name" is the name of the parameter, "value" is the
default value, and "type" is the type of the parameter. The
type is enclosed in square brackets.

Examples:
infile [file]
    Describes an automatic (queried) file-type parameter with
    no default value.

(plot = yes) [bool]
    Describes a hidden bool-type parameter named plot, whose
    default value is yes (true).
\endverbatim

    \subsection general gtpsearch General Parameters
\verbatim
algorithm = Chi2 [string]
    This parameter sets the type of search algortihm to use.
    Current valid choices are Chi2, Z2n, or H. The Rayleigh
    test may be obtained by choosing Z2n here and setting the
    number of bins used for each trial (numbins parameter) to 1.

evfile [file name]
    Name of input event file, FT1 format or equivalent.

psrname = ANY [string]
    The name of the pulsar, used to select only ephemerides
    valid for a particular pulsar.

ephstyle = DB [string]
    Specifies how the ephemerides will be supplied. If ephstyle
    is DB, a pulsar database file (GLAST D4 FITS format) will
    be used. If ephstyle is FREQ (PER), the user will supply values
    for the frequency (period) and its derivatives at the time
    given by the epoch parameter.

scanstep = 0.5 [double]
    Size of steps for trials, in units of the Fourier resolution.

ephepoch = 0. [string]
    The epoch, or time origin, for the ephemeris used in this search.

timeformat = FILE [string]
    String describing the representation used for the ephepoch.
    Valid choices are FILE, MJD and GLAST (MET). If FILE is chosen,
    the time format specified in the input event file header will be
    used.

timesys = FILE [string]
    String describing the time system used for the ephepoch.
    Valid choices are FILE, TAI, TDB, TT and UTC. If FILE is chosen,
    the time system specified in the input event file header (TIMESYS
    keyword) will be used.

numtrials = 100 [integer]
    The number of separate trials to perform. The larger this number,
    the wider the search around the central frequency.

numbins = 10 [integer]
    The number of bins in each trial. For the Rayleigh test, set
    algorithm to Z2n, and set numbins to 1.

timeorigin = MIDDLE [string]
    Selects the origin of time for the periodicity test. Valid
    choices are START (taken from the evfile), STOP, MIDDLE (mid-way
    between START and STOP) and USER (user will supply explicitly using
    usertime, userformat and usersys parameters.)

usertime = 0. [string]
    User-specified time origin for the periodicity test, used only
    if timeorigin parameter is USER.

userformat = FILE [string]
    String describing the representation used for the usertime.
    Valid choices are FILE, MJD and GLAST (MET). If FILE is chosen,
    the time format specified in the input event file header will be
    used. Used only if timeorigin parameter is USER.

usersys = FILE [string]
    String describing the time system used for the usertime.
    Valid choices are FILE, TAI, TDB, TT and UTC. If FILE is chosen,
    the time system specified in the input event file header (TIMESYS
    keyword) will be used. Used only if timeorigin parameter is USER.

(psrdbfile = DEFAULT) [file name]
    Name of pulsar ephemerides database file, in GLAST D4
    FITS format. If psrdbfile is DEFAULT, the canonical pulsar
    database file (master_pulsardb.fits), which is distributed
    with the pulsar tools, will be used.

(cancelpdot = no) [bool]
    Determines whether frequency derivatives will be used to
    account for known variations in the frequency.

(demodbin = AUTO) [string]
    A three-way switch. If demodbin is AUTO, binary demodulation
    will be applied only to pulsars for which orbital parameters
    exist in the pulsar database. If demodbin is YES, the computation
    is the same as for AUTO, but orbital parameters are required
    to be present and an error will be thrown if none are found.
    If demodbin is NO, binary demodulation will not be performed
    regardless of whether orbital ephemerides exist in the database.

(evtable = EVENTS) [string]
    The name of the table containing the event data. The default
    value is correct for the FT1 format.

(timefield = TIME) [string]
    This is the name of the field containing the time values for
    time binning. The default value is consistent with the FT1
    format.

(plot = yes) [bool]
    Display plot of results.

(title = DEFAULT) [string]
    Title for the graph. By default a title indicating the type of
    test and other pertinent information will be displayed.

(leapsecfile = DEFAULT) [file name]
    The file containing the name of the leap second table, in
    OGIP-compliant leap second table format. If leapsecfile is
    the string DEFAULT, the default leapsec file (leapsec.fits),
    which is distributed with the extFiles package, will be used.

\endverbatim

    \subsection gtpsearch_freq_par gtpsearch Frequency Parameters

\verbatim
f0 = 1. [double]
    The value of the frequency at the epoch. This is the central
    frequency for the search.

f1 = 0. [double]
    The value of the first time derivative of the frequency at the
    epoch. Only used if the cancelpdot parameter is yes.

f2 = 0. [double]
    The value of the second time derivative of the frequency at the
    epoch. Only used if the cancelpdot parameter is yes.
\endverbatim

    \subsection gtpsearch_per_par gtpsearch Period Parameters

\verbatim
p0 = 1. [double]
    The value of the period at the epoch. Only used if
    ephstyle is PER.

p1 = 0. [double]
    The value of the first time derivative of the period at the
    epoch. Only used if ephstyle is PER.

p2 = 0. [double]
    The value of the second time derivative of the period at the
    epoch. Only used if ephstyle is PER.
\endverbatim

*/
