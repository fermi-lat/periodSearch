/**
    \mainpage periodSearch package

    \author  Masaharu Hirayama hirayama@jca.umbc.edu
    \author  James Peachey James.Peachey-1@nasa.gov

    \section synopsis Synopsis

    This package contains a library and two applications, gtpsearch and gtpspec.
    The application gtpsearch searches for pulsations at frequencies
    near to a known, guessed or estimated reference frequency. It has
    several test statistics available, including Chi-squared (X2), Z-squared-n (Z2N),
    H, and Rayleigh (Z2N with N == 1). 
    The application gtpspec searches for pulsations in a much wider range of frequencies.
    It uses the Discrete Fast Fourier Transfer (FFT) technique to compute power spectrum 
    density. The application gtptest applies all the statistical tests available with
    gtpsearch application to a series of pulse phase values stored in given event files.
    For details about how each of these is computed, see
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

    \subsection gtpsearch_parameters gtpsearch Parameters
\verbatim
evfile [file name]
    Name of input event file, FT1 format or equivalent.

scfile [file name]
    Name of input spacecraft data file, FT2 format or equivalent.

psrdbfile [file name]
    Name of pulsar ephemerides database file, in GLAST D4
    FITS format.

psrname = ANY [string]
    Name of the pulsar, used to select only ephemerides valid for a
    particular pulsar.

outfile [file name]
    Name of output FITS file that contains a search result.  If
    outfile is NONE, no FITS output will be created.

algorithm = CHI2 [enumerated string (CHI2|Z2N|H)]
    Type of search algorithm to use.  Current valid choices are CHI2
    (for the chi-squared test), Z2N (the Z2n/Rayleigh test), or H (the
    H test). The Rayleigh test may be obtained by choosing Z2N here
    and setting the number of harmonics used for each trial (the
    numharm parameter) to 1.

numphase = 10 [integer]
    Number of phase bins in each trial for the chi-squared test.  This
    parameter only has effect if algorithm is CHI2.

numharm = 10 [integer]
    Number of harmonics in each trial for the Z2n/Rayleigh test.  This
    parameter only has effect if algorithm is Z2N.  For the Rayleigh
    test, set algorithm to Z2N, and set numharm to 1.

maxharm = 10 [integer]
    Maximum number of harmonics in each trial for the H test.  This
    parameter only has effect if algorithm is H.

scanstep = 0.5 [double]
    Size of steps for trials, in units of the Fourier resolution.

numtrials = 100 [integer]
    Number of separate trials to perform. The larger this number, the
    wider the search around the central frequency.

timeorigin = MIDDLE [enumerated string (START|STOP|MIDDLE|USER)]
    Origin of time for the periodicity test. If START or STOP is
    chosen, the start or stop time is taken from the input event
    file(s) and used as the time origin. If MIDDLE is chosen, the
    mid-time between START and STOP is used.  If USER is chosen, user
    will supply explicitly using usertime, userformat and usersys
    parameters.

usertime = 0. [string]
    User-specified time origin for the periodicity test, used only
    if timeorigin parameter is USER.

userformat = FILE [enumerated string (FILE|MJD|GLAST)]
    String describing the representation used for the usertime.  If
    FILE is chosen, the time format specified in the input event file
    header will be used. Used only if timeorigin parameter is USER.

usersys = FILE [enumerated string (FILE|TAI|TDB|TT|UTC)]
    String describing the time system used for the usertime. If FILE
    is chosen, the time system specified in the input event file
    header (TIMESYS keyword) will be used. Used only if timeorigin
    parameter is USER.

ephstyle = DB [enumerated string (DB|FREQ|PER)]
    Method to specify how the ephemerides will be supplied. If
    ephstyle is DB, a pulsar database file (GLAST D4 FITS format) will
    be used. If ephstyle is FREQ (PER), the user will supply values
    for the frequency (period) and its derivatives at the time given
    by the epoch parameter.

ephepoch = 0. [string]
    Reference epoch of the ephemeris that will be manually specified.
    This parameter only has effect if ephstyle is FREQ or PER.

timeformat = FILE [enumerated string (FILE|MJD|GLAST)]
    String describing the representation used for the ephepoch.  If
    FILE is chosen, the time format specified in the input event file
    header will be used.

timesys = FILE [enumerated string (FILE|TAI|TDB|TT|UTC)]
    String describing the time system used for the ephepoch.  If FILE
    is chosen, the time system specified in the input event file
    header (TIMESYS keyword) will be used.

ra [double]
    Right Ascension of point source in degrees for which to perform
    the barycentric correction.  This parameter only has effect if
    ephstyle is FREQ or PER.

dec [double]
    Declination of point source in degrees for which to perform the
    barycentric correction.  This parameter only has effect if
    ephstyle is FREQ or PER.

f0 = 1. [double]
    Value of the frequency at the time given by the epoch parameter.
    This parameter only has effect if ephstyle is FREQ.

f1 = 0. [double]
    Value of the first time derivative of the frequency at the time
    given by the epoch parameter.  This parameter only has effect if
    ephstyle is FREQ.

f2 = 0. [double]
    Value of the second time derivative of the frequency at the time
    given by the epoch parameter.  This parameter only has effect if
    ephstyle is FREQ.

p0 = 1. [double]
    Value of the period at the time given by the epoch parameter.
    This parameter only has effect if ephstyle is PER.

p1 = 0. [double]
    Value of the first time derivative of the period at the time given
    by the epoch parameter.  This parameter only has effect if
    ephstyle is PER.

p2 = 0. [double]
    Value of the second time derivative of the period at the time
    given by the epoch parameter.  This parameter only has effect if
    ephstyle is PER.

(tcorrect = AUTO) [enumerated string (NONE|AUTO|BARY|BIN|PDOT|ALL)]
    Set of arrival time corrections to apply. If tcorrect is NONE, no
    corrections will be applied. If tcorrect is BARY, only the
    barycentric correction will be applied. If tcorrect is BIN, an
    appropriate orbital ephemeris is searched for in the pulsar
    database, and if found, binary demodulation will be applied after
    the barycentric correction, and if not, an error will be
    generated.  If tcorrect is PDOT, an appropriate spin ephemeris is
    searched for in the pulsar database, and if found, pdot
    cancellation will be applied after the barycentric correction, and
    if not, an error will be generated.  If tcorrect is ALL, both
    actions for the BIN option and the PDOT option will be taken.  If
    tcorrect is AUTO, the barycentric correction will be applied, and
    the binary demodulation will be applied only when an orbital
    ephemeris is available in the pulsar database, then the pdot
    cancellation will applied only when a spin ephemeris is available
    in the pulsar database.

(solareph = JPL DE405) [enumerated string (JPL DE200|JPL DE405)]
    Solar system ephemeris for the barycentric correction.

(matchsolareph = ALL) [enumerated string (NONE|EVENT|PSRDB|ALL)]
    String that controls whether to use the name of the solar system
    ephemeris given by the solareph parameter to check the input event
    data file and to select ephemerides in the pulsar database.  If
    matchsolareph is EVENT, all the input event files are required to
    use the same solar system ephemeris for the barycentric correction
    as the one given by the solareph parameter, and an error will be
    generated if otherwise. Such an error may occur when an input
    event data file has already been applied the barycentric
    correction with a different solar system ephemeris. If
    matchsolareph is PSRDB, the string given by the solareph parameter
    is used to select the ephemerides.  If matchsolareph is ALL, both
    actions for the EVENT option and the PSRDB option will be taken.
    If matchsolareph is NONE, no selection will be performed.

(angtol = 1.e-8) [double]
    Angular tolerance in degrees in comparison of two source
    positions, one for which the barycentric correction is performed,
    and another given by RA_NOM and DEC_NOM header keyword of an event
    file to which the barycentric correction has already been applied.
    This parameter only has effect if the barycentric correction has
    been applied to an input event data file.  If the two source
    positions are separate from each other by this amount or less,
    then they will be considered to be the same position.  Otherwise
    an error will be generated.

(evtable = EVENTS) [string]
    Name of the FITS table containing the event data.

(timefield = TIME) [string]
    Name of the field containing the time values for temporal
    analysis.

(sctable = SC_DATA) [string]
    Name of the FITS table containing the spacecraft data.

(plot = yes) [bool]
    If plot is yes, the result will be displayed in a separate plot
    window, as well as the numerical results will be output in a text
    screen.  If plot is no, only the text output will be displayed,
    and no plot window will be displayed.

(title = DEFAULT) [string]
    Title for the graph. By default a title indicating the type of
    test and other pertinent information will be displayed.

(leapsecfile = DEFAULT) [file name]
    Name of the file containing the name of the leap second table, in
    OGIP-compliant leap second table format. If leapsecfile is the
    string DEFAULT, the default leap-second file (leapsec.fits), which
    is distributed with the extFiles package, will be used.
\endverbatim

    \subsection gtpspec_parameters gtpspec Parameters

\verbatim
evfile [file name]
    Name of input event file, FT1 format or equivalent.

scfile [file name]
    Name of input spacecraft data file, FT2 format or equivalent.

psrdbfile [file name]
    Name of pulsar ephemerides database file, in GLAST D4
    FITS format.

psrname = ANY [string]
    Name of the pulsar, used to select only ephemerides valid for a
    particular pulsar.

outfile [file name]
    Name of output FITS file that contains a search result.  If
    outfile is NONE, no FITS output will be created.

binwidth = 1.e-2 [double]
    Width of time bins in seconds to be used to internally create a
    binned light curve which will be passed to the Fast Fourier
    Transform (FFT) algorithm.  The product of binwidth and numbins
    will be the length of data to be transformed at once.

numbins = 1000000 [integer]
    Number of time bins to be used to internally create a binned light
    curve which will be passed to the Fast Fourier Transform (FFT)
    algorithm.  The product of binwidth and numbins will be the length
    of data to be transformed at once.

timeorigin = MIDDLE [enumerated string (START|STOP|MIDDLE|USER)]
    Origin of time for the periodicity test. If START or STOP is
    chosen, the start or stop time is taken from the input event
    file(s) and used as the time origin. If MIDDLE is chosen, the
    mid-time between START and STOP is used.  If USER is chosen, user
    will supply explicitly using usertime, userformat and usersys
    parameters.

usertime = 0. [string]
    User-specified time origin for the periodicity test, used only
    if timeorigin parameter is USER.

userformat = FILE [enumerated string (FILE|MJD|GLAST)]
    String describing the representation used for the usertime.  If
    FILE is chosen, the time format specified in the input event file
    header will be used. Used only if timeorigin parameter is USER.

usersys = FILE [enumerated string (FILE|TAI|TDB|TT|UTC)]
    String describing the time system used for the usertime. If FILE
    is chosen, the time system specified in the input event file
    header (TIMESYS keyword) will be used. Used only if timeorigin
    parameter is USER.

ra [double]
    Right Ascension of point source in degrees for which to perform
    the barycentric correction.

dec [double]
    Declination of point source in degrees for which to perform the
    barycentric correction.

ephstyle = FREQ [enumerated string (FREQ|PER)]
    Method to specify how the ephemeris for pdot cancellation will be
    supplied.  If ephstyle is FREQ, the user will supply values for
    the frequency and its derivatives at the time given by the epoch
    parameter.  If ephstyle if PER, the user will supply values for
    the period and its derivatives at the time given by the epoch
    parameter.

f1f0ratio = 0. [double]
    Ratio of frequency first time derivative to frequency at the time
    given by the timeorigin parameter.  This parameter only has effect
    if ephstyle is FREQ.

f2f0ratio = 0. [double]
    Ratio of frequency second time derivative to frequency at the time
    given by the timeorigin parameter.  This parameter only has effect
    if ephstyle is FREQ.

p1p0ratio = 0. [double]
    Ratio of period first time derivative to period at the time given
    by the timeorigin parameter.  This parameter only has effect if
    ephstyle is PER.

p2p0ratio = 0. [double]
    Ratio of period second time derivative to period at the time given
    by the timeorigin parameter.  This parameter only has effect if
    ephstyle is PER.

(tcorrect = AUTO) [enumerated string (NONE|AUTO|BARY|BIN|PDOT|ALL)]
    Set of arrival time corrections to apply. If tcorrect is NONE, no
    corrections will be applied. If tcorrect is BARY, only the
    barycentric correction will be applied. If tcorrect is BIN, an
    appropriate orbital ephemeris is searched for in the pulsar
    database, and if found, binary demodulation will be applied after
    the barycentric correction, and if not, an error will be
    generated.  If tcorrect is PDOT, an appropriate spin ephemeris is
    searched for in the pulsar database, and if found, pdot
    cancellation will be applied after the barycentric correction, and
    if not, an error will be generated.  If tcorrect is ALL, both
    actions for the BIN option and the PDOT option will be taken.  If
    tcorrect is AUTO, the barycentric correction will be applied, and
    the binary demodulation will be applied only when an orbital
    ephemeris is available in the pulsar database, then the pdot
    cancellation will applied only when a spin ephemeris is available
    in the pulsar database.

(solareph = JPL DE405) [enumerated string (JPL DE200|JPL DE405)]
    Solar system ephemeris for the barycentric correction.

(matchsolareph = ALL) [enumerated string (NONE|EVENT|PSRDB|ALL)]
    String that controls whether to use the name of the solar system
    ephemeris given by the solareph parameter to check the input event
    data file and to select ephemerides in the pulsar database.  If
    matchsolareph is EVENT, all the input event files are required to
    use the same solar system ephemeris for the barycentric correction
    as the one given by the solareph parameter, and an error will be
    generated if otherwise. Such an error may occur when an input
    event data file has already been applied the barycentric
    correction with a different solar system ephemeris. If
    matchsolareph is PSRDB, the string given by the solareph parameter
    is used to select the ephemerides.  If matchsolareph is ALL, both
    actions for the EVENT option and the PSRDB option will be taken.
    If matchsolareph is NONE, no selection will be performed.

(angtol = 1.e-8) [double]
    Angular tolerance in degrees in comparison of two source
    positions, one for which the barycentric correction is performed,
    and another given by RA_NOM and DEC_NOM header keyword of an event
    file to which the barycentric correction has already been applied.
    This parameter only has effect if the barycentric correction has
    been applied to an input event data file.  If the two source
    positions are separate from each other by this amount or less,
    then they will be considered to be the same position.  Otherwise
    an error will be generated.

(evtable = EVENTS) [string]
    Name of the FITS table containing the event data.

(timefield = TIME) [string]
    Name of the field containing the time values for temporal
    analysis.

(sctable = SC_DATA) [string]
    Name of the FITS table containing the spacecraft data.

(lowfcut = .01) [double]
    Low frequency cut-off for chance probability computation.  This
    parameter may be used to avoid the frequency range that is
    severely affected by artifacts such as data gaps due to the
    orbital motion of the spacecraft.

(plot = yes) [bool]
    If plot is yes, the result will be displayed in a separate plot
    window, as well as the numerical results will be output in a text
    screen.  If plot is no, only the text output will be displayed,
    and no plot window will be displayed.

(title = DEFAULT) [string]
    Title for the graph. By default a title indicating the type of
    test and other pertinent information will be displayed.

(leapsecfile = DEFAULT) [file name]
    Name of the file containing the name of the leap second table, in
    OGIP-compliant leap second table format. If leapsecfile is the
    string DEFAULT, the default leap-second file (leapsec.fits), which
    is distributed with the extFiles package, will be used.
\endverbatim

    \subsection gtptest_parameters gtptest Parameters
\verbatim
evfile [file name]
    Name of input event file, FT1 format or equivalent.

outfile [file name]
    Name of output FITS file that contains a search result.  If
    outfile is NONE, no FITS output will be created.

numphase = 10 [integer]
    Number of phase bins in each trial for the chi-squared test.  This
    parameter only has effect if algorithm is CHI2.

numharm = 10 [integer]
    Number of harmonics in each trial for the Z2n/Rayleigh test.  This
    parameter only has effect if algorithm is Z2N.  For the Rayleigh
    test, set algorithm to Z2N, and set numharm to 1.

maxharm = 10 [integer]
    Maximum number of harmonics in each trial for the H test.  This
    parameter only has effect if algorithm is H.

(evtable = EVENTS) [string]
    Name of the FITS table containing the event data.

(pphasefield = PULSE_PHASE) [string]
    Name of the output column to contain the assigned pulse phase.

(plot = yes) [bool]
    If plot is yes, the result will be displayed in a separate plot
    window, as well as the numerical results will be output in a text
    screen.  If plot is no, only the text output will be displayed,
    and no plot window will be displayed.

(title = DEFAULT) [string]
    Title for the graph. By default a title indicating the type of
    test and other pertinent information will be displayed.
\endverbatim
*/
