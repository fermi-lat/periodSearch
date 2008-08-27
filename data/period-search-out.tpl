SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                    8 / number of bits per data pixel
NAXIS   =                    0 / number of data axes
EXTEND  =                    T / FITS dataset may contain extensions
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
END
 
XTENSION= 'BINTABLE'           / binary table extension
BITPIX  =                    8 / 8-bit bytes
NAXIS   =                    2 / 2-dimensional binary table
NAXIS1  =                    0 / width of table in bytes
NAXIS2  =                    0 / number of rows in table
PCOUNT  =                    0 / size of special data area
GCOUNT  =                    1 / one data group (required keyword)
TFIELDS =                    5 / number of fields in each row
TTYPE#  = 'FREQUENCY'          / label for field
TFORM#  = 'D       '           / data format of field: 8-byte DOUBLE
TUNIT#  = 'Hz      '           / physical unit of field
TTYPE#  = 'STATISTIC'          / label for field
TFORM#  = 'E       '           / data format of field: 4-byte REAL
EXTNAME = 'POWER_SPECTRUM'     / name of this binary table extension
DATE    = '        '           / file creation date (YYYY-MM-DDThh:mm:ss UT)
CREATOR = '        '           / Name of program that created this file
ORIGIN  = 'LISOC'              / name of organization making file
OBJECT  = '        '
# The following commented-out keywords are used by the HEASOFT/Ftools/Xronos
# application powspec, but are not used by current applications in the
# periodSearch package.
#HDUCLASS= 'OGIP    '
#HDUCLAS1= 'TEMPORALDATA'
#HDUCLAS2= 'POWER SPECTRA'
#HDUCLAS3= 'RESULTS '
#CONTENT = 'GTPSPEC OUTPUT'
#TIMVERSN= 'OGIP/93-003'
#TSTARTI =                    0 / Start time for this extension
#TSTARTF =                   0. / Start time for this extension
#TSTOPI  =                    0 / Stop time for this extension
#TSTOPF  =                   0. / Stop time for this extension
#TIMEUNIT= '        '           / Units for header timing keywords
#TIMEZERI=                    0 / Zero-point offset for TIME column
#TIMEZERF=                   0. / Zero-point offset for TIME column
#AVRGE_1 =                   0. / Avg, count/s in interval
#FREXP_1 =                   0. / Avg. Frac exposure in frame
#VAROB_1 =                   0. / observed variance
#VAREX_1 =                   0. / expected variance
#THRDM_1 =                   0. / third moment
#MININ_1 =                   0. / Minimun Intensity
#MAXIN_1 =                   0. / Maximun Intensity
#EXVAR_1 =                   0. / excess of variance
#CHI2_1  =                   0. / Chi squared
#RMS_1   =                   0. /  RMS variability
#AVRGE_E1=                   0. / Error Avg count/s
#VAROB_E1=                   0. / Error observed variance
#VAREX_E1=                   0. / Error expected variance
#CHI2_E1 =                   0. / Error Chi2
#RMS_E1  =                   0. / RMS variability
END
