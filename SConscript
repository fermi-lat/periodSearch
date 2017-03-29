# -*- python -*-
# $Id: SConscript,v 1.32 2011/05/16 21:02:19 hirayama Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>, Joe Asercion <joseph.a.asercion@nasa.gov>
# Version: periodSearch-10-08-07
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

periodSearchLib = libEnv.StaticLibrary('periodSearch', listFiles(['src/*.cxx']))

progEnv.Tool('periodSearchLib')
gtpsearchBin = progEnv.Program('gtpsearch', listFiles(['src/gtpsearch/*.cxx']))
gtpspecBin = progEnv.Program('gtpspec', listFiles(['src/gtpspec/*.cxx']))
gtptestBin = progEnv.Program('gtptest', listFiles(['src/gtptest/*.cxx']))
test_periodSearchBin = progEnv.Program('test_periodSearch', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerTargets', package = 'periodSearch',
             staticLibraryCxts = [[periodSearchLib, libEnv]],
             binaryCxts = [[gtpsearchBin,progEnv], [gtpspecBin,progEnv],
                           [gtptestBin,progEnv]],
             testAppCxts = [[test_periodSearchBin, progEnv]],
             includes = listFiles(['periodSearch/*.h']),
             pfiles = listFiles(['pfiles/*.par']),
             data = listFiles(['data/*'], recursive = True))
