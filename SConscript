# -*- python -*-
# $Id: SConscript,v 1.24 2009/10/23 22:49:33 hirayama Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>
# Version: periodSearch-10-08-01
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('periodSearchLib', depsOnly = 1)
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
