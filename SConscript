# -*- python -*-
# $Id: SConscript,v 1.10 2008/09/30 22:30:35 glastrm Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>
# Version: periodSearch-10-02-01
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

progEnv.Tool('registerObjects', package = 'periodSearch', libraries = [periodSearchLib], binaries = [gtpsearchBin, gtpspecBin, gtptestBin],
             testApps = [test_periodSearchBin], includes = listFiles(['periodSearch/*.h']), pfiles = listFiles(['pfiles/*.par']),
             data = listFiles(['data/*'], recursive = True))
