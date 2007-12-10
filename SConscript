import glob,os,platform

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

periodSearchLib = libEnv.StaticLibrary('periodSearch', listFiles(['src/*.cxx']))

progEnv.Tool('periodSearchLib')
gtpsearchBin = progEnv.Program('gtpsearch', listFiles(['src/gtpsearch/*.cxx']))

progEnv.Tool('registerObjects', package = 'periodSearch', libraries = [periodSearchLib], binaries = [gtpsearchBin], includes = listFiles(['periodSearch/*.h']), pfiles = listFiles(['pfiles/*.par']))