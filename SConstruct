#!python
import subprocess
import os
import platform

# SConscript('color_SConscript')
# Import( 'env' )

SConscript('color_SConscript')
Import( 'env' )

SConscript( "modules/FemtoDstFormat/SConstruct" )
SConscript( "modules/FemtoPairFormat/SConstruct" )

env.Append(CPPPATH    = [ "modules/" ])
env.Append(CXXFLAGS   = ['-std=c++11', '-DJDB_LOG_LEVEL=0'])
JDB_LIB = os.environ.get("JDB_LIB", "" )
env.Append(LIBPATH    = [ JDB_LIB + "/lib/" ] )
env.Append(LIBS       = [ "TMVA", "libRooBarbCore.a", "libRooBarbConfig.a", "libRooBarbTasks.a", "libRooBarbRootAna.a", "libRooBarbUnitTest.a", "libRooBarbExtra.a" ] )

debug = ARGUMENTS.get('debug', 0)
if int(debug):
	env.Append(CXXFLAGS = '-g')

# # ROOT
ROOTCFLAGS      = subprocess.check_output( ['root-config',  '--cflags'] ).rstrip().decode('utf-8').split( " " )
ROOTLIBS        = subprocess.check_output( ["root-config",  "--libs"] ).decode('utf-8')
ROOT_SYS        = os.environ.get( "ROOTSYS", "" )

env.Append(CPPPATH   = [ ROOT_SYS + "/include/", JDB_LIB + "/include", "modules/StarClassLibrary" ])
env.Append(CXXFLAGS  = ROOTCFLAGS )
env.Append(LINKFLAGS = ROOTCFLAGS )

env[ "_LIBFLAGS" ] = env[ "_LIBFLAGS" ] + " " + ROOTLIBS + " "
# REMOVE "-pthread" from clang link options in OS X
# assuming darwin=clang but must be a better way...
if "Darwin" in platform.platform() :
	env[ "LINKFLAGS" ].remove( "-pthread" )

env.Program( target="bin/pico.app", source=[ "modules/FemtoPairFormat/DictionaryFemtoPair.cpp", "modules/FemtoDstFormat/DictionaryFemtoDst.cpp", "Engine.cpp"] )
