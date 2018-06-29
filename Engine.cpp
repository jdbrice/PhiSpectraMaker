

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "PhiHistogram.h"
#include "PhiMvaHistogram.h"
#include "PhiFitter.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"




int main( int argc, char* argv[] ) {
	// loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);

	loguru::init(argc, argv);
	Logger::setGlobalLogLevel( "none" );

	TaskFactory::registerTaskRunner<PhiHistogram>( "PhiHistogram" );
	TaskFactory::registerTaskRunner<PhiMvaHistogram>( "PhiMvaHistogram" );
	TaskFactory::registerTaskRunner<PhiFitter>( "PhiFitter" );

	TaskEngine engine( argc, argv, "PhiMvaHistogram" );


	return 0;
}
