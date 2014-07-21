#include <cppunit/extensions/HelperMacros.h>
#include "l1menu/TriggerMenu.h"

//
// Forward definitions
//
namespace l1menu
{
	class ISample;
	class TriggerRatePlot;
	class MenuRatePlots;
}

/** @brief A cppunit TestFixture to test TriggerRatePlot objects.
 *
 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
 * @date 04/Aug/2013
 */
class TriggerRatePlotUnitTestSuite : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE(TriggerRatePlotUnitTestSuite);
	CPPUNIT_TEST(testConstructingFromTH1);
	CPPUNIT_TEST(testFindThresholdError);
	CPPUNIT_TEST(testLowThresholdPlateau);
	CPPUNIT_TEST_SUITE_END();

protected:
	std::ostream* pVerboseOutput_;
	std::unique_ptr<l1menu::ISample> pSample_;
	std::unique_ptr<l1menu::TriggerMenu> pTriggerMenu_;
	std::unique_ptr<l1menu::MenuRatePlots> pRatePlotsFromDisk_;
public:
	TriggerRatePlotUnitTestSuite();
	void setUp();

protected:
	void testConstructingFromTH1();
	void testFindThresholdError();
	void testLowThresholdPlateau();
};





#include <cppunit/config/SourcePrefix.h>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include "l1menu/ISample.h"
#include "l1menu/TriggerTable.h"
#include "l1menu/ITrigger.h"
#include "l1menu/TriggerRatePlot.h"
#include "l1menu/MenuRatePlots.h"
#include "l1menu/tools/miscellaneous.h"
#include "l1menu/tools/fileIO.h"
#include "TestParameters.h"
#include <TFile.h>
#include <TH1.h>

CPPUNIT_TEST_SUITE_REGISTRATION(TriggerRatePlotUnitTestSuite);

TriggerRatePlotUnitTestSuite::TriggerRatePlotUnitTestSuite()
{
	pVerboseOutput_=nullptr;
	//pVerboseOutput_=&std::cout;
}

void TriggerRatePlotUnitTestSuite::setUp()
{
	std::string inputSampleFilename=TestParameters<std::string>::instance().getParameter( "TEST_SAMPLE_FILENAME" );
	std::string inputMenuFilename=TestParameters<std::string>::instance().getParameter( "TEST_MENU_FILENAME" );
	std::string ratePlotsFilename=TestParameters<std::string>::instance().getParameter( "TEST_RATEPLOT_FILENAME" );

	// Add a newline, because cppunit starts this function with half a line already written
	if( pVerboseOutput_!=nullptr ) *pVerboseOutput_ << "\n";

	if( pVerboseOutput_!=nullptr ) *pVerboseOutput_ << "Loading sample from file " << inputSampleFilename << std::endl;
	CPPUNIT_ASSERT_NO_THROW( pSample_=l1menu::tools::loadSample( inputSampleFilename ) );
	CPPUNIT_ASSERT_MESSAGE( "Trying to load sample from disk gave a null pointer", pSample_!=nullptr );

	if( pVerboseOutput_!=nullptr ) *pVerboseOutput_ << "Loading menu from file " << inputMenuFilename << std::endl;
	CPPUNIT_ASSERT_NO_THROW( pTriggerMenu_=l1menu::tools::loadMenu( inputMenuFilename ) );
	CPPUNIT_ASSERT_MESSAGE( "TriggerMenu supplied needs at least one trigger for the tests", pTriggerMenu_->numberOfTriggers()>=1 );

	if( pVerboseOutput_!=nullptr ) *pVerboseOutput_ << "Loading rate plots from file " << ratePlotsFilename << std::endl;
	std::unique_ptr<TFile> pRatePlotsRootFile( TFile::Open( ratePlotsFilename.c_str() ) );
	CPPUNIT_ASSERT_NO_THROW( pRatePlotsFromDisk_.reset( new l1menu::MenuRatePlots( pRatePlotsRootFile.get() ) ) );
	CPPUNIT_ASSERT_MESSAGE( "Trying to load rate plots from disk gave a null pointer", pRatePlotsFromDisk_!=nullptr );
}

void TriggerRatePlotUnitTestSuite::testConstructingFromTH1()
{
	const l1menu::TriggerTable& triggerTable=l1menu::TriggerTable::instance();

	//
	// First I need to create a TriggerRatePlot the normal way, so that I can
	// test creating a different one from its TH1 histogram. I'll take the first
	// trigger from the menu.
	//
	l1menu::ITrigger& trigger=pTriggerMenu_->getTrigger(1);
	std::vector<std::string> thresholdNames=l1menu::tools::getThresholdNames( trigger );
	CPPUNIT_ASSERT( !thresholdNames.empty() );
	std::string& mainThreshold=thresholdNames.front();

	unsigned int numberOfBins=100;
	float lowerEdge=0;
	float upperEdge=100;
	try
	{
		numberOfBins=triggerTable.getSuggestedNumberOfBins( trigger.name(), mainThreshold );
		lowerEdge=triggerTable.getSuggestedLowerEdge( trigger.name(), mainThreshold );
		upperEdge=triggerTable.getSuggestedUpperEdge( trigger.name(), mainThreshold );
	}
	catch( std::exception& error) { /* Do nothing. If no binning suggestions have been set for this trigger use the defaults I set above. */ }

	l1menu::TriggerRatePlot ratePlot( trigger, "testRatePlot", numberOfBins, lowerEdge, upperEdge, mainThreshold, thresholdNames );
	if( pVerboseOutput_!=nullptr ) *pVerboseOutput_ << "Adding test sample to rate plot. This could take a while." << std::endl;
	ratePlot.addSample( *pSample_ );

	//
	// Now I have a functional TriggerRatePlot, I'll try and create another one just from the
	// underlying TH1 and see if the two objects are the same.
	//
	l1menu::TriggerRatePlot duplicateRatePlot( ratePlot.getPlot() );

	// Check that the threshold plotted against is the same
	CPPUNIT_ASSERT_EQUAL( ratePlot.versusParameter(), duplicateRatePlot.versusParameter() );

	std::vector<std::string> parametersAlreadyChecked;
	parametersAlreadyChecked.push_back( ratePlot.versusParameter() );

	// Check that parameters that are scaled along with the main threshold are all the same.
	const auto expectedScaledParameters=ratePlot.otherScaledParameters();
	const auto actualScaledParameters=duplicateRatePlot.otherScaledParameters();

	CPPUNIT_ASSERT_EQUAL( expectedScaledParameters.size(), actualScaledParameters.size() );

	for( auto iExpected=expectedScaledParameters.begin(), iActual=actualScaledParameters.begin();
			iExpected!=expectedScaledParameters.end() && iActual!=actualScaledParameters.end();
			++iExpected, ++iActual )
	{
		CPPUNIT_ASSERT_EQUAL( iExpected->first, iActual->first );
		CPPUNIT_ASSERT_EQUAL( iExpected->second, iActual->second );
		parametersAlreadyChecked.push_back( iExpected->first );
	}

	// Now compare the triggers
	const l1menu::ITriggerDescription& expectedTrigger=ratePlot.getTrigger();
	const l1menu::ITriggerDescription& actualTrigger=duplicateRatePlot.getTrigger();
	CPPUNIT_ASSERT_EQUAL( expectedTrigger.name(), actualTrigger.name() );
	CPPUNIT_ASSERT_EQUAL( expectedTrigger.version(), actualTrigger.version() );

	// I'll assume that if the trigger name and version are the same (checked above) that
	// I don't need to test that the parameter names are all the same.
	for( const auto& parameterName : expectedTrigger.parameterNames() )
	{
		// Make sure this parameter isn't one of the ones that are scaled
		if( std::find(parametersAlreadyChecked.begin(),parametersAlreadyChecked.end(),parameterName)==parametersAlreadyChecked.end() )
		{
			CPPUNIT_ASSERT_DOUBLES_EQUAL( expectedTrigger.parameter(parameterName), actualTrigger.parameter(parameterName), std::pow(10,-7) );
		}
	}
}

void TriggerRatePlotUnitTestSuite::testFindThresholdError()
{
	//
	// Not much I can test here except that it runs. To test anything else I need
	// to know what the answers should be. Maybe later I'll create a histogram
	// with purely fake data so that I can determine the results.
	//
	l1menu::TriggerRatePlot& testTriggerRatePlot_=pRatePlotsFromDisk_->triggerRatePlots().front();
	std::pair<float,float> thresholdErrors;
	// Brian added an extra parameter to findThresholdError, which I don't fully understand. It doesn't look as
	// though it's used other than for a printout. For now I'll just pass -1.
	CPPUNIT_ASSERT_NO_THROW( thresholdErrors=testTriggerRatePlot_.findThresholdError( 50, -1 ) );
	// Although I can make sure it throws an exception if you get
	// the underflow and overflow bins.
	CPPUNIT_ASSERT_THROW( thresholdErrors=testTriggerRatePlot_.findThresholdError( -1, -1 ), std::runtime_error );
	CPPUNIT_ASSERT_THROW( thresholdErrors=testTriggerRatePlot_.findThresholdError( 99999999, -1 ), std::runtime_error );
}

void TriggerRatePlotUnitTestSuite::testLowThresholdPlateau()
{
	//
	// Loop over all of the trigger rate plots that were loaded from disk, and try and get a
	// threshold for a rate that is higher than than the rate for zero threshold.
	//
	for( const auto& triggerRatePlot : pRatePlotsFromDisk_->triggerRatePlots() )
	{
		// Find out what the rate is at minimum threshold, and ask for a rate
		// less than that. The plot should return zero. There have been some
		// cases where it returns minus infinity.
		const float maximumRate=triggerRatePlot.getPlot()->GetBinContent(1);
		CPPUNIT_ASSERT_EQUAL_MESSAGE( std::string("Plot that fails has the title ")+triggerRatePlot.getPlot()->GetTitle(), 0.0f, triggerRatePlot.findThreshold( maximumRate*2 ) );
	}
}
