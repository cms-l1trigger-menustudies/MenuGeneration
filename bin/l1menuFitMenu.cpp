#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>

#include "l1menu/ISample.h"
#include "l1menu/IMenuRate.h"
#include "l1menu/MenuFitter.h"
#include "l1menu/TriggerTable.h"
#include "l1menu/ITrigger.h"
#include "l1menu/TriggerMenu.h"
#include "l1menu/MenuRatePlots.h"
#include "l1menu/MenuRateMuonScaling.h"
#include "l1menu/MenuRateOfflineScaling.h"
#include "l1menu/IL1MenuFile.h"
#include "l1menu/tools/fileIO.h"
#include "l1menu/tools/stringManipulation.h"
#include "l1menu/tools/CommandLineParser.h"
#include <TFile.h>
#include <TH1.h>

void printUsage( const std::string& executableName, std::ostream& output=std::cout )
{
	output << "Usage:" << "\n"
			<< "\t" << executableName << " --totalrate <total rate in kHz> [--rateplots <rateplot filename>] [--output <output filename>] [--format <CSV | OLD | XML>] <sample filename> <menu filename> <totalRate1> [totalRate2 [totalRate3 [...] ] ]" << "\n"
			<< "\t" << "\t" << "Tries to fit the supplied menu using the sample provided. The optional \"rateplots\" option" << "\n"
			<< "\t" << "\t" << "allows you to reuse a valid file created by l1menuCreateRatePlots which will significantly" << "\n"
			<< "\t" << "\t" << "speed up execution. If the option \"outputprefix\" is supplied the results will be saved to" << "\n"
			<< "\t" << "\t" << "a file with the supplied name. If this option isn't provided the output will be printed to" << "\n"
			<< "\t" << "\t" << "standard output." << "\n"
			<< "\t" << "\t" << "The 'format' option allows you specify what format the output will be in. XML (the default)" << "\n"
			<< "\t" << "\t" << "is required to do the scaling with l1menuScaleMenuRates." << "\n"
			<< "\n"
			<< "\t" << executableName << " --help" << "\n"
			<< "\t" << "\t" << "prints this help message"
			<< "\n"
			<< std::endl;
}

int main( int argc, char* argv[] )
{
	std::string sampleFilename;
	std::string menuFilename;
	std::string ratePlotsFilename; // Filename for rateplots. This is optional and can be empty.
	std::string outputFilename; // The filename to save output to. This is optional and can be empty, in which case no files are produced.
	l1menu::IL1MenuFile::FileFormat fileFormat=l1menu::IL1MenuFile::FileFormat::XML;
	float totalTriggerRatekHz; // The rate if every single event passed
	std::vector<float> totalRates;

	l1menu::tools::CommandLineParser commandLineParser;
	try
	{
		commandLineParser.addOption( "help", l1menu::tools::CommandLineParser::NoArgument );
		commandLineParser.addOption( "totalrate", l1menu::tools::CommandLineParser::RequiredArgument );
		commandLineParser.addOption( "rateplots", l1menu::tools::CommandLineParser::RequiredArgument );
		commandLineParser.addOption( "output", l1menu::tools::CommandLineParser::RequiredArgument );
		commandLineParser.addOption( "format", l1menu::tools::CommandLineParser::RequiredArgument );
		commandLineParser.parse( argc, argv );

		if( commandLineParser.optionHasBeenSet( "help" ) )
		{
			printUsage( commandLineParser.executableName() );
			return 0;
		}

		if( commandLineParser.nonOptionArguments().size()<3 ) throw std::runtime_error( "Not enough command line arguments" );
		if( commandLineParser.optionHasBeenSet( "rateplots" ) ) ratePlotsFilename=commandLineParser.optionArguments("rateplots").back();
		if( commandLineParser.optionHasBeenSet( "format" ) )
		{
			std::string formatString=commandLineParser.optionArguments("format").back();
			if( formatString=="XML" ) fileFormat=l1menu::IL1MenuFile::FileFormat::XML;
			else if( formatString=="OLD" ) fileFormat=l1menu::IL1MenuFile::FileFormat::OLD;
			else if( formatString=="CSV" ) fileFormat=l1menu::IL1MenuFile::FileFormat::CSV;
			else throw std::runtime_error( "format must be one of 'XML', 'OLD', or 'CSV'" );
		}
		if( commandLineParser.optionHasBeenSet( "output" ) )
		{
			outputFilename=commandLineParser.optionArguments("output").back();
			std::string outputSuffix;
			if( fileFormat==l1menu::IL1MenuFile::FileFormat::XML ) outputSuffix=".xml";
			else if( fileFormat==l1menu::IL1MenuFile::FileFormat::OLD ) outputSuffix=".txt";
			else if( fileFormat==l1menu::IL1MenuFile::FileFormat::CSV ) outputSuffix=".csv";
			// Make sure the user hasn't already added the file extension
			if( outputFilename.substr(outputFilename.size()-4)!=outputSuffix ) outputFilename+=outputSuffix;
		}

		//
		// Code to work out what to scale to
		//
		if( commandLineParser.optionHasBeenSet("totalrate") )
		{
			totalTriggerRatekHz=l1menu::tools::convertStringToFloat( commandLineParser.optionArguments("totalrate").back() );
		}
		else
		{
			std::cerr << "An option to set the scaling has not been set (with 'totalrate'). Scaling to a default of 100 pileup at 4.4E34." << std::endl;
			// I'll leave all these numbers in to make clear how the final result is calculated.
			const float scaleToKiloHz=1.0/1000.0;
			const float orbitsPerSecond=11246;
			const float bunchSpacing=25;
			float numberOfBunches;
			if( bunchSpacing==50 ) numberOfBunches=1380;
			else if( bunchSpacing==25 ) numberOfBunches=2760;
			else throw std::logic_error( "The number of bunches has not been programmed for the bunch spacing selected" );

			totalTriggerRatekHz=orbitsPerSecond*numberOfBunches*scaleToKiloHz;
		}


		const std::vector<std::string>& arguments=commandLineParser.nonOptionArguments();
		sampleFilename=arguments[0];
		menuFilename=arguments[1];
		for( size_t index=2; index<arguments.size(); ++index )
		{
			try
			{
				totalRates.push_back( l1menu::tools::convertStringToFloat(arguments[index]) );
			}
			catch( std::runtime_error& exception )
			{
				std::cerr << "Problem converting argument " << index << " into a float: " << exception.what() << std::endl;
				return -1;
			}
		}
	} // end of try block
	catch( std::exception& error )
	{
		std::cerr << "Error parsing the command line: " << error.what() << std::endl;
		printUsage( commandLineParser.executableName(), std::cerr );
		return -1;
	}


	try
	{
		std::cout << "Loading sample from the file " << sampleFilename << std::endl;
		std::unique_ptr<l1menu::ISample> pSample=l1menu::tools::loadSample( sampleFilename );
		std::cout << "Loading menu from file " << menuFilename << std::endl;
		std::unique_ptr<l1menu::TriggerMenu> pMenu=l1menu::tools::loadMenu( menuFilename );
		pSample->setEventRate( totalTriggerRatekHz );

		std::unique_ptr<l1menu::MenuFitter> pMenuFitter;
		if( ratePlotsFilename.empty() ) pMenuFitter.reset( new l1menu::MenuFitter(*pSample,*pMenu) );
		else
		{
			// If the user has specified a rateplots file on the command line try and
			// load it up and use it to create the MenuFitter. This will save a significant
			// amount of time during the MenuFitter constructor.
			std::unique_ptr<TFile> pRatePlotsRootFile( TFile::Open( ratePlotsFilename.c_str() ) );
			l1menu::MenuRatePlots ratePlots( pRatePlotsRootFile.get() );

			pMenuFitter.reset( new l1menu::MenuFitter( *pSample, *pMenu, ratePlots ) );
		}


		std::unique_ptr<l1menu::IL1MenuFile> pOutputL1MenuFile;
		if( !outputFilename.empty() ) pOutputL1MenuFile=l1menu::IL1MenuFile::getOutputFile( fileFormat, outputFilename );
		// Otherwise dump the information to standard output
		else
		{
			std::cout << "outputprefix not specified so dumping results to standard output" << "\n";
			pOutputL1MenuFile=l1menu::IL1MenuFile::getOutputFile( fileFormat, std::cout );
		}

		for( const auto& totalRate : totalRates )
		{
			std::cout << "Fitting menu for a rate of " << totalRate << "kHz..." << std::endl; std::cout.flush();
			try
			{
				std::shared_ptr<const l1menu::IMenuRate> pRates=pMenuFitter->fit( totalRate, 5.0 ); //totalRate*0.05
				pOutputL1MenuFile->add( *pRates );
				//l1menu::tools::dumpTriggerRates( *pOutputStream, *pRates, fileFormat );
				std::cout << "done." << std::endl;
			}
			catch( std::exception& error )
			{
				std::cerr << "An exception occurred while trying to fit for " << totalRate << "kHz: " << error.what() << "\n";
				std::cout << "--------------------    Start of fit log    --------------------" << "\n"
						<< pMenuFitter->debugLog() << "\n"
						<< "--------------------     End of fit log     --------------------" << "\n";
			}
		}

	}
	catch( std::exception& error )
	{
		std::cerr << "Exception caught: " << error.what() << std::endl;
		return -1;
	}

	return 0;
}
