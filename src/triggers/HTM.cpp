#include "HTM.h"

#include <stdexcept>
#include <cmath>
#include "../implementation/RegisterTriggerMacro.h"
#include "l1menu/L1TriggerDPGEvent.h"
#include "UserCode/L1TriggerUpgrade/interface/L1AnalysisDataFormat.h"


namespace l1menu
{
	namespace triggers
	{

		/* The REGISTER_TRIGGER macro will make sure that the given trigger is registered in the
		 * l1menu::TriggerTable when the program starts. I also want to provide some suggested binning
		 * however. The REGISTER_TRIGGER_AND_CUSTOMISE macro does exactly the same but lets me pass
		 * a pointer to a function that will be called directly after the trigger has been registered
		 * at program startup. The function takes no parameters and returns void. In this case I'm
		 * giving it a lambda function.
		 */
		REGISTER_TRIGGER_AND_CUSTOMISE( HTM_v0,
			[]() // Use a lambda function to customise rather than creating a named function that never gets used again.
			{
				l1menu::TriggerTable& triggerTable=l1menu::TriggerTable::instance();
				HTM_v0 tempTriggerInstance;
				triggerTable.registerSuggestedBinning( tempTriggerInstance.name(), "threshold1", 150, 0, 300 );
			} // End of customisation lambda function
		) // End of REGISTER_TRIGGER_AND_CUSTOMISE macro call

		// No need to register suggested binning, it will use the binning for version 0
		REGISTER_TRIGGER( HTM_v1 )

	} // end of namespace triggers

} // end of namespace l1menu


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//---------------  Definitions below         ---------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


//
//  Version 1
//
l1menu::triggers::HTM_v1::HTM_v1()
	: regionCut_(4), ptCut_(0)
{
	// No operation besides the initialiser list
}

bool l1menu::triggers::HTM_v1::apply( const l1menu::L1TriggerDPGEvent& event ) const
{
	const L1Analysis::L1AnalysisDataFormat& analysisDataFormat=event.rawEvent();
	const bool* PhysicsBits=event.physicsBits();

	bool raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	double htmValueX=0.;
	double htmValueY=0.;

	// Calculate our own HT and HTM from the jets that survive the double jet removal.
	for( int i=0; i<analysisDataFormat.Njet; i++ )
	{
		if( analysisDataFormat.Bxjet.at( i )==0 && !analysisDataFormat.Taujet.at(i) )
		{
			if( analysisDataFormat.Etajet.at( i )>=regionCut_ and analysisDataFormat.Etajet.at( i )<=(21-regionCut_) )
			{
				if( analysisDataFormat.Etjet[i]>=ptCut_ )
				{
					//  Get the phi angle  towers are 0-17 (this is probably not real mapping but OK for just magnitude of HTM
					float phi=2*M_PI*(analysisDataFormat.Phijet.at( i )/18.);
					htmValueX+=std::cos( phi )*analysisDataFormat.Etjet.at( i );
					htmValueY+=std::sin( phi )*analysisDataFormat.Etjet.at( i );
				} // Jet is of minimum pT
			} //in proper eta range
		} //correct beam crossing
	} //loop over cleaned jets

	if( std::sqrt( htmValueX*htmValueX+htmValueY*htmValueY ) < threshold1_ ) return false;
	return true;
}

bool l1menu::triggers::HTM_v1::thresholdsAreCorrelated() const
{
	return false;
}

unsigned int l1menu::triggers::HTM_v1::version() const
{
	return 1;
}

const std::vector<std::string> l1menu::triggers::HTM_v1::parameterNames() const
{
	// First get the values from the base class, then add the extra entry
	std::vector<std::string> returnValue=HTM::parameterNames();
	returnValue.push_back("regionCut");
	returnValue.push_back("ptCut");
	return returnValue;
}

float& l1menu::triggers::HTM_v1::parameter( const std::string& parameterName )
{
	// Check if it's the parameter I've added, otherwise defer to the base class
	if( parameterName=="regionCut" ) return regionCut_;
	if( parameterName=="ptCut" ) return ptCut_;
	else return HTM::parameter(parameterName);
}

const float& l1menu::triggers::HTM_v1::parameter( const std::string& parameterName ) const
{
	// Check if it's the parameter I've added, otherwise defer to the base class
	if( parameterName=="regionCut" ) return regionCut_;
	if( parameterName=="ptCut" ) return ptCut_;
	else return HTM::parameter(parameterName);
}

//
//  Version 0
//
bool l1menu::triggers::HTM_v0::apply( const l1menu::L1TriggerDPGEvent& event ) const
{
	const L1Analysis::L1AnalysisDataFormat& analysisDataFormat=event.rawEvent();
	const bool* PhysicsBits=event.physicsBits();

	bool raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	float adc = analysisDataFormat.HTM ;
	float TheHTM = adc; // / 2. ;

	if (TheHTM < threshold1_) return false;
	return true;
}

bool l1menu::triggers::HTM_v0::thresholdsAreCorrelated() const
{
	return false;
}

unsigned int l1menu::triggers::HTM_v0::version() const
{
	return 0;
}

l1menu::triggers::HTM::HTM()
	: threshold1_(50)
{
	// No operation other than the initialiser list
}

const std::string l1menu::triggers::HTM::name() const
{
	return "L1_HTM";
}

const std::vector<std::string> l1menu::triggers::HTM::parameterNames() const
{
	std::vector<std::string> returnValue;
	returnValue.push_back("threshold1");
	return returnValue;
}

float& l1menu::triggers::HTM::parameter( const std::string& parameterName )
{
	if( parameterName=="threshold1" ) return threshold1_;
	else throw std::logic_error( "Not a valid parameter name" );
}

const float& l1menu::triggers::HTM::parameter( const std::string& parameterName ) const
{
	if( parameterName=="threshold1" ) return threshold1_;
	else throw std::logic_error( "Not a valid parameter name" );
}
