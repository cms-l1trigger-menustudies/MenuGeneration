#include "../implementation/RegisterTriggerMacro.h"
#include "l1menu/L1TriggerDPGEvent.h"

#include <stdexcept>
#include "UserCode/L1TriggerUpgrade/interface/L1AnalysisDataFormat.h"

#include "l1menu/ITrigger.h"

#include <string>
#include <vector>

namespace l1menu
{
	namespace triggers
	{
		/** @brief Base class for all versions of the HTT trigger.
		 *
		 * Note that this class is abstract because it doesn't implement the "version"
		 * and "apply" methods. That's left up to the implementations of the different
		 * versions.
		 *
		 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
		 * @date 02/Jun/2013
		 */
		class HTT : public l1menu::ITrigger
		{
		public:
			HTT();

			virtual const std::string name() const;
			virtual const std::vector<std::string> parameterNames() const;
			virtual float& parameter( const std::string& parameterName );
			virtual const float& parameter( const std::string& parameterName ) const;
		protected:
			float threshold1_;
		}; // end of the HTT base class

		/** @brief First version of the HTT trigger.
		 *
		 * @author probably Brian Winer
		 * @date sometime
		 */
		class HTT_v0 : public HTT
		{
		public:
			virtual unsigned int version() const;
			virtual bool apply( const l1menu::L1TriggerDPGEvent& event ) const;
			virtual bool thresholdsAreCorrelated() const;
		}; // end of version 0 class

		/** @brief HTT trigger where the jets are looped over here rather than in FullSample.
		 *
		 * Allows eta cuts to be applied here, rather than having a hard coded eta cut in FullSample.
		 * Also added a minimum pT value if the jet is to be used to calculate the sum.
		 *
		 * @author Mark Grimes, but just copied Brian's code from FullSample.cpp
		 * @date 28/Mar/2014
		 */
		class HTT_v1 : public HTT
		{
		public:
			HTT_v1(); // Need a constructor to initiate regionCut_
			virtual unsigned int version() const;
			virtual bool apply( const l1menu::L1TriggerDPGEvent& event ) const;
			virtual bool thresholdsAreCorrelated() const;
			// Also need to override these because I've added a parameter
			virtual const std::vector<std::string> parameterNames() const;
			virtual float& parameter( const std::string& parameterName );
			virtual const float& parameter( const std::string& parameterName ) const;
		protected:
			float regionCut_;
			float ptCut_;
		}; // end of version 1 class

		/* The REGISTER_TRIGGER macro will make sure that the given trigger is registered in the
		 * l1menu::TriggerTable when the program starts. I also want to provide some suggested binning
		 * however. The REGISTER_TRIGGER_AND_CUSTOMISE macro does exactly the same but lets me pass
		 * a pointer to a function that will be called directly after the trigger has been registered
		 * at program startup. The function takes no parameters and returns void. In this case I'm
		 * giving it a lambda function.
		 */
		REGISTER_TRIGGER_AND_CUSTOMISE( HTT_v0,
			[]() // Use a lambda function to customise rather than creating a named function that never gets used again.
			{
				l1menu::TriggerTable& triggerTable=l1menu::TriggerTable::instance();
				HTT_v0 tempTriggerInstance;
				triggerTable.registerSuggestedBinning( tempTriggerInstance.name(), "threshold1", 200, 0, 1000 );
			} // End of customisation lambda function
		) // End of REGISTER_TRIGGER_AND_CUSTOMISE macro call

		// No need to register suggested binning, it will use the binning for version 0
		REGISTER_TRIGGER( HTT_v1 )


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
l1menu::triggers::HTT_v1::HTT_v1()
	: regionCut_(4), ptCut_(0)
{
	// No operation besides the initialiser list
}

bool l1menu::triggers::HTT_v1::apply( const l1menu::L1TriggerDPGEvent& event ) const
{
	const L1Analysis::L1AnalysisDataFormat& analysisDataFormat=event.rawEvent();
	const bool* PhysicsBits=event.physicsBits();

	bool raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	double httValue=0.;

	// Calculate our own HT and HTM from the jets that survive the double jet removal.
	for( int i=0; i<analysisDataFormat.Njet; i++ )
	{
		if( analysisDataFormat.Bxjet.at( i )==0 && !analysisDataFormat.Taujet.at(i) )
		{
			if( analysisDataFormat.Etajet.at( i )>=regionCut_ and analysisDataFormat.Etajet.at( i )<=(21-regionCut_) )
			{
				if( analysisDataFormat.Etjet[i]>=ptCut_ ) httValue+=analysisDataFormat.Etjet.at( i );
			} //in proper eta range
		} //correct beam crossing
	} //loop over cleaned jets

	if (httValue < threshold1_) return false;
	return true;
}

bool l1menu::triggers::HTT_v1::thresholdsAreCorrelated() const
{
	return false;
}

unsigned int l1menu::triggers::HTT_v1::version() const
{
	return 1;
}

const std::vector<std::string> l1menu::triggers::HTT_v1::parameterNames() const
{
	// First get the values from the base class, then add the extra entries
	std::vector<std::string> returnValue=HTT::parameterNames();
	returnValue.push_back("regionCut");
	returnValue.push_back("ptCut");
	return returnValue;
}

float& l1menu::triggers::HTT_v1::parameter( const std::string& parameterName )
{
	// Check if it's the parameter I've added, otherwise defer to the base class
	if( parameterName=="regionCut" ) return regionCut_;
	if( parameterName=="ptCut" ) return ptCut_;
	else return HTT::parameter(parameterName);
}

const float& l1menu::triggers::HTT_v1::parameter( const std::string& parameterName ) const
{
	// Check if it's the parameter I've added, otherwise defer to the base class
	if( parameterName=="regionCut" ) return regionCut_;
	if( parameterName=="ptCut" ) return ptCut_;
	else return HTT::parameter(parameterName);
}

//
//  Version 0
//
bool l1menu::triggers::HTT_v0::apply( const l1menu::L1TriggerDPGEvent& event ) const
{
	const L1Analysis::L1AnalysisDataFormat& analysisDataFormat=event.rawEvent();
	const bool* PhysicsBits=event.physicsBits();

	bool raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	float adc = analysisDataFormat.HTT ;
	float TheHTT = adc; // / 2. ;

	if (TheHTT < threshold1_) return false;
	return true;
}

bool l1menu::triggers::HTT_v0::thresholdsAreCorrelated() const
{
	return false;
}

unsigned int l1menu::triggers::HTT_v0::version() const
{
	return 0;
}

l1menu::triggers::HTT::HTT()
	: threshold1_(100)
{
	// No operation other than the initialiser list
}

const std::string l1menu::triggers::HTT::name() const
{
	return "L1_HTT";
}

const std::vector<std::string> l1menu::triggers::HTT::parameterNames() const
{
	std::vector<std::string> returnValue;
	returnValue.push_back("threshold1");
	return returnValue;
}

float& l1menu::triggers::HTT::parameter( const std::string& parameterName )
{
	if( parameterName=="threshold1" ) return threshold1_;
	else throw std::logic_error( "Not a valid parameter name" );
}

const float& l1menu::triggers::HTT::parameter( const std::string& parameterName ) const
{
	if( parameterName=="threshold1" ) return threshold1_;
	else throw std::logic_error( "Not a valid parameter name" );
}
