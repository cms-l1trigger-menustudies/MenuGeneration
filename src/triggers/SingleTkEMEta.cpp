#include "SingleTkEMEta.h"

#include <stdexcept>
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
		REGISTER_TRIGGER_AND_CUSTOMISE( SingleTkEMEta_v0,
			[]() // Use a lambda function to customise rather than creating a named function that never gets used again.
			{
				l1menu::TriggerTable& triggerTable=l1menu::TriggerTable::instance();
				SingleTkEMEta_v0 tempTriggerInstance;
				triggerTable.registerSuggestedBinning( tempTriggerInstance.name(), "threshold1", 100, 0, 100 );
			} // End of customisation lambda function
		) // End of REGISTER_TRIGGER_AND_CUSTOMISE macro call


	} // end of namespace triggers

} // end of namespace l1menu


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//---------------  Definitions below         ---------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


bool l1menu::triggers::SingleTkEMEta_v0::apply( const l1menu::L1TriggerDPGEvent& event ) const
{
	const L1Analysis::L1AnalysisDataFormat& analysisDataFormat=event.rawEvent();
	const bool* PhysicsBits=event.physicsBits();

	bool raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	bool ok=false;
	int NTkem = analysisDataFormat.NTkem;
	for (int ue=0; ue < NTkem; ue++) {
		int bx = analysisDataFormat.BxTkem[ue];
		if (bx != 0) continue;
		float eta = analysisDataFormat.EtaTkem[ue];
		if (eta < regionCut_ || eta > 21.-regionCut_) continue;  // eta = 5 - 16
		float pt = analysisDataFormat.EtTkem[ue];    // the rank of the electron
		float tkI= analysisDataFormat.tIsoTkem[ue];
		if (pt >= threshold1_ && tkI<tkIsoCut_) ok = true;
	}  // end loop over EM objects

	return ok;
}

bool l1menu::triggers::SingleTkEMEta_v0::thresholdsAreCorrelated() const
{
	return false;
}

unsigned int l1menu::triggers::SingleTkEMEta_v0::version() const
{
	return 0;
}

l1menu::triggers::SingleTkEMEta::SingleTkEMEta()
	: threshold1_(20), regionCut_(4.5), tkIsoCut_(999.)
{
	// No operation other than the initialiser list
}

const std::string l1menu::triggers::SingleTkEMEta::name() const
{
	return "L1_SingleTkEM";
}

const std::vector<std::string> l1menu::triggers::SingleTkEMEta::parameterNames() const
{
	std::vector<std::string> returnValue;
	returnValue.push_back("threshold1");
	returnValue.push_back("regionCut");
	returnValue.push_back("trkIsolCut");
	return returnValue;
}

float& l1menu::triggers::SingleTkEMEta::parameter( const std::string& parameterName )
{
	if( parameterName=="threshold1" ) return threshold1_;
	else if( parameterName=="regionCut" ) return regionCut_;
	else if( parameterName=="trkIsolCut" ) return tkIsoCut_;
	else throw std::logic_error( "Not a valid parameter name" );
}

const float& l1menu::triggers::SingleTkEMEta::parameter( const std::string& parameterName ) const
{
	if( parameterName=="threshold1" ) return threshold1_;
	else if( parameterName=="regionCut" ) return regionCut_;
	else if( parameterName=="trkIsolCut" ) return tkIsoCut_;
	else throw std::logic_error( "Not a valid parameter name" );
}
