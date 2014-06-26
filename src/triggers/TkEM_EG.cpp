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
		/** @brief Base class for all versions of the TkEM_EG trigger.
		 *
		 * Note that this class is abstract because it doesn't implement the "version"
		 * and "apply" methods. That's left up to the implementations of the different
		 * versions.
		 *
		 * It would have been nicer to implement this as a derived class of CrossTrigger,
		 * but there is a check on whether the electron and jet have the same value for
		 * eta and phi, and so the two triggers are not independent.
		 *
		 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
		 * @date 02/Jun/2013
		 */
		class TkEM_EG : public l1menu::ITrigger
		{
		public:
			TkEM_EG();

			virtual const std::string name() const;
			virtual const std::vector<std::string> parameterNames() const;
			virtual float& parameter( const std::string& parameterName );
			virtual const float& parameter( const std::string& parameterName ) const;
		protected:
			float leg1threshold1_;
			float leg2threshold1_;
			float leg1regionCut_;
			float leg2regionCut_;
			float tkIsoCut_;
		}; // end of the TkEM_EG base class


		/** @brief First version of the TkEM_EG trigger.
		 *
		 * @author probably Brian Winer
		 * @date sometime
		 */
		class TkEM_EG_v0 : public TkEM_EG
		{
		public:
			virtual unsigned int version() const;
			virtual bool apply( const l1menu::L1TriggerDPGEvent& event ) const;
			virtual bool thresholdsAreCorrelated() const;
		}; // end of version 0 class


		/* The REGISTER_TRIGGER macro will make sure that the given trigger is registered in the
		 * l1menu::TriggerTable when the program starts. I also want to provide some suggested binning
		 * however. The REGISTER_TRIGGER_AND_CUSTOMISE macro does exactly the same but lets me pass
		 * a pointer to a function that will be called directly after the trigger has been registered
		 * at program startup. The function takes no parameters and returns void. In this case I'm
		 * giving it a lambda function.
		 */
		REGISTER_TRIGGER_AND_CUSTOMISE( TkEM_EG_v0,
			[]() // Use a lambda function to customise rather than creating a named function that never gets used again.
			{
				l1menu::TriggerTable& triggerTable=l1menu::TriggerTable::instance();
				TkEM_EG_v0 tempTriggerInstance;
				triggerTable.registerSuggestedBinning( tempTriggerInstance.name(), "leg1threshold1", 100, 0, 100 );
				triggerTable.registerSuggestedBinning( tempTriggerInstance.name(), "leg2threshold1", 100, 0, 100 );
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

bool l1menu::triggers::TkEM_EG_v0::apply( const l1menu::L1TriggerDPGEvent& event ) const
{
	const L1Analysis::L1AnalysisDataFormat& analysisDataFormat=event.rawEvent();
	const bool* PhysicsBits=event.physicsBits();

	bool raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	bool ok = false;

	int Neg = analysisDataFormat.Nele;

	int Nem = analysisDataFormat.NTkem;
	for (int ue=0; ue < Nem; ue++) {
		int bx = analysisDataFormat.BxTkem[ue];
		if (bx != 0 ) continue;
		float eta = analysisDataFormat.EtaTkem[ue];
		if (eta < leg1regionCut_ || eta > 21.-leg1regionCut_) continue;  // eta = 5 - 16
		float pt = analysisDataFormat.EtTkem[ue];    // the rank of the electron
		float tkI= analysisDataFormat.tIsoTkem[ue];
		if (pt >= leg1threshold1_ && tkI<tkIsoCut_){

                        
			for (int uj=0; uj < Neg; uj++) {
				int bxj = analysisDataFormat.Bxel[uj];
				if (bxj != 0) continue;
				float rankj = analysisDataFormat.Etel[uj];
				float ptj = rankj; //CorrectedL1JetPtByGCTregions(analysisDataFormat.Etajet[uj],rank*4.,theL1JetCorrection);

				if (analysisDataFormat.Etael[uj] < leg2regionCut_ || analysisDataFormat.Etael[uj] > 21.-leg2regionCut_) continue;
				if (ptj >= leg2threshold1_ &&
					!(analysisDataFormat.Etael[uj]==analysisDataFormat.EtaTkem[ue] &&
					  analysisDataFormat.Phiel[uj]==analysisDataFormat.PhiTkem[ue]) ) ok = true;
			}

		} // if good TkEle
	}  // end loop over TkEle objects

	return ok;
}



bool l1menu::triggers::TkEM_EG_v0::thresholdsAreCorrelated() const
{
	return true;
}

unsigned int l1menu::triggers::TkEM_EG_v0::version() const
{
	return 0;
}

l1menu::triggers::TkEM_EG::TkEM_EG()
	: leg1threshold1_(20), leg2threshold1_(20), leg1regionCut_(4.5), leg2regionCut_(4.5), tkIsoCut_(999.)
{
	// No operation other than the initialiser list
}

const std::string l1menu::triggers::TkEM_EG::name() const
{
	return "L1_TkEM_EG";
}

const std::vector<std::string> l1menu::triggers::TkEM_EG::parameterNames() const
{
	std::vector<std::string> returnValue;
	returnValue.push_back("leg1threshold1");
	returnValue.push_back("leg1regionCut");
	returnValue.push_back("leg2threshold1");
	returnValue.push_back("leg2regionCut");
	returnValue.push_back("trkIsolCut");
	return returnValue;
}

float& l1menu::triggers::TkEM_EG::parameter( const std::string& parameterName )
{
	if( parameterName=="leg1threshold1" ) return leg1threshold1_;
	else if( parameterName=="leg1regionCut" ) return leg1regionCut_;
	else if( parameterName=="leg2threshold1" ) return leg2threshold1_;
	else if( parameterName=="leg2regionCut" ) return leg2regionCut_;
	else if( parameterName=="trkIsolCut" ) return tkIsoCut_;
	else throw std::logic_error( "Not a valid parameter name" );
}

const float& l1menu::triggers::TkEM_EG::parameter( const std::string& parameterName ) const
{
	if( parameterName=="leg1threshold1" ) return leg1threshold1_;
	else if( parameterName=="leg1regionCut" ) return leg1regionCut_;
	else if( parameterName=="leg2threshold1" ) return leg2threshold1_;
	else if( parameterName=="leg2regionCut" ) return leg2regionCut_;
	else if( parameterName=="trkIsolCut" ) return tkIsoCut_;
	else throw std::logic_error( "Not a valid parameter name" );
}
