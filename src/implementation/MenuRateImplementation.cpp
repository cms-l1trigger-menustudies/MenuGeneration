#include "MenuRateImplementation.h"

#include <string>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include "l1menu/ITrigger.h"
#include "l1menu/ICachedTrigger.h"
#include "l1menu/ITriggerRate.h"
#include "l1menu/TriggerMenu.h"
#include "l1menu/ISample.h"
#include "l1menu/IEvent.h"
#include "l1menu/TriggerRatePlot.h"
#include "l1menu/MenuRatePlots.h"
#include "TriggerRateImplementation.h"
#include "TriggerDescriptionWithErrorsFromXML.h"
#include "l1menu/tools/XMLFile.h"
#include "l1menu/tools/XMLElement.h"
#include "l1menu/tools/fileIO.h"


void l1menu::implementation::MenuRateImplementation::commonConstruction( const l1menu::TriggerMenu& menu, const l1menu::ISample& sample )
{

	// The sum of event weights that pass each trigger
	std::vector<float> weightOfEventsPassed( menu.numberOfTriggers() );
	// The sume of weights squared that pass each trigger. Used to calculate the error.
	std::vector<float> weightSquaredOfEventsPassed( menu.numberOfTriggers() );

	// The number of events that only pass the given trigger
	std::vector<float> weightOfEventsPure( menu.numberOfTriggers() );
	std::vector<float> weightSquaredOfEventsPure( menu.numberOfTriggers() );

	float weightOfEventsPassingAnyTrigger=0;
	float weightSquaredOfEventsPassingAnyTrigger=0;
	float weightOfAllEvents=0;

	// Using cached triggers significantly increases speed for ReducedSample
	// because it cuts out expensive string comparisons when querying the trigger
	// parameters.
	std::vector< std::unique_ptr<l1menu::ICachedTrigger> > cachedTriggers;
	for( size_t triggerNumber=0; triggerNumber<menu.numberOfTriggers(); ++triggerNumber )
	{
		cachedTriggers.push_back( sample.createCachedTrigger( menu.getTrigger( triggerNumber ) ) );
	}

	size_t numberOfLastPassedTrigger=0; // This is just so I can work out the pure rate

	for( size_t eventNumber=0; eventNumber<sample.numberOfEvents(); ++eventNumber )
	{
		const l1menu::IEvent& event=sample.getEvent(eventNumber);
		float weight=event.weight();
		weightOfAllEvents+=weight;

		size_t numberOfTriggersPassed=0;

		for( size_t triggerNumber=0; triggerNumber<cachedTriggers.size(); ++triggerNumber )
		{
			if( cachedTriggers[triggerNumber]->apply(event) )
			{
				// If the event passes the trigger, increment the counters
				++numberOfTriggersPassed;
				weightOfEventsPassed[triggerNumber]+=weight;
				weightSquaredOfEventsPassed[triggerNumber]+=(weight*weight);
				numberOfLastPassedTrigger=triggerNumber; // If only one event passes, this is used to increment the pure counter
			}
		}

		// See if I should increment any of the pure or total counters
		if( numberOfTriggersPassed==1 )
		{
			weightOfEventsPure[numberOfLastPassedTrigger]+=weight;
			weightSquaredOfEventsPure[numberOfLastPassedTrigger]+=(weight*weight);
		}
		if( numberOfTriggersPassed>0 )
		{
			weightOfEventsPassingAnyTrigger+=weight;
			weightSquaredOfEventsPassingAnyTrigger+=(weight*weight);
		}
	}

	float scaling=sample.eventRate();

	for( size_t triggerNumber=0; triggerNumber<cachedTriggers.size(); ++triggerNumber )
	{
		float fraction=weightOfEventsPassed[triggerNumber]/weightOfAllEvents;
		float fractionError=std::sqrt(weightSquaredOfEventsPassed[triggerNumber])/weightOfAllEvents;
		float pureFraction=weightOfEventsPure[triggerNumber]/weightOfAllEvents;
		float pureFractionError=std::sqrt(weightSquaredOfEventsPure[triggerNumber])/weightOfAllEvents;
		triggerRates_.push_back( std::move(TriggerRateImplementation(menu.getTrigger(triggerNumber),fraction,fractionError,fraction*scaling,fractionError*scaling,pureFraction,pureFractionError,pureFraction*scaling,pureFractionError*scaling) ) );
		//triggerRates_.push_back( std::move(TriggerRateImplementation(menu.getTrigger(triggerNumber),weightOfEventsPassed[triggerNumber],weightSquaredOfEventsPassed[triggerNumber],weightOfEventsPure[triggerNumber],weightSquaredOfEventsPure[triggerNumber],*this)) );
	}

	//
	// Now I have everything I need to calculate all of the values required by the interface
	//
	totalFraction_=weightOfEventsPassingAnyTrigger/weightOfAllEvents;
	totalFractionError_=std::sqrt(weightSquaredOfEventsPassingAnyTrigger)/weightOfAllEvents;
	totalRate_=totalFraction_*scaling;
	totalRateError_=totalFractionError_*scaling;
}

l1menu::implementation::MenuRateImplementation::MenuRateImplementation( const l1menu::TriggerMenu& menu, const l1menu::ISample& sample )
{
	commonConstruction( menu, sample );
}

l1menu::implementation::MenuRateImplementation::MenuRateImplementation( const l1menu::TriggerMenu& menu, const l1menu::ISample& sample, const l1menu::MenuRatePlots& menuRatePlots )
{
	commonConstruction( menu, sample );
	// Loop over each of the trigger rates and try to set their threshold errors
	// from the information in the rate plots.
	for( auto& triggerRate : triggerRates_ )
	{
		// Search through the rate plots to find one which matches the trigger
		// for this rate.
		const l1menu::TriggerRatePlot* pRatePlot=nullptr;
		for( auto& triggerRatePlot : menuRatePlots.triggerRatePlots() )
		{
			if( triggerRatePlot.triggerMatches( triggerRate.trigger() ) )
			{
				pRatePlot=&triggerRatePlot;
				break;
			}
		}

		if( pRatePlot!=nullptr )
		{
			try
			{
				std::string mainThresholdName=pRatePlot->versusParameter();
				std::pair<float,float> primaryThresholdErrors=pRatePlot->findThresholdError( triggerRate.trigger().parameter(mainThresholdName), triggerRate.rate() );
				triggerRate.setParameterErrors( mainThresholdName, primaryThresholdErrors.first, primaryThresholdErrors.second );
			}
			catch( std::runtime_error& exception )
			{
				// If for some reason I couldn't get the treshold errors, then I can't set the
				// errors. I.e., just do nothing.

				// TODO - log this error somehow. Bad idea catching an exception and ignoring it, particularly
				// such a general one.
			}
		}
	}
}

l1menu::implementation::MenuRateImplementation::MenuRateImplementation( const l1menu::tools::XMLElement& xmlDescription )
{
	std::vector<l1menu::tools::XMLElement> parameterElements=xmlDescription.getChildren("totalFraction");
	if( parameterElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because the element did not have one and only one 'totalFraction' child." );
	totalFraction_=parameterElements.front().getFloatValue();

	parameterElements=xmlDescription.getChildren("totalFractionError");
	if( parameterElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because the element did not have one and only one 'totalFractionError' child." );
	totalFractionError_=parameterElements.front().getFloatValue();

	parameterElements=xmlDescription.getChildren("totalRate");
	if( parameterElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because the element did not have one and only one 'totalRate' child." );
	totalRate_=parameterElements.front().getFloatValue();

	parameterElements=xmlDescription.getChildren("totalRateError");
	if( parameterElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because the element did not have one and only one 'totalRateError' child." );
	totalRateError_=parameterElements.front().getFloatValue();

	parameterElements=xmlDescription.getChildren("TriggerRate");
	for( const auto& triggerRateElement : parameterElements )
	{
		std::vector<l1menu::tools::XMLElement> triggerElements=triggerRateElement.getChildren("Trigger");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'Trigger' child." );
		l1menu::implementation::TriggerDescriptionWithErrorsFromXML triggerDescription( triggerElements.front() );

		triggerElements=triggerRateElement.getChildren("fraction");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'fraction' child." );
		float fraction=triggerElements.front().getFloatValue();

		triggerElements=triggerRateElement.getChildren("fractionError");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'fractionError' child." );
		float fractionError=triggerElements.front().getFloatValue();

		triggerElements=triggerRateElement.getChildren("rate");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'rate' child." );
		float rate=triggerElements.front().getFloatValue();

		triggerElements=triggerRateElement.getChildren("rateError");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'rateError' child." );
		float rateError=triggerElements.front().getFloatValue();

		triggerElements=triggerRateElement.getChildren("pureFraction");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'pureFraction' child." );
		float pureFraction=triggerElements.front().getFloatValue();

		triggerElements=triggerRateElement.getChildren("pureFractionError");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'pureFractionError' child." );
		float pureFractionError=triggerElements.front().getFloatValue();

		triggerElements=triggerRateElement.getChildren("pureRate");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'pureRate' child." );
		float pureRate=triggerElements.front().getFloatValue();

		triggerElements=triggerRateElement.getChildren("pureRateError");
		if( triggerElements.size()!=1 ) throw std::runtime_error( "Failed to create IMenuRate from XML because one of the TriggerRate elements did not have one and only one 'pureRateError' child." );
		float pureRateError=triggerElements.front().getFloatValue();

		triggerRates_.push_back( std::move(TriggerRateImplementation(triggerDescription,fraction,fractionError,rate,rateError,pureFraction,pureFractionError,pureRate,pureRateError) ) );
	}
}

void l1menu::implementation::MenuRateImplementation::setTotalFraction( float totalFraction )
{
	totalFraction_=totalFraction;
}

void l1menu::implementation::MenuRateImplementation::setTotalFractionError( float totalFractionError )
{
	totalFractionError_=totalFractionError;
}

void l1menu::implementation::MenuRateImplementation::setTotalRate( float totalRate )
{
	totalRate_=totalRate;
}

void l1menu::implementation::MenuRateImplementation::setTotalRateError( float totalRateError )
{
	totalRateError_=totalRateError;
}

void l1menu::implementation::MenuRateImplementation::addTriggerRate( l1menu::implementation::TriggerRateImplementation&& triggerRate )
{
	triggerRates_.push_back( std::move(triggerRate) );
}

l1menu::implementation::MenuRateImplementation::MenuRateImplementation()
{
	// No operation.
}

float l1menu::implementation::MenuRateImplementation::totalFraction() const
{
	return totalFraction_;
}

float l1menu::implementation::MenuRateImplementation::totalFractionError() const
{
	return totalFractionError_;
}

float l1menu::implementation::MenuRateImplementation::totalRate() const
{
	return totalRate_;
}

float l1menu::implementation::MenuRateImplementation::totalRateError() const
{
	return totalRateError_;
}

const std::vector<const l1menu::ITriggerRate*>& l1menu::implementation::MenuRateImplementation::triggerRates() const
{
	// If the sizes are the same I'll assume nothing has changed and the references
	// are still valid. I don't expect this method to be called until the triggerRates_
	// vector is complete anyway.
	if( triggerRates_.size()!=baseClassPointers_.size() )
	{
		baseClassPointers_.clear();
		for( const auto& triggerRate : triggerRates_ )
		{
			baseClassPointers_.push_back( &triggerRate );
		}
	}

	return baseClassPointers_;
}
