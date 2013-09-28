#include "MenuRateImplementation.h"

#include <string>
#include <utility>
#include <cmath>
#include "l1menu/ITrigger.h"
#include "l1menu/ICachedTrigger.h"
#include "l1menu/ITriggerRate.h"
#include "l1menu/TriggerMenu.h"
#include "l1menu/ISample.h"
#include "l1menu/IEvent.h"
#include "TriggerRateImplementation.h"

l1menu::implementation::MenuRateImplementation::MenuRateImplementation( const l1menu::TriggerMenu& menu, const l1menu::ISample& sample )
{

	// The sum of event weights that pass each trigger
	std::vector<float> weightOfEventsPassed( menu.numberOfTriggers() );
	// The sume of weights squared that pass each trigger. Used to calculate the error.
	std::vector<float> weightSquaredOfEventsPassed( menu.numberOfTriggers() );

	// The number of events that only pass the given trigger
	std::vector<float> weightOfEventsPure( menu.numberOfTriggers() );
	std::vector<float> weightSquaredOfEventsPure( menu.numberOfTriggers() );

	weightOfEventsPassingAnyTrigger_=0;
	weightSquaredOfEventsPassingAnyTrigger_=0;
	weightOfAllEvents_=0;

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
		weightOfAllEvents_+=weight;

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
			weightOfEventsPassingAnyTrigger_+=weight;
			weightSquaredOfEventsPassingAnyTrigger_+=(weight*weight);
		}
	}

	scaling_=sample.eventRate();

	for( size_t triggerNumber=0; triggerNumber<cachedTriggers.size(); ++triggerNumber )
	{
		triggerRates_.push_back( std::move(TriggerRateImplementation(menu.getTrigger(triggerNumber),weightOfEventsPassed[triggerNumber],weightSquaredOfEventsPassed[triggerNumber],weightOfEventsPure[triggerNumber],weightSquaredOfEventsPure[triggerNumber],*this)) );
	}

}

float l1menu::implementation::MenuRateImplementation::weightOfAllEvents() const
{
	return weightOfAllEvents_;
}

float l1menu::implementation::MenuRateImplementation::weightOfAllEventsPassingAnyTrigger() const
{
	return weightOfEventsPassingAnyTrigger_;
}

float l1menu::implementation::MenuRateImplementation::weightSquaredOfAllEventsPassingAnyTrigger() const
{
	return weightSquaredOfEventsPassingAnyTrigger_;
}

float l1menu::implementation::MenuRateImplementation::scaling() const
{
	return scaling_;
}

float l1menu::implementation::MenuRateImplementation::totalFraction() const
{
	return weightOfEventsPassingAnyTrigger_/weightOfAllEvents_;
}

float l1menu::implementation::MenuRateImplementation::totalFractionError() const
{
	return std::sqrt(weightSquaredOfEventsPassingAnyTrigger_)/weightOfAllEvents_;
}

float l1menu::implementation::MenuRateImplementation::totalRate() const
{
	return totalFraction()*scaling_;
}

float l1menu::implementation::MenuRateImplementation::totalRateError() const
{
	return totalFractionError()*scaling_;
}

const std::vector<const l1menu::ITriggerRate*>& l1menu::implementation::MenuRateImplementation::triggerRates() const
{
	// If the sizes are the same I'll assume nothing has changed and the references
	// are still valid. I don't expect this method to be called until the triggerRates_
	// vector is complete anyway.
	if( triggerRates_.size()!=baseClassReferences_.size() )
	{
		baseClassReferences_.clear();
		for( const auto& triggerRate : triggerRates_ )
		{
			baseClassReferences_.push_back( &triggerRate );
		}
	}

	return baseClassReferences_;
}
