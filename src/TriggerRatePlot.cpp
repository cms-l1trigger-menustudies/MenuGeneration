#include "l1menu/TriggerRatePlot.h"

#include "l1menu/ITrigger.h"
#include "l1menu/ICachedTrigger.h"
#include "l1menu/L1TriggerDPGEvent.h"
#include "l1menu/IEvent.h"
#include "l1menu/ISample.h"
#include "l1menu/TriggerTable.h"
#include "l1menu/tools/miscellaneous.h"
#include "l1menu/tools/stringManipulation.h"
#include <TH1F.h>
#include <sstream>
#include <algorithm>
#include <stdexcept>

l1menu::TriggerRatePlot::TriggerRatePlot( const l1menu::ITriggerDescription& trigger, std::unique_ptr<TH1> pHistogram, const std::string& versusParameter, const std::vector<std::string> scaledParameters )
	: pHistogram_( std::move(pHistogram) ), versusParameter_(versusParameter), histogramOwnedByMe_(true)
{
	initiate( trigger, scaledParameters );
}

l1menu::TriggerRatePlot::TriggerRatePlot( const l1menu::ITriggerDescription& trigger, const std::string& name, size_t numberOfBins, float lowEdge, float highEdge, const std::string& versusParameter, const std::vector<std::string> scaledParameters )
	: pHistogram_( new TH1F( name.c_str(), "This title gets changed later", numberOfBins, lowEdge, highEdge ) ), versusParameter_(versusParameter), histogramOwnedByMe_(true)
{
	pHistogram_->SetDirectory( NULL ); // Hold in memory only
	initiate( trigger, scaledParameters );
}

l1menu::TriggerRatePlot::TriggerRatePlot( const TH1* pPreExisitingHistogram )
{
	// All of the information about the trigger is stored in the title, so examine
	// that to see what the trigger is, the version, and all of the parameters.
	// Format is "<name> rate versus <thresholdName> [v<version>,<param1>=<value>,<param2>=<value>,...]"
	// so first split by whitespace.
	std::vector<std::string> titleSplitByWhitespace=l1menu::tools::splitByWhitespace( pPreExisitingHistogram->GetTitle() );

	if( titleSplitByWhitespace.size()!=5 ) throw std::runtime_error( "TriggerRatePlot cannot be loaded from the supplied histogram because the title is not as expected" );
	std::string triggerName=titleSplitByWhitespace[0];
	versusParameter_=titleSplitByWhitespace[3];

	// Check that the final element begins with "[" and ends with "]"
	if( titleSplitByWhitespace[4].front()!='[' || titleSplitByWhitespace[4].back()!=']' ) throw std::runtime_error( "TriggerRatePlot cannot be loaded from the supplied histogram because the title is not as expected" );
	// If so, chop off those two parts
	std::string parameterList=titleSplitByWhitespace[4].substr( 1, titleSplitByWhitespace[4].size()-2 );
	// and split by commas
	std::vector<std::string> parameterListSplitByCommas=l1menu::tools::splitByDelimeters( parameterList, "," );
	// Need to at least have the version listed
	if( parameterListSplitByCommas.empty() ) throw std::runtime_error( "TriggerRatePlot cannot be loaded from the supplied histogram because the title is not as expected" );

	// Figure out the version number of the trigger. It should be prefixed with a "v".
	if( parameterListSplitByCommas[0].front()!='v' ) throw std::runtime_error( "TriggerRatePlot cannot be loaded from the supplied histogram because the title is not as expected" );
	std::string versionAsString=parameterListSplitByCommas[0].substr( 1 );
	unsigned int versionNumber=l1menu::tools::convertStringToFloat( versionAsString );

	// Now I have the trigger name and version number, I can create a copy from the TriggerTable.
	std::unique_ptr<l1menu::ITrigger> pTemporaryTriggerCopy=l1menu::TriggerTable::instance().getTrigger( triggerName, versionNumber );
	if( pTemporaryTriggerCopy==nullptr ) throw std::runtime_error( "TriggerRatePlot cannot be loaded from the supplied histogram because trigger is not in the TriggerTable" );

	// Set the main threshold to 1 just as an arbitrary starting point. I might need to scale thresholds
	// off that in a minute.
	pTemporaryTriggerCopy->parameter( versusParameter_ )=1;

	std::vector<std::string> scaledParameters; // This will build into a list of anything that should scale with the main threshold

	// Then loop through the other parameters and set them to whatever they're listed as.
	// Start at 1 because I've already done the version.
	for( size_t index=1; index<parameterListSplitByCommas.size(); ++index )
	{
		// Split by equals sign
		std::vector<std::string> nameAndValue=l1menu::tools::splitByDelimeters( parameterListSplitByCommas[index], "=" );
		if( nameAndValue.size()!=2 ) throw std::runtime_error( "TriggerRatePlot cannot be loaded from the supplied histogram because the title is not as expected" );

		// The value could be scaled against the main threshold, so I need to check for that.
		std::string valueAsString=nameAndValue[1];
		if( valueAsString.size()>2 )
		{
			if( valueAsString[0]=='x' && valueAsString[1]=='*' )
			{
				// Parameter is scaled with the main threshold. Since I've set that to 1
				// I can just use the number after the "x*" as is, so just chop off the
				// first two characters.
				valueAsString=valueAsString.substr( 2 );
				scaledParameters.push_back( nameAndValue[0] );
			}
		}

		pTemporaryTriggerCopy->parameter( nameAndValue[0] )=l1menu::tools::convertStringToFloat( valueAsString );
	}

	// Everything seems to have worked okay, so I can take a copy of the histogram
	// and delegate to the normal initialisation routine.
	pHistogram_.reset( static_cast<TH1*>( pPreExisitingHistogram->Clone() ) );
	pHistogram_->SetDirectory( NULL ); // Hold in memory only
	histogramOwnedByMe_=true;
	initiate( *pTemporaryTriggerCopy, scaledParameters );
}

l1menu::TriggerRatePlot::TriggerRatePlot( const l1menu::TriggerRatePlot& otherTriggerRatePlot )
	: pHistogram_( static_cast<TH1*>(otherTriggerRatePlot.pHistogram_->Clone()) ),
	  versusParameter_( otherTriggerRatePlot.versusParameter_ ),
	  histogramOwnedByMe_(true)
{
	// Make sure the cloned histogram doesn't think it belongs to a TDirectory
	pHistogram_->SetDirectory( nullptr );
	initiate( *otherTriggerRatePlot.pTrigger_, otherTriggerRatePlot.otherScaledParameters_ );
}

l1menu::TriggerRatePlot::TriggerRatePlot( l1menu::TriggerRatePlot&& otherTriggerRatePlot ) noexcept
	: pTrigger_( std::move(otherTriggerRatePlot.pTrigger_) ),
	  pHistogram_( std::move(otherTriggerRatePlot.pHistogram_) ),
	  versusParameter_( std::move(otherTriggerRatePlot.versusParameter_) ),
	  pParameter_(&pTrigger_->parameter(versusParameter_)),
	  otherScaledParameters_( std::move(otherTriggerRatePlot.otherScaledParameters_) ),
	  otherParameterScalings_( std::move(otherTriggerRatePlot.otherParameterScalings_) ),
	  histogramOwnedByMe_(otherTriggerRatePlot.histogramOwnedByMe_)
{
	// No operation besides the initaliser list
}

l1menu::TriggerRatePlot& l1menu::TriggerRatePlot::operator=( l1menu::TriggerRatePlot&& otherTriggerRatePlot ) noexcept
{
	// Free up whatever was there before
	pTrigger_.reset();
	if( histogramOwnedByMe_ ) pHistogram_.reset();
	else pHistogram_.release();

	// Then copy from the other object
	pTrigger_=std::move(otherTriggerRatePlot.pTrigger_);
	pHistogram_=std::move(otherTriggerRatePlot.pHistogram_);
	versusParameter_=std::move(otherTriggerRatePlot.versusParameter_);
	pParameter_=&pTrigger_->parameter(versusParameter_);
	otherScaledParameters_=std::move(otherTriggerRatePlot.otherScaledParameters_);
	otherParameterScalings_=std::move(otherTriggerRatePlot.otherParameterScalings_);
	histogramOwnedByMe_=otherTriggerRatePlot.histogramOwnedByMe_;

	return *this;
}

void l1menu::TriggerRatePlot::initiate( const l1menu::ITriggerDescription& trigger, const std::vector<std::string>& scaledParameters )
{
	// Before making any histograms make sure errors are done properly
	TH1::SetDefaultSumw2();

	// Take a copy of the trigger
	l1menu::TriggerTable& table=l1menu::TriggerTable::instance();
	pTrigger_=table.copyTrigger( trigger );

	// Make sure the versusParameter_ supplied is valid. If it's not then this call will
	// throw an exception. Take a pointer to the parameter so I don't need to keep performing
	// expensive string comparisons.
	pParameter_=&pTrigger_->parameter(versusParameter_);

	// If any parameters have been requested to be scaled along with the versusParameter, figure
	// out what the scaling should be and take a note of pointers.
	otherScaledParameters_.clear();
	otherParameterScalings_.clear();
	for( const auto& parameterToScale : scaledParameters )
	{
		if( parameterToScale!=versusParameter_ )
		{
			otherScaledParameters_.push_back( parameterToScale ); // Keep a record of the name
			otherParameterScalings_.push_back( std::make_pair( &pTrigger_->parameter(parameterToScale), pTrigger_->parameter(parameterToScale)/(*pParameter_) ) );
		}
	}
	// I want to make a note of the other parameters set for the trigger. As far as I know TH1
	// has no way of adding annotations so I'll tag it on the end of the title.
	std::stringstream description;
	description << pTrigger_->name() << " rate versus " << versusParameter_;

	// Loop over all of the parameters and add their values to the description
	std::vector<std::string> parameterNames=pTrigger_->parameterNames();
	description << " [v" << pTrigger_->version();
	if( parameterNames.size()>1 ) description << ",";
	for( std::vector<std::string>::const_iterator iName=parameterNames.begin(); iName!=parameterNames.end(); ++iName )
	{
		if( *iName==versusParameter_ ) continue; // Don't bother adding the parameter I'm plotting against

		// First check to see if this is one of the parameters that are being scaled
		if( std::find(scaledParameters.begin(),scaledParameters.end(),*iName)==scaledParameters.end() )
		{
			// This parameter isn't being scaled, so write the absoulte value in the description
			description << *iName << "=" << pTrigger_->parameter(*iName);
		}
		else
		{
			// This parameter is being scaled, so write what the scaling is in the description
			description << *iName << "=x*" << pTrigger_->parameter(*iName)/(*pParameter_);
		}

		if( iName+1!=parameterNames.end() ) description << ","; // Add delimeter between parameter names
	}
	description << "]";

	pHistogram_->SetTitle( description.str().c_str() );
}

l1menu::TriggerRatePlot::~TriggerRatePlot()
{
	// If the flag has been set telling me that this instance doesn't own the histogram,
	// then I need to tell the unique_ptr not to delete it.
	if( !histogramOwnedByMe_ )
	{
		pHistogram_.release();
	}
}

void l1menu::TriggerRatePlot::addEvent( const l1menu::IEvent& event )
{
	const l1menu::ISample& sample=event.sample();
	float weightPerEvent=sample.eventRate()/sample.sumOfWeights();

	// For some implementations of ISample, it is significantly faster to
	// create ICachedTriggers and then loop over those. The addEvent overload
	// I'm about to delegate to loops over the histogram bins, so even for
	// one event this trigger can be called multiple times.
	std::unique_ptr<l1menu::ICachedTrigger> pCachedTrigger=sample.createCachedTrigger( *pTrigger_ );

	addEvent( event, pCachedTrigger, weightPerEvent );
}

void l1menu::TriggerRatePlot::addSample( const l1menu::ISample& sample )
{
	float weightPerEvent=sample.eventRate()/sample.sumOfWeights();

	// Create a cached trigger, which depending on the concrete type of the ISample
	// may or may not significantly increase the speed at which this next loop happens.
	std::unique_ptr<l1menu::ICachedTrigger> pCachedTrigger=sample.createCachedTrigger( *pTrigger_ );

	for( size_t eventNumber=0; eventNumber<sample.numberOfEvents(); ++eventNumber )
	{
		addEvent( sample.getEvent(eventNumber), pCachedTrigger, weightPerEvent );
	} // end of loop over events

}

void l1menu::TriggerRatePlot::addEvent( const l1menu::IEvent& event, const std::unique_ptr<l1menu::ICachedTrigger>& pCachedTrigger, float weightPerEvent )
{
	//
	// Use bisection to find the bin that passes the trigger and the one
	// immediately after it that fails.
	//
	size_t lowBin=1;
	size_t highBin=pHistogram_->GetNbinsX();

	//
	// First need to perform a check that the first bin passes. If it doesn't then
	// the histogram doesn't need filling at all and I can return.
	//
	(*pParameter_)=pHistogram_->GetBinLowEdge(lowBin);
	// Scale accordingly any other parameters that should be scaled. Remember that
	// in parameterScalingPair, 'first' is a pointer to the threshold to be changed
	// and 'second' is the ratio of the first threshold it should be.
	for( const auto& parameterScalingPair : otherParameterScalings_ ) *(parameterScalingPair.first)=parameterScalingPair.second*(*pParameter_);
	if( !pCachedTrigger->apply(event) ) return;

	//
	// Also check the highest bin. If that passes then I just fill every bin,
	// otherwise I need to find the point at which the trigger fails.
	//
	(*pParameter_)=pHistogram_->GetBinLowEdge(highBin);
	for( const auto& parameterScalingPair : otherParameterScalings_ ) *(parameterScalingPair.first)=parameterScalingPair.second*(*pParameter_);

	if( pCachedTrigger->apply(event) ) lowBin=highBin;
	else
	{
		while( highBin-lowBin>1 ) // Loop until I find two bins next to each other
		{
			size_t middleBin=(highBin+lowBin)/2;

			(*pParameter_)=pHistogram_->GetBinLowEdge(middleBin);
			for( const auto& parameterScalingPair : otherParameterScalings_ ) *(parameterScalingPair.first)=parameterScalingPair.second*(*pParameter_);

			if( pCachedTrigger->apply(event) ) lowBin=middleBin;
			else highBin=middleBin;
		}
	}

	//
	// Now I know which bins need filling, loop over them and fill.
	//
	for( size_t binNumber=1; binNumber<=lowBin; ++binNumber )
	{
		pHistogram_->Fill( pHistogram_->GetBinCenter(binNumber), event.weight()*weightPerEvent );
	}

}

const l1menu::ITriggerDescription& l1menu::TriggerRatePlot::getTrigger() const
{
	return *pTrigger_;
}

const std::string& l1menu::TriggerRatePlot::versusParameter() const
{
	return versusParameter_;
}

std::vector<std::pair<std::string,float> > l1menu::TriggerRatePlot::otherScaledParameters() const
{
	std::vector<std::pair<std::string,float> > returnValue;

	auto iParamName=otherScaledParameters_.begin();
	auto iParamScaling=otherParameterScalings_.begin();

	for( ; iParamName!=otherScaledParameters_.end() && iParamScaling!=otherParameterScalings_.end(); ++iParamName, ++iParamScaling )
	{
		returnValue.push_back( std::make_pair(*iParamName,iParamScaling->second) );
	}

	return returnValue;
}

bool l1menu::TriggerRatePlot::triggerMatches( const l1menu::ITriggerDescription& trigger, bool matchTriggerVersion ) const
{
	if( trigger.name()!=pTrigger_->name() ) return false;
	if( matchTriggerVersion && trigger.version()!=pTrigger_->version() ) return false;

	// If control gets this far, then it's the same trigger but the parameters could still be different.
	// First check all of the parameters that aren't the versus parameter or a parameter that scales with
	// the versus parameter.
	for( const auto& parameterName : pTrigger_->parameterNames() )
	{
		// Skip over versus parameter and scaled parameters
		if( parameterName==versusParameter_ ) continue;
		if( std::find( otherScaledParameters_.begin(), otherScaledParameters_.end(), parameterName )!=otherScaledParameters_.end() ) continue;

		if( pTrigger_->parameter(parameterName)!=trigger.parameter(parameterName) ) return false;
	}

	// Now need to check that any scaled parameters are scaled the same
	float mainThreshold=trigger.parameter(versusParameter_);
	std::vector<std::string>::const_iterator iScaledParameterName=otherScaledParameters_.begin();
	// I already know the otherScaledParameters_ and otherParameterScalings_ vectors
	// have equal size.
	for( const auto& paramScalingPair : otherParameterScalings_ )
	{
		float thisTriggerParameterScaling=paramScalingPair.second;
		float parameterScaling=trigger.parameter(*iScaledParameterName)/mainThreshold;
		// Make sure they're equal. Float equality is a bit dodgy so I'll just check the
		// difference is less than an arbitrarily small amount.
		if( std::fabs(parameterScaling-thisTriggerParameterScaling) > std::pow( 10, -4 ) ) return false;

		++iScaledParameterName;
	}

	// If control has gotten this far, then all of the tests have passed.
	return true;
}

float l1menu::TriggerRatePlot::findThreshold( float targetRate ) const
{

        bool prt=false;
	if(strcmp(pHistogram_->GetName(),"L1_SingleMu_v_threshold1")==0 ) prt=true;
	if(prt) printf("FINDThreshold: Found trigger to debug %s  with target rate %f \n",pHistogram_->GetName(), targetRate);
	//
	// Loop over all of the bins in the plot and find the first one
	// that is less than the requested rate.
	//
	int binNumber;
	for( binNumber=1; binNumber<pHistogram_->GetNbinsX(); ++binNumber )
	{
		if( pHistogram_->GetBinContent(binNumber)<targetRate ) break;
	}

	// If the search got to the end, put it one back so that when I fill
	// the binNumberForLinearFit vector below it includes the last four
	// bins for the linear fit.
	if( binNumber==pHistogram_->GetNbinsX() ) binNumber=pHistogram_->GetNbinsX()-1;
	// Likewise make sure if I'm at the start that when I fill the vector
	// the bins it uses actually exist.
	if( binNumber<3 ) binNumber=3;

        if(prt) printf("Using bin number %i\n",binNumber);
	
	// I now have the bin number of the bin after the point I'm looking
	// for. Now do a linear fit to interpolate between the bins using the
	// two bins for the point and the two after.
	std::vector<int> binNumberForLinearFit;
	binNumberForLinearFit.push_back( binNumber-2 );
	binNumberForLinearFit.push_back( binNumber-1 );
	binNumberForLinearFit.push_back( binNumber );
	binNumberForLinearFit.push_back( binNumber+1 );

	std::vector< std::pair<double,double> > dataPoints;
	for( auto& number : binNumberForLinearFit )
	{
		// Make sure all of of the bin numbers are valid
		if( number<1 ) number=1;
		else if( number>pHistogram_->GetNbinsX() ) number=pHistogram_->GetNbinsX();
		dataPoints.push_back( std::make_pair( pHistogram_->GetBinLowEdge(number), pHistogram_->GetBinContent(number) ) );
	}

	// Now do a simple linear fit on the data points
	std::pair<double,double> slopeAndIntercept=l1menu::tools::simpleLinearFit( dataPoints );
	double slope=slopeAndIntercept.first;
	double intercept=slopeAndIntercept.second;

        if(prt) printf("slope and intercept %f %f \n",slope,intercept);
	  
	const double zeroEqualityWithTolerance=std::pow(10,-4);
	if( std::fabs(slope)<zeroEqualityWithTolerance ) // see if the slope==0 within a little tolerance
	{
		// Some of the triggers, e.g. muons, can only have the threshold set in discrete steps. So the rate plot
		// comes out looking as though it has steps in it. All of the bins for the linear fit could be on the same
		// step (I think this will only happen if it's the bins at the high threshold end, but anyway). In this
		// case I'll just come back along the plot and return the value for the lowest bin, i.e. the lowest
		// threshold for this particular "step".
		while( binNumber>1 && pHistogram_->GetBinContent(binNumber)==pHistogram_->GetBinContent(binNumber-1) ) --binNumber;
		if(prt) printf("Slope too small returning lower edge of bin %i \n",binNumber);
		return pHistogram_->GetBinLowEdge(binNumber);
	}
	else
	{
		double newThreshold=(targetRate-intercept)/slope;
		// Apply some sanity checks. If the bin for this new threshold is more than an arbitrary
		// number away from the two straddling the required threshold (i.e. binNumber-1 and binNumber)
		// then redo the fit with just those two bins in it, so that the result is guaranteed to be between
		// the two. This can happen at points of inflection. Also need to make sure the newThreshold is
		// on the histogram.
		int newThresholdBin=pHistogram_->FindBin(newThreshold);
		
		if(prt) printf("New Thr %f  New Bin %i \n",newThreshold,newThresholdBin);
		
		// Make sure it's not the under or overflow bin, and that it's arbitrarily close to the original bin.
		if( std::abs(binNumber-newThresholdBin)>-1 || newThresholdBin==0 || newThresholdBin==pHistogram_->GetNbinsX()+1 )
		{
		        if(prt) printf("Revaluation using fit to just two bins base binNumber %i\n",binNumber);
			dataPoints.clear();
			dataPoints.push_back( std::make_pair( pHistogram_->GetBinLowEdge(binNumber-1), pHistogram_->GetBinContent(binNumber-1) ) );
			dataPoints.push_back( std::make_pair( pHistogram_->GetBinLowEdge(binNumber), pHistogram_->GetBinContent(binNumber) ) );
			slopeAndIntercept=l1menu::tools::simpleLinearFit( dataPoints );
			// Need to check that the slope isn't zero (to within a little tolerance)
			if( std::fabs(slopeAndIntercept.first)<zeroEqualityWithTolerance )
			{
				// I want to return the lower edge of the bin, but the binNumber has been forced to greater
				// than 3, so I need to check I'm still on the first bin lower than the requested rate.
				while( binNumber>1 && pHistogram_->GetBinContent(binNumber)<targetRate ) --binNumber;
				newThreshold=pHistogram_->GetBinLowEdge(binNumber);
				if(prt) printf("Slope too small returning lower edge of bin %i \n",binNumber);
			}
			else newThreshold=(targetRate-slopeAndIntercept.second)/slopeAndIntercept.first;
			if(prt) printf("Revaluation: New Thr %f  \n",newThreshold);
		}
		if(prt) printf("Final thres %f \n",newThreshold); 
		//newThreshold = pHistogram_->GetBinLowEdge(pHistogram_->FindBin(newThreshold));
		if(prt) printf("Return Bin Edge: New Thr %f   \n",newThreshold);
		if( newThreshold<0 ) return 0; // Apply some sanity checks
		// Do I want to extrapolate out past where the histogram goes? Not sure. I won't for the time being.
		else if( newThreshold>pHistogram_->GetXaxis()->GetXmax() ) return pHistogram_->GetXaxis()->GetXmax();
		else return newThreshold;
	}
}

std::pair<float,float> l1menu::TriggerRatePlot::findThresholdError( float threshold, float rate ) const
{

        bool prt=false;
	if(strcmp(pHistogram_->GetName(),"L1_SingleMu_v_threshold1")==0 ) prt=true;
	if(prt) printf("FindError: Found trigger to debug %s with threshold %f and rate %f\n",pHistogram_->GetName(),threshold,rate);

	int binNumber=pHistogram_->FindFixBin(threshold);
	// If the bin is the under or overflow bin then there's not much I can do.
	if( binNumber==0 || binNumber==pHistogram_->GetNbinsX()+1 )
		throw std::runtime_error("TriggerRatePlot::findThresholdError was called with a threshold not on the histogram "\
				"(threshold="+std::to_string(threshold)+" title="+pHistogram_->GetTitle()+")" );

	//float rate=pHistogram_->GetBinContent(binNumber);
	
	// Estimate rate with this threshold
	float rl =pHistogram_->GetBinContent(binNumber);
	float tl =pHistogram_->GetBinLowEdge(binNumber);
	float rh =pHistogram_->GetBinContent(binNumber+1);
	float th =pHistogram_->GetBinLowEdge(binNumber+1);
	rate = (rl-rh)*(threshold-th)/(tl-th)+rh;
	if(prt) printf("New Rate estimate %f \n",rate);
	
	float sysRate = 0.3;
        if(prt) printf("----------------------------------------------\nTrigger: %s \n",pHistogram_->GetName());
	int lowestMatchingBin=binNumber+1;
	int binlowKeep = lowestMatchingBin;
	while( lowestMatchingBin>1 )
	{
		// The online documentation says TH1 has GetBinErrorLow but it's giving me
		// compiler errors. Must be a different version of Root.
		if( (1.-sysRate)*(pHistogram_->GetBinContent(lowestMatchingBin-1)-pHistogram_->GetBinError(lowestMatchingBin-1))>rate ) {
		  binlowKeep = lowestMatchingBin;
		  while(lowestMatchingBin>1) {
		     // Keep backing up for bins that are identical...dealing with muons
		     if(pHistogram_->GetBinContent(lowestMatchingBin-1) != pHistogram_->GetBinContent(binlowKeep-1)) break;
		     --lowestMatchingBin;
		  }
		   break;
		}   
		--lowestMatchingBin;
	}
	// Not sure what to do if I get the underflow bin. I'll just return
	// the lower edge of the histogram for now - i.e. the lower edge of
	// bin number one. The condition on the while loop makes sure I'll
	// get bin number one if the break is never hit.

	int highestMatchingBin=binNumber-1;
	while( highestMatchingBin<pHistogram_->GetNbinsX()+1 )
	{
		// Same thing for GetBinErrorUp - gives me compilation errors
		if( (1.+sysRate)*(pHistogram_->GetBinContent(highestMatchingBin)+pHistogram_->GetBinError(highestMatchingBin))<rate ) break;
		++highestMatchingBin;
	}
	// Again, don't know what to do if I get the overflow bin. I'll just
	// return the high edge of the histogram, i.e. the low edge of the
	// overflow bin. The condition of the while loop should make sure I
	// get the overflow bin if the break is never hit.

        //return std::make_pair( threshold-pHistogram_->GetBinLowEdge(lowestMatchingBin), pHistogram_->GetBinLowEdge(highestMatchingBin)-threshold );

        // Let's try something more sophisticated with the errors.  Interpolate between bins to get a better estimate of the error.
	float t1 = pHistogram_->GetBinLowEdge(binlowKeep);
	float r1 = (1.-sysRate)*(pHistogram_->GetBinContent(binlowKeep)-pHistogram_->GetBinError(binlowKeep));
	float t2 = pHistogram_->GetBinLowEdge(lowestMatchingBin-1);
	float r2 = (1.-sysRate)*(pHistogram_->GetBinContent(lowestMatchingBin-1)-pHistogram_->GetBinError(lowestMatchingBin-1));
	float lowThr = (rate-r1)*((t2-t1)/(r2-r1))+t1;
	if(prt) printf("lowBin %i t1,r1 %f, %f  t2,r2 %f, %f lowThr %f \n",lowestMatchingBin, t1,r1, t2,r2, lowThr);

	t1 = pHistogram_->GetBinLowEdge(highestMatchingBin);
	r1 = (1.+sysRate)*(pHistogram_->GetBinContent(highestMatchingBin)+pHistogram_->GetBinError(highestMatchingBin));
	t2 = pHistogram_->GetBinLowEdge(highestMatchingBin-1);
	r2 = (1.+sysRate)*(pHistogram_->GetBinContent(highestMatchingBin-1)+pHistogram_->GetBinError(highestMatchingBin-1));
	float highThr = (rate-r1)*((t2-t1)/(r2-r1))+t1;
	if(prt) printf("highBin %i t1,r1 %f, %f  t2,r2 %f, %f lowThr %f \n",highestMatchingBin, t1,r1, t2,r2, highThr);

	if(prt) printf("\n Rate %f  lowThrErr %f threshold %f  highThrErr %f\n",rate,threshold-lowThr,threshold,highThr-threshold);


	return std::make_pair( threshold-lowThr, highThr-threshold );
}

TH1* l1menu::TriggerRatePlot::getPlot()
{
	return pHistogram_.get();
}

const TH1* l1menu::TriggerRatePlot::getPlot() const
{
	return pHistogram_.get();
}

TH1* l1menu::TriggerRatePlot::relinquishOwnershipOfPlot()
{
	// Rather than call release() on the unique_ptr, I'll set a flag so that I know
	// to release() in the destructor instead. This way it's possible to relinquish
	// ownership but still perform operations on the histograms.
	histogramOwnedByMe_=false;
	return pHistogram_.get();
}

void l1menu::TriggerRatePlot::addSample( const l1menu::ISample& sample, std::vector<TriggerRatePlot>& ratePlots )
{
	float weightPerEvent=sample.eventRate()/sample.sumOfWeights();

	// Create cached triggers for each of the rate plots, which depending on the concrete type
	// of the ISample may or may not significantly increase the speed at which this next loop happens.
	std::vector< std::unique_ptr<l1menu::ICachedTrigger> > cachedTriggers;
	for( const auto& ratePlot : ratePlots ) cachedTriggers.push_back( sample.createCachedTrigger( *ratePlot.pTrigger_ ) );

	// Now instead of calling addSample() for each TriggerRatePlot individually, get each IEvent from the sample
	// and pass that to each rate plot. This is because (depending on the ISample concrete type) getting the
	// IEvent can be computationally expensive.
	std::vector< std::unique_ptr<l1menu::ICachedTrigger> >::const_iterator iTrigger;
	std::vector<TriggerRatePlot>::iterator iRatePlot;
	for( size_t eventNumber=0; eventNumber<sample.numberOfEvents(); ++eventNumber )
	{
		const l1menu::IEvent& event=sample.getEvent(eventNumber);

		for( iTrigger=cachedTriggers.begin(), iRatePlot=ratePlots.begin();
			iTrigger!=cachedTriggers.end() && iRatePlot!=ratePlots.end();
			++iTrigger, ++iRatePlot )
		{
			iRatePlot->addEvent( event, *iTrigger, weightPerEvent );
		}
	} // end of loop over events

}
