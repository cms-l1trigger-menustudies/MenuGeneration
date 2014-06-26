#include "l1menu/scalings/MCDataScaling.h"

#include <stdexcept>
#include "l1menu/MenuRatePlots.h"
#include "l1menu/ITrigger.h"
#include "l1menu/TriggerTable.h"
#include "../implementation/MenuRateImplementation.h"
#include <TFile.h>
#include <TH1.h>
#include <iostream>

namespace l1menu
{
	namespace scalings
	{
		/** @brief Private members for the MCDataScaling class, using the pimple idiom.
		 *
		 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
		 * @date 16/Oct/2013
		 */
		class MCDataScalingPrivateMembers
		{
		public:
			std::string detailedDescription_;
			std::unique_ptr<l1menu::MenuRatePlots> pMonteCarloRatePlots_;
			std::unique_ptr<l1menu::MenuRatePlots> pDataRatePlots_;
			std::unique_ptr<l1menu::MenuRatePlots> pUnscaledRatePlots_;

			/** @brief Effectively the same as the method in MCDataScaling to scale TriggerRatePlots except that it
			 * doesn't take a copy.
			 *
			 * Definition is at the bottom of the file.
			 * Required because scale(MenuRatePlots) needs to delegate to this, but needs to act directly on
			 * the histogram rather than a copy.
			 */
			void scaleTriggerRatePlot( l1menu::TriggerRatePlot& unscaledPlot );
			/** @brief Convenience method that just searches the IMenuRate for an ITriggerRate for the given trigger. */
			TH1* findRawPlot( l1menu::MenuRatePlots& ratePlotsToSearch, const l1menu::ITriggerDescription& trigger );
			const TH1* findRawPlot( const l1menu::MenuRatePlots& ratePlotsToSearch, const l1menu::ITriggerDescription& trigger );
		};

	} // end of namespace scalings
} // end of namespace l1menu

l1menu::scalings::MCDataScaling::MCDataScaling( const std::string& monteCarloRatesFilename, const std::string& dataRatesFilename )
	: pImple( new MCDataScalingPrivateMembers )
{
	pImple->detailedDescription_="Data rate filename: '"+dataRatesFilename+"', "
			+ "Monte Carlo rate filename: '"+monteCarloRatesFilename+"', "
			+ "no unscaled rate filename supplied";


	// Create copies of all of the rate plots
	std::unique_ptr<TFile> pRatePlotsRootFile( TFile::Open( monteCarloRatesFilename.c_str() ) );
	pImple->pMonteCarloRatePlots_.reset( new l1menu::MenuRatePlots( pRatePlotsRootFile.get() ) );

	pRatePlotsRootFile.reset( TFile::Open( dataRatesFilename.c_str() ) );
	pImple->pDataRatePlots_.reset( new l1menu::MenuRatePlots( pRatePlotsRootFile.get() ) );

	// If only plots are going to be scaled, then I don't need the unscaled rate plots. I'll
	// give people the option of creating the object without the unscaled rate plots, then if
	// if they later try and scale an IMenuRate I'll throw an exception.
}

l1menu::scalings::MCDataScaling::MCDataScaling( const std::string& monteCarloRatesFilename, const std::string& dataRatesFilename, const std::string& unscaledRatesFilename )
	: pImple( new MCDataScalingPrivateMembers )
{
	pImple->detailedDescription_="Data rate filename: '"+dataRatesFilename+"', "
			+ "Monte Carlo rate filename: '"+monteCarloRatesFilename+"', "
			+ "unscaled rate filename: '"+unscaledRatesFilename+"'";

	// Create copies of all of the rate plots
	std::unique_ptr<TFile> pRatePlotsRootFile( TFile::Open( monteCarloRatesFilename.c_str() ) );
	pImple->pMonteCarloRatePlots_.reset( new l1menu::MenuRatePlots( pRatePlotsRootFile.get() ) );

	pRatePlotsRootFile.reset( TFile::Open( dataRatesFilename.c_str() ) );
	pImple->pDataRatePlots_.reset( new l1menu::MenuRatePlots( pRatePlotsRootFile.get() ) );

	pRatePlotsRootFile.reset( TFile::Open( unscaledRatesFilename.c_str() ) );
	pImple->pUnscaledRatePlots_.reset( new l1menu::MenuRatePlots( pRatePlotsRootFile.get() ) );
}

l1menu::scalings::MCDataScaling::~MCDataScaling()
{
	// No operation
}

std::string l1menu::scalings::MCDataScaling::briefDescription()
{
	return "Data/MonteCarlo scaling";
}

std::string l1menu::scalings::MCDataScaling::detailedDescription()
{
	return pImple->detailedDescription_;
}

std::unique_ptr<l1menu::TriggerRatePlot> l1menu::scalings::MCDataScaling::scale( const l1menu::TriggerRatePlot& unscaledPlot )
{
	// Take a copy and work on that
	std::unique_ptr<l1menu::TriggerRatePlot> pReturnValue( new l1menu::TriggerRatePlot( unscaledPlot ) );

	TH1* pRawPlot=pReturnValue->getPlot();
	// First make sure it's held in memory. The user can change that later if they want.
	pRawPlot->SetDirectory( nullptr );

	// Delegate to the private method which acts directly on the object.
	pImple->scaleTriggerRatePlot( *pReturnValue );

	return pReturnValue;
}

std::unique_ptr<l1menu::MenuRatePlots> l1menu::scalings::MCDataScaling::scale( const l1menu::MenuRatePlots& unscaledPlots )
{
	// First take a copy of the whole lot, then run through them and scale
	std::unique_ptr<l1menu::MenuRatePlots> pReturnValue( new l1menu::MenuRatePlots(unscaledPlots) );

	for( auto& triggerRatePlot : pReturnValue->triggerRatePlots() )
	{
		// Delegate to the private method to do the work
		pImple->scaleTriggerRatePlot( triggerRatePlot );
	}

	return pReturnValue;
}

std::unique_ptr<l1menu::IMenuRate> l1menu::scalings::MCDataScaling::scale( const l1menu::IMenuRate& unscaledMenuRate )
{
	// If this instance was constructed without the unscaled rate plots then I can't scale the
	// IMenuRates, so I  have to throw an exception.
	if( pImple->pUnscaledRatePlots_==nullptr ) throw std::runtime_error( "MCDataScaling was asked to scale an IMenuRate, but the" \
			" unscaled IMenuRate was not provided in the constructor. You need to provided the unscaled IMenuRate." );

	// Create a new MenuRateImplementation. I need to access the extra setters that aren't in IMenuRate so
	// I'll keep hold of a raw pointer and operate on that.
	l1menu::implementation::MenuRateImplementation* pMenuRate=new l1menu::implementation::MenuRateImplementation;
	std::unique_ptr<l1menu::IMenuRate> pReturnValue( pMenuRate );

	// None of the rates change. What changes are the thresholds on the triggers, which are
	// modified so that they would give the same rate on the scaled rate plot as they do on
	// the unscaled plot. So I can copy over all of the rates directly.
	pMenuRate->setTotalFraction( unscaledMenuRate.totalFraction() );
	pMenuRate->setTotalFractionError( unscaledMenuRate.totalFractionError() );
	pMenuRate->setTotalRate( unscaledMenuRate.totalRate() );
	pMenuRate->setTotalRateError( unscaledMenuRate.totalRateError() );

	// Loop over each of the trigger rates. Again, I can copy over the rate directly but
	// I need to modify the thresholds on the trigger.
	for( const auto& pUnscaledTriggerRate : unscaledMenuRate.triggerRates() )
	{
		const l1menu::ITriggerDescription& trigger=pUnscaledTriggerRate->trigger();
		//
		// First I need to find the correct plots for data, Monte Carlo, and unscaled
		//
		const TH1* pDataRateRawPlot=pImple->findRawPlot( *pImple->pDataRatePlots_, trigger );
		if( pDataRateRawPlot==nullptr ) throw std::runtime_error( "Unable to scale for data/MC for rate plot for trigger "+trigger.name()+" because there is no matching data histogram." );

		// Now do the same for the Monte Carlo plots.
		const TH1* pMonteCarloRateRawPlot=pImple->findRawPlot( *pImple->pMonteCarloRatePlots_, trigger );
		if( pMonteCarloRateRawPlot==nullptr ) throw std::runtime_error( "Unable to scale for data/MC for rate plot for trigger "+trigger.name()+" because there is no matching Monte Carlo histogram." );

		// And finally for the unscaled plots
		const TH1* pUnscaledRateRawPlot=pImple->findRawPlot( *pImple->pUnscaledRatePlots_, trigger );
		if( pUnscaledRateRawPlot==nullptr ) throw std::runtime_error( "Unable to scale for data/MC for rate plot for trigger "+trigger.name()+" because there is no matching unscaled histogram." );

		// Make sure all of the binning is the same
		if( pUnscaledRateRawPlot->GetXaxis()->GetNbins()!=pDataRateRawPlot->GetXaxis()->GetNbins()
			|| pUnscaledRateRawPlot->GetXaxis()->GetXmin()!=pDataRateRawPlot->GetXaxis()->GetXmin()
			|| pUnscaledRateRawPlot->GetXaxis()->GetXmax()!=pDataRateRawPlot->GetXaxis()->GetXmax() )
				throw std::runtime_error( "Unable to scale for data/MC for trigger "+trigger.name()+" because the binning of the data histogram does not match." );
		if( pUnscaledRateRawPlot->GetXaxis()->GetNbins()!=pMonteCarloRateRawPlot->GetXaxis()->GetNbins()
			|| pUnscaledRateRawPlot->GetXaxis()->GetXmin()!=pMonteCarloRateRawPlot->GetXaxis()->GetXmin()
			|| pUnscaledRateRawPlot->GetXaxis()->GetXmax()!=pMonteCarloRateRawPlot->GetXaxis()->GetXmax() )
				throw std::runtime_error( "Unable to scale for data/MC for trigger "+trigger.name()+" because the binning of the Monte Carlo histogram does not match." );

		// It's quite handy to have the plot as a TriggerRatePlot object (instead of raw root histogram)
		// so that I can use the findThreshold method. This constructor takes a copy of the histogram so
		// I don't need to worry about the modifications I'm about to do.
		l1menu::TriggerRatePlot scaledTriggerRatePlot( pUnscaledRateRawPlot );
		// Now scale it for the data/MC differences
		scaledTriggerRatePlot.getPlot()->Multiply( pDataRateRawPlot );
		scaledTriggerRatePlot.getPlot()->Divide( pMonteCarloRateRawPlot );

		// Now that I have the scaled plot, I can read off what threshold gives the same rate
		// as the unscaled threshold.
		float scaledThreshold=scaledTriggerRatePlot.findThreshold( pUnscaledTriggerRate->rate() );

		// I also need to find out what the error on this threshold is
		std::map< std::string, std::pair<float,float> > thresholdErrors; // indexed on threshold name, pair.first is error low, second is error high
		std::pair<float,float> primaryThresholdErrors;
		try
		{
			primaryThresholdErrors=scaledTriggerRatePlot.findThresholdError( scaledThreshold, pUnscaledTriggerRate->rate() );
			thresholdErrors[scaledTriggerRatePlot.versusParameter()]=primaryThresholdErrors;
		}
		catch( std::runtime_error& exception )
		{
			/* Debugging why this wouldn't work */
			std::cerr << "Exception when trying to get " << scaledTriggerRatePlot.versusParameter() << " errors at " << scaledThreshold <<std::endl;
		}

		// Get a copy of the trigger
		std::unique_ptr<l1menu::ITrigger> pScaledTrigger=l1menu::TriggerTable::instance().copyTrigger( trigger );
		// Then modify it's thresholds to match the scaling
		pScaledTrigger->parameter( scaledTriggerRatePlot.versusParameter() )=scaledThreshold;
		for( const auto& nameScalingPair : scaledTriggerRatePlot.otherScaledParameters() )
		{
			pScaledTrigger->parameter( nameScalingPair.first )=nameScalingPair.second*scaledThreshold;
			if( !thresholdErrors.empty() )
			{
				thresholdErrors[nameScalingPair.first]=std::make_pair( nameScalingPair.second*primaryThresholdErrors.first, nameScalingPair.second*primaryThresholdErrors.second );
			}
		}

		l1menu::implementation::TriggerRateImplementation scaledTriggerRate( *pScaledTrigger,
				pUnscaledTriggerRate->fraction(), pUnscaledTriggerRate->fractionError(),
				pUnscaledTriggerRate->rate(), pUnscaledTriggerRate->rateError(),
				pUnscaledTriggerRate->pureFraction(), pUnscaledTriggerRate->pureFractionError(),
				pUnscaledTriggerRate->pureRate(), pUnscaledTriggerRate->pureRateError() );
		// Now that I've created the TriggerRate object I can set the threshold errors that I've calculated
		for( const auto& nameErrorsPair : thresholdErrors )
		{
			scaledTriggerRate.setParameterErrors( nameErrorsPair.first, nameErrorsPair.second.first, nameErrorsPair.second.second );
		}
		// And add it to the MenuRate
		pMenuRate->addTriggerRate( std::move(scaledTriggerRate) );
	}

	return pReturnValue;
}

void l1menu::scalings::MCDataScalingPrivateMembers::scaleTriggerRatePlot( l1menu::TriggerRatePlot& unscaledPlot )
{
	const l1menu::ITriggerDescription& trigger=unscaledPlot.getTrigger();

	// Loop over all of the data rate plots and see if one of them matches the plot
	// that needs to be scaled.
	const TH1* pDataRateRawPlot=findRawPlot( *pDataRatePlots_, trigger );
	if( pDataRateRawPlot==nullptr ) throw std::runtime_error( "Unable to scale for data/MC for rate plot for trigger "+trigger.name()+" because there is no matching data histogram." );

	// Now do the same for the Monte Carlo plots.
	const TH1* pMonteCarloRateRawPlot=findRawPlot( *pMonteCarloRatePlots_, trigger );
	if( pMonteCarloRateRawPlot==nullptr ) throw std::runtime_error( "Unable to scale for data/MC for rate plot for trigger "+trigger.name()+" because there is no matching Monte Carlo histogram." );

	//
	// Might as well check and make sure all of the binning matches
	//
	TH1* pRawPlotToScale=unscaledPlot.getPlot();
	if( pRawPlotToScale->GetXaxis()->GetNbins()!=pDataRateRawPlot->GetXaxis()->GetNbins()
		|| pRawPlotToScale->GetXaxis()->GetXmin()!=pDataRateRawPlot->GetXaxis()->GetXmin()
		|| pRawPlotToScale->GetXaxis()->GetXmax()!=pDataRateRawPlot->GetXaxis()->GetXmax() )
			throw std::runtime_error( "Unable to scale for data/MC for rate plot for trigger "+trigger.name()+" because the binning of the data histogram does not match." );
	if( pRawPlotToScale->GetXaxis()->GetNbins()!=pMonteCarloRateRawPlot->GetXaxis()->GetNbins()
		|| pRawPlotToScale->GetXaxis()->GetXmin()!=pMonteCarloRateRawPlot->GetXaxis()->GetXmin()
		|| pRawPlotToScale->GetXaxis()->GetXmax()!=pMonteCarloRateRawPlot->GetXaxis()->GetXmax() )
			throw std::runtime_error( "Unable to scale for data/MC for rate plot for trigger "+trigger.name()+" because the binning of the Monte Carlo histogram does not match." );

	//
	// Everything should be fine now, so I can just do the scaling
	//
	pRawPlotToScale->Multiply( pDataRateRawPlot );
	pRawPlotToScale->Divide( pMonteCarloRateRawPlot );
}

TH1* l1menu::scalings::MCDataScalingPrivateMembers::findRawPlot( l1menu::MenuRatePlots& ratePlotsToSearch, const l1menu::ITriggerDescription& trigger )
{
	for( auto& triggerRatePlot : ratePlotsToSearch.triggerRatePlots() )
	{
		if( triggerRatePlot.triggerMatches( trigger ) ) return triggerRatePlot.getPlot();
	}

	// If control gets to here a suitable trigger rate plot wasn't found
	return nullptr;
}

const TH1* l1menu::scalings::MCDataScalingPrivateMembers::findRawPlot( const l1menu::MenuRatePlots& ratePlotsToSearch, const l1menu::ITriggerDescription& trigger )
{
	for( const auto& triggerRatePlot : ratePlotsToSearch.triggerRatePlots() )
	{
		if( triggerRatePlot.triggerMatches( trigger ) ) return triggerRatePlot.getPlot();
	}

	// If control gets to here a suitable trigger rate plot wasn't found
	return nullptr;
}
