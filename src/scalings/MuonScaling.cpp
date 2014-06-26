#include "l1menu/scalings/MuonScaling.h"

#include <stdexcept>
#include "l1menu/TriggerTable.h"
#include "l1menu/ITrigger.h"
#include "l1menu/ITriggerDescriptionWithErrors.h"
#include "l1menu/TriggerRatePlot.h"
#include "l1menu/MenuRatePlots.h"
#include "l1menu/IMenuRate.h"
#include "../implementation/MenuRateImplementation.h"
#include <TH1.h>
#include <TFile.h>

namespace l1menu
{
	namespace scalings
	{
		/** @brief Private members for the MuonScaling class, using the pimple idiom.
		 *
		 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
		 * @date 16/Oct/2013
		 */
		class MuonScalingPrivateMembers
		{
		public:
			std::string detailedDescription_;
			std::unique_ptr<TH1> pMuonScalePtPlot;
			std::unique_ptr<TH1> pMuonScaleIsolationPlot;
			std::unique_ptr<l1menu::MenuRatePlots> pScaledRatePlots_;

			/** @brief Effectively the same as the method in MCDataScaling to scale TriggerRatePlots except that it
			 * doesn't take a copy.
			 *
			 * Definition is at the bottom of the file.
			 * Required because scale(MenuRatePlots) needs to delegate to this, but needs to act directly on
			 * the histogram rather than a copy.
			 */
			void scaleTriggerRatePlot( l1menu::TriggerRatePlot& unscaledPlot );
		};

	} // end of namespace scalings
} // end of namespace l1menu

l1menu::scalings::MuonScaling::MuonScaling( const std::string& muonScalingFilename )
	: pImple( new MuonScalingPrivateMembers )
{
	pImple->detailedDescription_="Muon scaling filename: '"+muonScalingFilename+"'";

	std::unique_ptr<TFile> pMuonScalingFile( TFile::Open(muonScalingFilename.c_str()) );
	if( pMuonScalingFile==nullptr ) throw std::runtime_error( "Couldn't create MenuRateMuonScaling because couldn't open TFile "+muonScalingFilename );

	pImple->pMuonScalePtPlot.reset( static_cast<TH1*>( pMuonScalingFile->Get("muonPtScale")->Clone() ) );
	pImple->pMuonScalePtPlot->SetDirectory( nullptr ); // Make sure it's only held in memory. Might be the default for clones but not sure.
	if( pImple->pMuonScalePtPlot==nullptr ) throw std::runtime_error( "Couldn't create MenuRateMuonScaling because couldn't get the TH1 \"muonPtScale\" from the  TFile" );

	pImple->pMuonScaleIsolationPlot.reset( static_cast<TH1*>( pMuonScalingFile->Get("muonIsoScale")->Clone() ) );
	pImple->pMuonScaleIsolationPlot->SetDirectory( nullptr ); // Make sure it's only held in memory.
	if( pImple->pMuonScaleIsolationPlot==nullptr ) throw std::runtime_error( "Couldn't create MenuRateMuonScaling because couldn't get the TH1 \"muonIsoScale\" from the  TFile" );

	pMuonScalingFile->Close(); // Taken copies of the muon plots so I can close the file now.
}

l1menu::scalings::MuonScaling::MuonScaling( const std::string& muonScalingFilename, const std::string& unscaledRatesFilename )
	: pImple( new MuonScalingPrivateMembers )
{
	pImple->detailedDescription_="Muon scaling filename: '"+muonScalingFilename+"', unscaled rate filename: '"+unscaledRatesFilename+"'";

	std::unique_ptr<TFile> pMuonScalingFile( TFile::Open(muonScalingFilename.c_str()) );
	if( pMuonScalingFile==nullptr ) throw std::runtime_error( "Couldn't create MenuRateMuonScaling because couldn't open TFile "+muonScalingFilename );

	pImple->pMuonScalePtPlot.reset( static_cast<TH1*>( pMuonScalingFile->Get("muonPtScale")->Clone() ) );
	pImple->pMuonScalePtPlot->SetDirectory( nullptr ); // Make sure it's only held in memory. Might be the default for clones but not sure.
	if( pImple->pMuonScalePtPlot==nullptr ) throw std::runtime_error( "Couldn't create MenuRateMuonScaling because couldn't get the TH1 \"muonPtScale\" from the  TFile" );

	pImple->pMuonScaleIsolationPlot.reset( static_cast<TH1*>( pMuonScalingFile->Get("muonIsoScale")->Clone() ) );
	pImple->pMuonScaleIsolationPlot->SetDirectory( nullptr ); // Make sure it's only held in memory.
	if( pImple->pMuonScaleIsolationPlot==nullptr ) throw std::runtime_error( "Couldn't create MenuRateMuonScaling because couldn't get the TH1 \"muonIsoScale\" from the  TFile" );

	pMuonScalingFile->Close(); // Taken copies of the muon plots so I can close the file now.


	//
	// Take a copy of all of the unscaled rate plots
	//
	std::unique_ptr<TFile> pRatePlotsRootFile( TFile::Open( unscaledRatesFilename.c_str() ) );
	pImple->pScaledRatePlots_.reset( new l1menu::MenuRatePlots( pRatePlotsRootFile.get() ) );
	// Loop over all of them and scale any of the muon plots by the required amount
	for( auto& triggerRatePlot : pImple->pScaledRatePlots_->triggerRatePlots() ) pImple->scaleTriggerRatePlot( triggerRatePlot );

}

l1menu::scalings::MuonScaling::~MuonScaling()
{
	// No operation
}

std::string l1menu::scalings::MuonScaling::briefDescription()
{
	return "Muon scaling for expected improvement in muon pT assignment and isolation";
}

std::string l1menu::scalings::MuonScaling::detailedDescription()
{
	return pImple->detailedDescription_;
}

std::unique_ptr<l1menu::TriggerRatePlot> l1menu::scalings::MuonScaling::scale( const l1menu::TriggerRatePlot& unscaledPlot )
{
	// Take a copy and work on that
	std::unique_ptr<l1menu::TriggerRatePlot> pReturnValue( new l1menu::TriggerRatePlot( unscaledPlot ) );

	// Delegate to a private method to do the work
	pImple->scaleTriggerRatePlot( *pReturnValue );

	return pReturnValue;
}

std::unique_ptr<l1menu::MenuRatePlots> l1menu::scalings::MuonScaling::scale( const l1menu::MenuRatePlots& unscaledPlots )
{
	// First take a copy of the input and modify that directly
	std::unique_ptr<l1menu::MenuRatePlots> pReturnValue( new l1menu::MenuRatePlots( unscaledPlots ) );

	// Loop over all of the plots and do the business
	for( auto& triggerRatePlot : pReturnValue->triggerRatePlots() )
	{
		pImple->scaleTriggerRatePlot( triggerRatePlot );
	}

	return pReturnValue;
}

std::unique_ptr<l1menu::IMenuRate> l1menu::scalings::MuonScaling::scale( const l1menu::IMenuRate& unscaledMenuRate )
{
	// If this instance was constructed without the unscaled rate plots then I can't scale the
	// IMenuRates, so I  have to throw an exception.
	if( pImple->pScaledRatePlots_==nullptr ) throw std::runtime_error( "MCDataScaling was asked to scale an IMenuRate, but the" \
			"unscaled IMenuRate was not provided in the constructor. You need to provided the unscaled IMenuRate." );

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

	for( const auto& pUnscaledTriggerRate : unscaledMenuRate.triggerRates() )
	{
		const l1menu::ITriggerDescriptionWithErrors& trigger=pUnscaledTriggerRate->trigger();

		//
		// First I'll get a copy of the trigger whether it's a muon trigger or not. If it is a muon
		// trigger I'll overwrite the threshold information in the next "if" block. If it isn't, I'll
		// just leave it as it is.
		//
		std::unique_ptr<l1menu::ITrigger> pScaledTrigger=l1menu::TriggerTable::instance().copyTrigger( trigger );

		//
		// I'll take a copy of all the available parameter errors from the unscaled TriggerRate.
		// Some of these might be overwritten in a minute if they're for muons. I'll later
		// write these into the scaled TriggerRate
		//
		std::map< std::string, std::pair<float,float> > parameterErrors;
		for( const auto& parameterName : trigger.parameterNames() )
		{
			if( trigger.parameterErrorsAreAvailable(parameterName) )
			{
				parameterErrors[parameterName]=std::make_pair( trigger.parameterErrorLow(parameterName), trigger.parameterErrorHigh(parameterName) );
			}
		}

		// Brian's original comment: "Scaling for Muon System Upgrades  (Problem: Can't easily scale
		// the EG_Mu trigger because x-axis is EG Et not muon Pt...need to fix)"
		//
		// No need to do a check for "IsoMu" or "isoMu" because the plot I'm going to use to create the new
		// threshold has had both the pT assignment and isolation scaling already applied.
		if( trigger.name().find("Mu")!=std::string::npos && trigger.name().find("EG_Mu")==std::string::npos )
		{
			const l1menu::TriggerRatePlot* pScaledTriggerRatePlot=nullptr;
			for( const auto& triggerRatePlot : pImple->pScaledRatePlots_->triggerRatePlots() )
			{
				if( triggerRatePlot.triggerMatches( trigger ) )
				{
					pScaledTriggerRatePlot=&triggerRatePlot;
					break;
				}
			}
			if( pScaledTriggerRatePlot==nullptr ) throw std::runtime_error( "Can't scale MenuRate because the trigger "+trigger.name()+" has no scaled TriggerRatePlot." );

			// Now that I have the scaled plot, I can read off what threshold gives the same rate
			// as the unscaled threshold.
			float scaledThreshold=pScaledTriggerRatePlot->findThreshold( pUnscaledTriggerRate->rate() );
			// Also add the threshold errors.
			std::pair<float,float> primaryThresholdErrors=pScaledTriggerRatePlot->findThresholdError( scaledThreshold, pUnscaledTriggerRate->rate()  );
			parameterErrors[pScaledTriggerRatePlot->versusParameter()]=primaryThresholdErrors;

			// Then modify the thresholds on the trigger
			pScaledTrigger->parameter( pScaledTriggerRatePlot->versusParameter() )=scaledThreshold;
			for( const auto& nameScalingPair : pScaledTriggerRatePlot->otherScaledParameters() )
			{
				pScaledTrigger->parameter( nameScalingPair.first )=nameScalingPair.second*scaledThreshold;
				parameterErrors[nameScalingPair.first]=std::make_pair( nameScalingPair.second*primaryThresholdErrors.first, nameScalingPair.second*primaryThresholdErrors.second );
			}

		}

		// Now I can add the TriggerRate to the MenuRate. It it was a muon trigger the thresholds will have been
		// scaled, if it wasn't then it will be exactly the same. I can always copy the rates straight over because
		// it's the thresholds that are scaled.
		l1menu::implementation::TriggerRateImplementation scaledTriggerRate( *pScaledTrigger,
				pUnscaledTriggerRate->fraction(), pUnscaledTriggerRate->fractionError(),
				pUnscaledTriggerRate->rate(), pUnscaledTriggerRate->rateError(),
				pUnscaledTriggerRate->pureFraction(), pUnscaledTriggerRate->pureFractionError(),
				pUnscaledTriggerRate->pureRate(), pUnscaledTriggerRate->pureRateError() );
		// Now that I've created the TriggerRate object I can set the threshold errors that I've calculated
		for( const auto& nameErrorsPair : parameterErrors )
		{
			scaledTriggerRate.setParameterErrors( nameErrorsPair.first, nameErrorsPair.second.first, nameErrorsPair.second.second );
		}
		pMenuRate->addTriggerRate( std::move(scaledTriggerRate) );

	} // end of loop over unscaled TriggerRates

	return pReturnValue;
}

void l1menu::scalings::MuonScalingPrivateMembers::scaleTriggerRatePlot( l1menu::TriggerRatePlot& unscaledPlot )
{
	const std::string triggerName=unscaledPlot.getTrigger().name();

	// Brian's original comment: Scaling for Muon System Upgrades  (Problem: Can't easily scale
	// the EG_Mu trigger because x-axis is EG Et not muon Pt...need to fix)
	if( triggerName.find("Mu")!=std::string::npos && triggerName.find("EG_Mu")==std::string::npos )
	{
		// I changed the binning on a lot of the histograms, but haven't done so for the muon scaling
		// plots. As a quick thing to get results for the TPSWG meeting I'll multiply each bin by hand.
		// The root TH1::Rebin method doesn't work very well unless the bin edges line up.
		// The pMuonScalePtPlot plot is quite discretised anyway, and the binning is almost similar. So
		// it won't make much difference. Mark Grimes 18/Oct/2013.
		TH1* pRawPlot=unscaledPlot.getPlot();
		for( int bin=1; bin<=pRawPlot->GetNbinsX(); ++bin )
		{
			// Using the bin low edge because that's the value I use to see if a bin has
			// passed a trigger.
			int scaleBin=pMuonScalePtPlot->FindBin( pRawPlot->GetBinLowEdge( bin ) );
			float scaleFactor=pMuonScalePtPlot->GetBinContent(scaleBin);
			pRawPlot->SetBinContent( bin, pRawPlot->GetBinContent(bin)*scaleFactor );
		}
	}

	if( triggerName.find("IsoMu")!=std::string::npos || triggerName.find("isoMu")!=std::string::npos )
	{
		// See note above for pt assignment scale factor
		TH1* pRawPlot=unscaledPlot.getPlot();
		for( int bin=1; bin<=pRawPlot->GetNbinsX(); ++bin )
		{
			int scaleBin=pMuonScaleIsolationPlot->FindBin( pRawPlot->GetBinLowEdge( bin ) );
			float scaleFactor=pMuonScaleIsolationPlot->GetBinContent(scaleBin);
			pRawPlot->SetBinContent( bin, pRawPlot->GetBinContent(bin)*scaleFactor );
		}
	}

}
