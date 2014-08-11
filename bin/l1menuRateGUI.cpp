#include <stdexcept>
#include <iostream>

#include <QtGui>
#include "l1menu/TriggerMenu.h"
#include "l1menu/ITrigger.h"
#include "l1menu/ITriggerDescription.h"
#include "l1menu/ITriggerDescriptionWithErrors.h"
#include "l1menu/IMenuRate.h"
#include "l1menu/ITriggerRate.h"
#include "l1menu/TriggerConstraint.h"
#include "l1menu/ReducedSample.h"
#include "l1menu/IL1MenuFile.h"
#include "l1menu/tools/miscellaneous.h"

/*                                     --------------------------
 *                                     -----   IMPORTANT   ------
 *                                     --------------------------
 *
 * If you modify this file by adding more slots you will need to regenerate the l1menuRateGUI.moc file.
 * This can be done with the command:
 * "moc $CMSSW_BASE/src/L1Trigger/MenuGeneration/bin/l1menuRateGUI.cpp > $CMSSW_BASE/src/L1Trigger/MenuGeneration/bin/l1menuRateGUI.moc"
 *
 */

//
// Use a new namespace for things only used in this file
//
namespace menuwidgets
{
	/** @brief Widget representation of a trigger, that provides GUI elements to modify thresholds
	 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
	 * @date 26/May/2014
	 */
	class TriggerWidget : public QWidget
	{
		Q_OBJECT
	public:
		TriggerWidget( l1menu::ITrigger& trigger );
		bool isEnabled() const;
		const l1menu::ITriggerDescription& trigger() const;
		void setRate( const l1menu::ITriggerRate* pRate, double totRate );
	private slots:
		void stateChanged( int state );
		void thresholdChanged( double value );
	private:
		bool enabled_;
		l1menu::ITrigger& trigger_;
		std::vector<QDoubleSpinBox*> freeThresholds_;
		std::vector< std::pair<QDoubleSpinBox*,double> > scaledThresholds_;
		std::unique_ptr<QLabel> pRateLabel_;
		std::unique_ptr<QLabel> pPureRateLabel_;
	};

	/** @brief Widget for the main window, which comprises everything else
	 *
	 * Note that the l1menu::ReducedSample passed in the constructor must exist for
	 * the entire time this object does.
	 *
	 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
	 * @date 27/May/2014
	 */
	class MainWidget : public QWidget
	{
		Q_OBJECT
	public:
		MainWidget( l1menu::ReducedSample& sample );
	private slots:
		void calculateRates();
		void saveMenu();
	private:
		l1menu::ReducedSample& sample_;
		l1menu::TriggerMenu menu_;
		std::vector<TriggerWidget*> triggerWidgets_;
		std::unique_ptr<QDoubleSpinBox> pCollisionRate_;
		std::unique_ptr<QPushButton> pCalculateRateButton_;
		std::unique_ptr<QLabel> pTotalRateLabel_;
	};
} // end of the menuwidgets namespace

int main( int argc, char **argv )
{
	if( argc!=2 )
	{
		std::cerr << "The only argument should be the filename of the ReducedSample" << std::endl;
		return 0;
	}

	QApplication app( argc, argv );

	l1menu::ReducedSample sample( argv[1] );

	menuwidgets::MainWidget myMainWidget( sample );
	myMainWidget.show();

	return app.exec();
}

//------------------------------------------------------------------
// Implementation of the utility classes
//------------------------------------------------------------------

menuwidgets::TriggerWidget::TriggerWidget( l1menu::ITrigger& trigger ) : enabled_(true), trigger_(trigger)
{
	std::unique_ptr<QHBoxLayout> pLayout( new QHBoxLayout );

	std::unique_ptr<QCheckBox> pCheckBox( new QCheckBox( (trigger.name()+" v"+std::to_string(trigger.version())).c_str() ) );
	pCheckBox->setTristate( false );
	pCheckBox->setCheckState( Qt::Checked );
	connect( pCheckBox.get(), SIGNAL(stateChanged(int)), this, SLOT(stateChanged(int)) );
	pLayout->addWidget( pCheckBox.release() );

	// Use unique_ptrs in case something goes wrong before they're added to the layout.
	// Once added the layout takes ownership so I can release the unique_ptrs.
	QDoubleSpinBox* pPrimaryThreshold=nullptr;
	for( const auto& thresholdName : l1menu::tools::getThresholdNames(trigger) )
	{
		std::unique_ptr<QDoubleSpinBox> pButton( new QDoubleSpinBox );
		pButton->setMinimum(0);
		pButton->setMaximum(1000);
		pButton->setValue( trigger_.parameter(thresholdName) );
		pButton->setToolTip( thresholdName.c_str() );
		connect( pButton.get(), SIGNAL(valueChanged(double)), this, SLOT(thresholdChanged(double)) );

		if( pPrimaryThreshold==nullptr )
		{
			pPrimaryThreshold=pButton.get();
			freeThresholds_.push_back( pButton.get() );
		}
		else
		{
			if( trigger.thresholdsAreCorrelated() )
			{
				pButton->setDisabled(true);
				double scaledValue=pButton->value()/pPrimaryThreshold->value();
				scaledThresholds_.push_back( std::make_pair(pButton.get(),scaledValue) );
			}
			else freeThresholds_.push_back( pButton.get() );
		}

		pLayout->addWidget( pButton.release() );
	}

	pLayout->addStretch();

	pRateLabel_.reset( new QLabel("? kHz") );
	pLayout->addWidget( pRateLabel_.get() );

	pLayout->addSpacing(5);

	pPureRateLabel_.reset( new QLabel("? kHz") );
	pLayout->addWidget( pPureRateLabel_.get() );

	this->setLayout( pLayout.release() );
}

bool menuwidgets::TriggerWidget::isEnabled() const
{
	return enabled_;
}

const l1menu::ITriggerDescription& menuwidgets::TriggerWidget::trigger() const
{
	return trigger_;
}

void menuwidgets::TriggerWidget::setRate( const l1menu::ITriggerRate* pRate, double totRate )
{
	//
	// First perform a check to make sure this rate matches this trigger
	if( trigger_.name()!=pRate->trigger().name() ) throw std::runtime_error("TriggerWidget::setRate was called with an incompatible trigger");
	if( trigger_.version()!=pRate->trigger().version() ) throw std::runtime_error("TriggerWidget::setRate was called with an incompatible trigger version");

	pRateLabel_->setText( (std::to_string(pRate->rate())/*+" +/-"+std::to_string(pRate->rateError())*/+" kHz").c_str() );
	pPureRateLabel_->setText( (std::to_string(pRate->pureRate())/*+" +/-"+std::to_string(pRate->pureRateError())*/+" kHz").c_str() );

	printf("%20s   %8.1f (%8.6f)",pRate->trigger().name().c_str(),pRate->rate(),pRate->rate()/totRate);
	for( const auto& thresholdName : l1menu::tools::getThresholdNames(pRate->trigger()) )
	{
	   printf(" %5.1f ",pRate->trigger().parameter(thresholdName));
        }
	printf("\n");

}

void menuwidgets::TriggerWidget::stateChanged( int state )
{
	if( state==Qt::Checked )
	{
		enabled_=true;
		for( auto& pThresholdWidget : freeThresholds_ ) pThresholdWidget->setDisabled(false);
		pRateLabel_->setDisabled(false);
		pPureRateLabel_->setDisabled(false);
	}
	else if( state==Qt::Unchecked )
	{
		enabled_=false;
		for( auto& pThresholdWidget : freeThresholds_ ) pThresholdWidget->setDisabled(true);
		pRateLabel_->setDisabled(true);
		pPureRateLabel_->setDisabled(true);
	}
}

void menuwidgets::TriggerWidget::thresholdChanged( double value )
{
	// The "value" parameter isn't actually used, because I can't tell which threshold
	// it's for. I just loop over all of them and read the values from the boxes.

	// Loop over the free thresholds and propagate the values in the boxes
	// to the trigger.
	for( auto& pWidget : freeThresholds_ )
	{
		// Parameter name is in the tool tip, but need to convert QString to C
		// string with toLocal8Bit().data() methods.
		trigger_.parameter( pWidget->toolTip().toLocal8Bit().data() )=pWidget->value();
	}

	// If there are scaled thresholds, then I need to modify those two.
	double primaryThresholdValue=freeThresholds_.front()->value();
	for( auto& widgetScalePair : scaledThresholds_ )
	{
		double newValue=primaryThresholdValue*widgetScalePair.second;
		widgetScalePair.first->setValue( newValue );
		// Also propagate the change to the trigger. The name of the threshold
		// parameter is stored in the tool tip for the widget (convert QString
		// to C string with toLocal8Bit().data()).
		trigger_.parameter( widgetScalePair.first->toolTip().toLocal8Bit().data() )=newValue;
	}

}

menuwidgets::MainWidget::MainWidget( l1menu::ReducedSample& sample ) : sample_(sample), menu_(sample.getTriggerMenu())
{
	std::unique_ptr<QVBoxLayout> pTriggerListLayout( new QVBoxLayout );

	for( size_t index=0; index<menu_.numberOfTriggers(); ++index )
	{
		std::unique_ptr<menuwidgets::TriggerWidget> pTriggerWidget( new menuwidgets::TriggerWidget( menu_.getTrigger(index) ) );
		triggerWidgets_.push_back( pTriggerWidget.get() );
		pTriggerListLayout->addWidget( pTriggerWidget.release() );
	}

	// For some reason the QScrollArea isn't accepting the QLayout, so get around that
	// I'm creating a dummy parent widget, setting the layout of that to pTriggerListLayout
	// then setting the dummy parent as the widget for the QScrollArea.
	std::unique_ptr<QWidget> pParentToAllTriggers( new QWidget );
	pParentToAllTriggers->setLayout(pTriggerListLayout.release());

	std::unique_ptr<QScrollArea> pScrollArea( new QScrollArea );
	pScrollArea->setWidget(pParentToAllTriggers.release());
	pScrollArea->setWidgetResizable( true );

	std::unique_ptr<QVBoxLayout> pMainLayout( new QVBoxLayout );
	pMainLayout->addWidget( pScrollArea.release() );

	pCollisionRate_.reset( new QDoubleSpinBox() );
	pCollisionRate_->setMinimum(0);
	pCollisionRate_->setMaximum(100000);
	{ // Block to limit the scope of temporary variables
		const float scaleToKiloHz=1.0/1000.0;
		const float orbitsPerSecond=11246;
		//const float numberOfBunches=1380; // for 50 ns bunch spacing
		const float numberOfBunches=2760; // for 25 ns bunch spacing
		pCollisionRate_->setValue( orbitsPerSecond*numberOfBunches*scaleToKiloHz );
	}
	pCalculateRateButton_.reset( new QPushButton( "&Calculate rates" ) );
	QObject::connect( pCalculateRateButton_.get(), SIGNAL( clicked() ), this, SLOT( calculateRates() ) );

	std::unique_ptr<QHBoxLayout> pCollisionRateLineLayout( new QHBoxLayout );
	pCollisionRateLineLayout->addWidget( new QLabel("Collision rate (kHz)=") );
	pCollisionRateLineLayout->addWidget( pCollisionRate_.get() );
	pCollisionRateLineLayout->addWidget( pCalculateRateButton_.get() );
	std::unique_ptr<QWidget> pCollisionRateLine( new QWidget );
	pCollisionRateLine->setLayout( pCollisionRateLineLayout.release() );

	pMainLayout->addWidget( pCollisionRateLine.release() );

	pTotalRateLabel_.reset( new QLabel("Total rate= ? kHz") );
	pMainLayout->addWidget( pTotalRateLabel_.get() );
	std::unique_ptr<QPushButton> pSaveMenuButton( new QPushButton( "&Save menu" ) );
	QObject::connect( pSaveMenuButton.get(), SIGNAL( clicked() ), this, SLOT( saveMenu() ) );
	pMainLayout->addWidget( pSaveMenuButton.release() );

	this->setLayout( pMainLayout.release() );
}

void menuwidgets::MainWidget::calculateRates()
{
	// Disable the button while doing the calculation, so that the user can't
	// click on it multiple times.
	pCalculateRateButton_->setDisabled(true);

	// Create a menu by copying just the triggers that are active
	l1menu::TriggerMenu menuForCalculation;
	for( const auto& pTriggerWidget : triggerWidgets_ )
	{
		if( pTriggerWidget->isEnabled() ) menuForCalculation.addTrigger( pTriggerWidget->trigger() );
	}

	sample_.setEventRate( pCollisionRate_->value() );
	std::shared_ptr<const l1menu::IMenuRate> pMenuRate=sample_.rate(menuForCalculation);
	pTotalRateLabel_->setText( ("Total rate= "+std::to_string(pMenuRate->totalRate())+" +/- "+std::to_string(pMenuRate->totalRateError())+" kHz").c_str() );

	// The trigger rates should be in the same order as the active trigger widgets
	printf("\n===========================================================\n");
	std::vector<const l1menu::ITriggerRate*>::const_iterator iTriggerRate=pMenuRate->triggerRates().begin();
	for( const auto& pTriggerWidget : triggerWidgets_ )
	{
		if( pTriggerWidget->isEnabled() )
		{
			pTriggerWidget->setRate( *iTriggerRate, pMenuRate->totalRate() );
			++iTriggerRate;
			if( iTriggerRate==pMenuRate->triggerRates().begin() ) break;
		}
	}

	// Re-enable the button
	printf(" Total Rate        %6.1f\n",pMenuRate->totalRate() );
	pCalculateRateButton_->setDisabled(false);
}

void menuwidgets::MainWidget::saveMenu()
{
	QString filename=QFileDialog::getSaveFileName(this, tr("Save File"), "", tr("XML Files (*.xml)"));

	if( filename!="" )
	{
		std::unique_ptr<l1menu::IL1MenuFile> pOutputL1MenuFile=l1menu::IL1MenuFile::getOutputFile( l1menu::IL1MenuFile::FileFormat::XML, filename.toLocal8Bit().data() );

		// Create a menu by copying just the triggers that are active
		l1menu::TriggerMenu menuToSave;
		for( const auto& pTriggerWidget : triggerWidgets_ )
		{
			if( pTriggerWidget->isEnabled() ) menuToSave.addTrigger( pTriggerWidget->trigger() );
		}

		// I also want to set the fraction of bandwidth constraint, which is only
		// used if the menu is scaled for a particular bandwidth. To work out what
		// this is I'll fit the menu again. I could read the rates from the text
		// boxes but the user could have edited the menu since the last fit.
		sample_.setEventRate( pCollisionRate_->value() );
		std::shared_ptr<const l1menu::IMenuRate> pMenuRate=sample_.rate(menuToSave);

		for( size_t index=0; index<menuToSave.numberOfTriggers(); ++index )
		{
			float fraction=pMenuRate->triggerRates()[index]->rate()/pMenuRate->totalRate();
			l1menu::TriggerConstraint& newConstraint=menuToSave.getTriggerConstraint(index);
			newConstraint.type( l1menu::TriggerConstraint::Type::FRACTION_OF_BANDWIDTH );
			newConstraint.value( fraction );
		}

		pOutputL1MenuFile->add( menuToSave );
	}
}

#include "l1menuRateGUI.moc"
