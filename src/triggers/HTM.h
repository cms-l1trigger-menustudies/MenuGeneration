#ifndef l1menu_triggers_HTM_h
#define l1menu_triggers_HTM_h

#include <string>
#include <vector>
#include "l1menu/ITrigger.h"

//
// Forward declarations
//
namespace l1menu
{
	class L1TriggerDPGEvent;
}

namespace l1menu
{
	namespace triggers
	{
		/** @brief Base class for all versions of the HTM trigger.
		 *
		 * Note that this class is abstract because it doesn't implement the "version"
		 * and "apply" methods. That's left up to the implementations of the different
		 * versions.
		 *
		 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
		 * @date 02/Jun/2013
		 */
		class HTM : public l1menu::ITrigger
		{
		public:
			HTM();

			virtual const std::string name() const;
			virtual const std::vector<std::string> parameterNames() const;
			virtual float& parameter( const std::string& parameterName );
			virtual const float& parameter( const std::string& parameterName ) const;
		protected:
			float threshold1_;
		}; // end of the HTM base class

		/** @brief First version of the HTM trigger.
		 *
		 * @author probably Brian Winer
		 * @date sometime
		 */
		class HTM_v0 : public HTM
		{
		public:
			virtual unsigned int version() const;
			virtual bool apply( const l1menu::L1TriggerDPGEvent& event ) const;
			virtual bool thresholdsAreCorrelated() const;
		}; // end of version 0 class


		/** @brief HTM trigger where the jets are looped over here rather than in FullSample.
		 *
		 * Allows eta cuts to be applied here, rather than having a hard coded eta cut in FullSample.
		 * Also added a minimum pT value if the jet is to be used to calculate the vector sum.
		 *
		 * @author Mark Grimes, but just copied Brian's code from FullSample.cpp
		 * @date 28/Mar/2014
		 */
		class HTM_v1 : public HTM
		{
		public:
			HTM_v1(); // Need a constructor to initiate regionCut_
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

	} // end of namespace triggers

} // end of namespace l1menu

#endif
