#ifndef TestParameters_h
#define TestParameters_h

#include <map>
#include <stdexcept>


/** @brief Class for test suites to query global parameters for the tests.
 *
 * For example, a series of test suites might want to have the filename of an input
 * file to run the tests on. The main function would set the filename from the calling
 * arguments and then each of the test suites could query it from this class.
 *
 * This class only allows retrieving of the parameters. The parameters are set by a
 * subclass that is implemented in the same file as main() (currently unitTestsMain.cpp).
 *
 * It follows the Meyer's singleton pattern. I've since made the parameters collection
 * and the getter static so there's not a lot of point getting an instance.
 *
 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
 * @date 09/Aug/2013
 */
template<class T>
class TestParameters
{
public:
	/** @brief Not actually a lot of use for this now I've made the getter static */
	static TestParameters& instance();
	static T getParameter( const std::string& parameterName );
	static const std::map<std::string,T>& parameters();
protected:
	TestParameters();
	virtual ~TestParameters();
	TestParameters( const TestParameters& otherParameters ) = delete;
	TestParameters( const TestParameters&& otherParameters ) = delete;
	TestParameters& operator=( const TestParameters& otherParameters ) = delete;
	TestParameters& operator=( const TestParameters&& otherParameters ) = delete;
protected:
	static std::map<std::string,T> parameters_;
};

template<class T>
std::map<std::string,T> TestParameters<T>::parameters_;

template<class T>
TestParameters<T>::TestParameters()
{
	// No operation
}

template<class T>
TestParameters<T>::~TestParameters()
{
	// No operation
}

template<class T>
TestParameters<T>& TestParameters<T>::instance()
{
	static TestParameters<T> onlyInstance;
	return onlyInstance;
}

template<class T>
T TestParameters<T>::getParameter( const std::string& parameterName )
{
	typename std::map<std::string,T>::const_iterator iFindResult=parameters_.find( parameterName );
	if( iFindResult!=parameters_.end() ) return iFindResult->second;
	throw std::invalid_argument( "TestParameters::getParameter(\""+parameterName+"\") - Invalid parameter name" );
}

template<class T>
const std::map<std::string,T>& TestParameters<T>::parameters()
{
	return parameters_;
}

#endif
