#include "CaloObject/Pad.h"

#include <TEfficiency.h>

namespace caloobject
{

Pad::Pad(int _id)
	: id(_id) ,
	  position()
{
	thresholds.push_back(0) ;
	nDetected.push_back(0) ;
	efficiencies.push_back(0.0) ;
	multiSumVec.push_back(0.0) ;
	multiSquareSumVec.push_back(0.0) ;
	multiplicities.push_back(0.0) ;
}

void Pad::reset()
{
	nTracks = 0 ;
	nDetected.clear() ;

	thresholds.clear() ;
	efficiencies.clear() ;
	multiSumVec.clear() ;
	multiSquareSumVec.clear() ;
	multiplicities.clear() ;
}

void Pad::setThresholds(const std::vector<double>& thr)
{
	reset() ;

	thresholds = thr ;

	nDetected = std::vector<int>(thresholds.size() , 0) ;
	efficiencies = std::vector<double>(thresholds.size() , 0.0) ;

	multiSumVec = std::vector<double>(thresholds.size() , 0.0) ;
	multiSquareSumVec = std::vector<double>(thresholds.size() , 0.0) ;
	multiplicities = std::vector<double>(thresholds.size() , 0.0) ;
}

std::vector<double> Pad::getEfficienciesError() const
{
	std::vector<double> toReturn ;

	for ( unsigned int i = 0 ; i < nDetected.size() ; ++i )
	{
		double eff , effErr ;

		if ( nDetected.at(i) == nTracks )
			eff = (nTracks-1.0)/nTracks ;
		else if ( nDetected.at(i) == 0 )
			eff = 1.0/nTracks ;
		else
			eff = static_cast<double>( nDetected.at(i) )/nTracks ;


		effErr = sqrt( eff*(1-eff)/nTracks ) ;
		toReturn.push_back(effErr) ;
	}
	return toReturn ;
}

std::array< std::vector<double> , 2> Pad::getEfficienciesBound() const
{
	constexpr double level = 0.683 ;

	std::vector<double> lowerBoundVec = {} ;
	std::vector<double> upperBoundVec = {} ;

	lowerBoundVec.reserve( nDetected.size() ) ;
	upperBoundVec.reserve( nDetected.size() ) ;

	for ( unsigned int i = 0 ; i < nDetected.size() ; ++i )
	{
		double a = double( nDetected.at(i) ) + 1 ;
		double b = double( nTracks - nDetected.at(i) ) + 1 ;

		double lowerBound = 0 ;
		double upperBound = 0 ;
		TEfficiency::BetaShortestInterval( level , a , b , lowerBound , upperBound ) ;

		lowerBoundVec.push_back( lowerBound ) ;
		upperBoundVec.push_back( upperBound ) ;
	}

	std::array< std::vector<double> , 2> toReturn = {{ lowerBoundVec , upperBoundVec }} ;
	return toReturn ;
}

std::vector<double> Pad::getMultiplicitiesError() const
{
	std::vector<double> toReturn ;

	for ( unsigned int i = 0 ; i < nDetected.size() ; ++i )
	{
		double mulErr ;
		if ( !nDetected.at(0) )
		{
			mulErr = 0.0 ;
			toReturn.push_back(mulErr) ;
			continue ;
		}

		double var = multiSquareSumVec.at(i)/nDetected.at(i) - (multiSumVec.at(i)/nDetected.at(i))*(multiSumVec.at(i)/nDetected.at(i)) ;

		if ( var < std::numeric_limits<double>::epsilon() )
			var = 1.0/( std::sqrt(12*nDetected.at(i)) ) ;

		if ( nDetected.at(i) < 2 )
			mulErr = multiSumVec.at(i)/nDetected.at(i) ;
		else
			mulErr = sqrt( var/(nDetected.at(i)-1.0) ) ;

		toReturn.push_back(mulErr) ;
	}
	return toReturn ;
}

void Pad::update(CaloCluster2D* cluster)
{
	nTracks++ ;

	if (cluster)
	{
		for ( unsigned int i = 0 ; i < thresholds.size() ; ++i )
		{
			if ( cluster->getMaxEnergy() >= thresholds.at(i) )
			{
				nDetected.at(i) ++ ;

				int nHitsOfThisThreshold = 0 ;
				for ( const auto hit : cluster->getHits() )
				{
					if ( hit->getEnergy() >= thresholds.at(i) )
						nHitsOfThisThreshold ++ ;
				}

				multiSumVec.at(i) += nHitsOfThisThreshold ;
				multiSquareSumVec.at(i) += nHitsOfThisThreshold*nHitsOfThisThreshold ;
			}
			else
				break ;
		}
	}
}

void Pad::updateEfficiencies()
{
	for ( unsigned int i = 0 ; i < efficiencies.size() ; ++i )
		efficiencies.at(i) = 1.0*nDetected.at(i)/nTracks ;
}

void Pad::updateMultiplicities()
{
	for ( unsigned int i = 0 ; i < multiplicities.size() ; ++i )
	{
		if ( nDetected.at(i) )
			multiplicities.at(i) = multiSumVec.at(i)/nDetected.at(i) ;
		else
			multiplicities.at(i) = 0.0 ;
	}
}

//===============
//======SDHCALPad
//===============

SDHCALPad::SDHCALPad(int _id)
	: Pad(_id)
{
}


} //namespace caloobject
