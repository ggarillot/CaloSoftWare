#include <CaloObject/Layer.h>

#include <TEfficiency.h>

namespace caloobject
{

Layer::Layer(int _id)
	: id(_id) ,
	  position() ,
	  asics()
{
	thresholds.push_back(0) ;
	nDetected.push_back(0) ;
	efficiencies.push_back(0.0) ;
	multiSumVec.push_back(0.0) ;
	multiSquareSumVec.push_back(0.0) ;
	multiplicities.push_back(0.0) ;
}

void Layer::reset()
{
	nTracks = 0 ;
	nDetected.clear() ;

	thresholds.clear() ;
	efficiencies.clear() ;
	multiSumVec.clear() ;
	multiSquareSumVec.clear() ;
	multiplicities.clear() ;

	for ( AsicMap::const_iterator it = asics.begin() ; it != asics.end() ; ++it )
		it->second->reset() ;
}

void Layer::setThresholds(const std::vector<double>& thr)
{
	reset() ;

	thresholds = thr ;

	nDetected = std::vector<int>(thresholds.size() , 0) ;
	efficiencies = std::vector<double>(thresholds.size() , 0.0) ;

	multiSumVec = std::vector<double>(thresholds.size() , 0.0) ;
	multiSquareSumVec = std::vector<double>(thresholds.size() , 0.0) ;
	multiplicities = std::vector<double>(thresholds.size() , 0.0) ;

	for ( AsicMap::const_iterator it = asics.begin() ; it != asics.end() ; ++it )
		it->second->setThresholds(thresholds) ;
}

std::vector<double> Layer::getEfficienciesError() const
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

std::array< std::vector<double> , 2> Layer::getEfficienciesBound() const
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

std::vector<double> Layer::getMultiplicitiesError() const
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

void Layer::update(const CLHEP::Hep3Vector& impactPos , CaloCluster2D* cluster)
{
	nTracks++ ;

	Asic* asic = findAsic(impactPos) ;
	asic->update(impactPos , cluster) ;

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

void Layer::updateEfficiencies()
{
	for ( unsigned int i = 0 ; i < efficiencies.size() ; ++i )
		efficiencies.at(i) = 1.0*nDetected.at(i)/nTracks ;
}

void Layer::updateMultiplicities()
{
	for ( unsigned int i = 0 ; i < multiplicities.size() ; ++i )
	{
		if ( nDetected.at(i) )
			multiplicities.at(i) = multiSumVec.at(i)/nDetected.at(i) ;
		else
			multiplicities.at(i) = 0.0 ;
	}
}


SDHCALLayer::SDHCALLayer(int _id, int difL, int difC, int difR)
	: Layer(_id)
{
	difID[0] = difL ;
	difID[1] = difC ;
	difID[2] = difR ;
}

void SDHCALLayer::buildAsics()
{
	for ( int j = 0 ; j < 3 ; j++ )
	{
		for ( int i = 1 ; i < 49 ; i++ )
		{
			SDHCALAsic* asic = new SDHCALAsic(i , difID[j]) ;
			asic->setLayer(this) ;
			CLHEP::Hep3Vector asicPos = position + 10.408*( (j*32 + jAsicTab[i-1]*8)*CLHEP::Hep3Vector(0 , 1 , 0) + iAsicTab[i-1]*8*CLHEP::Hep3Vector(1 , 0 , 0) ) ;
			asic->setPosition(asicPos) ;
			asic->buildPads() ;
			asics.insert( std::make_pair(difID[j]*100+i , asic) ) ;
		}
	}
}

Asic* SDHCALLayer::findAsic(const CLHEP::Hep3Vector& pos) const
{
	CLHEP::Hep3Vector posInLayer = pos - position ;

	int i = static_cast<int>( posInLayer.x()/83.264 ) ; //8*10.408 = 83.264
	int j = static_cast<int>( posInLayer.y()/83.264 ) ;

	if ( j/4 > 2 )
	{
		std::cout << "Error in SDHCALLayer::findAsic : dif not possible" << std::endl ;
		std::cout << "layer : " << id << std::endl ;
		std::cout << "pos : " << pos << std::endl ;
		std::cout << "posInLayer : " << posInLayer << std::endl ;
		return NULL ;
	}

	int difj = difID[ j/4 ] ;

	AsicMap::const_iterator it = asics.find( 100*difj + asicTab[i][j%4] ) ;
	if ( it == asics.end() )
	{
		std::cout << "Error in SDHCALLayer::findAsic : non existing asic : "  << std::endl ;
		std::cout << "dif : " << difj << " , iAsic : " << i << std::endl ;
		return NULL ;
	}

	return it->second ;
}


const int SDHCALLayer::asicTab[12][4] =
{
	{ 4, 3, 2, 1},
	{ 5, 6, 7, 8},
	{12,11,10, 9},
	{13,14,15,16},
	{20,19,18,17},
	{21,22,23,24},
	{28,27,26,25},
	{29,30,31,32},
	{36,35,34,33},
	{37,38,39,40},
	{44,43,42,41},
	{45,46,47,48},
} ;

const int SDHCALLayer::iAsicTab[48] =
{
	0, 0, 0, 0, 1, 1, 1, 1, 2,
	2, 2, 2, 3, 3, 3, 3, 4, 4, 4,
	4, 5, 5, 5, 5, 6, 6, 6, 6, 7,
	7, 7, 7, 8, 8, 8, 8, 9, 9, 9,
	9,10,10,10,10,11,11,11,11
} ;
const int SDHCALLayer::jAsicTab[48] =
{
	3, 2, 1, 0, 0, 1, 2, 3, 3,
	2, 1, 0, 0, 1, 2, 3, 3, 2, 1,
	0, 0, 1, 2, 3, 3, 2, 1, 0, 0,
	1, 2, 3, 3, 2, 1, 0, 0, 1, 2,
	3, 3, 2, 1, 0, 0, 1, 2, 3
} ;

} //namespace caloobject
