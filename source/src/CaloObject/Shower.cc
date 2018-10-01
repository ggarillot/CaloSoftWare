#include "CaloObject/Shower.h"

#include "Algorithm/PCA.h"
#include "Algorithm/LinearFit3D.h"
#include "Algorithm/Distance.h"
#include "Algorithm/ManipMap.h"

namespace caloobject
{

Shower::Shower()
	: interactionSettings() ,
	  geometrySettings()
{
}

Shower::Shower(const HitVec& vec)
	: interactionSettings() ,
	  geometrySettings() ,
	  hits(vec)
{
}

Shower::Shower(const Cluster2DVec& vec)
	: interactionSettings() ,
	  geometrySettings() ,
	  clusters(vec)
{
	std::map<int , std::tuple<int,double,double,double,double>> dMap ; //for interactingLayers

	int lastLayer = -10 ;
	for ( const auto& cluster : clusters )
	{
		for ( const auto& hit : cluster->getHits() )
		{
			hits.push_back(hit) ;
			hitTimes.push_back( hit->getTime() ) ;

			std::get<0>( dMap[hit->getCellID()[2]] ) ++ ;
			std::get<1>( dMap[hit->getCellID()[2]] ) += hit->getPosition().x() ;
			std::get<2>( dMap[hit->getCellID()[2]] ) += hit->getPosition().x()*hit->getPosition().x() ;
			std::get<3>( dMap[hit->getCellID()[2]] ) += hit->getPosition().y() ;
			std::get<4>( dMap[hit->getCellID()[2]] ) += hit->getPosition().y()*hit->getPosition().y() ;
		}

		lastLayer = std::max(cluster->getLayerID() , lastLayer) ;

		firedLayers.insert( cluster->getLayerID() ) ;
	}

	// compute interacting layers
	for ( const auto& it : dMap )
	{
		const int& n = std::get<0>(it.second) ;
		const double& x = std::get<1>(it.second) ;
		const double& x2 = std::get<2>(it.second) ;
		const double& y = std::get<3>(it.second) ;
		const double& y2 = std::get<4>(it.second) ;

		if ( n <= 5 )
			continue ;

		double varLay = x2/n - x*x/(n*n) + y2/n - y*y/(n*n) ;

		if ( varLay > 100 )
			interactingLayers.insert(it.first) ;

	}

	lastClusterLayer = lastLayer ;
}


Shower::Shower( const caloobject::CaloCluster3D* cluster )
	: interactionSettings() ,
	  geometrySettings() ,
	  hits( cluster->getHits() )
{}

Shower::~Shower()
{}

AnalogShower::~AnalogShower()
{}

AnalogShower::AnalogShower(const HitVec& vec)
	: Shower(vec)
{}

AnalogShower::AnalogShower(const Cluster2DVec& vec)
	: Shower(vec)
{
	edep = 0.0f ;
	rmsEdep = 0.0f ;
	for ( const auto& it : clusters )
		clustersEnergy.push_back( it->getEnergy() ) ;

	for ( const auto& hit : hits )
	{
		edep += hit->getEnergy() ;
		rmsEdep += hit->getEnergy()*hit->getEnergy() ;
		edepPerCell.push_back( hit->getEnergy() ) ;
	}

	meanEdep = edep/getHits().size() ;
	rmsEdep = std::sqrt( rmsEdep/getHits().size() - meanEdep*meanEdep ) ;
}

AnalogShower::AnalogShower(const caloobject::CaloCluster3D* cluster)
	: Shower(cluster)
{}

DigitalShower::DigitalShower(const Cluster2DVec& vec)
	: Shower(vec)
{

	for (const auto& hit : hits)
	{
		nHits.at(0) ++ ;
		unsigned int e = static_cast<unsigned int>( hit->getEnergy() ) ;
		if (e < 4) //1,2 or 3
			nHits.at( e ) ++ ;
	}
}

double AnalogShower::valueOfHitForProfile(const caloobject::CaloHit* hit) const
{
	return hit->getEnergy() ;
}
double DigitalShower::valueOfHitForProfile(const caloobject::CaloHit*) const
{
	return 1.0 ;
}


void Shower::computePCA()
{
	std::vector<double> x , y , z ;
	x.reserve( hits.size() ) ;
	y.reserve( hits.size() ) ;
	z.reserve( hits.size() ) ;

	for ( const auto& hit : hits )
	{
		x.push_back( hit->getPosition().x() ) ;
		y.push_back( hit->getPosition().y() ) ;
		z.push_back( hit->getPosition().z() ) ;
	}

	VectorList vecList ;
	vecList.push_back(x) ;
	vecList.push_back(y) ;
	vecList.push_back(z) ;

	algorithm::PCA* algo_PCA = new algorithm::PCA(vecList) ;
	transverseRatio = static_cast<float>( std::sqrt( algo_PCA->eigenValues()[0]*algo_PCA->eigenValues()[0] + algo_PCA->eigenValues()[1]*algo_PCA->eigenValues()[1] )/algo_PCA->eigenValues()[2] ) ;
	delete algo_PCA ;

	return ;
}

void Shower::computeThrust()
{
	isThrustComputed = false ;
	std::vector<CLHEP::Hep3Vector> hitPos ;
	hitPos.reserve( hits.size() ) ;
	std::vector<int> clSize(hits.size() , 1) ;

	for ( const auto& hit : hits )
		hitPos.push_back( hit->getPosition() ) ;

	algorithm::LinearFit3D* algo_LinearFit3D = new algorithm::LinearFit3D(hitPos,clSize) ;
	thrust.clear() ;
	for(int i = 0 ; i < 4 ; i++)
		thrust.push_back(algo_LinearFit3D->getFitParameters()[i]);

	CLHEP::Hep3Vector orientation = CLHEP::Hep3Vector(-1,0,thrust.at(1)).cross(CLHEP::Hep3Vector(0,-1,thrust.at(3))) ;
	reconstructedCosTheta = static_cast<float>( orientation.cosTheta() ) ;
	eta = static_cast<float>( orientation.eta() ) ;
	phi = static_cast<float>( orientation.phi() ) ;

	delete algo_LinearFit3D ;
	isThrustComputed = true ;
	return ;
}

void Shower::computeInteraction()
{
	if ( !isThrustComputed )
		computeThrust() ;

	firstIntCluster = nullptr ;
	algorithm::InteractionFinder* algo = new algorithm::InteractionFinder() ;
	algo->SetInteractionFinderParameterSetting(interactionSettings) ;

	algo->Run(clusters , thrust) ;
	if( algo->FindInteraction() )
	{
		firstIntCluster = algo->getFirstInteractionCluster() ;
		startingPosition = algo->getFirstInteractionCluster()->getPosition() ;
	}
	isInteractionComputed = true ;
	return ;
}

void Shower::computeProfile()
{
	if ( !isThrustComputed )
		computeThrust() ;

	longitudinalProfile.clear() ;
	transverseProfile.clear() ;
	distanceToAxis.clear() ;

	longitudinalProfile = std::vector<double>( static_cast<unsigned int>(geometrySettings.nLayers) , 0 ) ;
	transverseProfile = std::vector<double>(static_cast<unsigned int>(1000.0f/geometrySettings.pixelSize) , 0 ) ; //1000 is actually a parameter : maximumRadius (need to change)

	algorithm::Distance<caloobject::CaloHit,CLHEP::Hep3Vector> dist ;

	for ( const auto& hit : hits )
	{
		if( hit->getCellID()[2] >= geometrySettings.nLayers )
		{
			std::cout << "Problem in Shower::computeProfile() \n\t => Hit found at " << hit->getCellID()[2] << " while geometrySettings.nLayers = " << geometrySettings.nLayers << std::endl;
			continue ;
		}
		longitudinalProfile.at( hit->getCellID()[2] ) += valueOfHitForProfile(hit) ;

		CLHEP::Hep3Vector vec(thrust.at(0) + thrust.at(1)*hit->getPosition().z() ,
							  thrust.at(2) + thrust.at(3)*hit->getPosition().z() ,
							  hit->getPosition().z() ) ;
		distanceToAxis.push_back( dist.getDistance(hit , vec) ) ;
		int ring = (int)(dist.getDistance( hit,vec )/geometrySettings.pixelSize );
		if( ring > (int)(1000.0f/geometrySettings.pixelSize-1) || ring<0 )
			continue ;
		else
			transverseProfile.at(ring) += valueOfHitForProfile(hit) ;

	}
}

} //namespace caloobject
