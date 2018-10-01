#include "Algorithm/InteractionFinder.h"

#include <algorithm>
#include <set>

namespace algorithm
{
void InteractionFinder::Run(const std::vector<caloobject::CaloCluster2D*>& clusters, const std::vector<float> params)
{
	findInteraction = false ;
	firstCluster = nullptr ;
	Distance<caloobject::CaloCluster2D, CLHEP::Hep3Vector> dist ;
	CLHEP::Hep3Vector vec;
	for(std::vector<caloobject::CaloCluster2D*>::const_iterator it=clusters.begin(); it!=clusters.end(); ++it){

		vec = CLHEP::Hep3Vector(params[0]+params[1]*(*it)->getPosition().z(),
				params[2]+params[3]*(*it)->getPosition().z(),
				(*it)->getPosition().z());
		if( settings.printDebug )
			std::cout << "size = " << (*it)->getHits().size() << "\t"
					  << "layer = " << (*it)->getLayerID() << "\t"
					  << "distance = " << dist.getDistance( (*it),vec ) << "\t"
					  << "energy = " << (*it)->getEnergy()*1000 << std::endl;

		if( settings.useAnalogEnergy==false )
		{
			if( (*it)->getHits().size()<(unsigned int)settings.minSize || dist.getDistance( (*it),vec ) > settings.maxRadius )
				continue ;
		}
		else
		{
			if( (*it)->getEnergy()<settings.minEnergy || dist.getDistance( (*it),vec ) > settings.maxRadius )
				continue ;
		}
		int counter=0;

		for(std::vector<caloobject::CaloCluster2D*>::const_iterator jt=it+1; jt!=clusters.end(); ++jt){
			vec = CLHEP::Hep3Vector(params[0]+params[1]*(*jt)->getPosition().z(),
					params[2]+params[3]*(*jt)->getPosition().z(),
					(*jt)->getPosition().z());

			if( settings.useAnalogEnergy==false )
			{
				if( (*jt)->getHits().size() >= (unsigned int)settings.minSize &&
					(*jt)->getPosition().z()-(*it)->getPosition().z() > std::numeric_limits<float>::epsilon() &&
					(*jt)->getLayerID()-(*it)->getLayerID() < settings.maxDepth&&
					dist.getDistance( (*jt),vec ) < settings.maxRadius )
					counter++;
			}
			else
			{
				if( (*jt)->getEnergy() >= settings.minEnergy &&
					(*jt)->getPosition().z()-(*it)->getPosition().z() > std::numeric_limits<float>::epsilon() &&
					(*jt)->getLayerID()-(*it)->getLayerID() < settings.maxDepth&&
					dist.getDistance( (*jt),vec ) < settings.maxRadius )
					counter++;
			}

			if( settings.printDebug )
				std::cout << "\t size = " << (*jt)->getHits().size() << "\t"
						  << "\t layer = " << (*jt)->getLayerID() << "\t"
						  << "\t distance = " << dist.getDistance( (*jt),vec ) << "\t"
						  << "\t energy = " << (*jt)->getEnergy()*1000 << "\t"
						  << "\t counter = " << counter << std::endl;
		}
		if(counter>=settings.minNumberOfCluster)
		{
			firstCluster = (*it);
			findInteraction = true;
			break;
		}
	}
}

void InteractionFinder::Run(const std::vector<caloobject::CaloCluster2D*>& clusters, std::vector<caloobject::CaloTrack*> tracks, const std::vector<float> params)
{
	Run(clusters , params) ;

	if ( firstCluster )
		return ;

	int end = -10 ;
	std::map<int , HitVec> layers ;

	Distance<caloobject::CaloCluster2D, CLHEP::Hep3Vector> dist ;
	Distance<caloobject::CaloHit, CLHEP::Hep3Vector> disth ;
	for (const auto& cluster : clusters)
	{
		CLHEP::Hep3Vector layerThrustPos(params[0] + params[1]*cluster->getPosition().z() ,
				params[2] + params[3]*cluster->getPosition().z() ,
				cluster->getPosition().z()) ;

		if ( dist.getDistance(cluster , layerThrustPos) < settings.maxRadius )
		{
			HitVec& hitV = layers[cluster->getLayerID()] ;
			hitV.insert(hitV.end() , cluster->getHits().begin() , cluster->getHits().end() ) ;
		}
		end = std::max( end , cluster->getLayerID() ) ;
	}

	//to avoid dealing with muons
	if ( !(end < settings.maxLastFiredLayer || layers.size() < settings.maxFiredLayers) )
		return ;

	auto sortTrackByFirstClusterLayer = [](std::vector<caloobject::CaloTrack*>::value_type a , std::vector<caloobject::CaloTrack*>::value_type b)
			-> bool { return a->getTrackStartingCluster()->getLayerID() < b->getTrackStartingCluster()->getLayerID() ; } ;
	std::sort(tracks.begin() , tracks.end() , sortTrackByFirstClusterLayer) ;

	int intLayer = -10 ;
	if ( !tracks.empty() )
	{
		firstCluster = tracks.at(0)->getTrackLastCluster() ;
		intLayer = firstCluster->getLayerID() ;
		findInteraction = true ;

		if ( tracks.size() > 1 )
		{
			if ( tracks.at(1)->getTrackStartingCluster()->getLayerID() < intLayer )
			{
				firstCluster = tracks.at(1)->getTrackStartingCluster() ;
				intLayer = firstCluster->getLayerID() ;
				findInteraction = true ;
			}
		}
	}

	for ( std::map<int , HitVec>::const_iterator it = layers.begin() ; it != layers.end() ; ++it )
	{
		if ( it->first >= intLayer )
			break ;

		if ( it->second.size() >= 5 )
		{
			intLayer = it->first ;
			double minDist = std::numeric_limits<double>::max() ;
			for ( const auto& hit : it->second )
			{
				CLHEP::Hep3Vector layerThrustPos(params[0] + params[1]*hit->getPosition().z() ,
						params[2] + params[3]*hit->getPosition().z() ,
						hit->getPosition().z()) ;

				if ( disth.getDistance(hit,layerThrustPos) < minDist )
				{
					firstCluster = dynamic_cast<caloobject::CaloCluster2D*>( hit->getCluster() ) ;
					findInteraction = true ;
				}
			}
			break ;
		}
	}

	if ( findInteraction )
		return ;

	int lastClusterZ = -10 ;
	for ( const auto& cluster : clusters )
	{
		if ( cluster->getLayerID() > layers.size() )
			continue ;

		if ( cluster->getLayerID() > lastClusterZ )
		{
			lastClusterZ = cluster->getLayerID() ;
			firstCluster = cluster ;
			findInteraction = true ;
		}
	}
}

}//namespace algorithm
