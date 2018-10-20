#include <Algorithm/Hough.h>

#include <ctime>

const float PI = 3.1415927f ;

namespace algorithm
{

void Hough::createHoughObjects(std::vector<caloobject::CaloCluster2D*> &clusters)
{
	std::vector<caloobject::CaloCluster2D*> mipCandidate ;
	selectNonDensePart( clusters,mipCandidate ) ;
	for( std::vector<caloobject::CaloCluster2D*>::iterator it = mipCandidate.begin() ; it != mipCandidate.end() ; ++it )
	{
		HoughObject* obj = new HoughObject() ;
		obj->cluster = (*it) ;
		obj->tag = fMip ;
		if( settings.printDebug )
			std::cout << (*it)->getPosition() << std::endl;
		for( unsigned int theta = 0 ; theta < settings.thetaSteps ; theta++ )
		{
			obj->thetas.push_back(theta) ;
			obj->rhoXVec.push_back( static_cast<float>( (*it)->getPosition().z()*std::cos(theta*PI/settings.thetaSteps) + (*it)->getPosition().x()*std::sin(theta*PI/settings.thetaSteps) ) ) ;
			obj->rhoYVec.push_back( static_cast<float>( (*it)->getPosition().z()*std::cos(theta*PI/settings.thetaSteps) + (*it)->getPosition().y()*std::sin(theta*PI/settings.thetaSteps) ) ) ;
		}
		houghObjects.push_back(obj) ;
	}
}

void Hough::selectNonDensePart( std::vector<caloobject::CaloCluster2D*> &clusters, std::vector<caloobject::CaloCluster2D*> &mipCandidate )
{
	for ( const auto& cluster : clusters )
	{
		if( cluster->getHits().size() > settings.maximumClusterSizeForMip )
			continue ;

		int neighbour=0;
		int big_neighbour=0;

		float X = static_cast<float>( cluster->getPosition().x() ) ;
		float Y = static_cast<float>( cluster->getPosition().y() ) ;
		int Z = cluster->getLayerID() ;

		for ( const auto& cluster2 : clusters )
		{
			if( cluster ==cluster2 )
				continue ;

			float x = static_cast<float>( cluster2->getPosition().x() ) ;
			float y = static_cast<float>( cluster2->getPosition().y() ) ;
			int z = cluster2->getLayerID() ;

			if ( std::abs(X-x) < 50 && std::abs(Y-y) < 50 && Z==z )
				neighbour++ ;
			if( std::abs(X-x) < 50 && std::abs(Y-y) < 50 && Z==z && cluster2->getHits().size() > 5 )
				big_neighbour++ ;
		}

		if( neighbour <= settings.maximumNumberOfNeighboursForMip || big_neighbour <= settings.maximumNumberOfCoreNeighboursForMip )
			mipCandidate.push_back(cluster);
	}
}


std::vector< HoughBin > Hough::getHoughBinsFromZX()
{
	std::vector< HoughBin > outputBins;
	for( std::vector<HoughObject*>::iterator it = houghObjects.begin() ; it != houghObjects.end() ; ++it )
	{
		for( unsigned int theta = 0 ; theta < settings.thetaSteps ; theta++ )
		{
			std::vector< HoughBin >::iterator jt ;
			for( jt = outputBins.begin() ; jt != outputBins.end() ; ++jt )
			{
				if( (*it)->thetas.at(theta) == (*jt).theta &&
						fabs( (*it)->rhoXVec.at(theta)-(*jt).rho ) < settings.geometry.pixelSize + settings.geometry.pixelSize*settings.rhoTolerance )
				{
					(*jt).houghObjects.push_back(*it) ;
					break ;
				}
			}
			if( jt != outputBins.end() )
				continue ;
			else
			{
				HoughBin bin ;
				bin.theta = theta ;
				bin.rho = (*it)->rhoXVec.at(theta) ;
				bin.houghObjects.push_back(*it) ;
				outputBins.push_back(bin) ;
			}
		}
	}
	for( std::vector< HoughBin >::iterator it=outputBins.begin(); it!=outputBins.end(); ++it)
	{
		(*it).rmTag = false ;
		if( TestHoughBinSize(*it) )
			(*it).rmTag = true ;
	}
	outputBins.erase(std::remove_if(outputBins.begin(),outputBins.end(),HasToBeDeleted),outputBins.end());
	std::sort(outputBins.begin(), outputBins.end(),SortHoughBinBySize::Sort);

	return outputBins;
}

HoughBin Hough::getBestHoughBinFromZY( HoughBin& inputBin )
{
	std::vector< HoughBin > outputBins;
	for( std::vector< HoughObject*  >::iterator it=inputBin.houghObjects.begin(); it!=inputBin.houghObjects.end(); ++it ){
		for( unsigned int theta=0; theta<settings.thetaSteps; theta++ ){
			std::vector< HoughBin >::iterator jt;
			for( jt=outputBins.begin(); jt!=outputBins.end(); ++jt ){
				if( (*it)->thetas.at(theta)==(*jt).theta &&
						fabs( (*it)->rhoYVec.at(theta)-(*jt).rho ) < settings.geometry.pixelSize + settings.geometry.pixelSize*settings.rhoTolerance )
				{
					(*jt).houghObjects.push_back(*it);
					break;
				}
			}
			if( jt!=outputBins.end() ) continue;
			else{
				HoughBin bin;
				bin.theta=theta;
				bin.rho=(*it)->rhoYVec.at(theta);
				bin.houghObjects.push_back(*it);
				outputBins.push_back(bin);
			}
		}
	}
	std::sort(outputBins.begin(), outputBins.end(),SortHoughBinBySize::Sort);
	return ( *outputBins.begin() );
}

void Hough::RemoveIsolatedClusterInHoughBin(HoughBin &a)
{
	//Distance<caloobject::CaloCluster2D,caloobject::CaloCluster2D> dist;
	for( std::vector<HoughObject*>::iterator it=a.houghObjects.begin(); it!=a.houghObjects.end(); ++it ){
		std::vector<HoughObject*>::iterator jt;
		for( jt=a.houghObjects.begin(); jt!=a.houghObjects.end(); ++jt ){
			if( it!=jt &&
					std::abs( (*it)->cluster->getLayerID()-(*jt)->cluster->getLayerID() ) <= settings.isolationDistance )
				break;
		}
		if( jt==a.houghObjects.end() ){
			if( settings.printDebug )
				std::cout << "Hough::RemoveIsolatedClusterInHoughBin << DEBUG : Find one isloated cluster at : " << (*it)->cluster->getPosition() << " " << (*it)->cluster << std::endl;
			a.houghObjects.erase(it);
			--it;
		}
	}
}

void Hough::RemoveTrackedObjects(std::vector<HoughBin> &houghBins)
{
	for(std::vector<HoughBin>::iterator it=houghBins.begin(); it!=houghBins.end(); ++it)
		(*it).houghObjects.erase( std::remove_if( (*it).houghObjects.begin(),(*it).houghObjects.end(), RemoveTrackedObject::rm ),(*it).houghObjects.end() );

}

void Hough::runHough(std::vector<caloobject::CaloCluster2D*> &clusters, std::vector<caloobject::CaloTrack*> &tracks, algorithm::Tracking *algo_Tracking)
{

	if( NULL==algo_Tracking )
	{
		std::cout << "ERROR : an algorithm::Tracking must be initialised before calling Hough::runHough => return" << std::endl;
		return;
	}

	createHoughObjects(clusters) ;

	std::vector< HoughBin > houghBins = getHoughBinsFromZX() ;

	while(1)
	{
		if (houghBins.empty() )
			break ;
		std::vector< HoughBin >::iterator it = houghBins.begin() ;
		HoughBin bestBin = getBestHoughBinFromZY( *it ) ;

		RemoveIsolatedClusterInHoughBin( bestBin ) ;
//		if( TestHoughBinSize( bestBin ) )
		if( bestBin.houghObjects.size() < 7 )
		{
			houghBins.erase(it) ;
			continue ;
		}
		caloobject::CaloTrack* track = nullptr ;
		std::vector<caloobject::CaloCluster2D*> temp ;
		for( std::vector<HoughObject*>::iterator jt = bestBin.houghObjects.begin() ; jt != bestBin.houghObjects.end() ; ++jt )
			temp.push_back( (*jt)->cluster ) ;

		algo_Tracking->Run( temp,track ) ;
		if ( track != nullptr )
		{
			for(std::vector<HoughObject*>::iterator jt=houghObjects.begin(); jt!=houghObjects.end(); ++jt)
			{
				if( std::find( track->getClusters().begin(), track->getClusters().end(), (*jt)->cluster ) != track->getClusters().end() )
					continue;
				algo_Tracking->TryToAddAClusterInTrack((*jt)->cluster, track) ;
			}

			if( track->getClusters().size() <= 4 )
			{
				delete track;
				houghBins.erase(it);
				continue;
			}
			algo_Tracking->splitTrack(track) ;
			tracks.push_back(track) ;
			for(std::vector<caloobject::CaloCluster2D*>::const_iterator jt = track->getClusters().begin() ; jt != track->getClusters().end() ; ++jt)
			{
				for( std::vector< HoughObject* >::iterator kt = houghObjects.begin() ; kt != houghObjects.end() ; ++kt)
				{
					if( (*kt)->cluster == (*jt) )
					{
						(*kt)->tag = fTrack ;
						break ;
					}
				}
			}
			RemoveTrackedObjects( houghBins ) ;
		}
		else
			houghBins.erase(it) ;

		for( std::vector< HoughBin >::iterator jt = houghBins.begin() ; jt != houghBins.end() ; ++jt )
		{
			if( TestHoughBinSize(*jt) )
			{
				houghBins.erase(jt) ;
				jt-- ;
			}
		}

		if( houghBins.size() == 0 )
			break ;

		std::sort(houghBins.begin(), houghBins.end(),SortHoughBinBySize::Sort) ;
	}

	houghBins.clear() ;

	for(std::vector< HoughObject* >::iterator it = houghObjects.begin() ; it != houghObjects.end() ; ++it)
		delete (*it) ;

	houghObjects.clear() ;

	for ( auto it = tracks.begin() ; it != tracks.end() ; ++it)
	{
		if( (*it)->getClusters().size()<4 )
		{
			tracks.erase(it);
			it--;
		}
	}

}

} //namespace algorithm
