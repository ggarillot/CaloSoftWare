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
	Distance<caloobject::CaloCluster2D,caloobject::CaloCluster2D> dist;
	for( std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it )
	{
		if( settings.useAnalogEnergy == false )
		{
			// use it for sdhcal like detector
			if( (*it)->getHits().size() > settings.maximumClusterSizeForMip )
			{
				if( settings.printDebug )
					std::cout << "cluster at " << (*it)->getPosition() << " with too many hits : " << (*it)->getHits().size() << std::endl;

				continue ;
			}
		}
		else{
			// use it for hgcal like detector
			if( (*it)->getEnergy() > settings.maxEnergy )
				continue;
		}
		int nNeighbours=0;
		int nCoreNeighbours=0;
		for( std::vector<caloobject::CaloCluster2D*>::iterator jt=clusters.begin(); jt!=clusters.end(); ++jt ){
			if( (*it)==(*jt) || fabs( (*it)->getPosition().z()-(*jt)->getPosition().z() ) > std::numeric_limits<float>::epsilon() ) continue;
			if( dist.getDistance( (*it),(*jt) ) < settings.transversalDistance ){
				if( settings.printDebug )
					std::cout << "Distance --->>>" << (*it)->getPosition() << "\t" << (*jt)->getPosition() << "\t" << dist.getDistance( (*it),(*jt) ) << std::endl;
				nNeighbours++;
				if( (*jt)->getHits().size() > settings.maximumClusterSizeForMip && settings.useAnalogEnergy==false )
					nCoreNeighbours++;
				else if( (*jt)->getEnergy() > settings.maxEnergy && settings.useAnalogEnergy==true )
					nCoreNeighbours++;
			}
		}
		if( nNeighbours > settings.maximumNumberOfNeighboursForMip &&
			nCoreNeighbours > settings.maximumNumberOfCoreNeighboursForMip )
			continue;
		else
			mipCandidate.push_back(*it);
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
					fabs( (*it)->rhoYVec.at(theta)-(*jt).rho ) < settings.geometry.pixelSize + settings.geometry.pixelSize/10.0 ){
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

		if( TestHoughBinSize(bestBin) )
			break ;

		RemoveIsolatedClusterInHoughBin( bestBin ) ;
		if( TestHoughBinSize( bestBin ) )
		{
			if( settings.printDebug )
				std::cout << "Hough::runHough << DEBUG : Hough bin size is no longer big enough to create a track after removing isloated clusters" << std::endl;
			houghBins.erase(it) ;
			continue ;
		}
		caloobject::CaloTrack* track = NULL ;
		std::vector<caloobject::CaloCluster2D*> temp ;
		for( std::vector<HoughObject*>::iterator jt = bestBin.houghObjects.begin() ; jt != bestBin.houghObjects.end() ; ++jt )
			temp.push_back( (*jt)->cluster ) ;

		algo_Tracking->Run( temp,track ) ;
		if ( track != NULL )
		{
			algo_Tracking->splitTrack(track);
			for(std::vector<HoughObject*>::iterator jt=houghObjects.begin(); jt!=houghObjects.end(); ++jt)
			{
				if( std::find( track->getClusters().begin(), track->getClusters().end(), (*jt)->cluster ) != track->getClusters().end() )
					continue;
				algo_Tracking->TryToAddAClusterInTrack((*jt)->cluster, track) ;
			}
			//			if( track->getClusters().size()<settings.minimumNBins ){
			if( track->getClusters().size() < 4 )
			{
				delete track;
				houghBins.erase(it);
				continue;
			}

			tracks.push_back(track);
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

}

} //namespace algorithm
