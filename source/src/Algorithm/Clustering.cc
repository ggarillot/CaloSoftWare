#include "Algorithm/Clustering.h"

namespace algorithm
{

void Clustering::Run(HitVec& hits , std::vector<caloobject::CaloCluster2D*>& outClusterVec)
{
	HitVec temp ;
	for( caloobject::CaloHit* it : hits )
	{
		if ( std::find(temp.begin() , temp.end() , it) != temp.end() )
			continue ;
		_clusterHitList.push_back(it) ;
		temp.push_back(it) ;
		BuildCluster(temp,hits,it) ;
		caloobject::CaloCluster2D* calocluster = new caloobject::CaloCluster2D(_clusterHitList, settings.useEnergyToWeightPosition ) ;
		outClusterVec.push_back(calocluster) ;
		_clusterHitList.clear() ;
	}
}

void Clustering::BuildCluster(HitVec& temp, HitVec& hits, caloobject::CaloHit*& hit)
{
	for( caloobject::CaloHit* it : hits )
	{
		if ( std::find( temp.begin() , temp.end() , it ) != temp.end() )
			continue ;

		if ( !settings.useDistanceInsteadCellID )
		{
			if( std::abs( (hit)->getCellID()[0]-it->getCellID()[0] ) + std::abs( (hit)->getCellID()[1]-it->getCellID()[1] ) <= settings.maxTransversal &&
				std::abs( (hit)->getCellID()[2]-it->getCellID()[2] ) <= settings.maxLongitudinal )
			{
				_clusterHitList.push_back(it) ;
				temp.push_back(it);
				this->BuildCluster(temp,hits,it);
			}
		}
		else
		{
			double deltaX = std::abs( (hit)->getPosition().x()-it->getPosition().x() ) ;
			double deltaY = std::abs( (hit)->getPosition().y()-it->getPosition().y() ) ;

			if( deltaX*deltaX + deltaY*deltaY <= settings.maxTransversalDistance*settings.maxTransversalDistance &&
				std::abs( (hit)->getPosition().z()-it->getPosition().z() ) <= settings.maxLongitudinalDistance )
			{
				_clusterHitList.push_back(it);
				temp.push_back(it);
				this->BuildCluster(temp,hits, it);
			}
		}
	}
}

}
