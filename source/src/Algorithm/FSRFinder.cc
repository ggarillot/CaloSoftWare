#include "Algorithm/FSRFinder.h"

#include "Algorithm/ClusteringHelper.h"
#include <algorithm>

namespace algorithm
{

void FSRFinder::Run(std::vector<caloobject::CaloHit*> hits, caloobject::CaloTrack* track, caloobject::CaloCluster3D* gammaCluster)
{
	rmTrackedHitsFromList(hits,track);

	algorithm::Cluster3D *algo_Cluster3D=new algorithm::Cluster3D();
	std::vector<caloobject::CaloCluster3D*> outClusterVec;
	algo_Cluster3D->SetCluster3DParameterSetting( settings.clustering );
	algo_Cluster3D->Run( hits,
						 outClusterVec,
						 track->expectedTrackProjection( settings.geometry.firstLayerZ ),
						 track->orientationVector() );
	delete algo_Cluster3D;
	if( outClusterVec.size()==0 )
		return;
	algorithm::SortClusterBySize<caloobject::CaloCluster3D,caloobject::CaloCluster3D> sorter;
	std::sort( outClusterVec.begin(), outClusterVec.end(), sorter.sort );
	caloobject::CaloCluster3D* fsrCandidate=(*outClusterVec.begin());
	if( settings.printDebug ) std::cout << "gammaCluster " << gammaCluster << std::endl;
	if( outClusterVec.size()>1 )
		for( std::vector<caloobject::CaloCluster3D*>::iterator it=outClusterVec.begin()+1; it!=outClusterVec.end(); ++it){
			if( settings.printDebug ) std::cout << "(*it) " << (*it) << std::endl;
			delete (*it);
		}
	if( isGammaCluster( fsrCandidate ) )
		gammaCluster=fsrCandidate;
}

void FSRFinder::rmTrackedHitsFromList(HitVec& hits, caloobject::CaloTrack* track)
{
	std::vector<caloobject::CaloHit*> trackedHits;
	for(std::vector<caloobject::CaloCluster2D*>::const_iterator it=track->getClusters().begin(); it!=track->getClusters().end(); ++it)
		trackedHits.insert( trackedHits.begin(), (*it)->getHits().begin(), (*it)->getHits().end() );

	for(std::vector<caloobject::CaloHit*>::iterator it=hits.begin(); it!=hits.end(); ++it)
		if( std::find( trackedHits.begin(), trackedHits.end(), (*it) )!=trackedHits.end() ){
			hits.erase(it);
			it--;
		}
}

bool FSRFinder::isGammaCluster( caloobject::CaloCluster3D* fsrCandidate )
{
	return false ;
}

}
