#include "Algorithm/Tracking.h"
#include "set"
namespace algorithm
{
void Tracking::Run(std::vector<caloobject::CaloCluster2D*> &vec , caloobject::CaloTrack* &track)
{
	track = nullptr ;
	if( vec.size() == 0 )
		return ;
	_transverseRatio=0.0;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	std::vector<CLHEP::Hep3Vector> clusterPos;
	std::vector<int> clusterSize;
	std::set<int> nLayerSet;
	for(std::vector<caloobject::CaloCluster2D*>::iterator it=vec.begin(); it!=vec.end(); ++it)
	{
		clusterPos.push_back( (*it)->getPosition() );
		clusterSize.push_back( (*it)->getHits().size() );
		x.push_back( (*it)->getPosition().x() );
		y.push_back( (*it)->getPosition().y() );
		z.push_back( (*it)->getPosition().z() );
		nLayerSet.insert( (*it)->getLayerID() );
	}

	if( nLayerSet.size() < settings.minNumberOfFiredLayers )
		return ;

	std::vector< std::vector<double> > pcaVec;
	pcaVec.push_back(x);
	pcaVec.push_back(y);
	pcaVec.push_back(z);
	algorithm::PCA* pca = new algorithm::PCA(pcaVec);
	std::vector<double> eigenValues = pca->eigenValues();
	_transverseRatio = sqrt(eigenValues.at(0)*eigenValues.at(0) + eigenValues.at(1)*eigenValues.at(1)) / eigenValues.at(2);
	delete pca ;
	if( _transverseRatio > settings.maxTransverseRatio ){
		if( settings.printDebug ) std::cout << "Tracking::Run << DEBUG : _transverseRatio = " << _transverseRatio << std::endl;
		return;
	}
	algorithm::LinearFit3D* fit = new algorithm::LinearFit3D(clusterPos,clusterSize);
	if( fit->getChi2() > settings.chiSquareLimit ){
		if( settings.printDebug ) std::cout << "Tracking::Run << DEBUG : chiSquare = " << fit->getChi2() << std::endl;
		delete fit ;
		return ;
	}
	if( fit->getFitParameters()[0]!=fit->getFitParameters()[0] ||
		fit->getFitParameters()[1]!=fit->getFitParameters()[1] ||
		fit->getFitParameters()[2]!=fit->getFitParameters()[2] ||
		fit->getFitParameters()[3]!=fit->getFitParameters()[3] ){
		if( settings.printDebug )
			std::cout << "Tracking::Run << DEBUG : Problem with fit : fit->getFitParameters() = "
					  << fit->getFitParameters()[0] << "\t"
					  << fit->getFitParameters()[1] << "\t"
					  << fit->getFitParameters()[2] << "\t"
					  << fit->getFitParameters()[3] << std::endl;
		delete fit;
		return;
	}
	CLHEP::Hep3Vector nx(-1.,0.,fit->getFitParameters()[1]);
	CLHEP::Hep3Vector ny(0.,-1.,fit->getFitParameters()[3]);
	if( (nx.cross(ny)).cosTheta() < settings.cosThetaLimit ){
		if( settings.printDebug )
			std::cout << "Tracking::Run << DEBUG : Track angle is too high : cos(theta) = "
					  << (nx.cross(ny)).cosTheta()
					  << std::endl;
		delete fit;
		return;
	}

	track=new caloobject::CaloTrack(vec);
	track->setChi2(fit->getChi2());
	track->setTrackParameters(fit->getFitParameters());
	std::sort( track->getClusters().begin(), track->getClusters().end(), algorithm::ClusteringHelper::SortClusterByLayer );
	if( settings.printDebug )
		std::cout << "Tracking::Run << DEBUG : a track has been created at  = "
				  << track
				  << std::endl;
	delete fit;
}

void Tracking::TryToAddAClusterInTrack(caloobject::CaloCluster2D* cluster, caloobject::CaloTrack* &track)
{
	if( track->getTrackStartingCluster()->getLayerID()-cluster->getLayerID() > settings.maxDiffBetweenLayer ||
		cluster->getLayerID()-track->getTrackLastCluster()->getLayerID() > settings.maxDiffBetweenLayer )
		return;

	Distance<caloobject::CaloCluster2D,caloobject::CaloTrack> dist;
	if( dist.getDistanceInLayer(cluster,track) > settings.maxDistance )
		return;

	std::vector<caloobject::CaloCluster2D*> clusters=track->getClusters();
	if( clusters.size()==0 ) return;
	clusters.push_back(cluster);
	_transverseRatio=0.0;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	std::vector<CLHEP::Hep3Vector> clusterPos;
	std::vector<int> clusterSize;
	for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
		clusterPos.push_back( (*it)->getPosition() );
		clusterSize.push_back( (*it)->getHits().size() );
		x.push_back( (*it)->getPosition().x() );
		y.push_back( (*it)->getPosition().y() );
		z.push_back( (*it)->getPosition().z() );
	}
	std::vector< std::vector<double> > pcaVec;
	pcaVec.push_back(x);
	pcaVec.push_back(y);
	pcaVec.push_back(z);
	algorithm::PCA* pca=new algorithm::PCA(pcaVec);
	std::vector<double> eigenValues=pca->eigenValues();
	_transverseRatio = sqrt(eigenValues.at(0)*eigenValues.at(0) +
							eigenValues.at(1)*eigenValues.at(1)) / eigenValues.at(2);
	delete pca;
	if( _transverseRatio > settings.maxTransverseRatio ){
		if( settings.printDebug ) std::cout << "Tracking::Run << DEBUG : _transverseRatio = " << _transverseRatio << std::endl;
		return;
	}
	algorithm::LinearFit3D *fit=new algorithm::LinearFit3D(clusterPos,clusterSize);
	if( fit->getChi2() > settings.chiSquareLimit ){
		if( settings.printDebug ) std::cout << "Tracking::Run << DEBUG : chiSquare = " << fit->getChi2() << std::endl;
		delete fit;
		return;
	}
	if( fit->getFitParameters()[0]!=fit->getFitParameters()[0] ||
		fit->getFitParameters()[1]!=fit->getFitParameters()[1] ||
		fit->getFitParameters()[2]!=fit->getFitParameters()[2] ||
		fit->getFitParameters()[3]!=fit->getFitParameters()[3] ){
		if( settings.printDebug )
			std::cout << "Tracking::Run << DEBUG : Problem with fit : fit->getFitParameters() = "
					  << fit->getFitParameters()[0] << "\t"
					  << fit->getFitParameters()[1] << "\t"
					  << fit->getFitParameters()[2] << "\t"
					  << fit->getFitParameters()[3] << std::endl;
		delete fit;
		return;
	}
	CLHEP::Hep3Vector nx(-1.,0.,fit->getFitParameters()[1]);
	CLHEP::Hep3Vector ny(0.,-1.,fit->getFitParameters()[3]);
	if( (nx.cross(ny)).cosTheta() < settings.cosThetaLimit ){
		if( settings.printDebug )
			std::cout << "Tracking::Run << DEBUG : Track angle is too high : cos(theta) = "
					  << (nx.cross(ny)).cosTheta()
					  << std::endl;
		delete fit;
		return;
	}
	track->addCluster(cluster);
	track->setChi2(fit->getChi2());
	track->setTrackParameters(fit->getFitParameters());
	std::sort( track->getClusters().begin(), track->getClusters().end(), algorithm::ClusteringHelper::SortClusterByLayer );
	delete fit;
}

void Tracking::splitTrack(caloobject::CaloTrack* track)
{
	int ref=track->getTrackStartingCluster()->getLayerID();
	for(std::vector<caloobject::CaloCluster2D*>::const_iterator it=track->getClusters().begin(); it!=track->getClusters().end(); ++it){
		if( (*it)->getLayerID()-ref>settings.maxDiffBetweenLayer )
			break;
		else
			ref=(*it)->getLayerID();
	}

	if( ref != track->getTrackLastCluster()->getLayerID() ){
		if( settings.printDebug )
			std::cout << "Tracking::Run << DEBUG : I split the track, hum, I split the track, hum... at layer " << track->getTrackLastCluster()->getLayerID() << std::endl;
		track->getClusters().erase( std::remove_if( track->getClusters().begin(),
													track->getClusters().end(),
													removeClusterFromTrackIfLayerBiggerThanValue(ref) ),
									track->getClusters().end() );
	}
}
}
