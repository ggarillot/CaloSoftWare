#ifndef TRACKING_HH
#define TRACKING_HH

#include <cmath>
#include <algorithm>

#include "CLHEP/Vector/ThreeVector.h"

#include "CaloObject/CaloCluster.h"
#include "CaloObject/CaloTrack.h"
#include "Algorithm/PCA.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/LinearFit3D.h"
#include "Algorithm/Distance.h"

namespace algorithm
{ 

struct TrackingParameterSetting
{
		float maxTransverseRatio = 1.0f ;
		float chiSquareLimit = 100.0f ;
		float cosThetaLimit = 0.0f ;
		float maxDistance = 20.0f ;
		int maxDiffBetweenLayer = 3 ;
		unsigned int minNumberOfFiredLayers = 4 ;
		bool printDebug = false ;
} ;

class Tracking
{
	public :
		Tracking()
			: settings()
		{;}
		~Tracking(){;}

		void Run(std::vector<caloobject::CaloCluster2D*>& vec , caloobject::CaloTrack*& track) ;
		void TryToAddAClusterInTrack(caloobject::CaloCluster2D* cluster , caloobject::CaloTrack*& track) ;
		void splitTrack(caloobject::CaloTrack*) ;
		inline void SetTrackingParameterSetting(TrackingParameterSetting params) { settings = params ; }
		inline float getTransverseRatio() { return _transverseRatio ; }

	private :
		TrackingParameterSetting settings ;
		float _transverseRatio = 0.0f ;
} ;


class removeClusterFromTrackIfLayerBiggerThanValue
{
	public:
		//use :
		//track->getClusters().erase( std::remove_if(track->getClusters().begin(), track->getClusters().end(), removeClusterFromTrackIfLayerBiggerThanValue(val) ) , track->getClusters().end() )
		removeClusterFromTrackIfLayerBiggerThanValue(int val)
			: _value(val)
		{}
		bool operator()(caloobject::CaloCluster2D* cluster) { return cluster->getLayerID() > _value ; }
	private:
		int _value ;
} ;

}//namespace algorithm

#endif
