#ifndef CLUSTER_HH
#define CLUSTER_HH

#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloCluster.h"
#include <cmath>
#include <algorithm>

namespace algorithm
{ 

struct clusterParameterSetting
{
		int maxLongitudinal = 0 ;
		int maxTransversal = 1 ;
		bool useDistanceInsteadCellID = false ;
		float maxLongitudinalDistance = 0.0 ;
		float maxTransversalDistance = 11.0 ;
		bool useEnergyToWeightPosition = false ;
} ;

class Clustering
{
	public :
		Clustering()
			: settings()
		{;}
		~Clustering(){;}

		void Run(HitVec& hits, std::vector<caloobject::CaloCluster2D*>& outClusterVec) ;
		inline void SetClusterParameterSetting(clusterParameterSetting val) { settings = val ; }

	private :
		void BuildCluster(HitVec& temp , HitVec& calohit , caloobject::CaloHit*& hit) ;

		HitVec _clusterHitList {} ;

		clusterParameterSetting settings ;
} ;

} //namespace algorithm

#endif
