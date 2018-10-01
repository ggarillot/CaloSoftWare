#ifndef CLUSTER3D_HH
#define CLUSTER3D_HH

#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloCluster.h"
#include "CaloObject/CaloGeom.h"
#include <cmath>
#include <algorithm>

namespace algorithm
{ 

struct cluster3DParameterSetting
{
		int maxLongitudinal = 3 ;
		float maxTransverseDistance = 30.0f ;
		bool useEnergyToWeightPosition = false ;
		bool printDebug = false ;
		caloobject::GeomParameterSetting geometry ;

		cluster3DParameterSetting()
			: geometry()
		{;}
} ;

class Cluster3D
{
	public:
		Cluster3D()
			: _clusterHitList() , settings()
		{;}
		~Cluster3D(){;}

		void Run(HitVec &hits, std::vector<caloobject::CaloCluster3D*> &outClusterVec) ;
		void Run(HitVec &hits, std::vector<caloobject::CaloCluster3D*> &outClusterVec, CLHEP::Hep3Vector seed, CLHEP::Hep3Vector direction) ;
		inline void SetCluster3DParameterSetting(cluster3DParameterSetting val) { settings = val ; }

	private:
		void BuildCluster(HitVec &temp,
						  HitVec &calohit,
						  caloobject::CaloHit* &hit,
						  CLHEP::Hep3Vector seed,
						  CLHEP::Hep3Vector direction) ;

		HitVec _clusterHitList ;

		cluster3DParameterSetting settings ;
} ;

} //namespace algorithm

#endif
