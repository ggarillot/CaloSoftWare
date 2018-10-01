#ifndef CALOCLUSTER_HH
#define CALOCLUSTER_HH

#include "iostream"
#include "vector"
#include <map>

#include "CLHEP/Vector/ThreeVector.h"

#include "CaloObject/CaloHit.h"

namespace caloobject
{

class Cluster
{
	public :
		Cluster(std::vector<caloobject::CaloHit*>& vec, bool useWeight = 0) ;
		virtual ~Cluster() ;

		inline const std::vector<caloobject::CaloHit*>& getHits() const { return hits ; }
		inline CLHEP::Hep3Vector getPosition() const { return clusterPosition ; }
		inline float getEnergy() const { return energy ; }
		inline float getMaxEnergy() const { return maxEnergy ; }

	protected :
		HitVec hits ;
		CLHEP::Hep3Vector clusterPosition ;
		float energy = 0.f ;
		float maxEnergy = 0.f ;
} ;

class CaloCluster2D : public Cluster
{
	public :
		CaloCluster2D(std::vector<caloobject::CaloHit*>& vec, bool useWeight = 0) ;
		virtual ~CaloCluster2D() ;
		inline int getLayerID() const { return layerID ; }

	private :
		int layerID ;
} ;

class CaloCluster3D : public Cluster
{
	public :
		CaloCluster3D(std::vector<caloobject::CaloHit*>& vec, bool useWeight = 0) ;
		virtual ~CaloCluster3D() ;
		inline std::map<int, CLHEP::Hep3Vector> getBarycenter() const { return barycenter ; }
		void printInfo() ;

	private :
		std::map<int, CLHEP::Hep3Vector> barycenter ;
} ;

}//namespace caloobject

typedef std::vector<caloobject::Cluster*> ClusterVec ;
typedef std::vector<caloobject::CaloCluster2D*> Cluster2DVec ;

#endif
