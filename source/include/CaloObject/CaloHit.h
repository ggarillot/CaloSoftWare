#ifndef CALOHIT_HH
#define CALOHIT_HH

#include <iostream>
#include <vector>
#include "CLHEP/Vector/ThreeVector.h"


namespace caloobject
{ 

class Cluster ;
class CaloHit
{
	public :
		CaloHit(int* cellID, CLHEP::Hep3Vector vec, float energy, float time, CLHEP::Hep3Vector posShift) ;
		// algorithms assume that the zero is located at the middle of the first layer.

		~CaloHit(){ ; }

		inline const CLHEP::Hep3Vector getPosition() const { return _hitPosition ; }
		inline const int* getCellID() const { return _cellID ; }
		inline float getEnergy() const { return _energy ; }
		inline float getTime() const { return _time ; }

		Cluster* getCluster() const { return cluster ; }
		void setCluster(Cluster* _cluster) { cluster = _cluster ; }

		CaloHit(const CaloHit &toCopy) = delete ;
		void operator=(const CaloHit &toCopy) = delete ;

	private :
		int _cellID[3] ;
		CLHEP::Hep3Vector _hitPosition ;
		Cluster* cluster = nullptr ;
		float _energy ;
		float _time ;
} ;

} //namespace caloobject

typedef std::vector<caloobject::CaloHit*> HitVec ;

#endif
