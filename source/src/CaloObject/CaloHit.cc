#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloCluster.h"

namespace caloobject
{

CaloHit::CaloHit(int* cellID, CLHEP::Hep3Vector vec, float energy, float time, CLHEP::Hep3Vector posShift)
	: _cellID() ,
	  _hitPosition(vec-posShift) ,
	  _energy(energy) ,
	  _time(time)
{
	for(int i = 0 ; i < 3 ; i++)
		_cellID[i] = cellID[i] ;
}

}
