#include "Algorithm/Density.h"
#include "Algorithm/Distance.h"

namespace algorithm
{

Density::Density()
	: densityFunc(new DensityDistanceFunction3By3) ,
	  densityPerHit()
{

}

Density::~Density()
{
	delete densityFunc ;
}

float Density::compute(HitVec& hitVec)
{
	densityPerHit.clear() ;

	for( HitVec::const_iterator it = hitVec.begin() ; it != hitVec.end() ; ++it )
	{
		for( HitVec::const_iterator jt = it ; jt != hitVec.end() ; ++jt )
		{
			float d = densityFunc->eval(*it , *jt) ;
			densityPerHit[*it] += d ;

			if ( it != jt )
				densityPerHit[*jt] += d ;
		}

	}

	float density = 0 ;
	for ( std::map<caloobject::CaloHit* , float>::const_iterator it = densityPerHit.begin() ; it != densityPerHit.end() ; ++it )
		density += it->second ;

	density /= densityPerHit.size() ;
	return density ;
}


DensityDistanceFunction::DensityDistanceFunction()
{}

DensityDistanceFunction::~DensityDistanceFunction()
{}

DensityDistanceFunction3By3::DensityDistanceFunction3By3()
	: DensityDistanceFunction()
{}

DensityDistanceFunction3By3::~DensityDistanceFunction3By3()
{}

float DensityDistanceFunction3By3::eval(caloobject::CaloHit* hitA , caloobject::CaloHit* hitB) const
{
	if ( std::abs(hitA->getCellID()[0] - hitB->getCellID()[0]) <= 1 && std::abs(hitA->getCellID()[1] - hitB->getCellID()[1]) <= 1 && hitA->getCellID()[2] == hitB->getCellID()[2]
		 && Distance<caloobject::CaloHit,caloobject::CaloHit>().getDistance(hitA,hitB)<20 )
		return 1.0f ;
	else
		return 0.0f ;
}




} //namespace algorithm


