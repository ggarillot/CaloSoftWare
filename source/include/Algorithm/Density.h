#ifndef Density_h
#define Density_h

#include "CaloObject/CaloHit.h"

#include <map>
#include <cmath>
#include <algorithm>

namespace algorithm
{

class DensityDistanceFunction
{
	public :
		DensityDistanceFunction() ;
		virtual ~DensityDistanceFunction() ;

		virtual float eval(caloobject::CaloHit* hitA , caloobject::CaloHit* hitB) const = 0 ;

	protected :

} ;

class DensityDistanceFunction3By3 : public DensityDistanceFunction
{
	public :
		DensityDistanceFunction3By3() ;
		virtual ~DensityDistanceFunction3By3() ;

		virtual float eval(caloobject::CaloHit* hitA , caloobject::CaloHit* hitB) const ;

	protected :

} ;


class Density
{
	public :
		Density() ;
		~Density() ;

		float compute(HitVec& hitVec) ;

		const std::map<caloobject::CaloHit* , float>& getDensityPerHit() const { return densityPerHit ; }

		Density(const Density&) = delete ;
		void operator=(const Density&) = delete ;


	protected :
		DensityDistanceFunction* densityFunc ;

		std::map<caloobject::CaloHit* , float> densityPerHit ;
} ;


} //namespace algorithm

#endif //Density_h
