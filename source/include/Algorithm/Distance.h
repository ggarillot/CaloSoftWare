#ifndef DISTANCE_HH
#define DISTANCE_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloCluster.h"
#include "CaloObject/CaloTrack.h"

#include "CLHEP/Vector/ThreeVector.h"

namespace algorithm
{

template <typename T, typename S>
class Distance
{
	public :
		Distance(){;}
		~Distance(){;}

		float getDistance(T* t, S* s){ return (t->getPosition()-s->getPosition()).mag(); }
};

template <typename T>
class Distance<T, CLHEP::Hep3Vector>
{
	public :
		Distance(){;}
		~Distance(){;}

		float getDistance(T* t, CLHEP::Hep3Vector vec){ return (t->getPosition()-vec).mag(); }
};

template <typename S>
class Distance<S, caloobject::CaloTrack>
{
	public :
		Distance(){;}
		~Distance(){;}

		float getDistance(S* s,caloobject::CaloTrack *t)
		{
			/*
	  point H(x,y,z) (S* s)
	  track T : x_t = trackparam[0] + trackparam[1]*z_t => plan equation; normal vector Nx(-1,0,trackparam[1])
	  track T : y_t = trackparam[2] + trackparam[3]*z_t => plan equation; normal vector Ny(0,-1,trackparam[3])
	  track T orientation vector u  : u = Nx * Ny
	  d(H,T) = || vec(BH) * u || / || u || where B is a point from the track
	*/

			//Nx,Ny : plans containing the track
			CLHEP::Hep3Vector Nx(-1., 0., t->getTrackParameters()[1]);
			CLHEP::Hep3Vector Ny(0., -1., t->getTrackParameters()[3]);
			//u : track orientation vector
			CLHEP::Hep3Vector u=Nx.cross(Ny);

			//B : a track point
			CLHEP::Hep3Vector B(t->getTrackParameters()[0],t->getTrackParameters()[2],0.);

			CLHEP::Hep3Vector H=s->getPosition();
			if(u.mag()>0)
				return ((B-H).cross(u)).mag()/u.mag();
			else{
				std::cout << "ORIENTATION VECTOR u IS NULL => return -10 ( float getDistance(S* s,caloobject::CaloTrack *t) )" << std::endl;
				return -10;
			}
		}
		float getDistanceInLayer(S* s,caloobject::CaloTrack *t)
		{
			CLHEP::Hep3Vector impact=t->expectedTrackProjection(s->getPosition().z());
			return (impact-s->getPosition()).mag();
		}
};

template <typename S>
class Distance<S, float*>
{
	public :
		Distance(){;}
		~Distance(){;}

		float getDistance(S s,float* trackParams)
		{
			/*
	  point s(x,y,z)
	  track T : x_t = trackparam[0] + trackparam[1]*z_t => plan equation; normal vector Nx(-1,0,trackparam[1])
	  track T : y_t = trackparam[2] + trackparam[3]*z_t => plan equation; normal vector Ny(0,-1,trackparam[3])
	  track T orientation vector u  : u = Nx * Ny
	  d(s,T) = || vec(Bs) * u || / || u || where B is a point from the track
	*/

			//Nx,Ny : plans containing the track
			CLHEP::Hep3Vector Nx(-1., 0., trackParams[1]);
			CLHEP::Hep3Vector Ny(0., -1., trackParams[3]);
			//u : track orientation vector
			CLHEP::Hep3Vector u=Nx.cross(Ny);

			//B : a track point
			CLHEP::Hep3Vector B(trackParams[0],trackParams[2],0.);

			if(u.mag()>0)
				return ((B-s).cross(u)).mag()/u.mag();
			else{
				std::cout << "ORIENTATION VECTOR u IS NULL => return -10  ( float getDistance(S s,float* trackParams) ) " << std::endl;
				return -10;
			}
		}
};
}

#endif
