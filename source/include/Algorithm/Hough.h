#ifndef HOUGH_HH
#define HOUGH_HH

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <limits>
#include <algorithm>

#include "CLHEP/Vector/ThreeVector.h"

#include "CaloObject/CaloCluster.h"
#include "CaloObject/CaloTrack.h"
#include "CaloObject/CaloGeom.h"
#include "Algorithm/Distance.h"
#include "Algorithm/Tracking.h"

namespace algorithm
{

struct HoughParameterSetting
{
		unsigned int thetaSteps = 100 ;
		unsigned int minimumNBins = 6 ;
		unsigned int maximumClusterSizeForMip = 4 ;
		int maximumNumberOfNeighboursForMip = 1 ;
		int maximumNumberOfCoreNeighboursForMip = 0 ;
		float transversalDistance = 85.0f ;
		float rhoTolerance = 0.1f ; // 10% (of the pixel size) tolerance for rho to be in the same (theta,rho) HoughBin
		int isolationDistance = 3 ;
		bool printDebug = false ;
		bool useAnalogEnergy = false ;
		float maxEnergy = 0.001f ; // 1 MeV (MIP arround 85 keV)
		caloobject::GeomParameterSetting geometry ;
		HoughParameterSetting()
			: geometry()
		{;}
} ;

enum HoughTag
{
	fMip ,
	fTrack ,
	fIsolated
} ;

struct HoughObject
{
		std::vector<unsigned int> thetas ;
		std::vector<float> rhoXVec ;
		std::vector<float> rhoYVec ;
		HoughTag tag = fMip ;
		caloobject::CaloCluster2D* cluster = nullptr ;

		HoughObject()
			: thetas() , rhoXVec() , rhoYVec()
		{ ; }

		HoughObject(const HoughObject&) = delete ;
		void operator=(const HoughObject&) = delete ;

} ;

struct HoughBin
{
		std::vector<HoughObject*> houghObjects ;
		bool rmTag ;
		unsigned int theta ;
		float rho ;

		HoughBin()
			: houghObjects() , rmTag() , theta() , rho()
		{ ; }
} ;

class Hough
{
	public :
		Hough()
			: settings() , houghObjects()
		{;}
		~Hough(){;}

		inline void SetHoughParameterSetting(HoughParameterSetting params) { settings = params ; }
		void runHough(std::vector<caloobject::CaloCluster2D*>& clusters , std::vector<caloobject::CaloTrack*>& tracks, algorithm::Tracking* algo_Tracking) ;

		void selectNonDensePart( std::vector<caloobject::CaloCluster2D*>& clusters, std::vector<caloobject::CaloCluster2D*>& mipCandidate ) ;

	private :
		void createHoughObjects( std::vector<caloobject::CaloCluster2D*>& clusters ) ;

		std::vector<HoughBin> getHoughBinsFromZX() ;
		HoughBin getBestHoughBinFromZY( HoughBin& inputBin ) ;
		void RemoveIsolatedClusterInHoughBin(HoughBin& a) ;
		void RemoveTrackedObjects(std::vector<HoughBin>& houghBins) ;

		HoughParameterSetting settings ;
		std::vector<HoughObject*> houghObjects ;

		inline bool TestHoughBinSize( HoughBin a ){ return a.houghObjects.size() < settings.minimumNBins ; }
		static bool HasToBeDeleted( HoughBin a ){ return a.rmTag; }
} ;

class SortHoughBinBySize
{
	public:
		SortHoughBinBySize(){;}
		static bool Sort(HoughBin a, HoughBin b){return a.houghObjects.size() > b.houghObjects.size() ; }
} ;

class RemoveTrackedObject
{
	public:
		RemoveTrackedObject(){;}
		static bool rm(HoughObject* a) {return a->tag == fTrack ; }
} ;

} //namespace algorithm

#endif
