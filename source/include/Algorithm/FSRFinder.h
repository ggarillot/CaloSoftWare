#ifndef FSRFINDER_HH
#define FSRFINDER_HH

#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloTrack.h"
#include "CaloObject/CaloGeom.h"
#include "CaloObject/CaloCluster.h"
#include "Algorithm/Cluster3D.h"
#include "Algorithm/ShowerAnalyser.h"

namespace algorithm
{ 

struct FSRFinderParameterSetting
{
		float minTransverseRatio = 0.02f ;
		float minMeanEdepPerCell = 0.001f ;
		float minRMSEdepPerCell = 0.001f ;
		float minEdepPerX0 = 0.1f ;
		float minSurroundedEdepFraction = 0.8f ;
		bool printDebug = false ;
		std::vector<float> minEnergyFraction ;
		caloobject::GeomParameterSetting geometry ;
		algorithm::cluster3DParameterSetting clustering ;
//		algorithm::ShowerAnalyserParameterSetting showerAnalyser ;

		FSRFinderParameterSetting()
			: minEnergyFraction() , geometry() , clustering() //, showerAnalyser()
		{;}
} ;

class FSRFinder
{
	public :
		FSRFinder()
			: settings()
		{;}
		~FSRFinder(){;}

		void Run(std::vector<caloobject::CaloHit*> hits, caloobject::CaloTrack* track, caloobject::CaloCluster3D* gammaCluster) ;
		inline void SetFSRFinderParameterSetting(FSRFinderParameterSetting val) {settings = val ; }

	private :
		void rmTrackedHitsFromList(HitVec& hits, caloobject::CaloTrack* track) ;
		bool isGammaCluster( caloobject::CaloCluster3D* fsrCandidate ) ;
		FSRFinderParameterSetting settings ;
} ;

} //namespace algorithm
#endif
