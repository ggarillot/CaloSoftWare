#ifndef SHOWER_HH
#define SHOWER_HH

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <cstring>

#include <Algorithm/InteractionFinder.h>
#include <CaloObject/CaloGeom.h>
#include <CaloObject/CaloHit.h>
#include <CaloObject/CaloCluster.h>
#include <CLHEP/Vector/ThreeVector.h>

namespace algorithm 
{ 

class ShowerAnalyser ;

}  //namespace algorithm

namespace caloobject
{

class Shower
{
	public :
		Shower(const HitVec& vec) ;
		Shower(const Cluster2DVec& vec) ;
		Shower(const caloobject::CaloCluster3D* cluster) ;
		virtual ~Shower() ;


		void computePCA() ;
		void computeThrust() ;
		void computeInteraction() ;

		void computeProfile() ;


		inline const std::vector<caloobject::CaloHit*>& getHits() const { return hits ; }
		inline const std::vector<caloobject::CaloCluster2D*>& getClusters() const { return clusters ; }

		inline unsigned int getNhit() const { return static_cast<unsigned int>( hits.size() ) ; }

		inline const std::set<int>& getFiredLayers() const { return firedLayers ; }
		inline const std::set<int>& getInteractingLayers() const { return interactingLayers ; }

		inline caloobject::CaloCluster2D* getFirstIntCluster() const { return firstIntCluster ; }
		inline CLHEP::Hep3Vector getStartingPosition() const { return startingPosition ; }

		inline int getLastClusterLayer() const { return lastClusterLayer ; }

		inline std::vector<float> getThrust() const { return thrust ; }
		inline float getReconstructedCosTheta() const { return reconstructedCosTheta ; }
		inline float getTransverseRatio() const { return transverseRatio ; }
		inline float getEta() const { return eta ; }
		inline float getPhi() const { return phi ; }

//		inline float getShowerMax() const { return showerMax ; }

		inline const std::vector<double>& getLongitudinalProfile() const { return longitudinalProfile ; }
		inline const std::vector<double>& getTransverseProfile() const { return transverseProfile ; }
		inline const std::vector<double>& getDistanceToAxis() const { return distanceToAxis ; }
		inline const std::vector<double>& getHitTimes() const { return hitTimes ; }




		inline void setInteractionSettings(algorithm::InteractionFinderParameterSetting s) { interactionSettings = s ; }
		inline void setGeometrySettings(caloobject::GeomParameterSetting s) { geometrySettings = s ; }

		Shower(const Shower& toCopy) = delete ;
		void operator=(const Shower& toCopy) = delete ;

	protected :
		Shower() ;

		virtual double valueOfHitForProfile(const caloobject::CaloHit* hit) const = 0 ;

		algorithm::InteractionFinderParameterSetting interactionSettings ;
		caloobject::GeomParameterSetting geometrySettings ;

		HitVec hits = {} ;
		Cluster2DVec clusters = {} ;
		caloobject::CaloCluster2D* firstIntCluster = nullptr ;
		CLHEP::Hep3Vector startingPosition = CLHEP::Hep3Vector(0,0,0) ;

		int lastClusterLayer = 0 ;

		std::vector<double> hitTimes = {} ;
		std::set<int> firedLayers = {} ;
		std::set<int> interactingLayers = {} ;

		std::vector<float> thrust = {} ;
		float reconstructedCosTheta = 0.0f ;
		float transverseRatio = 0.0f ;
		float eta = 0.0f ;
		float phi = 0.0f ;

		std::vector<double> longitudinalProfile = {} ;
		std::vector<double> transverseProfile = {} ;
		std::vector<double> distanceToAxis = {} ;


		bool isThrustComputed = false ;
		bool isInteractionComputed = false ;

} ;

class AnalogShower : public Shower
{

	public :
		AnalogShower(const HitVec& vec) ;
		AnalogShower(const Cluster2DVec& vec) ;
		AnalogShower(const caloobject::CaloCluster3D* cluster) ;
		virtual ~AnalogShower() ;


		inline float getEnergy() const { return energy ; }
		inline float getEdep() const { return edep ; }
		inline float getMeanEdep() const { return meanEdep ; }
		inline float getRMSEdep() const { return rmsEdep ; }


//		inline float getEdepAtMax(){return edepAtMax;}
		inline const std::vector<double>& getEdepPerCell() const { return edepPerCell ; }

		inline const std::vector<double>& getClustersEnergy() const { return clustersEnergy ; }


	protected :

		virtual double valueOfHitForProfile(const caloobject::CaloHit* hit) const ;

		float energy = 0.0f ;
		float edep = 0.0f ;
		std::vector<double> edepPerCell = {} ;

		float meanEdep = 0.0f ;
		float rmsEdep = 0.0f ;
//		float edepAtMax = 0.0f ;

		std::vector<double> clustersEnergy = {} ;
} ;

class DigitalShower : public Shower
{

	public :
		DigitalShower(const Cluster2DVec& vec) ;
		~DigitalShower() { ; }

		inline std::vector<int> getNHits() const { return nHits ; }

	protected :

		virtual double valueOfHitForProfile(const caloobject::CaloHit*) const ;

		std::vector<int> nHits = {0,0,0,0} ;
} ;

} //namespace caloobject

#endif
