#ifndef Asic_h
#define Asic_h

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>
#include <CLHEP/Vector/ThreeVector.h>

#include "CaloObject/Pad.h"

namespace caloobject
{

class Layer ;

class Asic
{
	public :
		Asic(int _id , int _difID = -1) ;
		virtual ~Asic() = default ;

		void setLayer(caloobject::Layer* _layer) { layer = _layer ; }
		void setPosition(CLHEP::Hep3Vector pos) { position = pos ; }

		void setThresholds(const std::vector<double>& thr) ;

		int getNTrack() const { return nTracks ; }
		double getEfficiency() const { return efficiencies.at(0) ; }
		virtual const std::vector<double>& getEfficiencies() { updateEfficiencies() ; return efficiencies ; }
		virtual std::vector<double> getEfficienciesError() const ;
		virtual std::array< std::vector<double> , 2> getEfficienciesBound() const ;

		virtual const std::vector<double>& getMultiplicities() { updateMultiplicities() ; return multiplicities ; }
		virtual std::vector<double> getMultiplicitiesError() const ;


		virtual Pad* findPad(const CLHEP::Hep3Vector& pos) const = 0 ;

		inline int getID() const { return id ; }
		inline int getDifID() const { return difID ; }
		inline const CLHEP::Hep3Vector &getPosition() const { return position ; }
		inline const PadMap& getPads() const { return pads ; }

		virtual void update(const CLHEP::Hep3Vector& impactPos , CaloCluster2D* cluster = nullptr) ;

		virtual void buildPads() = 0 ;


		void reset() ;


		Asic(const Asic &toCopy) = delete ;
		void operator=(const Asic &toCopy) = delete ;

	protected :
		virtual void updateEfficiencies() ;
		virtual void updateMultiplicities() ;

	protected :

		int id ;
		int difID ;

		caloobject::Layer* layer = nullptr ;

		int nTracks = 0 ;
		std::vector<int> nDetected = {} ;

		std::vector<double> thresholds = {} ;
		std::vector<double> efficiencies = {} ;


		std::vector<double> multiSumVec = {} ;
		std::vector<double> multiSquareSumVec = {} ;

		std::vector<double> multiplicities = {} ;

		CLHEP::Hep3Vector position ; //top left corner position

		PadMap pads ;
} ;

class SDHCALAsic : public Asic
{
	public :
		SDHCALAsic(int _id , int _difID) ;
		virtual ~SDHCALAsic() = default ;

		static const int padTab[8][8] ;
		static const int iPadTab[64] ;
		static const int jPadTab[64] ;

		virtual Pad* findPad(const CLHEP::Hep3Vector& pos) const ;

		virtual void buildPads() ;


	protected :

} ;


} //namespace caloobject

typedef std::map<int , caloobject::Asic*> AsicMap ;

#endif //Asic_h
