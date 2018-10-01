#ifndef Layer_h
#define Layer_h

#include <iostream>
#include <limits>
#include <vector>
#include <array>

#include <CaloObject/Asic.h>

namespace caloobject
{

class Layer
{
	public :
		Layer(int _id) ;
		virtual ~Layer() = default ;

		void setPosition(CLHEP::Hep3Vector pos) { position = pos ; }

		void setThresholds(const std::vector<double>& thr) ;

		int getNTrack() const { return nTracks ; }
		double getEfficiency() const { return efficiencies.at(0) ; }
		virtual const std::vector<double>& getEfficiencies() { updateEfficiencies() ; return efficiencies ; }
		virtual std::vector<double> getEfficienciesError() const ;
		virtual std::array< std::vector<double> , 2> getEfficienciesBound() const ;


		virtual const std::vector<double>& getMultiplicities() { updateMultiplicities() ; return multiplicities ; }
		virtual std::vector<double> getMultiplicitiesError() const ;

		virtual Asic* findAsic(const CLHEP::Hep3Vector& pos) const = 0 ;


		inline int getID() const { return id ; }
		inline const CLHEP::Hep3Vector& getPosition() const { return position ; }
		inline const AsicMap& getAsics() const { return asics ; }


		virtual void update(const CLHEP::Hep3Vector& impactPos , CaloCluster2D* cluster = nullptr) ;

		virtual void buildAsics() = 0 ;

		void reset() ;

		Layer(const Layer &toCopy) = delete ;
		void operator=(const Layer &toCopy) = delete ;

	protected :
		virtual void updateEfficiencies() ;
		virtual void updateMultiplicities() ;

	protected :

		int id ;

		int nTracks = 0 ;
		std::vector<int> nDetected = {} ;

		std::vector<double> thresholds = {} ;
		std::vector<double> efficiencies = {} ;


		std::vector<double> multiSumVec = {} ;
		std::vector<double> multiSquareSumVec = {} ;

		std::vector<double> multiplicities = {} ;


		CLHEP::Hep3Vector position ; //top left corner position

		AsicMap asics ;
} ;

class SDHCALLayer : public Layer
{
	public :
		SDHCALLayer(int _id , int difL , int difC , int difR) ;
		virtual ~SDHCALLayer() = default ;

		static const int asicTab[12][4] ;
		static const int iAsicTab[48] ;
		static const int jAsicTab[48] ;

		virtual Asic* findAsic(const CLHEP::Hep3Vector& pos) const ;

		virtual void buildAsics() ;

	protected :

		int difID[3] ;
} ;

} //namespace caloobject


#endif //Layer_h
