#ifndef Pad_h
#define Pad_h

#include <limits>
#include <map>
#include <vector>
#include <array>

#include <CLHEP/Vector/ThreeVector.h>

#include <CaloObject/CaloCluster.h>


namespace caloobject
{

class Asic ;

class Pad
{
	public :
		Pad(int _id) ;
		virtual ~Pad() = default ;

		void setAsic(caloobject::Asic* _asic) { asic = _asic ; }
		void setPosition(CLHEP::Hep3Vector pos) { position = pos ; }

		void setThresholds(const std::vector<double>& thr) ;

		int getNTrack() const { return nTracks ; }
		double getEfficiency() const { return efficiencies.at(0) ; }
		virtual const std::vector<double>& getEfficiencies() { updateEfficiencies() ; return efficiencies ; }
		virtual std::vector<double> getEfficienciesError() const ;
		virtual std::array< std::vector<double> , 2> getEfficienciesBound() const ;

		virtual const std::vector<double>& getMultiplicities() { updateMultiplicities() ; return multiplicities ; }
		virtual std::vector<double> getMultiplicitiesError() const ;

		inline int getID() const { return id ; }
		inline const CLHEP::Hep3Vector &getPosition() const { return position ; }


		virtual void update(CaloCluster2D* cluster = nullptr) ;

		void reset() ;

		Pad(const Pad &toCopy) = delete ;
		void operator=(const Pad &toCopy) = delete ;

	protected :
		virtual void updateEfficiencies() ;
		virtual void updateMultiplicities() ;

	protected :
		int id ;
		caloobject::Asic* asic = nullptr ;

		int nTracks = 0 ;
		std::vector<int> nDetected = {} ;

		std::vector<double> thresholds = {} ;
		std::vector<double> efficiencies = {} ;


		std::vector<double> multiSumVec = {} ;
		std::vector<double> multiSquareSumVec = {} ;

		std::vector<double> multiplicities = {} ;

		CLHEP::Hep3Vector position ;

} ;

class SDHCALPad : public Pad
{
	public :
		SDHCALPad(int _id) ;
		virtual ~SDHCALPad() = default ;

	protected :

} ;

} //namespace caloobject

typedef std::map<int , caloobject::Pad*> PadMap ;


#endif //Pad_h
