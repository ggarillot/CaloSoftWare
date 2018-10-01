#ifndef ASICKEYFINDER_HH
#define ASICKEYFINDER_HH

#include <iostream>
#include <vector>

#include <CLHEP/Vector/ThreeVector.h>

namespace algorithm
{

enum TestBeamSetUp
{
	sdhcal ,
	siwecal_nov2015
} ;

struct AsicKeyFinderParameterSetting
{
		std::vector<int> keyFactors = std::vector<int>{1,12,100} ;
		int nPadX = 96 ; // number of pads in x direction per layer
		int nPadY = 96 ; // number of pads in y direction per layer
		int asicNPad = 8 ; // number of pads in x or y direction per asic (assuming a square)
		float layerGap = 26.131f ; // gap size (in mm) between 2 layers
		float padSize = 10.408f ;
		TestBeamSetUp setup = sdhcal ;
		float ecalGap = 0.1f ; //gap between ecal sensors
		bool printDebug = false ;
} ;


class AsicKeyFinder
{
	public :
		AsicKeyFinder()
			:	settings()
		{;}
		~AsicKeyFinder(){;}

		inline void SetAsicKeyFinderParameterSetting(AsicKeyFinderParameterSetting params){ settings = params ; }
		int BuildAsicKey(int i, int j, int k) ;
		int FindAsicKey(CLHEP::Hep3Vector vec) ;

	private :
		AsicKeyFinderParameterSetting settings ;
} ;

} //namespace algorithm

#endif
