#include "Algorithm/AsicKeyFinder.h"

namespace algorithm
{

int AsicKeyFinder::BuildAsicKey(int i, int j, int k)
{
	if( settings.keyFactors.size() != 3 )
	{
		std::cout << "Asic key factors vector has a wrong size : it should be 3 for AsicKeyFinder::BuildAsicKey(int i, int j, int k) >> THROW" << std::endl ;
		throw ;
	}

	return i*settings.keyFactors.at(0) + j*settings.keyFactors.at(1) + k*settings.keyFactors.at(2) ;
}

int AsicKeyFinder::FindAsicKey(CLHEP::Hep3Vector vec)
{
	if( settings.setup == sdhcal )
	{
		int I = static_cast<int>( round( vec.x()/settings.padSize ) ) ;
		int J = static_cast<int>( round( vec.y()/settings.padSize ) ) ;
		int K = static_cast<int>( round( vec.z()/settings.layerGap ) - 1 ) ;
		if( I>settings.nPadX || J>settings.nPadY || I<0 || J<0 )
		{
			if( settings.printDebug )
				std::cout << "Impossible to find the Asic key with AsicKeyFinder::FindAsicKey(CLHEP::Hep3Vector vec) with vec = " << vec << "\t return -1" << std::endl;
			return -1 ;
		}
		return
				K*settings.keyFactors.at(2) +
				((J-1)/settings.asicNPad)*settings.keyFactors.at(1) +
				((I-1)/settings.asicNPad)*settings.keyFactors.at(0) ;
	}
	else if(settings.setup==siwecal_nov2015)
	{
		int I = ( vec.x() > 0. ) ? static_cast<int>( round( (vec.x()-settings.ecalGap-settings.padSize)/settings.padSize ) ) : static_cast<int>( round( (vec.x()+settings.ecalGap+settings.padSize)/settings.padSize ) ) ;
		int J = ( vec.y() > 0. ) ? static_cast<int>( round( (vec.y()-settings.ecalGap-settings.padSize)/settings.padSize ) ) : static_cast<int>( round( (vec.y()+settings.ecalGap+settings.padSize)/settings.padSize ) ) ;

		int K = static_cast<int>( round( vec.z()/settings.layerGap ) ) ;
		if(I>settings.nPadX || J>settings.nPadY || I<0 || J<0)
		{
			if( settings.printDebug )
				std::cout << "Impossible to find the Asic key with AsicKeyFinder::FindAsicKey(CLHEP::Hep3Vector vec) with vec = " << vec << "\t return -1" << std::endl;
			return -1 ;
		}
		return
				K*settings.keyFactors.at(2) +
				((J-1)/settings.asicNPad)*settings.keyFactors.at(1) +
				((I-1)/settings.asicNPad)*settings.keyFactors.at(0) ;
	}
	else
		return 0 ;
}

}
