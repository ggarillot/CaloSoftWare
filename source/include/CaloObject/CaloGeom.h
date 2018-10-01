#ifndef CALOGEOM_HH
#define CALOGEOM_HH

#include <iostream>

namespace caloobject
{ 

struct GeomParameterSetting
{
		int nLayers = 48 ;
		int nPixelsPerLayer = 96 ;
		int nPixelsX = 96 ;
		int nPixelsY = 96 ;
		float xmin = -500.0f ;
		float xmax = 500.0f ;
		float ymin = -500.0f ;
		float ymax = 500.0f ;
		float pixelSize = 10.408f ;
		float layerGap = 26.131f ;
		float firstLayerZ = 10.0f ;
} ;

}// namespace caloobject

#endif
