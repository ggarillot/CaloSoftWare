#include "CaloObject/CaloTrack.h"

namespace caloobject
{

CaloTrack::CaloTrack(Cluster2DVec& vec)
	: clusters(vec) ,
	  _params()
{

}

} //namespace caloobject
