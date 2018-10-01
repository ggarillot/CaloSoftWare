#ifndef CLUSTERINGHELPER_HH
#define CLUSTERINGHELPER_HH

#include "CaloObject/CaloCluster.h"
#include <cmath>
#include <algorithm>

namespace algorithm
{

struct ClusteringHelperParameterSetting
{
		float longitudinalDistance = 100.0f ;
		float transversalDistance = 200.0f ;
} ;

class ClusteringHelper
{
	public:
		ClusteringHelper()
			: settings()
		{;}
		~ClusteringHelper(){;}
		inline void SetClusteringHelperParameterSetting(ClusteringHelperParameterSetting params) {settings = params ;}
		static bool SortClusterByLayer(caloobject::CaloCluster2D* a , caloobject::CaloCluster2D* b) { return a->getPosition().z() < b->getPosition().z() ; }
		inline bool IsIsolatedCluster(caloobject::CaloCluster2D* cl , std::vector<caloobject::CaloCluster2D*>& list)
		{
			for(std::vector<caloobject::CaloCluster2D*>::iterator it=list.begin(); it!=list.end(); ++it)
			{
				if( cl == (*it) )
					continue ;
				if( fabs(cl->getPosition().z()-(*it)->getPosition().z())<settings.longitudinalDistance &&
					fabs(cl->getPosition().y()-(*it)->getPosition().y())<settings.transversalDistance &&
					fabs(cl->getPosition().x()-(*it)->getPosition().x())<settings.transversalDistance )
				{
					return false ;
				}
			}
			return true ;
		}
		ClusteringHelperParameterSetting settings ;

} ;

template <typename T, typename S>
class SortClusterBySize
{
	public :
		SortClusterBySize(){;}
		~SortClusterBySize(){;}

		static bool sort(T* t, S* s){ return t->getHits().size() < s->getHits().size() ; }
} ;


} //namespace algorithm
#endif
