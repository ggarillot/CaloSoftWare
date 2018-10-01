#include "Algorithm/LinearFit3D.h"

namespace algorithm
{

LinearFit3D::LinearFit3D( std::vector<CLHEP::Hep3Vector>& pos , std::vector<int>& clustersSize)
	: _positions(pos) , _clusterSize(clustersSize)
{
	float xsum = 0.0;
	float ysum = 0.0;
	float zsum = 0.0;
	float zzsum = 0.0;
	float xzsum = 0.0;
	float yzsum = 0.0;

	for ( unsigned int i=0 ; i<_positions.size() ; i++ )
	{

		//for equation 1
		zsum = zsum + _positions.at(i).z();
		xsum = xsum + _positions.at(i).x();
		zzsum = zzsum + (_positions.at(i).z()*_positions.at(i).z());
		xzsum = xzsum + _positions.at(i).x()*_positions.at(i).z();

		//for equation 2
		ysum = ysum + _positions.at(i).y();
		yzsum = yzsum + _positions.at(i).y()*_positions.at(i).z();

	}

	float A1 = zsum;
	float B1 = _positions.size();
	float C1 = xsum;
	float D1 = zzsum;
	float E1 = xzsum;

	float C2 = ysum;
	float E2 = yzsum;

	_params[0] = (D1*C1-E1*A1)/(B1*D1-A1*A1);
	_params[1] = (E1*B1-C1*A1)/(B1*D1-A1*A1);
	_params[2] = (D1*C2-E2*A1)/(B1*D1-A1*A1);
	_params[3] = (E2*B1-C2*A1)/(B1*D1-A1*A1);

	_paramsError[0] = std::sqrt(D1/(B1*D1-A1*A1));
	_paramsError[1] = std::sqrt(B1/(B1*D1-A1*A1));
	_paramsError[2] = std::sqrt(D1/(B1*D1-A1*A1));
	_paramsError[3] = std::sqrt(B1/(B1*D1-A1*A1));

	ComputeChi2();
}

void LinearFit3D::ComputeChi2()
{

	_chi2 = 0;
	algorithm::Distance<CLHEP::Hep3Vector,float *> dist;
	for( unsigned int i=0 ; i<_positions.size() ; i++ )
	{
		float d=dist.getDistance(_positions.at(i),_params);
		float mult = _clusterSize.at(i);
		float err = mult/sqrt(12/100.);
		_chi2 += (d/err)*(d/err);
	}
	_chi2=_chi2/(_positions.size()-1);
}

} //namespace algorithm
