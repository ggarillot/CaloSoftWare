#include "Algorithm/PCA.h"
#include "numeric"

namespace algorithm
{

PCA::PCA(VectorList &vec)
	: _eigenValues() , _eigenVectors()
{
	Eigen::MatrixXd mat(vec.at(0).size(),vec.size()) ;

	unsigned int mat_column;
	unsigned int mat_line;
	double mean;

	for(VectorList::iterator it=vec.begin(); it!=vec.end(); ++it)
	{
		if( (*it).size()!=vec.at(0).size() )
		{
			std::cout << "ERROR IN PCA : VECTORS SHOULD HAVE THE SAME SIZE" << std::endl;
			abort();
		}
		mean=0;
		for(std::vector<double>::iterator jt=(*it).begin(); jt!=(*it).end(); ++jt)
			mean+= (*jt);
		mean/=(*it).size();

		mat_column=std::distance(vec.begin(),it);
		for(std::vector<double>::iterator jt=(*it).begin(); jt!=(*it).end(); ++jt)
		{
			mat_line=std::distance((*it).begin(),jt);
			mat(mat_line,mat_column)=((*jt)-mean);
		}
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver( mat.transpose()*mat );
	_eigenVectors=eigensolver.eigenvectors();
	for(unsigned int i=0; i<eigensolver.eigenvalues().size(); i++)
		_eigenValues.push_back(eigensolver.eigenvalues()[i]);
}

}//namespace algorithm
