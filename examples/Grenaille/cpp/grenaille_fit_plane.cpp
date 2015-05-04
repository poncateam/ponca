/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
\file examples/Grenaille/basic_cpu.h
\brief Basic use of Grenaille with plane fit

\author: Nicolas Mellado, Gautier Ciaudo
*/
#include <cmath>
#include <algorithm>
#include <iostream>

#include "Patate/grenaille.h"
#include "Eigen/Eigen"

#include <vector>

using namespace std;
using namespace Grenaille;


// This class defines the input data format
class MyPoint
{
public:
	enum {Dim = 3};
	typedef double Scalar;
	typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
    typedef Eigen::Matrix<Scalar, Dim+1, 1>	HVectorType;
	typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

	MULTIARCH inline MyPoint(   const VectorType& _pos    = VectorType::Zero(), 
		                        const VectorType& _normal = VectorType::Zero()
                            )
		: m_pos(_pos), m_normal(_normal) {}

	MULTIARCH inline const VectorType& pos()    const { return m_pos; }  
	MULTIARCH inline const VectorType& normal() const { return m_normal; }

	MULTIARCH inline VectorType& pos()    { return m_pos; }  
	MULTIARCH inline VectorType& normal() { return m_normal; }

	static inline MyPoint Random()
    {
		VectorType n = VectorType::Random().normalized();
		VectorType p = n * Eigen::internal::random<Scalar>(0.9,1.1);
		return MyPoint (p, (n + VectorType::Random()*0.1).normalized());
	};

private:
	VectorType m_pos, m_normal;
};

typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;

// Define related structure
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc;
typedef Basket<MyPoint, WeightFunc, CompactPlane, CovariancePlaneFit> CovPlaneFit;

template<typename Fit>
void test_fit(Fit& _fit, vector<MyPoint>& _vecs, const VectorType& _p)
{
	Scalar tmax = 100.0;

	// Set a weighting function instance
	_fit.setWeightFunc(WeightFunc(tmax));  

	// Set the evaluation position
	_fit.init(_p);

	// Iterate over samples and _fit the primitive
	for(vector<MyPoint>::iterator it = _vecs.begin(); it != _vecs.end(); it++)
    {
		_fit.addNeighbor(*it);  
    }

	// The plane fitting is a multipass method
	if(_fit.finalize() == NEED_OTHER_PASS)
    {
        // Iterate over samples and _fit the primitive
        for(vector<MyPoint>::iterator it = _vecs.begin(); it != _vecs.end(); it++)
        {
            _fit.addNeighbor(*it);  
        }

        //finalize fitting
        _fit.finalize();

        //Test if the fitting ended without errors
        if(_fit.isStable())
        {
            cout << "Value of the scalar field at the initial point: " 
                << _p.transpose() 
                << " is equal to " << _fit.potential(_p)
                << endl;

            cout << "It's gradient at this place is equal to: "
                << _fit.primitiveGradient(_p).transpose()
                << endl;

            cout << "The initial point " << _p.transpose()              << endl
                << "Is projected at   " << _fit.project(_p).transpose() << endl;

            cout << "Value of the surface variation: "
                << _fit.surfaceVariation()
                << endl;
        }
    }
}

int main()
{
	// set evaluation point and scale
	VectorType p = VectorType::Random();

	// init input data
	int n = 10000;
	vector<MyPoint> vecs (n);

	for(int k=0; k<n; ++k)
    {
		vecs[k] = MyPoint::Random();
    }

	p = vecs.at(0).pos();

	std::cout << "====================\nCovariancePlaneFit:\n";
	CovPlaneFit fit;
	test_fit(fit, vecs, p);
}
