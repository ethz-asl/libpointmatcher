// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2012,
François Pomerleau and Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the authors at <f dot pomerleau at gmail dot com> and
<stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <iostream>
#include <pointmatcher/PointMatcher.h>

#include "Eigen/QR"
#include "Eigen/Eigenvalues"

#include "ErrorMinimizersImpl.h"
#include "PointMatcherPrivate.h"
#include "Functions.h"

using namespace Eigen;
using namespace std;
using namespace PointMatcherSupport;

typedef PointMatcherSupport::Parametrizable Parametrizable;
typedef PointMatcherSupport::Parametrizable P;
typedef Parametrizable::Parameters Parameters;
typedef Parametrizable::ParameterDoc ParameterDoc;
typedef Parametrizable::ParametersDoc ParametersDoc;

template<typename T>
PointToGaussianErrorMinimizer<T>::PointToGaussianErrorMinimizer(const Parameters& params):
	PointToPlaneErrorMinimizer<T>(PointToGaussianErrorMinimizer::availableParameters(), params),
	//confidenceInPenalties(Parametrizable::get<T>("confidenceInPenalties")),
	force2D(Parametrizable::get<bool>("force2D")),
    force4DOF(Parametrizable::get<T>("force4DOF")),
    noiseSensor(Parametrizable::get<float>("noiseSensor"))
{
	if(force2D)
	{
		if (force4DOF)
		{
			throw PointMatcherSupport::ConfigurationError("Force 2D cannot be used together with force4DOF.");
		}
		else
		{
			LOG_INFO_STREAM("PointMatcher::PointToGaussianErrorMinimizer - minimization will be in 2D.");
		}
	}
	else if(force4DOF)
	{
		LOG_INFO_STREAM("PointMatcher::PointToGaussianErrorMinimizer - minimization will be in 4-DOF (yaw,x,y,z).");
	}
	else
	{
		LOG_INFO_STREAM("PointMatcher::PointToGaussianErrorMinimizer - minimization will be in full 6DOF.");
	}
}

template<typename T>
typename PointMatcher<T>::TransformationParameters PointToGaussianErrorMinimizer<T>::compute(const ErrorElements& mPts_const)
{
	auto errElements = convertCovariancesToNormals(mPts_const);
	typename PointMatcher<T>::TransformationParameters out = PointToPlaneErrorMinimizer<T>::compute_in_place(errElements);

	return out;
}


template<typename T>
typename PointToGaussianErrorMinimizer<T>::ErrorElements PointToGaussianErrorMinimizer<T>::convertCovariancesToNormals(const ErrorElements& mPts_const) const
{
	typedef typename DataPoints::View View;
	typedef typename DataPoints::ConstView ConstView;
	typedef typename PointMatcher<T>::ConvergenceError ConvergenceError;

	ErrorElements mPts = mPts_const;
	const size_t dim(mPts_const.reference.features.rows() - 1);
	const size_t newSize(dim * (mPts.weights.cols())); // + nbPenalty

	mPts.reference.conservativeResize(newSize);
	mPts.reading.conservativeResize(newSize);
	mPts.weights.conservativeResize(Eigen::NoChange, newSize);

	if (!mPts.reference.descriptorExists("normals")) {
		Labels cloudLabels;
		cloudLabels.push_back(Label("normals", dim));
		// Reserve memory
		mPts.reference.allocateDescriptors(cloudLabels);
//		mPts.reference.getDescriptorViewByName("normals").setZero();
	}
	View normals = mPts.reference.getDescriptorViewByName("normals");
	ConstView eigVectors = mPts_const.reference.getDescriptorViewByName("eigVectors");
	ConstView eigValues = mPts_const.reference.getDescriptorViewByName("eigValues");

    if ((eigValues.array() < -noiseSensor*noiseSensor).any()) {
        throw ConvergenceError("PointToGaussian(): Some of the eigen values are negative.");
    }

    for (long i = 0; i < mPts_const.reference.features.cols(); ++i) {
        mPts.reading.features.block(0, dim * i, dim + 1, dim) = mPts_const.reading.features.col(i).replicate(1, dim);
        mPts.reference.features.block(0, dim * i, dim + 1, dim) = mPts_const.reference.features.col(i).replicate(1, dim);

        // Convert a vector into a matrix
        const Matrix matEigVectors = Eigen::Map<const Matrix>(eigVectors.col(i).data(), dim, dim);

        normals.block(0, dim * i, dim, dim) = matEigVectors.transpose();

        // A eigen value of zero will have infinite weight, that's why we add a noise constant for the weight
        // Same solution as Bosse did with his Zebedee platefrom
        mPts.weights.block(0, dim * i, 1, dim) = (pow(eigValues.col(i).array()+noiseSensor*noiseSensor,1/2)).inverse().transpose();
    }

	return mPts;
}

template<typename T>
T PointToGaussianErrorMinimizer<T>::getResidualError(
				const DataPoints& filteredReading,
				const DataPoints& filteredReference,
				const OutlierWeights& outlierWeights,
				const Matches& matches) const
{
	assert(matches.ids.rows() > 0);

	// Fetch paired points
	const typename ErrorMinimizer::ErrorElements mPts(filteredReading, filteredReference, outlierWeights, matches);
	auto mPtsNormals = convertCovariancesToNormals(mPts);

	return PointToGaussianErrorMinimizer::computeResidualError(mPtsNormals, force2D);
}

template struct PointToGaussianErrorMinimizer<float>;
template struct PointToGaussianErrorMinimizer<double>;
