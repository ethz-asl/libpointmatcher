// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2018,
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
#include "OctreeGrid.h"

//Helper function
//FIXME: to put in DataPoints implementation
template<typename T>
void swapCols(typename PointMatcher<T>::DataPoints& self, 
	typename PointMatcher<T>::DataPoints::Index iCol, 
	typename PointMatcher<T>::DataPoints::Index jCol)
{
	//Switch columns j and i
	self.features.col(iCol).swap(self.features.col(jCol));
	if (self.descriptors.cols() > 0)
		self.descriptors.col(iCol).swap(self.descriptors.col(jCol));
	if (self.times.cols() > 0)
		self.times.col(iCol).swap(self.times.col(jCol));
}

//Define Visitor classes to apply processing
template<typename T>
OctreeGridDataPointsFilter<T>::FirstPtsSampler::FirstPtsSampler(DataPoints& dp):
	idx{0}, pts(dp) 
{
}

template <typename T>
template<template<typename> typename Tree>
bool OctreeGridDataPointsFilter<T>::FirstPtsSampler::operator()(Tree<T>& oc)
{
	if(oc.isLeaf() and not oc.isEmpty())
	{			
		auto* data = oc.getData();	
		const auto& d = (*data)[0];
		
		std::size_t j = d;
		
		//retrieve index from lookup table if sampling in already switched element
		if(std::size_t(d)<idx)
			j = mapidx[d];
			
		//Switch columns j and idx
		swapCols<T>(pts, idx, j);
				
		//Maintain new index position	
		mapidx[idx] = j;
		//Update index
		++idx;		
	}
	
	return true;
}
 
template <typename T>
bool OctreeGridDataPointsFilter<T>::FirstPtsSampler::finalize()
{
	//Resize point cloud
	pts.conservativeResize(idx);
	//Reset param
	idx=0;
	return true;
}


template<typename T>
OctreeGridDataPointsFilter<T>::RandomPtsSampler::RandomPtsSampler(DataPoints& dp) 
	: OctreeGridDataPointsFilter<T>::FirstPtsSampler{dp}, seed{1}
{
	std::srand(seed);
}
template<typename T>
OctreeGridDataPointsFilter<T>::RandomPtsSampler::RandomPtsSampler(DataPoints& dp, const std::size_t seed_): 
	OctreeGridDataPointsFilter<T>::FirstPtsSampler{dp}, seed{seed_}
{
	std::srand(seed);
}
template<typename T>
template<template<typename> typename Tree>
bool OctreeGridDataPointsFilter<T>::RandomPtsSampler::operator()(Tree<T>& oc)
{
	if(oc.isLeaf() and not oc.isEmpty())
	{			
		auto* data = oc.getData();
		const std::size_t nbData = (*data).size() - 1;
		const std::size_t randId = 
			static_cast<std::size_t>( nbData * 
				(static_cast<float>(std::rand()/static_cast<float>(RAND_MAX))));
				
		const auto& d = (*data)[randId];
		
		std::size_t j = d;
		
		//retrieve index from lookup table if sampling in already switched element
		if(std::size_t(d)<idx)
			j = mapidx[d];
			
		//Switch columns j and idx
		swapCols<T>(pts, idx, j);	
		
		//Maintain new index position	
		mapidx[idx] = j;
		//Update index
		++idx;		
	}
	
	return true;
}
	
template<typename T>
bool OctreeGridDataPointsFilter<T>::RandomPtsSampler::finalize()
{
	bool ret = FirstPtsSampler::finalize();
	//Reset seed
	std::srand(seed);
	
	return ret;			
}

template<typename T>
OctreeGridDataPointsFilter<T>::CentroidSampler::CentroidSampler(DataPoints& dp):
	OctreeGridDataPointsFilter<T>::FirstPtsSampler{dp}
{
}
	
template<typename T>
template<template<typename> typename Tree>
bool OctreeGridDataPointsFilter<T>::CentroidSampler::operator()(Tree<T>& oc)
{
	if(oc.isLeaf() and not oc.isEmpty())
	{			
		const int featDim(pts.features.rows());
		const int descDim(pts.descriptors.rows());
		const int timeDim(pts.times.rows());
		
		auto* data = oc.getData();
		const std::size_t nbData = (*data).size();
			
		const auto& d = (*data)[0]; //get first data
		std::size_t j = d; //j contains real index of first point
		
		//retrieve index from lookup table if sampling in already switched element
		if(std::size_t(d)<idx)
			j = mapidx[d];
			
		//We sum all the data in the first data
		for(std::size_t id=1;id<nbData;++id)
		{
			//get current idx
			const auto& curId = (*data)[id];
			std::size_t i = curId; //i contains real index
			
			//retrieve index from lookup table if sampling in already switched element
			if(std::size_t(curId)<idx)
				i = mapidx[curId];
			
			for (int f = 0; f < (featDim - 1); ++f)
				pts.features(f,j) += pts.features(f,i);
			
			if (pts.descriptors.cols() > 0)
				for (int d = 0; d < descDim; ++d)
					pts.descriptors(d,j) += pts.descriptors(d,i);
			
			if (pts.times.cols() > 0)
				for (int t = 0; t < timeDim; ++t)
					pts.times(t,j) += pts.times(t,i);	
		}
		
		// Normalize sums to get centroid (average)
		for (int f = 0; f < (featDim - 1); ++f)
			pts.features(f,j) /= T(nbData);
		
		if (pts.descriptors.cols() > 0)
			for (int d = 0; d < descDim; ++d)
				pts.descriptors(d,j) /= T(nbData);
		
		if (pts.times.cols() > 0)
			for (int t = 0; t < timeDim; ++t)
				pts.times(t,j) /= T(nbData);	
								
		//Switch columns j and idx
		swapCols<T>(pts, idx, j);
		
		//Maintain new index position	
		mapidx[idx] = j;
		//Update index
		++idx;		
	}
	
	return true;
}

// OctreeGridDataPointsFilter
template <typename T>
OctreeGridDataPointsFilter<T>::OctreeGridDataPointsFilter(const Parameters& params):
	PointMatcher<T>::DataPointsFilter("OctreeGridDataPointsFilter", 
		OctreeGridDataPointsFilter::availableParameters(), params),
	parallel_build{Parametrizable::get<bool>("buildParallel")},
	maxPointByNode{Parametrizable::get<std::size_t>("maxPointByNode")},
	maxSizeByNode{Parametrizable::get<T>("maxSizeByNode")}
{
	try 
	{
		const int bm = Parametrizable::get<int>("buildMethod");
		buildMethod = BuildMethod(bm);
		const int sm = Parametrizable::get<int>("samplingMethod");
		samplingMethod = SamplingMethod(sm);
	}
	catch (const InvalidParameter& e) 
	{
		buildMethod = BuildMethod::MAX_POINT;
		samplingMethod = SamplingMethod::FIRST_PTS;
	}
}

template <typename T>
typename PointMatcher<T>::DataPoints
OctreeGridDataPointsFilter<T>::filter(const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

template <typename T>
void OctreeGridDataPointsFilter<T>::inPlaceFilter(DataPoints& cloud)
{
	assert(cloud.features.rows() == 4 || cloud.features.rows() == 3 );
	
	if(cloud.features.rows() == 4) //3D case
		applySampler<Octree>(cloud);
	else if (cloud.features.rows() == 3) //2D case
		applySampler<Quadtree>(cloud);
}

template <typename T>
template<template<typename> typename Tree>
void OctreeGridDataPointsFilter<T>::applySampler(DataPoints& cloud){
	Tree<T> oc{};
	
	switch(buildMethod) 
	{
		case BuildMethod::MAX_POINT:
		{
			oc.build(cloud, maxPointByNode, parallel_build);
			break;
		}
		case BuildMethod::MAX_SIZE:
		{
			oc.build(cloud, maxSizeByNode, parallel_build);
			break;
		}
	}
	
	switch(samplingMethod)
	{
		case SamplingMethod::FIRST_PTS:
		{
			FirstPtsSampler sampler(cloud);
			oc.visit(sampler);
			sampler.finalize();
			break;
		}
		case SamplingMethod::RAND_PTS:
		{
			RandomPtsSampler sampler(cloud); //FIXME: add seed parameter
			oc.visit(sampler);
			sampler.finalize();
			break;
		}
		case SamplingMethod::CENTROID:
		{
			CentroidSampler sampler(cloud);
			oc.visit(sampler);
			sampler.finalize();
			break;
		}
	}	
}

template struct OctreeGridDataPointsFilter<float>;
template struct OctreeGridDataPointsFilter<double>;
